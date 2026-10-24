"""Fermionic -> qubit Hamiltonians, Pauli statistics, 1-norms and Z2 tapering.

Spin-orbital convention (interleaved): qubit 2p = orbital p alpha, qubit 2p+1 = orbital p beta.
This matches OpenFermion's `spinorb_from_spatial` and PennyLane's `hf_state`/`excitations`.
"""
from __future__ import annotations

import numpy as np
from openfermion import (InteractionOperator, QubitOperator, bravyi_kitaev, count_qubits,
                         jordan_wigner)
from openfermion.chem.molecular_data import spinorb_from_spatial

from .integrals import ActiveSpace

MAX_EXPLICIT_QUBITS = 24   # hard constraint: no explicit Pauli expansion above this


def interaction_operator(asp: ActiveSpace) -> InteractionOperator:
    """Chemist (pq|rs) -> OpenFermion physicist two-body tensor with the 1/2 factor included."""
    two_body_phys = 0.5 * asp.eri.transpose(0, 2, 3, 1)
    h1s, h2s = spinorb_from_spatial(asp.h1, two_body_phys)
    return InteractionOperator(asp.ecore, h1s, h2s)


def to_qubit_operator(asp: ActiveSpace, mapping: str = "jw") -> QubitOperator:
    n = asp.n_spin_orbitals
    if n > MAX_EXPLICIT_QUBITS:
        raise ValueError(f"{n} qubits > MAX_EXPLICIT_QUBITS={MAX_EXPLICIT_QUBITS}; use lambda_from_integrals")
    op = interaction_operator(asp)
    if mapping == "jw":
        q = jordan_wigner(op)
    elif mapping == "bk":
        q = bravyi_kitaev(op, n_qubits=n)
    else:
        raise ValueError(mapping)
    q.compress(1e-12)
    return q


def pauli_stats(q: QubitOperator, n_qubits: int | None = None) -> dict:
    """Term count (excluding identity), exact Pauli 1-norm (excluding identity), weight histogram."""
    n_qubits = n_qubits or count_qubits(q)
    lam = 0.0
    n_terms = 0
    max_weight = 0
    weights = {}
    const = 0.0
    for term, c in q.terms.items():
        if term == ():
            const = float(np.real(c))
            continue
        n_terms += 1
        lam += abs(c)
        w = len(term)
        max_weight = max(max_weight, w)
        weights[w] = weights.get(w, 0) + 1
    return {"n_qubits": int(n_qubits), "n_pauli_terms": n_terms, "lambda_pauli": float(lam),
            "constant": const, "max_pauli_weight": max_weight,
            "weight_histogram": {str(k): v for k, v in sorted(weights.items())}}


def lambda_from_integrals(h1: np.ndarray, eri: np.ndarray) -> dict:
    """Jordan-Wigner Pauli 1-norm computed directly from spatial-orbital integrals (no Pauli expansion).

    Formula of Koridon et al., Phys. Rev. Research 3, 033127 (2021), eqs. (7)-(9), for
    H = sum h_pq a+a + 1/2 sum (pq|rs) a+_p a+_r a_s a_q with real orbitals:
        lambda_T = sum_pq | h_pq + sum_r [ (pq|rr) - 1/2 (pr|rq) ] |      ... effective one-body part
        lambda_V = 1/2 sum_{p>r, s>q} | (pq|rs) - (ps|rq) |  +  1/4 sum_pqrs | (pq|rs) |
    This is validated against the explicit JW expansion in tests (should agree to ~1e-10 for real
    orbitals); for complex (two-component) integrals it is only a proxy.
    """
    h1 = np.asarray(h1)
    g = np.asarray(eri)
    n = h1.shape[0]
    # effective one-body: h_pq + sum_r (pq|rr) - 1/2 sum_r (pr|rq)
    t_eff = h1 + np.einsum("pqrr->pq", g) - 0.5 * np.einsum("prrq->pq", g)
    lam_t = float(np.abs(t_eff).sum())
    lam_v2 = 0.25 * float(np.abs(g).sum())
    # 1/2 sum_{p>r, s>q} |g_pqrs - g_psrq|
    diff = g - g.transpose(0, 3, 2, 1)           # (pq|rs) - (ps|rq)
    pr_mask = np.tril(np.ones((n, n), bool), -1)  # p>r
    sq_mask = np.tril(np.ones((n, n), bool), -1)  # s>q
    # diff indices [p,q,r,s]; need p>r and s>q
    m = pr_mask[:, None, :, None] & sq_mask.T[None, :, None, :]  # sq_mask.T[q,s] true when s>q
    lam_v1 = 0.5 * float(np.abs(diff[m]).sum())
    return {"lambda_T": lam_t, "lambda_V": lam_v1 + lam_v2, "lambda_total": lam_t + lam_v1 + lam_v2,
            "n_spatial_orbitals": int(n), "n_qubits": int(2 * n),
            "n_nonzero_eri": int(np.count_nonzero(np.abs(g) > 1e-12)),
            "n_nonzero_h1": int(np.count_nonzero(np.abs(h1) > 1e-12))}


def pauli_term_count_upper_bound(n_spatial: int) -> int:
    """Generic upper bound on distinct JW Pauli terms for a real two-body Hamiltonian in N spatial
    orbitals (used only for extrapolation plots; actual counts are lower due to symmetry/sparsity)."""
    n = 2 * n_spatial
    return int(n * (n - 1) // 2 * 4 + n * n + n * (n - 1) * (n - 2) * (n - 3) // 24 * 8)


# ---------------------------------------------------------------- PennyLane interop and tapering
def to_pennylane(q: QubitOperator):
    import pennylane as qml
    return qml.from_openfermion(q)


def taper_hamiltonian(q: QubitOperator, n_electrons: int, n_qubits: int):
    """Z2-symmetry tapering with PennyLane. Returns (H_tapered, info)."""
    import pennylane as qml
    H = to_pennylane(q)
    gens = qml.qchem.symmetry_generators(H)
    paulixops = qml.qchem.paulix_ops(gens, n_qubits)
    sector = qml.qchem.optimal_sector(H, gens, n_electrons)
    Ht = qml.qchem.taper(H, gens, paulixops, sector)
    Ht = Ht.simplify()
    wires = sorted(set(int(w) for w in Ht.wires))
    coeffs, ops = Ht.terms()
    lam = float(sum(abs(c) for c, o in zip(coeffs, ops) if not _is_identity(o)))
    nterms = int(sum(1 for o in ops if not _is_identity(o)))
    info = {"n_generators": len(gens), "sector": [int(s) for s in sector],
            "n_qubits_tapered": len(wires), "n_pauli_terms_tapered": nterms, "lambda_pauli_tapered": lam}
    return Ht, info


def _is_identity(op) -> bool:
    import pennylane as qml
    if isinstance(op, qml.Identity):
        return True
    name = getattr(op, "name", "")
    if name == "Identity":
        return True
    try:
        pw = op.pauli_rep
        if pw is not None and len(pw) == 1:
            (word, _), = pw.items()
            return len(word) == 0
    except Exception:
        pass
    return False
