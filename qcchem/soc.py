"""Two-component (spin-orbit) extension of an active-space Hamiltonian.

Model: scalar-relativistic (sfx2c1e) active space from `make_active_space` (real spatial orbitals,
spin-free two-electron integrals), plus the one-electron spin-orbit coupling taken from the full
two-component X2C one-electron Hamiltonian in the same orbital basis:

    h1_so[P,Q] = h1eff_scalar[p,q] delta_{s s'}  +  ( <p s| h_X2C(2c) |q s'>  -  <p| h_sfX2C |q> delta_{s s'} )

with spin-orbitals P=(p,s) in the interleaved ordering (2p = alpha, 2p+1 = beta).  The two-electron part
stays spin-free.  Limitations (stated in the report): one-electron SOC only (no two-electron / mean-field
screening, so 4f splittings are overestimated, typically by 20-40 %), no spin-dependent two-electron terms,
core orbitals frozen at the scalar level.

What this buys: the resulting qubit Hamiltonian has complex coefficients, no longer conserves Sz, and its
Pauli structure (term count, 1-norm, surviving Z2 symmetries) can be compared with the scalar one on an
equal footing (same orbitals, same integrals except the SOC block).  For 1- and 2-electron systems the
N-restricted sparse ED gives the spin-orbit multiplet levels (e.g. Ce3+ 2F5/2 / 2F7/2).
"""
from __future__ import annotations

import numpy as np
from openfermion import InteractionOperator, jordan_wigner
from openfermion.chem.molecular_data import spinorb_from_spatial
from pyscf import gto
from pyscf.x2c import sfx2c1e, x2c as x2c_mod

from .integrals import ActiveSpace


def x2c_oneelectron_in_spin_ao_basis(mol: gto.Mole) -> np.ndarray:
    """Full two-component X2C one-electron Hamiltonian transformed from the spinor basis to the
    (phi_mu alpha, phi_mu beta) real-spherical-AO x spin basis.  Returns a (2*nao, 2*nao) complex array
    in block order [alpha AO block ; beta AO block]."""
    xobj = x2c_mod.X2C(mol)
    h_sp = xobj.get_hcore(mol)                       # (nao_2c, nao_2c) complex, spinor basis
    ca, cb = mol.sph2spinor_coeff()                  # each (nao, nao_2c): spinor_i = sum_mu ca[mu,i] phi_mu a + cb[mu,i] phi_mu b
    U = np.vstack([ca, cb])                          # (2nao, nao_2c), columns = spinors in (AOa, AOb) basis
    # h in the (AOa, AOb) basis: h_sp = U^+ h_spin U  (U is square and invertible: same span)
    Uinv = np.linalg.inv(U)
    h_spin = Uinv.conj().T @ h_sp @ Uinv
    return h_spin


def scalar_x2c_hcore(mol: gto.Mole) -> np.ndarray:
    xobj = sfx2c1e.SpinFreeX2C(mol)
    return xobj.get_hcore(mol)


def soc_active_space(asp: ActiveSpace, mol: gto.Mole, mo_act: np.ndarray) -> dict:
    """Build the two-component active-space one-body matrix (interleaved spin-orbital order) from the
    scalar active space `asp` and its active orbital coefficients `mo_act` (nao, ncas)."""
    nao, ncas = mo_act.shape
    h_spin = x2c_oneelectron_in_spin_ao_basis(mol)
    h_sf = scalar_x2c_hcore(mol)
    C = np.zeros((2 * nao, 2 * ncas), dtype=complex)
    C[:nao, 0::2] = mo_act          # alpha spin-orbitals on even indices
    C[nao:, 1::2] = mo_act          # beta on odd
    h1_2c_mo = C.conj().T @ h_spin @ C                       # (2ncas, 2ncas) complex
    h_sf_mo = mo_act.T @ h_sf @ mo_act                       # (ncas, ncas) real
    h_sf_so = np.kron(h_sf_mo, np.eye(2))                    # interleaved
    delta_soc = h1_2c_mo - h_sf_so
    h1eff_so = np.kron(asp.h1, np.eye(2)).astype(complex) + delta_soc
    # hermiticity check and SOC magnitude diagnostics
    herm = float(np.max(np.abs(h1eff_so - h1eff_so.conj().T)))
    h1eff_so = 0.5 * (h1eff_so + h1eff_so.conj().T)
    diag_scalar_check = float(np.max(np.abs(0.5 * (h1_2c_mo[0::2, 0::2] + h1_2c_mo[1::2, 1::2]) - h_sf_mo)))
    return {"h1eff_so": h1eff_so, "delta_soc": delta_soc, "hermiticity_violation": herm,
            "soc_frobenius": float(np.linalg.norm(delta_soc)),
            "soc_max_abs": float(np.max(np.abs(delta_soc))),
            "spin_offdiag_max_abs": float(np.max(np.abs(delta_soc[0::2, 1::2]))),
            "scalar_part_vs_sfx2c_max_abs": diag_scalar_check}


def soc_interaction_operator(asp: ActiveSpace, h1eff_so: np.ndarray) -> InteractionOperator:
    two_body_phys = 0.5 * asp.eri.transpose(0, 2, 3, 1)
    _, h2s = spinorb_from_spatial(asp.h1, two_body_phys)     # spin-free two-body in interleaved order
    return InteractionOperator(asp.ecore, h1eff_so, h2s.astype(complex))


def soc_qubit_operator(asp: ActiveSpace, h1eff_so: np.ndarray):
    q = jordan_wigner(soc_interaction_operator(asp, h1eff_so))
    q.compress(1e-12)
    return q


def pauli_stats_complex(q) -> dict:
    lam = 0.0; n = 0; n_complex = 0; n_sz_breaking = 0
    for term, c in q.terms.items():
        if term == ():
            continue
        n += 1
        lam += abs(c)
        if abs(np.imag(c)) > 1e-10:
            n_complex += 1
        # a Pauli string changes Sz parity if it has an odd number of X/Y on alpha qubits (even index)
        # relative to beta qubits; simple diagnostic: count strings with X/Y on both an alpha and a beta qubit
        xy = [(w % 2) for (w, p) in term if p in ("X", "Y")]
        if xy and (sum(xy) % 2 == 1 or len(xy) % 2 == 1):
            n_sz_breaking += 1
    return {"n_pauli_terms": n, "lambda_pauli": float(lam), "n_terms_complex_coeff": n_complex,
            "n_terms_sz_nonconserving": n_sz_breaking}
