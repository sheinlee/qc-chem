"""Fermionic-pool ADAPT-VQE on an exact statevector (Grimsley et al., Nat. Commun. 10, 3007 (2019)).

Pool: anti-Hermitian generators A = T - T^+ of all spin-conserving single and double excitations from the
(n_alpha, n_beta) reference (the same operators that UCCSD uses, JW-mapped with OpenFermion), as sparse
matrices.  Ansatz |psi> = exp(theta_m A_m) ... exp(theta_1 A_1) |ref>.  Each iteration adds the pool operator
with the largest energy gradient |<psi|[H, A]|psi>| = 2 |Re <psi|H A|psi>| and re-optimises all angles (BFGS
with analytic gradients from a reverse sweep).  Exact for systems up to ~16 qubits; no sampling noise.

Circuit cost is reported in two gate sets: (i) PennyLane's fermionic (Jordan-Wigner string) excitation
gates, the same compilation as the UCCSD numbers in this repository; (ii) Givens-type particle-conserving
gates (SingleExcitation = 2 CNOT, DoubleExcitation = 14 CNOT) for a hardware-oriented lower count.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from openfermion import FermionOperator, get_sparse_operator, jordan_wigner
from scipy.optimize import minimize

from .vqe import uccsd_setup


def _ref_state(n_qubits: int, occ: np.ndarray) -> np.ndarray:
    idx = 0
    for q, o in enumerate(occ):           # OpenFermion/JW: qubit 0 is the most significant bit
        if o:
            idx |= 1 << (n_qubits - 1 - q)
    v = np.zeros(2 ** n_qubits, dtype=complex)
    v[idx] = 1.0
    return v


def build_pool(nelecas, n_qubits: int):
    singles, doubles, _, _, hf = uccsd_setup(nelecas, n_qubits)
    ops, labels = [], []
    for i, a in singles:
        T = FermionOperator(((a, 1), (i, 0)))
        ops.append(T - FermionOperator(((i, 1), (a, 0))))
        labels.append(("S", (i, a)))
    for i, j, a, b in doubles:
        T = FermionOperator(((b, 1), (a, 1), (j, 0), (i, 0)))
        ops.append(T - FermionOperator(((i, 1), (j, 1), (a, 0), (b, 0))))
        labels.append(("D", (i, j, a, b)))
    mats = [get_sparse_operator(jordan_wigner(op), n_qubits=n_qubits).tocsc() for op in ops]
    return mats, labels, hf


def adapt_vqe(H_sparse, nelecas, n_qubits: int, *, e_ref: float | None = None, max_ops: int = 60,
              grad_tol: float = 1e-3, energy_tol: float | None = 1e-6, bfgs_gtol: float = 1e-7,
              verbose: bool = False) -> dict:
    t0 = time.time()
    H = H_sparse.tocsc()
    pool, labels, hf = build_pool(nelecas, n_qubits)
    psi0 = _ref_state(n_qubits, hf)
    chosen: list[int] = []
    theta = np.zeros(0)
    history = []

    def apply_ansatz(th, return_states=False):
        states = [psi0]
        psi = psi0
        for k, idx in enumerate(chosen):
            psi = spla.expm_multiply(th[k] * pool[idx], psi)
            if return_states:
                states.append(psi)
        return (psi, states) if return_states else psi

    def energy_and_grad(th):
        psi, states = apply_ansatz(th, return_states=True)
        hpsi = H @ psi
        e = float(np.real(np.vdot(psi, hpsi)))
        # reverse sweep: lam_k = U_{k+1}^+ ... U_n^+ H psi ; dE/dtheta_k = 2 Re <lam_k| A_k |psi_k>
        g = np.zeros(len(th))
        lam = hpsi
        for k in range(len(chosen) - 1, -1, -1):
            g[k] = 2.0 * float(np.real(np.vdot(lam, pool[chosen[k]] @ states[k + 1])))
            lam = spla.expm_multiply(-th[k] * pool[chosen[k]], lam)
        return e, g

    e_cur = float(np.real(np.vdot(psi0, H @ psi0)))
    e_hf = e_cur
    n_energy_evals = 0
    for it in range(max_ops):
        psi = apply_ansatz(theta)
        hpsi = H @ psi
        grads = np.array([2.0 * np.real(np.vdot(hpsi, A @ psi)) for A in pool])
        k = int(np.argmax(np.abs(grads)))
        gmax = float(abs(grads[k]))
        if gmax < grad_tol:
            history.append({"iter": it, "stop": "gradient", "max_grad": gmax})
            break
        chosen.append(k)
        theta = np.append(theta, 0.0)
        res = minimize(lambda x: energy_and_grad(x), theta, jac=True, method="BFGS",
                       options={"gtol": bfgs_gtol, "maxiter": 500})
        n_energy_evals += res.nfev
        theta = res.x
        e_cur = float(res.fun)
        entry = {"iter": it + 1, "op": labels[k], "max_grad": gmax, "energy": e_cur, "n_params": len(theta),
                 "bfgs_nit": int(res.nit)}
        if e_ref is not None:
            entry["error"] = e_cur - e_ref
        history.append(entry)
        if verbose:
            print(f"  ADAPT it {it+1}: op {labels[k]} |g|={gmax:.2e} E={e_cur:.10f}"
                  + (f" err={e_cur - e_ref:.2e}" if e_ref is not None else ""), flush=True)
        if e_ref is not None and energy_tol is not None and abs(e_cur - e_ref) < energy_tol:
            history.append({"iter": it + 1, "stop": "energy", "error": e_cur - e_ref})
            break
    n_s = sum(1 for i in chosen if labels[i][0] == "S")
    n_d = len(chosen) - n_s
    return {"energy": e_cur, "e_hf": e_hf, "n_ops": len(chosen), "n_singles": n_s, "n_doubles": n_d,
            "ops": [labels[i] for i in chosen], "theta": [float(x) for x in theta], "history": history,
            "pool_size": len(pool), "n_energy_evals_total": int(n_energy_evals),
            "cnot_givens": 2 * n_s + 14 * n_d, "t_seconds": round(time.time() - t0, 1),
            "stop_reason": history[-1].get("stop", "max_ops") if history else "none"}


def fermionic_gate_counts(ops, n_qubits: int) -> dict:
    """CNOT count and depth of the chosen operators compiled with PennyLane's fermionic excitation gates
    (same gate set as the UCCSD numbers in this repository)."""
    import pennylane as qml
    from .vqe import GATE_SET

    def circuit(params):
        for p, (kind, idx) in zip(params, ops):
            if kind == "S":
                i, a = idx
                qml.FermionicSingleExcitation(p, wires=list(range(i, a + 1)))
            else:
                i, j, a, b = idx
                qml.FermionicDoubleExcitation(p, wires1=list(range(i, j + 1)), wires2=list(range(a, b + 1)))
        return qml.state()

    dev = qml.device("default.qubit", wires=n_qubits)
    qnode = qml.QNode(circuit, dev)
    tape = qml.workflow.construct_tape(qnode)(np.full(len(ops), 0.1))
    (tape_dec,), _ = qml.transforms.decompose(tape, gate_set=GATE_SET)
    res = tape_dec.specs["resources"]
    return {"cnot_fermionic": int(res.gate_types.get("CNOT", 0)), "depth_fermionic": int(res.depth),
            "n_gates_fermionic": int(res.num_gates)}
