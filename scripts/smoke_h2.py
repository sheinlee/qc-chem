"""Smoke test of the whole pipeline on H2 and LiH (run: scripts/qpy.sh scripts/smoke_h2.py)."""
import time
from qcchem.integrals import build_mol, make_active_space
from qcchem.hamiltonian import to_qubit_operator, pauli_stats, lambda_from_integrals, taper_hamiltonian, to_pennylane
from qcchem.benchmark import sparse_ground_state
from qcchem.vqe import run_uccsd_vqe

cases = [
    ("H2", "H 0 0 0; H 0 0 0.74", "sto-3g", 2, 2, 0),
    ("LiH_cas22", "Li 0 0 0; H 0 0 1.6", "sto-3g", 2, 2, 0),
    ("LiH_cas44", "Li 0 0 0; H 0 0 1.6", "sto-3g", 4, 4, 0),
    ("N2_cas66", "N 0 0 0; N 0 0 1.1", "sto-3g", 6, 6, 0),
]
for name, atom, basis, ncas, nelec, spin in cases:
    mol = build_mol(atom, basis, spin=spin)
    asp = make_active_space(name, mol, ncas=ncas, nelecas=nelec)
    nq = asp.n_spin_orbitals
    qjw = to_qubit_operator(asp, "jw"); qbk = to_qubit_operator(asp, "bk")
    st = pauli_stats(qjw, nq); stbk = pauli_stats(qbk, nq)
    lam = lambda_from_integrals(asp.h1, asp.eri)
    ed = sparse_ground_state(qjw, nq, n_electrons=asp.n_electrons, sz=0)
    edbk = sparse_ground_state(qbk, nq)   # BK: no simple number restriction; full spectrum min
    Ht, tinfo = taper_hamiltonian(qjw, asp.n_electrons, nq)
    import numpy as np, scipy.sparse.linalg as spla
    Hm = Ht.sparse_matrix()
    et = float(np.min(np.linalg.eigvalsh(Hm.toarray()))) if Hm.shape[0] <= 4096 else float(spla.eigsh(Hm, k=1, which="SA")[0][0])
    print(f"\n== {name}: {nq} qubits, nelec={asp.nelecas}, HF={asp.e_hf:.8f} CASCI={asp.e_casci:.10f}")
    print(f"  JW terms={st['n_pauli_terms']} lambda={st['lambda_pauli']:.8f} | BK terms={stbk['n_pauli_terms']} lambda={stbk['lambda_pauli']:.8f}")
    print(f"  lambda_from_integrals={lam['lambda_total']:.8f}  diff vs JW={lam['lambda_total']-st['lambda_pauli']:.2e}")
    print(f"  ED(N,Sz) e0={ed['e0']:.10f} diff={ed['e0']-asp.e_casci:.2e} dim={ed['dim']} | BK full-min={edbk['e0']:.10f} diff={edbk['e0']-asp.e_casci:.2e}")
    print(f"  taper: {tinfo}  E_taper={et:.10f} diff={et-asp.e_casci:.2e}")
    if nq <= 12:
        t = time.time()
        v = run_uccsd_vqe(to_pennylane(qjw), asp.n_electrons, nq)
        print(f"  VQE-UCCSD E={v['energy']:.10f} err={v['energy']-asp.e_casci:.2e} params={v['n_params']} iters={v['n_iter']} nfev={v['n_energy_evals']} depth={v['circuit_depth']} 2q={v['n_two_qubit_gates']} t={v['t_seconds']}s")
