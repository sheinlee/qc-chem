"""Check N2 CAS(6,6) UCCSD convergence and decomposed circuit resources."""
from qcchem.integrals import build_mol, make_active_space
from qcchem.hamiltonian import to_qubit_operator, to_pennylane
from qcchem.vqe import run_uccsd_vqe

mol = build_mol("N 0 0 0; N 0 0 1.1", "sto-3g")
asp = make_active_space("N2", mol, ncas=6, nelecas=6)
q = to_qubit_operator(asp)
v = run_uccsd_vqe(to_pennylane(q), 6, 12, verbose=True)
for k in ["energy", "success", "message", "n_iter", "final_grad_maxabs", "polished_lbfgsb", "n_energy_evals",
          "circuit_depth", "n_gates", "n_two_qubit_gates", "n_cnot", "n_template_ops", "t_seconds"]:
    print(f"  {k}: {v[k]}")
print("  err vs CASCI:", v["energy"] - asp.e_casci)
