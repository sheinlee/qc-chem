"""ADAPT-VQE vs UCCSD on stored active spaces.

usage: scripts/qpy.sh scripts/run_adapt.py <name> <tag> [<name> <tag> ...]   (e.g. n2_cas66 R=2.0)
Reads results/integrals/<name>_<tag>.npz, exact-diagonalises the (N,Sz) sector for the reference energy,
runs ADAPT until |E - E_exact| < 1e-6 Ha or the gradient falls below 1e-3, and records the UCCSD result for
the same point from results/<name>.json when available.  Writes results/adapt_<name>_<tag>.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from openfermion import get_sparse_operator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace  # noqa: E402
from qcchem.hamiltonian import to_qubit_operator  # noqa: E402
from qcchem.benchmark import sparse_ground_state  # noqa: E402
from qcchem.adapt import adapt_vqe, fermionic_gate_counts  # noqa: E402


import os
GRAD_TOL = float(os.environ.get("ADAPT_GRAD_TOL", "1e-3"))
MAX_OPS = int(os.environ.get("ADAPT_MAX_OPS", "60"))
SUFFIX = os.environ.get("ADAPT_SUFFIX", "")
POOL = os.environ.get("ADAPT_POOL", "fermionic")


def run(name, tag):
    asp = ActiveSpace.load(str(ROOT / "results" / "integrals" / f"{name}_{tag}.npz"))
    nq = asp.n_spin_orbitals
    q = to_qubit_operator(asp, "jw")
    H = get_sparse_operator(q, n_qubits=nq)
    sz = (asp.nelecas[0] - asp.nelecas[1]) / 2
    e_exact = sparse_ground_state(q, nq, n_electrons=asp.n_electrons, sz=sz)["e0"]
    print(f"[{name} {tag}] {nq} qubits, nelec {asp.nelecas}, E_exact {e_exact:.10f}", flush=True)
    res = adapt_vqe(H, asp.nelecas, nq, e_ref=e_exact, verbose=True, grad_tol=GRAD_TOL, max_ops=MAX_OPS, pool_type=POOL)
    res["grad_tol"] = GRAD_TOL; res["max_ops"] = MAX_OPS
    if POOL == "fermionic":
        res.update(fermionic_gate_counts(res["ops"], nq))
    elif POOL == "qubit":
        res["cnot_pauli"] = int(sum(2 * (len(lab[1]) - 1) for lab in res["ops"]))   # exp(i t P): 2(w-1) CNOTs
    res["error_vs_exact"] = res["energy"] - e_exact
    res["e_exact"] = e_exact
    res["e_casci"] = asp.e_casci
    # UCCSD comparison from the stored run (same geometry/tag)
    ucc = None
    try:
        d = json.load(open(ROOT / "results" / f"{name}.json", encoding="utf-8"))
        for p in d["points"]:
            if p.get("tag") == tag and "vqe_uccsd" in p:
                v = p["vqe_uccsd"]
                ucc = {"energy": v["energy"], "error_vs_casci": v["error_vs_casci"], "n_params": v["n_params"],
                       "n_singles": v["n_singles"], "n_doubles": v["n_doubles"], "cnot_fermionic": v.get("n_cnot"),
                       "depth_fermionic": v.get("circuit_depth"), "n_energy_evals": v.get("n_energy_evals"),
                       "cnot_givens_equivalent": 2 * v["n_singles"] + 14 * v["n_doubles"]}
    except Exception:
        pass
    res["uccsd"] = ucc
    out = {"name": name, "tag": tag, "n_qubits": nq, "nelecas": list(asp.nelecas), "adapt": res}
    json.dump(out, open(ROOT / "results" / f"adapt_{name}_{tag}{SUFFIX}.json", "w"), indent=1, default=str)
    print(f"  ADAPT: {res['n_ops']} ops ({res['n_singles']}S+{res['n_doubles']}D), err {res['error_vs_exact']:.2e}, "
          f"CNOT fermionic {res.get('cnot_fermionic')} / Givens {res.get('cnot_givens')} / Pauli {res.get('cnot_pauli')}, evals {res['n_energy_evals_total']}, "
          f"{res['t_seconds']} s, stop={res['stop_reason']}"
          + (f" | UCCSD: {ucc['n_params']} params, err {ucc['error_vs_casci']:.2e}, CNOT {ucc['cnot_fermionic']} / Givens-eq {ucc['cnot_givens_equivalent']}" if ucc else ""), flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    for i in range(0, len(args), 2):
        run(args[i], args[i + 1])
