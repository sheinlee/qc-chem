"""Two-component (one-electron X2C spin-orbit) extension of a scalar active space.

usage: scripts/qpy.sh scripts/run_soc.py configs/ce3_4f_sa.yaml
Reads results/integrals/<name>_single.npz (needs mo_act), rebuilds the molecule, adds the SOC block,
compares scalar vs two-component qubit Hamiltonians, and (<= 16 qubits) does N-restricted sparse ED.
Writes results/soc_<name>.json.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla
import yaml
from openfermion import get_sparse_operator, jw_number_restrict_operator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace, build_mol  # noqa: E402
from qcchem.hamiltonian import pauli_stats, to_qubit_operator  # noqa: E402
from qcchem.soc import pauli_stats_complex, soc_active_space, soc_qubit_operator  # noqa: E402

HA_TO_CM = 219474.63


def ed_levels(q, nq, nelec, k):
    H = get_sparse_operator(q, n_qubits=nq)
    H = jw_number_restrict_operator(H, nelec, n_qubits=nq)
    dim = H.shape[0]
    if dim <= 400:
        vals = np.linalg.eigvalsh(H.toarray())[:k]
    else:
        vals = np.sort(spla.eigsh(H, k=min(k, dim - 2), which="SA", tol=1e-10)[0])
    return [float(v) for v in vals], int(dim)


def main(cfg_path: str):
    cfg = yaml.safe_load(open(cfg_path, encoding="utf-8"))
    name = cfg["name"]
    npz = ROOT / "results" / "integrals" / f"{name}_single.npz"
    asp = ActiveSpace.load(str(npz))
    assert asp.mo_act is not None, "npz has no mo_act; rerun run_system.py for this config"
    mol = build_mol(cfg["atom"], cfg["basis"], charge=cfg.get("charge", 0), spin=cfg.get("spin", 0))
    t0 = time.time()
    soc = soc_active_space(asp, mol, asp.mo_act)
    nq = asp.n_spin_orbitals
    nelec = asp.n_electrons
    out = {"name": name, "n_qubits": nq, "n_electrons": nelec, "ncas": asp.ncas,
           "soc_diagnostics": {k: v for k, v in soc.items() if not isinstance(v, np.ndarray)},
           "model": "scalar sfx2c1e active space + one-electron X2C spin-orbit block (no 2e SOC screening)"}
    # SOC one-body block statistics
    d = soc["delta_soc"]
    out["soc_block"] = {"n_nonzero": int(np.count_nonzero(np.abs(d) > 1e-8)),
                        "max_abs_ha": float(np.max(np.abs(d))), "max_abs_cm": float(np.max(np.abs(d)) * HA_TO_CM)}
    if nq <= 24:
        q_sc = to_qubit_operator(asp, "jw")
        q_so = soc_qubit_operator(asp, soc["h1eff_so"])
        st_sc = pauli_stats(q_sc, nq)
        st_so = pauli_stats_complex(q_so)
        out["scalar"] = st_sc
        out["two_component"] = st_so
        out["ratio_terms"] = st_so["n_pauli_terms"] / st_sc["n_pauli_terms"]
        out["ratio_lambda"] = st_so["lambda_pauli"] / st_sc["lambda_pauli"]
        # Z2 symmetries surviving (PennyLane): scalar vs two-component
        try:
            import pennylane as qml
            from qcchem.hamiltonian import to_pennylane
            g_sc = qml.qchem.symmetry_generators(to_pennylane(q_sc))
            g_so = qml.qchem.symmetry_generators(to_pennylane(q_so))
            out["z2_generators_scalar"] = len(g_sc)
            out["z2_generators_two_component"] = len(g_so)
        except Exception as e:
            out["z2_error"] = repr(e)
        if nq <= 16:
            k = int(cfg.get("methods", {}).get("ed_k", 7))
            e_sc, dim = ed_levels(q_sc, nq, nelec, k)
            e_so, _ = ed_levels(q_so, nq, nelec, k)
            out["ed_scalar_levels_mha"] = [round((e - e_sc[0]) * 1e3, 4) for e in e_sc]
            out["ed_two_component_levels_mha"] = [round((e - e_so[0]) * 1e3, 4) for e in e_so]
            out["ed_two_component_levels_cm"] = [round((e - e_so[0]) * HA_TO_CM, 1) for e in e_so]
            out["ed_scalar_e0"] = e_sc[0]; out["ed_two_component_e0"] = e_so[0]
            out["soc_ground_state_shift_mha"] = (e_so[0] - e_sc[0]) * 1e3
            out["ed_dim_N_sector"] = dim
    out["t_seconds"] = round(time.time() - t0, 1)
    outp = ROOT / "results" / f"soc_{name}.json"
    json.dump(out, open(outp, "w", encoding="utf-8"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("scalar", "two_component")}, indent=1, default=str))
    print("wrote", outp)


if __name__ == "__main__":
    main(sys.argv[1])
