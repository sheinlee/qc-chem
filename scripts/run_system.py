"""Config-driven driver: YAML -> active space -> qubit Hamiltonian -> statistics / ED / VQE -> JSON.

usage: scripts/qpy.sh scripts/run_system.py configs/n2_cas66.yaml [--only-stats] [--points 3]
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from qcchem.integrals import build_mol, make_active_space, run_scf  # noqa: E402
from qcchem.hamiltonian import (MAX_EXPLICIT_QUBITS, lambda_from_integrals, pauli_stats,  # noqa: E402
                                taper_hamiltonian, to_pennylane, to_qubit_operator)
from qcchem.benchmark import MAX_ED_QUBITS, sparse_ground_state  # noqa: E402
from qcchem.resources import qpe_proxy  # noqa: E402
from qcchem.vqe import MAX_VQE_QUBITS, run_uccsd_vqe  # noqa: E402


def versions() -> dict:
    import openfermion, pennylane, pyscf, scipy
    return {"pyscf": pyscf.__version__, "openfermion": openfermion.__version__,
            "pennylane": pennylane.__version__, "scipy": scipy.__version__, "numpy": np.__version__,
            "python": platform.python_version(), "machine": platform.node(), "os": platform.platform(),
            "git_sha": _git_sha()}


def _git_sha() -> str:
    import subprocess
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "unknown"


def run_point(cfg: dict, tag: str, atom: str, methods: dict, outdir: Path) -> dict:
    t0 = time.time()
    mol = build_mol(atom, cfg["basis"], charge=cfg.get("charge", 0), spin=cfg.get("spin", 0))
    cas = cfg.get("cas", {})
    asp = make_active_space(
        f"{cfg['name']}_{tag}", mol, ncas=cas.get("ncas"), nelecas=cas.get("nelecas"),
        x2c=cfg.get("x2c", False), avas_labels=cas.get("avas_labels"),
        avas_threshold=cas.get("avas_threshold", 0.2), avas_minao=cas.get("avas_minao", "minao"),
        avas_openshell_option=cas.get("avas_openshell_option", 2), relax_reference=cas.get("relax_reference"),
        casscf=cas.get("casscf", False), solve_ci=cas.get("solve_ci", True),
        ci_max_memory_mb=cas.get("ci_max_memory_mb", 16000), frac_occ=cfg.get("frac_occ", False),
        state_average=cas.get("state_average", 0),
        drop_active=cas.get("drop_active"), casscf_max_cycle=cas.get("casscf_max_cycle", 50),
        casscf_conv_tol=cas.get("casscf_conv_tol", 1e-7),
        chkfile=(str(ROOT / "results" / "chk" / f"{cfg['name']}.chk") if cas.get("chkfile") else None),
        restart_from_chk=cas.get("restart_from_chk", False), freeze_core=cas.get("freeze_core", False),
        fix_spin=cas.get("fix_spin"), composition_labels=cfg.get("composition_labels"))
    nq = asp.n_spin_orbitals
    rec = {"tag": tag, "atom": atom, "n_active_electrons": asp.n_electrons, "nelecas": list(asp.nelecas),
           "n_active_orbitals": asp.ncas, "n_spin_orbitals": nq, "n_qubits": nq,
           "e_hf": asp.e_hf, "e_casci": asp.e_casci, "e_casscf": asp.e_casscf, "meta": asp.meta}
    (outdir / "chk").mkdir(parents=True, exist_ok=True)
    npz = outdir / "integrals" / f"{cfg['name']}_{tag}.npz"
    npz.parent.mkdir(parents=True, exist_ok=True)
    asp.save(str(npz))
    rec["integrals_file"] = str(npz.relative_to(ROOT))

    # exact lambda from integrals: always available (no Pauli expansion)
    rec["lambda_from_integrals"] = lambda_from_integrals(asp.h1, asp.eri)
    rec["lambda_from_integrals"]["kind"] = ("exact JW Pauli 1-norm (Koridon identity, real orbitals)"
                                            if np.isrealobj(asp.h1) else "proxy (complex integrals)")

    if nq <= MAX_EXPLICIT_QUBITS and methods.get("jw", True):
        qjw = to_qubit_operator(asp, "jw")
        rec["jw"] = pauli_stats(qjw, nq)
        rec["jw"]["lambda_identity_check"] = rec["jw"]["lambda_pauli"] - rec["lambda_from_integrals"]["lambda_total"]
        rec["qpe_proxy_jw"] = qpe_proxy(nq, rec["jw"]["n_pauli_terms"], rec["jw"]["lambda_pauli"])
        if methods.get("bk", False):
            rec["bk"] = pauli_stats(to_qubit_operator(asp, "bk"), nq)
        if methods.get("taper", True):
            try:
                _, tinfo = taper_hamiltonian(qjw, asp.n_electrons, nq)
                rec["taper"] = tinfo
                rec["qpe_proxy_tapered"] = qpe_proxy(tinfo["n_qubits_tapered"], tinfo["n_pauli_terms_tapered"],
                                                     tinfo["lambda_pauli_tapered"])
            except Exception as e:  # tapering is diagnostic; never abort the point
                rec["taper"] = {"error": repr(e)}
        if methods.get("ed", True) and nq <= MAX_ED_QUBITS:
            sz = (asp.nelecas[0] - asp.nelecas[1]) / 2
            ed = sparse_ground_state(qjw, nq, n_electrons=asp.n_electrons, sz=sz, k=methods.get("ed_k", 1))
            ed["error_vs_casci"] = ed["e0"] - asp.e_casci
            rec["ed"] = ed
        vqe_limit = 14 if methods.get("vqe_allow_14", False) else MAX_VQE_QUBITS
        if methods.get("vqe", True) and nq <= vqe_limit:
            Hpl = to_pennylane(qjw)
            v = run_uccsd_vqe(Hpl, asp.nelecas, nq, allow_14=methods.get("vqe_allow_14", False),
                              gtol=methods.get("vqe_gtol", 1e-6), maxiter=methods.get("vqe_maxiter", 500))
            v["error_vs_casci"] = v["energy"] - asp.e_casci
            # restart strategy (benchmark luxury: the reference energy is known): if the first run is more
            # than 1e-5 Ha above CASCI, retry from small random angles and keep the best. Count is recorded.
            rng = np.random.default_rng(0)
            restarts = []
            for _ in range(int(methods.get("vqe_restarts", 2))):
                if abs(v["error_vs_casci"]) <= 1e-5:
                    break
                th0 = rng.normal(0.0, 0.2, size=v["n_params"])
                v2 = run_uccsd_vqe(Hpl, asp.nelecas, nq, allow_14=methods.get("vqe_allow_14", False),
                                   gtol=methods.get("vqe_gtol", 1e-6), maxiter=methods.get("vqe_maxiter", 500),
                                   theta0=th0, with_resources=False)
                restarts.append(v2["energy"])
                if v2["energy"] < v["energy"]:
                    for k in ("energy", "theta", "n_iter", "n_energy_evals", "n_gradient_evals", "final_grad_maxabs",
                              "success", "message", "polished_lbfgsb"):
                        v[k] = v2[k]
                    v["t_seconds"] = round(v["t_seconds"] + v2["t_seconds"], 2)
                    v["error_vs_casci"] = v["energy"] - asp.e_casci
            v["restart_energies"] = restarts
            v["n_restarts_used"] = len(restarts)
            v["chemical_accuracy_reached"] = abs(v["error_vs_casci"]) < 1.6e-3
            v.pop("trace", None)
            rec["vqe_uccsd"] = v
    else:
        rec["qpe_proxy_from_integrals"] = qpe_proxy(
            nq, rec["lambda_from_integrals"]["n_nonzero_eri"] + rec["lambda_from_integrals"]["n_nonzero_h1"],
            rec["lambda_from_integrals"]["lambda_total"])
        rec["qpe_proxy_from_integrals"]["note"] += ("; term count proxy = nonzero elements of the full (pq|rs) "
                                                     "tensor + nonzero h_pq (no Pauli expansion; order of magnitude)")
    rec["t_seconds"] = round(time.time() - t0, 2)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--only-stats", action="store_true", help="skip ED and VQE")
    ap.add_argument("--points", type=int, default=None, help="use only the first N scan points")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config, encoding="utf-8"))
    methods = dict(cfg.get("methods", {}))
    if args.only_stats:
        methods.update(ed=False, vqe=False)
    outdir = ROOT / "results"
    outdir.mkdir(exist_ok=True)
    out = Path(args.out) if args.out else outdir / f"{cfg['name']}.json"

    scan = cfg.get("scan")
    points = []
    if scan:
        vals = scan["values"][: args.points] if args.points else scan["values"]
        for v in vals:
            atom = cfg["atom"].replace("{" + scan["param"] + "}", str(v))
            if cfg.get("geometry") == "symmetric_stretch":   # water-like: R = O-H, angle 104.5 deg
                half = np.deg2rad(104.5 / 2)
                atom = atom.replace("{HY}", f"{v * np.sin(half):.6f}").replace("{HZ}", f"{v * np.cos(half):.6f}")
            points.append((f"{scan['param']}={v}", atom))
    else:
        points.append(("single", cfg["atom"]))

    results = {"name": cfg["name"], "config": cfg, "versions": versions(), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
               "points": []}
    for tag, atom in points:
        print(f"[{cfg['name']}] {tag} ...", flush=True)
        try:
            rec = run_point(cfg, tag, atom, methods, outdir)
        except Exception as e:
            import traceback
            rec = {"tag": tag, "atom": atom, "error": repr(e), "traceback": traceback.format_exc()}
            print(f"  FAILED: {e!r}", flush=True)
        results["points"].append(rec)
        summary = {k: rec.get(k) for k in ("n_qubits", "e_hf", "e_casci")}
        if "jw" in rec:
            summary["jw_terms"] = rec["jw"]["n_pauli_terms"]; summary["lambda"] = round(rec["jw"]["lambda_pauli"], 6)
        if "ed" in rec:
            summary["ed_err"] = f"{rec['ed']['error_vs_casci']:.1e}"
        if "vqe_uccsd" in rec:
            summary["vqe_err"] = f"{rec['vqe_uccsd']['error_vs_casci']:.2e}"; summary["vqe_t"] = rec["vqe_uccsd"]["t_seconds"]
        print("  ", summary, flush=True)
        # write after every point so partial results survive interruption
        json.dump(results, open(out, "w", encoding="utf-8"), indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    results["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    json.dump(results, open(out, "w", encoding="utf-8"), indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    print("wrote", out)


if __name__ == "__main__":
    main()
