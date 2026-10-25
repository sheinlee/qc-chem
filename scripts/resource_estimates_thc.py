"""THC (tensor hypercontraction) resource estimates without pybtas: direct L-BFGS least-squares fit of the
ERI tensor with OpenFermion's jax objective (`lbfgsb_opt_thc`, random start, best of several seeds), lambda
per Lee et al. 2021 (eqs. 13, 18-20) and the THC costing of `openfermion.resource_estimates.thc`.

usage: scripts/qpy.sh scripts/resource_estimates_thc.py <name> [--tag single] [--nthc 68 102 136] [--seeds 2]
       [--casci-error] [--maxiter 20000]
Writes results/thc_<name>.json.  Note: the fit quality of a random-start L-BFGS THC is worse than the CP3 +
BFGS route of the original paper; the numbers here are upper bounds on what THC can do.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace  # noqa: E402

from openfermion.resource_estimates.thc.utils.thc_factorization import lbfgsb_opt_thc  # noqa: E402
from openfermion.resource_estimates.thc.compute_cost_thc import compute_cost as thc_compute_cost  # noqa: E402

DE, CHI, BETA = 1e-3, 10, 20


def thc_reconstruct(etaPp, MPQ):
    return np.einsum("Pp,Pr,Qq,Qs,PQ->prqs", etaPp, etaPp, etaPp, etaPp, MPQ, optimize=True)


def lambda_thc(h1, eri, etaPp, MPQ):
    SPQ = etaPp @ etaPp.T
    cP = np.diag(np.diag(SPQ))
    lam_z = 0.5 * np.sum(np.abs(cP @ MPQ @ cP))
    T = h1 - 0.5 * np.einsum("illj->ij", eri) + np.einsum("llij->ij", eri)
    lam_T = float(np.sum(np.abs(np.linalg.eigvalsh(T))))
    return float(lam_z + lam_T), float(lam_z), lam_T


def casci_energy(h1, eri, ecore, ncas, nelec):
    from pyscf import fci
    cis = fci.direct_spin1.FCI(); cis.conv_tol = 1e-12
    return float(cis.kernel(h1, eri, ncas, nelec, ecore=ecore)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name"); ap.add_argument("--tag", default="single")
    ap.add_argument("--nthc", nargs="*", type=int, default=None)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--maxiter", type=int, default=20000)
    ap.add_argument("--casci-error", action="store_true")
    a = ap.parse_args()
    asp = ActiveSpace.load(str(ROOT / "results" / "integrals" / f"{a.name}_{a.tag}.npz"))
    h1, eri, n = asp.h1, asp.eri, asp.ncas
    n_so = 2 * n
    nthc_list = a.nthc or [4 * n, 6 * n, 8 * n]
    e_exact = casci_energy(h1, eri, asp.ecore, n, asp.nelecas) if a.casci_error else None
    out = {"name": a.name, "ncas": n, "n_spin_orbitals": n_so, "dE": DE, "chi": CHI, "beta": BETA, "e_casci_exact": e_exact, "thc": []}
    scratch = Path(os.environ.get("TMP", "/tmp"))
    for nthc in nthc_list:
        best = None
        t0 = time.time()
        for seed in range(a.seeds):
            chk = str(scratch / f"thc_{a.name}_{nthc}_{seed}.h5")
            x = lbfgsb_opt_thc(eri, nthc, chkfile_name=chk, random_seed=seed, maxiter=a.maxiter)
            etaPp = x[: n * nthc].reshape(nthc, n); MPQ = x[n * nthc:].reshape(nthc, nthc)
            MPQ = 0.5 * (MPQ + MPQ.T)
            eri_thc = thc_reconstruct(etaPp, MPQ)
            err = float(np.linalg.norm(eri_thc - eri))
            if best is None or err < best[0]:
                best = (err, etaPp, MPQ, eri_thc, seed)
            try:
                os.remove(chk)
            except OSError:
                pass
        err, etaPp, MPQ, eri_thc, seed = best
        lam, lam_z, lam_T = lambda_thc(h1, eri, etaPp, MPQ)
        rec = {"nthc": nthc, "eri_error_fro": err, "lambda": lam, "lambda_z": lam_z, "lambda_T": lam_T,
               "best_seed": seed, "t_seconds": round(time.time() - t0, 1)}
        try:
            stps1 = thc_compute_cost(n_so, lam, DE, chi=CHI, beta=BETA, M=nthc, stps=20000)[0]
            step, total, qubits = thc_compute_cost(n_so, lam, DE, chi=CHI, beta=BETA, M=nthc, stps=stps1)
            rec.update(toffoli_per_step=int(step), toffoli_total=int(total), logical_qubits=int(qubits))
        except BaseException as e:
            rec["costing_error"] = str(e)[:120]
        if a.casci_error:
            rec["casci_error_mha"] = (casci_energy(h1, eri_thc, asp.ecore, n, asp.nelecas) - e_exact) * 1e3
        out["thc"].append(rec)
        print(f"THC nthc={nthc}: |dERI|={err:.2e} lam={lam:.2f} (z {lam_z:.2f}, T {lam_T:.2f}) Toffoli={rec.get('toffoli_total', float('nan')):.3e} qubits={rec.get('logical_qubits')}"
              + (f" dE_CASCI={rec.get('casci_error_mha', float('nan')):.3f} mHa" if a.casci_error else "") + f" [{rec['t_seconds']} s]", flush=True)
    out["note"] = "random-start L-BFGS THC fit (no CP3 initialisation); lambda per Lee et al. 2021; costs exclude error correction"
    json.dump(out, open(ROOT / "results" / f"thc_{a.name}.json", "w"), indent=1)
    print("wrote", ROOT / "results" / f"thc_{a.name}.json")


if __name__ == "__main__":
    main()
