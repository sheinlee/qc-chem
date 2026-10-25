"""lambda-minimising orbital rotations within a fixed active subspace (Koridon et al., PRR 3, 033127, 2021).

Minimise the JW Pauli 1-norm lambda(U) = lambda_from_integrals(U^T h1 U, (U x U x U x U) eri) over real
orthogonal rotations U = exp(K), K antisymmetric, of the active orbitals only. The active subspace, the CI
space and therefore the CASCI energy are invariant; this is checked by exact CASCI (pyscf FCI) before and
after. The non-smooth |.| is handled by a smoothed sqrt(x^2 + eps^2) objective with decreasing eps, L-BFGS,
and a final evaluation of the exact lambda.

usage: scripts/qpy.sh scripts/lambda_orbital_opt.py <name> [--tag single] [--restarts 3]
Writes results/lamopt_<name>.json (and the rotated integrals to results/integrals/<name>_<tag>_lamopt.npz).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace  # noqa: E402
from qcchem.hamiltonian import lambda_from_integrals, pauli_stats, to_qubit_operator  # noqa: E402


def rotate(h1, eri, U):
    h1r = U.T @ h1 @ U
    erir = np.einsum("pqrs,pi,qj,rk,sl->ijkl", eri, U, U, U, U, optimize=True)
    return h1r, erir


def smooth_lambda(h1, g, eps):
    """Smoothed version of lambda_from_integrals (same decomposition, |x| -> sqrt(x^2+eps^2))."""
    n = h1.shape[0]
    t_eff = h1 + np.einsum("pqrr->pq", g) - 0.5 * np.einsum("prrq->pq", g)
    sa = lambda x: np.sqrt(x * x + eps * eps)
    lam_t = sa(t_eff).sum()
    lam_v2 = 0.25 * sa(g).sum()
    diff = g - g.transpose(0, 3, 2, 1)
    pr = np.tril(np.ones((n, n), bool), -1)
    m = pr[:, None, :, None] & pr.T[None, :, None, :]
    lam_v1 = 0.5 * sa(diff[m]).sum()
    return lam_t + lam_v1 + lam_v2


def casci_energy(h1, eri, ecore, ncas, nelec):
    from pyscf import fci
    cis = fci.direct_spin1.FCI(); cis.conv_tol = 1e-12
    return float(cis.kernel(h1, eri, ncas, nelec, ecore=ecore)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name"); ap.add_argument("--tag", default="single"); ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--maxiter", type=int, default=300)
    ap.add_argument("--no-casci", action="store_true", help="skip the exact CASCI invariance check (CI space too large)")
    a = ap.parse_args()
    asp = ActiveSpace.load(str(ROOT / "results" / "integrals" / f"{a.name}_{a.tag}.npz"))
    h1, eri, n = asp.h1, asp.eri, asp.ncas
    iu = np.triu_indices(n, 1)
    nparam = len(iu[0])

    def unpack(x):
        K = np.zeros((n, n)); K[iu] = x; K = K - K.T
        return expm(K)

    def obj(x, eps):
        U = unpack(x); h1r, g = rotate(h1, eri, U)
        return smooth_lambda(h1r, g, eps)

    lam0 = lambda_from_integrals(h1, eri)["lambda_total"]
    e0 = None if a.no_casci else casci_energy(h1, eri, asp.ecore, n, asp.nelecas)
    rng = np.random.default_rng(0)
    best = (lam0, np.zeros(nparam))
    t0 = time.time()
    history = []
    for r in range(a.restarts):
        x = np.zeros(nparam) if r == 0 else rng.normal(0, 0.3, nparam)
        for eps in (1e-2, 1e-3, 1e-4):
            res = minimize(obj, x, args=(eps,), method="L-BFGS-B", options={"maxiter": a.maxiter})
            x = res.x
        U = unpack(x); h1r, g = rotate(h1, eri, U)
        lam = lambda_from_integrals(h1r, g)["lambda_total"]
        history.append({"restart": r, "lambda": lam, "nit": int(res.nit)})
        print(f"restart {r}: lambda {lam0:.4f} -> {lam:.4f}  ({lam/lam0:.3f})", flush=True)
        if lam < best[0]:
            best = (lam, x.copy())
    U = unpack(best[1]); h1r, g = rotate(h1, eri, U)
    e1 = None if a.no_casci else casci_energy(h1r, g, asp.ecore, n, asp.nelecas)
    out = {"name": a.name, "ncas": n, "lambda_initial": lam0, "lambda_optimised": best[0], "ratio": best[0] / lam0,
           "casci_energy_initial": e0, "casci_energy_rotated": e1,
           "casci_invariance_ha": None if a.no_casci else e1 - e0,
           "casci_check": "skipped (CI space too large; invariance holds by unitarity)" if a.no_casci else "exact CASCI before/after",
           "restarts": history, "t_seconds": round(time.time() - t0, 1)}
    if 2 * n <= 24:
        from openfermion import QubitOperator
        asp_r = ActiveSpace(name=a.name + "_lamopt", h1=h1r, eri=g, ecore=asp.ecore, ncas=n, nelecas=asp.nelecas,
                            e_hf=asp.e_hf, e_casci=asp.e_casci, mo_act=None if asp.mo_act is None else asp.mo_act @ U)
        st0 = pauli_stats(to_qubit_operator(asp, "jw"), 2 * n)
        st1 = pauli_stats(to_qubit_operator(asp_r, "jw"), 2 * n)
        out["pauli_terms_initial"] = st0["n_pauli_terms"]; out["pauli_terms_optimised"] = st1["n_pauli_terms"]
        out["lambda_jw_check"] = [st0["lambda_pauli"], st1["lambda_pauli"]]
        asp_r.save(str(ROOT / "results" / "integrals" / f"{a.name}_{a.tag}_lamopt.npz"))
    json.dump(out, open(ROOT / "results" / f"lamopt_{a.name}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "restarts"}, indent=1))


if __name__ == "__main__":
    main()
