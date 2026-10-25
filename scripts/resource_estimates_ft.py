"""Fault-tolerant resource estimates (sparse / single-factorised / double-factorised qubitization) for an
active-space Hamiltonian stored in results/integrals/<name>_<tag>.npz, using the OpenFermion
`resource_estimates` module (Lee et al. 2021 / von Burg et al. 2021 / Berry et al. 2019 costings).

usage: scripts/qpy.sh scripts/resource_estimates_ft.py <name> [--tag single] [--casci-error] [--df-thresh ...]
Writes results/ft_<name>.json. Error metric per factorisation: ||ERI - ERI_rr|| and, with --casci-error,
the exact CASCI energy change when the factorised integrals replace the exact ones (pyscf FCI).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace  # noqa: E402
from qcchem.hamiltonian import lambda_from_integrals  # noqa: E402

from openfermion.resource_estimates.sparse.costing_sparse import cost_sparse  # noqa: E402
from openfermion.resource_estimates.sf.factorize_sf import factorize as sf_factorize  # noqa: E402
from openfermion.resource_estimates.sf.compute_cost_sf import compute_cost as sf_compute_cost  # noqa: E402
from openfermion.resource_estimates.df.factorize_df import factorize as df_factorize  # noqa: E402
from openfermion.resource_estimates.df.compute_cost_df import compute_cost as df_compute_cost  # noqa: E402

DE = 1e-3     # target phase-estimation error (Ha)
CHI = 10      # bits for coefficient representation
BETA = 20     # bits for rotations (DF)


def lambda_sf(h1, eri, LR):
    T = h1 - 0.5 * np.einsum("pqqs->ps", eri) + np.einsum("pqrr->pq", eri)
    return float(np.sum(np.abs(T)) + 0.25 * np.einsum("ijP,klP->", np.abs(LR), np.abs(LR)))


def lambda_df(h1, eri, LR):
    T = h1 - 0.5 * np.einsum("illj->ij", eri) + np.einsum("llij->ij", eri)
    lam_T = float(np.sum(np.abs(np.linalg.eigvalsh(T))))
    lam_F = sum(0.25 * np.sum(np.abs(np.linalg.eigvalsh(LR[:, :, v]))) ** 2 for v in range(LR.shape[2]))
    return lam_T + float(lam_F)


def casci_energy(h1, eri, ecore, ncas, nelec):
    from pyscf import fci
    cis = fci.direct_spin1.FCI()
    cis.conv_tol = 1e-12
    e, _ = cis.kernel(h1, eri, ncas, nelec, ecore=ecore)
    return float(e)


def unique_nonzero_integrals(h1, eri, thresh):
    """Number of symmetry-unique nonzero one- and two-body integrals above `thresh` (sparse-method 'd')."""
    n = h1.shape[0]
    d1 = int(np.sum(np.abs(np.triu(h1)) > thresh))
    cnt = 0
    for p in range(n):
        for q in range(p + 1):
            pq = p * (p + 1) // 2 + q
            for r in range(n):
                for s in range(r + 1):
                    rs = r * (r + 1) // 2 + s
                    if rs <= pq and abs(eri[p, q, r, s]) > thresh:
                        cnt += 1
    return d1 + cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--tag", default="single")
    ap.add_argument("--casci-error", action="store_true")
    ap.add_argument("--df-thresh", nargs="*", type=float, default=[1e-2, 1e-3, 1e-4, 1e-5])
    ap.add_argument("--sf-ranks", nargs="*", type=int, default=None)
    ap.add_argument("--sparse-thresh", nargs="*", type=float, default=[0.0, 1e-5, 1e-4, 1e-3])
    a = ap.parse_args()
    asp = ActiveSpace.load(str(ROOT / "results" / "integrals" / f"{a.name}_{a.tag}.npz"))
    h1, eri, ncas = asp.h1, asp.eri, asp.ncas
    n_so = 2 * ncas
    nelec = asp.nelecas
    out = {"name": a.name, "ncas": ncas, "n_spin_orbitals": n_so, "nelecas": list(nelec), "dE": DE, "chi": CHI, "beta": BETA}
    e_exact = casci_energy(h1, eri, asp.ecore, ncas, nelec) if a.casci_error else None
    out["e_casci_exact"] = e_exact

    # ---- sparse (Lee et al.): lambda = exact JW Pauli 1-norm (Koridon identity), d = unique nonzeros
    lam_sp = lambda_from_integrals(h1, eri)["lambda_total"]
    out["sparse"] = []
    for th in a.sparse_thresh:
        d = unique_nonzero_integrals(h1, eri, th)
        stps1 = cost_sparse(n_so, lam_sp, d, DE, CHI, stps=20000)[0]
        step, total, qubits = cost_sparse(n_so, lam_sp, d, DE, CHI, stps=stps1)
        rec = {"threshold": th, "d_unique_nonzero": d, "lambda": lam_sp, "toffoli_per_step": int(step),
               "toffoli_total": int(total), "logical_qubits": int(qubits)}
        if a.casci_error and th > 0:
            eri_t = np.where(np.abs(eri) > th, eri, 0.0)
            rec["casci_error_mha"] = (casci_energy(h1, eri_t, asp.ecore, ncas, nelec) - e_exact) * 1e3
        out["sparse"].append(rec)
        print(f"sparse thr={th:g}: d={d} lam={lam_sp:.2f} Toffoli={total:.3e} qubits={qubits}" + (f" dE_CASCI={rec.get('casci_error_mha', float('nan')):.3f} mHa" if a.casci_error else ""), flush=True)

    # ---- single factorisation (Berry et al.)
    ranks = a.sf_ranks or sorted({max(1, ncas), ncas * 2, ncas * 4, ncas * (ncas + 1) // 2})
    out["sf"] = []
    for rank in ranks:
        eri_rr, LR = sf_factorize(eri, rank)
        lam = lambda_sf(h1, eri_rr, LR)
        rec = {"rank": int(rank), "eri_error_fro": float(np.linalg.norm(eri_rr - eri)), "lambda": lam}
        try:
            stps1 = sf_compute_cost(n_so, lam, DE, L=rank, chi=CHI, stps=20000)[0]
            step, total, qubits = sf_compute_cost(n_so, lam, DE, L=rank, chi=CHI, stps=stps1)
            rec.update(toffoli_per_step=int(step), toffoli_total=int(total), logical_qubits=int(qubits))
        except BaseException as e:   # (SystemExit from the QROM helper) the costing's QROM layout needs L > M; small spaces/ranks are not costable
            rec["costing_error"] = str(e)[:120]; total = float("nan"); qubits = -1
        if a.casci_error:
            rec["casci_error_mha"] = (casci_energy(h1, eri_rr, asp.ecore, ncas, nelec) - e_exact) * 1e3
        out["sf"].append(rec)
        print(f"SF rank={rank}: |dERI|={rec['eri_error_fro']:.2e} lam={lam:.2f} Toffoli={float(total):.3e} qubits={qubits}" + (f" dE_CASCI={rec.get('casci_error_mha', float('nan')):.3f} mHa" if a.casci_error else ""), flush=True)

    # ---- double factorisation (von Burg et al.)
    out["df"] = []
    for th in a.df_thresh:
        eri_rr, LR, L, Lxi = df_factorize(eri, thresh=th)
        lam = lambda_df(h1, eri_rr, LR)
        rec = {"threshold": th, "rank_L": int(L), "n_eigenvectors_Lxi": int(Lxi),
               "eri_error_fro": float(np.linalg.norm(eri_rr - eri)), "lambda": lam}
        try:
            stps1 = df_compute_cost(n_so, lam, DE, L=L, Lxi=Lxi, chi=CHI, beta=BETA, stps=20000)[0]
            step, total, qubits = df_compute_cost(n_so, lam, DE, L=L, Lxi=Lxi, chi=CHI, beta=BETA, stps=stps1)
            rec.update(toffoli_per_step=int(step), toffoli_total=int(total), logical_qubits=int(qubits))
        except BaseException as e:
            rec["costing_error"] = str(e)[:120]; total = float("nan"); qubits = -1
        if a.casci_error:
            rec["casci_error_mha"] = (casci_energy(h1, eri_rr, asp.ecore, ncas, nelec) - e_exact) * 1e3
        out["df"].append(rec)
        print(f"DF thr={th:g}: L={L} Lxi={Lxi} |dERI|={rec['eri_error_fro']:.2e} lam={lam:.2f} Toffoli={float(total):.3e} qubits={qubits}" + (f" dE_CASCI={rec.get('casci_error_mha', float('nan')):.3f} mHa" if a.casci_error else ""), flush=True)

    out["note"] = ("OpenFermion resource_estimates costings: sparse (Lee et al. 2021), SF (Berry et al. 2019), "
                   "DF (von Burg et al. 2021); dE = 1 mHa, chi = 10, beta = 20; Toffoli counts and logical qubits "
                   "exclude surface-code overhead. Sparse lambda is the exact JW Pauli 1-norm.")
    json.dump(out, open(ROOT / "results" / f"ft_{a.name}.json", "w"), indent=1)
    print("wrote", ROOT / "results" / f"ft_{a.name}.json")


if __name__ == "__main__":
    main()
