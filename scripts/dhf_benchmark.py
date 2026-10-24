"""Two-electron spin-orbit screening benchmark for the Ce 4f shell.

Compare, for the closed-shell Ce4+ (4f0) and Ce3+ ions, the 4f5/2 - 4f7/2 spinor orbital-energy splitting from
  (a) two-component X2C-HF (one-electron X2C spin-orbit only, Coulomb mean field without 2e-SOC),
  (b) four-component Dirac-Coulomb HF (2e spin-same-orbit screening included implicitly),
  (c) Dirac-Coulomb-Gaunt HF (adds spin-other-orbit) when affordable.
The ratio (b)/(a) is the mean-field screening factor of the 4f one-electron SOC. Closed-shell Ce4+ gives
a clean Koopmans comparison (both 4f spinors virtual); Ce3+ 4f1 is also run for reference.
usage: scripts/qpy.sh scripts/dhf_benchmark.py [basis]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import gto, scf
from pyscf.x2c import x2c

HA_TO_CM = 219474.63
ROOT = Path(__file__).resolve().parents[1]


def f_spinor_mask(mol):
    """Boolean mask over spinor AO labels selecting 4f (first f shell) functions."""
    return np.array([("4f" in l) for l in mol.spinor_labels()])


def analyse_4f(mol, mo_energy, mo_coeff, n_positive_offset=0, tag=""):
    """Find the 14 spinor MOs with the largest 4f character and group their energies into j=5/2 / 7/2 sets."""
    s = mol.intor("int1e_ovlp_spinor")
    mask = f_spinor_mask(mol)
    C = mo_coeff[:, n_positive_offset:]
    E = mo_energy[n_positive_offset:]
    # 4f weight of each spinor MO: C_f^+ S_ff C_f (block-restricted Mulliken-type projection)
    pop = np.real(np.einsum("pi,pq,qi->i", C[mask].conj(), s[mask][:, mask], C[mask]))
    idx = np.argsort(-pop)[:14]
    e = np.sort(E[idx])
    w = np.sort(pop[idx])
    # 4f5/2 (6 lowest) vs 4f7/2 (8 highest)
    e52, e72 = e[:6], e[6:]
    out = {"tag": tag, "e_4f_sorted": [float(x) for x in e], "min_4f_weight": float(w[0]),
           "e52_mean": float(e52.mean()), "e72_mean": float(e72.mean()),
           "spread52": float(e52.max() - e52.min()), "spread72": float(e72.max() - e72.min()),
           "split_ha": float(e72.mean() - e52.mean())}
    out["split_cm"] = out["split_ha"] * HA_TO_CM
    return out


def run(charge, basis, do_gaunt):
    mol = gto.M(atom="Ce 0 0 0", basis=basis, charge=charge, spin=0 if charge == 4 else 1, verbose=0)
    res = {"charge": charge, "basis": basis, "nao": mol.nao, "nao_2c": mol.nao_2c()}
    # (a) 2c X2C-HF: one-electron SOC + nonrelativistic Coulomb/exchange
    t = time.time()
    mf2c = x2c.UHF(mol)
    mf2c.conv_tol = 1e-9
    mf2c.kernel()
    res["x2c_hf"] = analyse_4f(mol, mf2c.mo_energy, mf2c.mo_coeff, 0, "X2C-HF (1e SOC)")
    res["x2c_hf"].update(e_tot=float(mf2c.e_tot), converged=bool(mf2c.converged), t_s=round(time.time() - t, 1))
    # (b) 4c Dirac-Coulomb HF
    t = time.time()
    mf4c = scf.DHF(mol)
    mf4c.conv_tol = 1e-9
    mf4c.kernel()
    n2c = mol.nao_2c()
    res["dc_hf"] = analyse_4f(mol, mf4c.mo_energy, mf4c.mo_coeff[:n2c], n2c, "Dirac-Coulomb HF (4c)")
    res["dc_hf"].update(e_tot=float(mf4c.e_tot), converged=bool(mf4c.converged), t_s=round(time.time() - t, 1))
    if do_gaunt:
        t = time.time()
        mfg = scf.DHF(mol)
        mfg.with_gaunt = True
        mfg.conv_tol = 1e-9
        mfg.kernel()
        res["dcg_hf"] = analyse_4f(mol, mfg.mo_energy, mfg.mo_coeff[:n2c], n2c, "Dirac-Coulomb-Gaunt HF (4c)")
        res["dcg_hf"].update(e_tot=float(mfg.e_tot), converged=bool(mfg.converged), t_s=round(time.time() - t, 1))
    res["screening_factor_dc_over_x2c"] = res["dc_hf"]["split_ha"] / res["x2c_hf"]["split_ha"]
    if "dcg_hf" in res:
        res["screening_factor_dcg_over_x2c"] = res["dcg_hf"]["split_ha"] / res["x2c_hf"]["split_ha"]
    return res


if __name__ == "__main__":
    basis = sys.argv[1] if len(sys.argv) > 1 else "cc-pvdz-dk"
    out = {"experiment_ce3_2F_split_cm": 2253.0, "results": []}
    for charge, gaunt in [(4, True), (3, False)]:
        r = run(charge, basis, gaunt)
        out["results"].append(r)
        for k in ("x2c_hf", "dc_hf", "dcg_hf"):
            if k in r:
                v = r[k]
                print(f"Ce{charge}+ {v['tag']:28s} split = {v['split_ha']*1e3:8.3f} mHa = {v['split_cm']:8.1f} cm-1  "
                      f"(min 4f weight {v['min_4f_weight']:.2f}, spreads {v['spread52']*1e3:.2f}/{v['spread72']*1e3:.2f} mHa, "
                      f"{v['t_s']} s, conv={v['converged']})", flush=True)
        print(f"   screening factor DC/X2C = {r['screening_factor_dc_over_x2c']:.3f}",
              f"DCG/X2C = {r.get('screening_factor_dcg_over_x2c', float('nan')):.3f}", flush=True)
    json.dump(out, open(ROOT / "results" / "soc_dhf_benchmark.json", "w"), indent=1)
    print("wrote results/soc_dhf_benchmark.json")
