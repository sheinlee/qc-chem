"""Nestedness and baseline checks for the CeCl3 ladder (run after the ladder JSON/npz files exist).

For each pair (A ⊂? B) report how much of A's active subspace lies inside B's active subspace, separately for
the occupied-like (ligand) and the metal orbitals, using the S-metric projector  P_B = C_B C_B^T S.
Also report E_CASCI - E_ROHF per level with one common ROHF baseline and the CI convergence settings.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace, build_mol  # noqa: E402

GEOM = "Ce 0 0 0; Cl 2.57 0 0; Cl -1.285 2.22568 0; Cl -1.285 -2.22568 0"
mol = build_mol(GEOM, "cc-pvdz-dk", charge=0, spin=1)
S = mol.intor_symmetric("int1e_ovlp")
levels = ["cecl3_L1_4f", "cecl3_L2_4f5d", "cecl3_C3_4f_sigma", "cecl3_L3_4f5d_sigma", "cecl3_L4_4f_cl3p_ci", "cecl3_L5_full"]
asp = {}
rec = {}
for n in levels:
    f = ROOT / "results" / "integrals" / f"{n}_single.npz"
    if f.exists():
        asp[n] = ActiveSpace.load(str(f))
        rec[n] = json.load(open(ROOT / "results" / f"{n}.json"))["points"][0]


def subspace_fraction(CA, CB):
    """Mean over columns of A of the squared norm of their projection onto span(B) (S-metric)."""
    M = CA.T @ S @ CB                     # overlaps between orthonormal sets
    return float(np.mean(np.sum(M * M, axis=1))), [round(float(x), 4) for x in np.sum(M * M, axis=1)]


e_rohf = None
print("level                  CAS        E_CASCI - E_ROHF / mHa   CI conv_tol   natocc(min occ, max virt)")
for n in levels:
    if n not in rec:
        print(f"{n:22s} missing"); continue
    p = rec[n]
    if e_rohf is None:
        e_rohf = p["e_hf"]
    de = (p["e_casci"] - e_rohf) * 1e3 if p["e_casci"] == p["e_casci"] else float("nan")
    no = p["meta"].get("ci_natural_occupations")
    occ = [x for x in no if x > 1.5] if no else []
    vir = [x for x in no if x < 0.5] if no else []
    print(f"{n:22s} ({p['n_active_electrons']},{p['n_active_orbitals']})   {de:9.4f}   {p['meta'].get('ci_conv_tol')}   "
          f"{min(occ) if occ else '-'}, {max(vir) if vir else '-'}   dE_HF vs common baseline {(p['e_hf']-e_rohf)*1e9:.1f} nHa")

print("\nsubspace containment (fraction of A inside B):")
pairs = [("cecl3_C3_4f_sigma", "cecl3_L3_4f5d_sigma"), ("cecl3_C3_4f_sigma", "cecl3_L4_4f_cl3p_ci"),
         ("cecl3_L1_4f", "cecl3_C3_4f_sigma"), ("cecl3_L1_4f", "cecl3_L4_4f_cl3p_ci"), ("cecl3_L4_4f_cl3p_ci", "cecl3_L5_full"),
         ("cecl3_L3_4f5d_sigma", "cecl3_L5_full"), ("cecl3_L2_4f5d", "cecl3_L3_4f5d_sigma")]
for a, b in pairs:
    if a in asp and b in asp:
        fr, per = subspace_fraction(asp[a].mo_act, asp[b].mo_act)
        print(f"  {a:22s} in {b:22s}: mean {fr:.4f}   per-orbital {per}")
