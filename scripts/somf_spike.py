"""Feasibility spike for SOMF screening.
1. Convention check: Breit-Pauli one-electron 2p splitting of hydrogenic Ne9+/Ar17+ vs Dirac.
2. Ce3+ (state-averaged 4f orbitals from results/integrals/ce3_4f_sa_single.npz): 2F7/2-2F5/2 from
   X2C-1e (current model), BP-1e, BP-1e + SOMF-2e, BP-1e + AMFI-2e, and X2C-1e + SOMF-2e.
usage: scripts/qpy.sh scripts/somf_spike.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sl
import yaml
from pyscf import gto, scf
from pyscf.data import nist

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace, build_mol, run_scf  # noqa: E402
from qcchem.soc import scalar_x2c_hcore, soc_active_space  # noqa: E402
from qcchem.somf import bp_soc_1e_ao, bp_soc_2e_ao, soc_block_spinorb  # noqa: E402

HA_TO_CM = 219474.63
c = nist.LIGHT_SPEED


def dirac(Z, n, j):
    za = Z / c; k = j + 0.5; d = k - np.sqrt(k * k - za * za)
    return c * c / np.sqrt(1.0 + (za / (n - d)) ** 2) - c * c


out = {}
# ---------------------------------------------------------------- 1. hydrogenic convention check
for Z, sym in [(10, "Ne"), (18, "Ar")]:
    basis = {sym: gto.etbs([(0, 26, 0.02, 2.4), (1, 26, 0.02, 2.4)])}
    mol = gto.M(atom=f"{sym} 0 0 0", basis=basis, charge=Z - 1, spin=1, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    # nonrelativistic 2p orbitals from the NR one-electron Hamiltonian (BP is a perturbation on NR)
    h_nr = mol.intor_symmetric("int1e_kin") + mol.intor_symmetric("int1e_nuc")
    es, Cs = sl.eigh(h_nr, s)
    p_mask = np.array([("p" in l.split()[2]) for l in mol.ao_labels()])
    pop_p = np.array([(Cs[:, i] * (s @ Cs[:, i]))[p_mask].sum() for i in range(mol.nao)])
    idx = [i for i in range(mol.nao) if pop_p[i] > 0.99][:3]
    C = Cs[:, idx]
    hso = soc_block_spinorb(bp_soc_1e_ao(mol), C)
    e = np.linalg.eigvalsh(hso)                      # SOC-only eigenvalues: 2p1/2 (x2) and 2p3/2 (x4)
    split_bp = e[2:].mean() - e[:2].mean()
    exact = dirac(Z, 2, 1.5) - dirac(Z, 2, 0.5)
    out[f"hydrogenic_{sym}"] = {"bp_1e_split_ha": float(split_bp), "dirac_split_ha": float(exact), "ratio": float(split_bp / exact)}
    print(f"hydrogenic {sym}{Z-1}+: BP-1e 2p split {split_bp*1e3:.4f} mHa vs Dirac {exact*1e3:.4f} mHa  ratio {split_bp/exact:.4f}")

# ---------------------------------------------------------------- 2. Ce3+ 4f
cfg = yaml.safe_load(open(ROOT / "configs/ce3_4f_sa.yaml", encoding="utf-8"))
asp = ActiveSpace.load(str(ROOT / "results/integrals/ce3_4f_sa_single.npz"))
mol = build_mol(cfg["atom"], cfg["basis"], charge=cfg["charge"], spin=cfg["spin"])
mf = run_scf(mol, x2c=True)                              # scalar reference for the SOMF density
dm = mf.make_rdm1(); dm = dm[0] + dm[1] if np.ndim(dm) == 3 else dm
mo = asp.mo_act
x2c = soc_active_space(asp, mol, mo)["delta_soc"]        # X2C one-electron SOC block (current model)
bp1 = soc_block_spinorb(bp_soc_1e_ao(mol), mo)
somf2 = soc_block_spinorb(bp_soc_2e_ao(mol, dm, amfi=False), mo)
amfi2 = soc_block_spinorb(bp_soc_2e_ao(mol, dm, amfi=True), mo)


def split(h):
    e = np.linalg.eigvalsh(h)
    return float(e[6:].mean() - e[:6].mean()), float(e[6:].max() - e[6:].min()), float(e[:6].max() - e[:6].min())


res = {}
for name, h in [("X2C-1e", x2c), ("BP-1e", bp1), ("BP-1e+SOMF-2e", bp1 + somf2), ("BP-1e+AMFI-2e", bp1 + amfi2),
                ("X2C-1e+SOMF-2e", x2c + somf2), ("X2C-1e+AMFI-2e", x2c + amfi2)]:
    sp, w72, w52 = split(h)
    res[name] = {"split_ha": sp, "split_cm": sp * HA_TO_CM, "spread72_mha": w72 * 1e3, "spread52_mha": w52 * 1e3}
    print(f"Ce3+ 4f  {name:16s} 2F7/2-2F5/2 = {sp*1e3:8.3f} mHa = {sp*HA_TO_CM:8.1f} cm-1   (spreads {w52*1e3:.3f}/{w72*1e3:.3f} mHa)")
res["screening_factor_somf"] = res["BP-1e+SOMF-2e"]["split_ha"] / res["BP-1e"]["split_ha"]
res["screening_factor_amfi"] = res["BP-1e+AMFI-2e"]["split_ha"] / res["BP-1e"]["split_ha"]
res["experiment_cm"] = 2253.0
print(f"screening factor SOMF {res['screening_factor_somf']:.3f}  AMFI {res['screening_factor_amfi']:.3f}; experiment 2253 cm-1")
out["ce3_4f"] = res
json.dump(out, open(ROOT / "results/soc_somf_spike.json", "w"), indent=1)
print("wrote results/soc_somf_spike.json")
