"""Probe: free Ce3+ (4f1) and CeCl3 SCF + AVAS behaviour with scalar X2C and cc-pVDZ-DK."""
import time, sys
import numpy as np
from pyscf.mcscf import avas
from qcchem.integrals import build_mol, run_scf, orbital_composition

which = sys.argv[1] if len(sys.argv) > 1 else "ion"
if which == "ion":
    cases = [("Ce3+", "Ce 0 0 0", 3, 1), ("Pr3+", "Pr 0 0 0", 3, 2), ("Eu3+", "Eu 0 0 0", 3, 6), ("Gd3+", "Gd 0 0 0", 3, 7)]
    basis = "cc-pvdz-dk"
else:
    cases = [("CeCl3", "Ce 0 0 0; Cl 2.57 0 0; Cl -1.285 2.22568 0; Cl -1.285 -2.22568 0", 0, 1)]
    basis = "cc-pvdz-dk"
for name, atom, charge, spin in cases:
    t = time.time()
    mol = build_mol(atom, basis, charge=charge, spin=spin)
    mf = run_scf(mol, x2c=True)
    print(f"\n== {name} nao={mol.nao} E_ROHF(sfx2c)={mf.e_tot:.6f} conv={mf.converged} t={time.time()-t:.1f}s")
    el = atom.split()[0]
    for labels in ([f"{el} 4f"], [f"{el} 4f", f"{el} 5d"]) + ((["Ce 4f", "Cl 3p"], ["Ce 4f", "Ce 5d", "Cl 3p"]) if which != "ion" else ()):
        for thr in (0.2,):
            ncas, nelec, mo = avas.avas(mf, labels, threshold=thr, minao="ano", canonicalize=True)
            ncore = (mol.nelectron - nelec) // 2
            comp = orbital_composition(mol, mo, range(ncore, ncore + ncas), [f"{el} 4f", f"{el} 5d", f"{el} 6s", f"{el} 5p", "Cl 3p", "Cl 3s"])
            print(f"  AVAS {labels} thr={thr}: ncas={ncas} nelec={nelec} -> {2*ncas} qubits")
            for i, c in comp.items():
                top = sorted(c.items(), key=lambda kv: -kv[1])[:3]
                print(f"     mo{i}: " + ", ".join(f"{k}={v:.2f}" for k, v in top))
