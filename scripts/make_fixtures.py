"""Generate tests/fixtures/h2_uccsd_reference.json (reference angles + energy for the regression test)."""
import json
from pathlib import Path
from qcchem.integrals import build_mol, make_active_space
from qcchem.hamiltonian import to_qubit_operator, to_pennylane
from qcchem.vqe import run_uccsd_vqe

asp = make_active_space("H2", build_mol("H 0 0 0; H 0 0 0.74", "sto-3g"), ncas=2, nelecas=2)
v = run_uccsd_vqe(to_pennylane(to_qubit_operator(asp)), 2, 4, gtol=1e-10, with_resources=False)
out = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "h2_uccsd_reference.json"
out.parent.mkdir(exist_ok=True)
json.dump({"theta": v["theta"], "energy": v["energy"], "e_casci": asp.e_casci, "geometry": "H 0 0 0; H 0 0 0.74",
           "basis": "sto-3g"}, open(out, "w"), indent=1)
print("wrote", out, "err", v["energy"] - asp.e_casci)
