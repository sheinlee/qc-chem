"""Noise resilience at fixed optimal angles (PLAN hard limits: depolarising with default.mixed <= 8 qubits,
shot noise <= 12 qubits; no re-optimisation under noise).

For each stored UCCSD result (results/<name>.json, point <tag>): rebuild the UCCSD circuit with the stored
theta*, then (a) depolarising channel of strength p after every gate of the decomposed circuit, energy vs p;
(b) finite-shot estimates of the energy on lightning.qubit, mean and standard deviation over repetitions.
usage: scripts/qpy.sh scripts/run_noise.py h2_sto3g R=0.74 lih_cas44 R=1.6 n2_cas66 R=1.2
Writes results/noise_<name>_<tag>.json.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from qcchem.integrals import ActiveSpace  # noqa: E402
from qcchem.hamiltonian import to_pennylane, to_qubit_operator  # noqa: E402
from qcchem.vqe import GATE_SET, uccsd_setup  # noqa: E402

P_GRID = [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
SHOTS = [1000, 10000, 100000]
REPS = 20


def run(name, tag):
    d = json.load(open(ROOT / "results" / f"{name}.json", encoding="utf-8"))
    pt = next(p for p in d["points"] if p["tag"] == tag)
    theta = pnp.array(pt["vqe_uccsd"]["theta"], requires_grad=False)
    asp = ActiveSpace.load(str(ROOT / "results" / "integrals" / f"{name}_{tag}.npz"))
    nq = asp.n_spin_orbitals
    H = to_pennylane(to_qubit_operator(asp, "jw"))
    s, dd, sw, dw, hf = uccsd_setup(asp.nelecas, nq)
    e_ref = pt["e_casci"]
    out = {"name": name, "tag": tag, "n_qubits": nq, "e_casci": e_ref, "e_uccsd_noiseless": pt["vqe_uccsd"]["energy"]}

    def ansatz(th):
        qml.UCCSD(th, wires=range(nq), s_wires=sw, d_wires=dw, init_state=hf)

    # (a) depolarising noise on the decomposed circuit (density matrix, only for small systems)
    if nq <= 8:
        dev0 = qml.device("default.qubit", wires=nq)

        @qml.qnode(dev0)
        def clean(th):
            ansatz(th); return qml.expval(H)

        tape = qml.workflow.construct_tape(clean)(theta)
        (tape_dec,), _ = qml.transforms.decompose(tape, gate_set=GATE_SET)
        ops = tape_dec.operations
        n_gates = len(ops)
        devm = qml.device("default.mixed", wires=nq)

        def noisy_energy(p):
            @qml.qnode(devm)
            def circ():
                for op in ops:
                    qml.apply(op)
                    if p > 0:
                        for w in op.wires:
                            qml.DepolarizingChannel(p, wires=w)
                return qml.expval(H)
            return float(circ())

        out["depolarizing"] = [{"p": p, "energy": noisy_energy(p), "error_vs_casci": noisy_energy(p) - e_ref} for p in P_GRID]
        out["n_gates_decomposed"] = n_gates
        print(f"[{name} {tag}] {nq} qubits, {n_gates} gates: " + ", ".join(f"p={r['p']:g}: {r['error_vs_casci']*1e3:+.2f} mHa" for r in out["depolarizing"]), flush=True)

    # (b) shot noise on the statevector simulator (no gate noise)
    shot_res = []
    for shots in SHOTS:
        dev = qml.device("lightning.qubit", wires=nq, shots=shots)

        @qml.qnode(dev)
        def sampled(th):
            ansatz(th); return qml.expval(H)

        vals = np.array([float(sampled(theta)) for _ in range(REPS)])
        shot_res.append({"shots": shots, "mean": float(vals.mean()), "std": float(vals.std(ddof=1)),
                         "mean_error_vs_casci": float(vals.mean() - e_ref), "reps": REPS})
        print(f"   shots={shots}: mean err {(vals.mean()-e_ref)*1e3:+.3f} mHa, std {vals.std(ddof=1)*1e3:.3f} mHa", flush=True)
    out["shot_noise"] = shot_res
    out["note"] = "fixed theta* from the noiseless optimisation; depolarising channel on every wire after every gate of the CNOT/rotation-decomposed circuit; shot noise = term-by-term sampling of H without grouping"
    json.dump(out, open(ROOT / "results" / f"noise_{name}_{tag}.json", "w"), indent=1)


if __name__ == "__main__":
    args = sys.argv[1:]
    for i in range(0, len(args), 2):
        t = time.time(); run(args[i], args[i + 1]); print(f"   ({time.time()-t:.0f} s)", flush=True)
