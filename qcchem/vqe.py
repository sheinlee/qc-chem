"""UCCSD-VQE on PennyLane lightning.qubit with adjoint gradients and scipy BFGS."""
from __future__ import annotations

import time

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
from scipy.optimize import minimize

MAX_VQE_QUBITS = 12          # hard constraint; 14 only with allow_14=True (tiny excitation spaces)

# Hardware-style gate set used for depth / two-qubit counts (CNOT + single-qubit rotations).
GATE_SET = {"CNOT", "RX", "RY", "RZ", "Hadamard", "PauliX", "PauliY", "PauliZ", "S", "T",
            "PhaseShift", "Rot", "SX", "Adjoint(S)"}


def _nelec_tuple(n_electrons) -> tuple[int, int]:
    if isinstance(n_electrons, (tuple, list)):
        return int(n_electrons[0]), int(n_electrons[1])
    n = int(n_electrons)
    return (n + 1) // 2, n // 2


def uccsd_setup(n_electrons, n_qubits: int, delta_sz: int = 0):
    """Spin-conserving singles/doubles from an (n_alpha, n_beta) reference in interleaved ordering.

    n_electrons may be an int (closed shell / standard PennyLane convention) or an (na, nb) tuple for
    open-shell high-spin references (e.g. Ce3+ 4f1 = (1, 0), Pr3+ 4f2 = (2, 0)).
    """
    na, nb = _nelec_tuple(n_electrons)
    occ = sorted([2 * p for p in range(na)] + [2 * p + 1 for p in range(nb)])
    virt = [q for q in range(n_qubits) if q not in occ]
    spin = lambda q: q % 2  # noqa: E731  (0 = alpha, 1 = beta)
    singles = [[i, a] for i in occ for a in virt if spin(i) == spin(a)]
    doubles = []
    for ii, i in enumerate(occ):
        for j in occ[ii + 1:]:
            for aa, a in enumerate(virt):
                for b in virt[aa + 1:]:
                    if spin(i) + spin(j) == spin(a) + spin(b):
                        doubles.append([i, j, a, b])
    s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)
    hf = np.zeros(n_qubits, dtype=int)
    hf[occ] = 1
    return singles, doubles, s_wires, d_wires, hf


def circuit_resources(qnode, *args) -> dict:
    """Gate counts and depth after decomposing to GATE_SET (CNOT + 1-qubit gates)."""
    tape = qml.workflow.construct_tape(qnode)(*args)
    (tape_dec,), _ = qml.transforms.decompose(tape, gate_set=GATE_SET)
    specs = tape_dec.specs
    res = specs["resources"]
    two_q = int(sum(v for k, v in res.gate_sizes.items() if k >= 2))
    return {"circuit_depth": int(res.depth), "n_gates": int(res.num_gates), "n_two_qubit_gates": two_q,
            "n_cnot": int(res.gate_types.get("CNOT", 0)),
            "gate_types": {str(k): int(v) for k, v in res.gate_types.items()},
            "n_template_ops": int(len(tape.operations))}


def run_uccsd_vqe(H, n_electrons, n_qubits: int, *, allow_14: bool = False,
                  method: str = "BFGS", gtol: float = 1e-6, maxiter: int = 500,
                  theta0: np.ndarray | None = None, verbose: bool = False,
                  with_resources: bool = True) -> dict:
    """Returns dict with energy, parameters, evaluations, circuit specs and timings.

    H: PennyLane operator (e.g. from qml.from_openfermion). Interleaved spin-orbital ordering.
    """
    limit = 14 if allow_14 else MAX_VQE_QUBITS
    if n_qubits > limit:
        raise ValueError(f"{n_qubits} qubits > VQE limit {limit}")
    singles, doubles, s_wires, d_wires, hf = uccsd_setup(n_electrons, n_qubits)
    n_params = len(singles) + len(doubles)
    dev = qml.device("lightning.qubit", wires=n_qubits)

    @qml.qnode(dev, diff_method="adjoint")
    def energy(theta):
        qml.UCCSD(theta, wires=range(n_qubits), s_wires=s_wires, d_wires=d_wires, init_state=hf)
        return qml.expval(H)

    grad = qml.grad(energy)
    n_eval = {"f": 0, "g": 0}
    trace = []

    def f(x):
        n_eval["f"] += 1
        e = float(energy(pnp.array(x, requires_grad=False)))
        trace.append(e)
        return e

    def g(x):
        n_eval["g"] += 1
        return np.asarray(grad(pnp.array(x, requires_grad=True)), dtype=float)

    x0 = np.zeros(n_params) if theta0 is None else np.asarray(theta0, float)
    t0 = time.time()
    e_hf = f(x0)
    res = minimize(f, x0, jac=g, method=method, options={"gtol": gtol, "maxiter": maxiter})
    # BFGS sometimes stops on "precision loss" before gtol; polish once with L-BFGS-B if needed
    gnorm = float(np.max(np.abs(g(res.x))))
    polished = False
    if gnorm > gtol and res.nit < maxiter:
        res2 = minimize(f, res.x, jac=g, method="L-BFGS-B",
                        options={"gtol": gtol, "ftol": 1e-14, "maxiter": maxiter})
        if res2.fun <= res.fun:
            res, polished = res2, True
        gnorm = float(np.max(np.abs(g(res.x))))
    t_opt = time.time() - t0

    out = {
        "ansatz": "UCCSD", "n_qubits": n_qubits, "n_electrons": list(_nelec_tuple(n_electrons)),
        "n_params": int(n_params), "n_singles": len(singles), "n_doubles": len(doubles),
        "e_hf_circuit": e_hf, "energy": float(res.fun), "success": bool(res.success),
        "message": str(res.message), "n_iter": int(res.nit), "final_grad_maxabs": gnorm,
        "polished_lbfgsb": polished,
        "n_energy_evals": n_eval["f"], "n_gradient_evals": n_eval["g"],
        "t_seconds": round(t_opt, 2), "theta": [float(x) for x in res.x],
        "optimizer": method, "gtol": gtol, "trace": [float(e) for e in trace],
    }
    if with_resources:
        out.update(circuit_resources(energy, pnp.array(res.x, requires_grad=False)))
    if verbose:
        print(f"UCCSD-VQE: E={res.fun:.10f} iters={res.nit} nfev={n_eval['f']} |g|={gnorm:.1e} t={t_opt:.1f}s")
    return out
