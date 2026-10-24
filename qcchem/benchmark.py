"""Classical references in qubit space: sparse exact diagonalisation (<= 16 qubits), sector-restricted."""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse.linalg as spla
from openfermion import (QubitOperator, count_qubits, get_sparse_operator,
                         jw_number_restrict_operator, jw_sz_restrict_operator)

MAX_ED_QUBITS = 16   # hard constraint: sparse eigsh only, never dense


def sparse_ground_state(q: QubitOperator, n_qubits: int, n_electrons: int | None = None,
                        sz: float | None = None, k: int = 1) -> dict:
    """Lowest k eigenvalues of the JW qubit Hamiltonian, optionally restricted to a particle-number
    (and Sz) sector.  Sz restriction assumes the interleaved alpha/beta ordering and even n_qubits."""
    if n_qubits > MAX_ED_QUBITS:
        raise ValueError(f"{n_qubits} qubits > MAX_ED_QUBITS={MAX_ED_QUBITS}")
    t0 = time.time()
    H = get_sparse_operator(q, n_qubits=n_qubits)
    sector = "full"
    if n_electrons is not None:
        if sz is not None:
            H = jw_sz_restrict_operator(H, sz, n_electrons=n_electrons, n_qubits=n_qubits)
            sector = f"N={n_electrons},Sz={sz}"
        else:
            H = jw_number_restrict_operator(H, n_electrons, n_qubits=n_qubits)
            sector = f"N={n_electrons}"
    dim = H.shape[0]
    if dim <= 64:
        vals = np.linalg.eigvalsh(H.toarray())[:k]
    else:
        vals = np.sort(spla.eigsh(H, k=min(k, dim - 2), which="SA", tol=1e-10)[0])
    return {"energies": [float(v) for v in vals], "e0": float(vals[0]), "sector": sector,
            "dim": int(dim), "t_seconds": round(time.time() - t0, 3)}
