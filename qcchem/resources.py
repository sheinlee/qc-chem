"""Proxy fault-tolerant resource estimates from (n_qubits, n_terms, lambda).

These are *proxies* based on the Pauli 1-norm, not optimised resource counts (no double
factorisation / THC, no surface-code assumptions).  Reported quantities:
  - queries_qpe = ceil(pi * lambda / (2 * eps)): number of walk-operator (qubitised) queries for
    phase estimation to precision eps (Heisenberg-limited, standard qubitization scaling).
  - toffoli_proxy_sparse = queries * ceil(log2(n_terms)) * n_qubits  (rough per-query cost of a
    sparse LCU SELECT/PREPARE over n_terms Pauli strings on n_qubits).
  - logical_qubits_proxy = n_qubits + ceil(log2(n_terms)) + ceil(log2(queries)) (system + LCU
    register + phase register).
"""
from __future__ import annotations

import math

CHEMICAL_ACCURACY_HA = 1.6e-3


def qpe_proxy(n_qubits: int, n_terms: int, lam: float, eps: float = 1e-3) -> dict:
    queries = math.ceil(math.pi * lam / (2 * eps))
    log_terms = max(1, math.ceil(math.log2(max(n_terms, 2))))
    return {
        "eps_ha": eps,
        "lambda": lam,
        "queries_qpe": queries,
        "toffoli_proxy_sparse": queries * log_terms * n_qubits,
        "logical_qubits_proxy": n_qubits + log_terms + math.ceil(math.log2(max(queries, 2))),
        "note": "Pauli-1-norm proxy; not an optimised fault-tolerant resource count",
    }
