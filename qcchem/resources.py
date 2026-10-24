"""Proxy fault-tolerant resource estimates from (n_qubits, n_terms, lambda).

These are *proxies* based on the Pauli 1-norm, not optimised resource counts (no double
factorisation / THC, no surface-code assumptions).  Reported quantities:
  - queries_qpe = ceil(pi * lambda / (2 * eps)): number of walk-operator (qubitised) queries for
    phase estimation to precision eps (Heisenberg-limited, standard qubitization scaling).
  - toffoli_proxy_sparse = queries * n_terms: a sparse-LCU walk operator costs O(L) Toffolis per query
    (SELECT by unary iteration over L = n_terms Pauli strings; PREPARE via QROM is of the same order).
    Double factorisation / THC reduce both L and lambda; they are not evaluated here.
  - lambda is the Fock-space 1-norm of the full qubit Hamiltonian. For spaces with very few electrons
    (e.g. CAS(1,7)) most two-body terms vanish in the physical N-sector, so lambda and the query count
    overstate the difficulty of that particular problem; they characterise the Hamiltonian, not the sector.
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
        "toffoli_proxy_sparse": queries * n_terms,
        "logical_qubits_proxy": n_qubits + log_terms + math.ceil(math.log2(max(queries, 2))),
        "note": "Fock-space Pauli-1-norm proxy (sparse LCU, SELECT ~ L Toffolis/query); not an optimised "
                "fault-tolerant resource count (no DF/THC, no error correction)",
    }
