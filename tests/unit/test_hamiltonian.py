"""Deterministic checks: integral folding, mappings, 1-norm identity, tapering, fixed-parameter VQE energy."""
import json
from pathlib import Path

import numpy as np
import pytest
from pennylane import numpy as pnp
import pennylane as qml

from qcchem.benchmark import sparse_ground_state
from qcchem.hamiltonian import (lambda_from_integrals, pauli_stats, taper_hamiltonian, to_pennylane,
                                to_qubit_operator)
from qcchem.vqe import uccsd_setup

FIX = Path(__file__).resolve().parents[1] / "fixtures"


@pytest.mark.parametrize("fx", ["h2", "lih44", "n2_stretched"])
def test_folding_ed_equals_casci(fx, request):
    """OpenFermion sparse ED in the (N, Sz) sector must reproduce PySCF CASCI e_tot (< 1e-8 Ha)."""
    asp = request.getfixturevalue(fx)
    q = to_qubit_operator(asp, "jw")
    sz = (asp.nelecas[0] - asp.nelecas[1]) / 2
    ed = sparse_ground_state(q, asp.n_spin_orbitals, n_electrons=asp.n_electrons, sz=sz)
    assert abs(ed["e0"] - asp.e_casci) < 1e-8


@pytest.mark.parametrize("fx", ["h2", "lih44", "n2_stretched"])
def test_lambda_identity(fx, request):
    """Koridon 1-norm from integrals equals the explicit JW Pauli 1-norm for real orbitals."""
    asp = request.getfixturevalue(fx)
    st = pauli_stats(to_qubit_operator(asp, "jw"), asp.n_spin_orbitals)
    lam = lambda_from_integrals(asp.h1, asp.eri)
    assert abs(st["lambda_pauli"] - lam["lambda_total"]) < 1e-9 * max(1.0, st["lambda_pauli"])


@pytest.mark.parametrize("fx", ["h2", "lih44"])
def test_jw_bk_same_spectrum(fx, request):
    asp = request.getfixturevalue(fx)
    nq = asp.n_spin_orbitals
    e_jw = sparse_ground_state(to_qubit_operator(asp, "jw"), nq, k=4)["energies"]
    e_bk = sparse_ground_state(to_qubit_operator(asp, "bk"), nq, k=4)["energies"]
    assert np.allclose(e_jw, e_bk, atol=1e-9)


@pytest.mark.parametrize("fx", ["h2", "lih44"])
def test_tapering_preserves_ground_state(fx, request):
    asp = request.getfixturevalue(fx)
    nq = asp.n_spin_orbitals
    Ht, info = taper_hamiltonian(to_qubit_operator(asp, "jw"), asp.n_electrons, nq)
    assert info["n_qubits_tapered"] < nq
    e = np.min(np.linalg.eigvalsh(Ht.sparse_matrix().toarray()))
    assert abs(e - asp.e_casci) < 1e-8


def test_uccsd_theta_zero_is_hf(h2):
    """With all angles zero the UCCSD circuit is the HF determinant: <H> = E_HF (deterministic)."""
    nq, ne = h2.n_spin_orbitals, h2.n_electrons
    H = to_pennylane(to_qubit_operator(h2, "jw"))
    s, d, sw, dw, hf = uccsd_setup(ne, nq)
    dev = qml.device("lightning.qubit", wires=nq)

    @qml.qnode(dev)
    def e(theta):
        qml.UCCSD(theta, wires=range(nq), s_wires=sw, d_wires=dw, init_state=hf)
        return qml.expval(H)

    assert abs(float(e(pnp.zeros(len(s) + len(d)))) - h2.e_hf) < 1e-8


def test_uccsd_reference_parameters(h2):
    """Energy at stored reference angles reproduces the stored energy (regression, no optimizer)."""
    ref = json.load(open(FIX / "h2_uccsd_reference.json"))
    nq, ne = h2.n_spin_orbitals, h2.n_electrons
    H = to_pennylane(to_qubit_operator(h2, "jw"))
    s, d, sw, dw, hf = uccsd_setup(ne, nq)
    dev = qml.device("lightning.qubit", wires=nq)

    @qml.qnode(dev)
    def e(theta):
        qml.UCCSD(theta, wires=range(nq), s_wires=sw, d_wires=dw, init_state=hf)
        return qml.expval(H)

    val = float(e(pnp.array(ref["theta"])))
    assert abs(val - ref["energy"]) < 1e-8
    assert abs(val - h2.e_casci) < 1e-8     # H2/STO-3G UCCSD is exact
