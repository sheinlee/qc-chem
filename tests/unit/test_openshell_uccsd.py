"""Open-shell UCCSD setup: reference determinant and spin-conserving excitation counts."""
import numpy as np
import pennylane as qml

from qcchem.vqe import uccsd_setup


def test_closed_shell_matches_pennylane():
    s, d, _, _, hf = uccsd_setup(4, 8)
    s_ref, d_ref = qml.qchem.excitations(4, 8)
    assert len(s) == len(s_ref) and len(d) == len(d_ref)
    assert np.array_equal(hf, qml.qchem.hf_state(4, 8))


def test_high_spin_reference_and_counts():
    # Pr3+ 4f2 high-spin: 2 alpha electrons in 7 spatial orbitals (14 qubits), interleaved ordering
    s, d, _, _, hf = uccsd_setup((2, 0), 14)
    assert hf.tolist() == [1, 0, 1, 0] + [0] * 10
    assert len(s) == 2 * 5          # each alpha electron -> 5 alpha virtuals
    assert len(d) == 10             # alpha-alpha pairs into C(5,2) virtual pairs
    for i, a in s:
        assert i % 2 == a % 2       # spin conserved
    for i, j, a, b in d:
        assert (i % 2 + j % 2) == (a % 2 + b % 2)


def test_single_electron_has_no_doubles():
    s, d, _, _, hf = uccsd_setup((1, 0), 14)
    assert hf.tolist() == [1] + [0] * 13
    assert len(s) == 6 and len(d) == 0
