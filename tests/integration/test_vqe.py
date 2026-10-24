"""Optimizer actually runs. Tolerances are loose enough that optimizer noise cannot fail them."""
import pytest

from qcchem.hamiltonian import to_pennylane, to_qubit_operator
from qcchem.vqe import run_uccsd_vqe

CHEM_ACC = 1.6e-3


@pytest.mark.parametrize("fx,tol", [("h2", 1e-6), ("lih44", 1e-5)])
def test_uccsd_vqe_reaches_fci(fx, tol, request):
    asp = request.getfixturevalue(fx)
    v = run_uccsd_vqe(to_pennylane(to_qubit_operator(asp, "jw")), asp.n_electrons, asp.n_spin_orbitals,
                      with_resources=False)
    err = v["energy"] - asp.e_casci
    assert err > -1e-9, "variational bound violated"
    assert abs(err) < tol
    assert abs(err) < CHEM_ACC
