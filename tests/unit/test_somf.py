"""Breit-Pauli one-electron SOC convention: hydrogenic 2p fine structure vs Dirac (BP is O(alpha^2))."""
import numpy as np
import pytest
import scipy.linalg as sl
from pyscf import gto
from pyscf.data import nist

from qcchem.somf import bp_soc_1e_ao, soc_block_spinorb


def dirac(Z, n, j, c):
    za = Z / c; k = j + 0.5; d = k - np.sqrt(k * k - za * za)
    return c * c / np.sqrt(1.0 + (za / (n - d)) ** 2) - c * c


@pytest.mark.parametrize("Z,sym", [(10, "Ne"), (18, "Ar")])
def test_bp_1e_2p_splitting(Z, sym):
    c = nist.LIGHT_SPEED
    mol = gto.M(atom=f"{sym} 0 0 0", basis={sym: gto.etbs([(0, 26, 0.02, 2.4), (1, 26, 0.02, 2.4)])},
                charge=Z - 1, spin=1, verbose=0)
    s = mol.intor_symmetric("int1e_ovlp")
    es, Cs = sl.eigh(mol.intor_symmetric("int1e_kin") + mol.intor_symmetric("int1e_nuc"), s)
    p_mask = np.array([("p" in l.split()[2]) for l in mol.ao_labels()])
    pop = np.array([(Cs[:, i] * (s @ Cs[:, i]))[p_mask].sum() for i in range(mol.nao)])
    C = Cs[:, [i for i in range(mol.nao) if pop[i] > 0.99][:3]]
    e = np.linalg.eigvalsh(soc_block_spinorb(bp_soc_1e_ao(mol), C))
    split = e[2:].mean() - e[:2].mean()
    exact = dirac(Z, 2, 1.5, c) - dirac(Z, 2, 0.5, c)
    assert abs(split / exact - 1.0) < 0.02          # BP vs Dirac: O((Z alpha)^2) difference, ~1 % for Ar17+
