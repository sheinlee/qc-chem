"""Validate the one-electron X2C spin-orbit block against the analytic Dirac fine structure of a
hydrogen-like ion.  For one electron there is no two-electron screening, so the 2p1/2 - 2p3/2 splitting
obtained from our (scalar h ⊗ 1 + ΔSOC) one-body matrix must match Dirac up to basis-set incompleteness."""
import numpy as np
import pytest
from pyscf import gto
from pyscf.data import nist

from qcchem.soc import x2c_oneelectron_in_spin_ao_basis, scalar_x2c_hcore


def dirac_energy(Z: float, n: int, j: float, c: float) -> float:
    za = Z / c
    k = j + 0.5
    delta = k - np.sqrt(k * k - za * za)
    return c * c / np.sqrt(1.0 + (za / (n - delta)) ** 2) - c * c


@pytest.mark.parametrize("Z,symbol", [(10, "Ne"), (18, "Ar")])
def test_2p_fine_structure_matches_dirac(Z, symbol):
    c = nist.LIGHT_SPEED
    # large even-tempered s/p basis: a neutral-atom contracted basis misdescribes the compact 2p of a
    # highly charged hydrogenic ion (2p energies off by ~1 Ha), so the fine structure cannot be compared
    basis = {symbol: gto.etbs([(0, 26, 0.02, 2.4), (1, 26, 0.02, 2.4)])}
    mol = gto.M(atom=f"{symbol} 0 0 0", basis=basis, charge=Z - 1, spin=1, verbose=0)
    # scalar 2p orbitals = the three lowest p-type eigenvectors of the spin-free one-electron Hamiltonian
    # (for one electron h alone is the Hamiltonian; the ROHF Fock operator would add the 1s electron's field)
    import scipy.linalg as sl
    s = mol.intor_symmetric("int1e_ovlp")
    hsf0 = scalar_x2c_hcore(mol)
    es, Cs = sl.eigh(hsf0, s)
    p_mask = np.array([("p" in l.split()[2]) for l in mol.ao_labels()])
    pop_p = np.array([(Cs[:, i] * (s @ Cs[:, i]))[p_mask].sum() for i in range(mol.nao)])
    idx = [i for i in range(mol.nao) if pop_p[i] > 0.99][:3]
    C = Cs[:, idx]
    h2c = x2c_oneelectron_in_spin_ao_basis(mol)
    hsf = scalar_x2c_hcore(mol)
    nao = mol.nao
    Cso = np.zeros((2 * nao, 6), dtype=complex)
    Cso[:nao, 0::2] = C
    Cso[nao:, 1::2] = C
    h_so = Cso.conj().T @ h2c @ Cso                      # 6x6 one-electron Hamiltonian incl. SOC
    e = np.linalg.eigvalsh(h_so)
    split = e[2:].mean() - e[:2].mean()                  # 2p3/2 (x4) above 2p1/2 (x2)
    assert np.allclose(e[:2], e[0], atol=1e-6) and np.allclose(e[2:], e[2], atol=1e-6)
    exact = dirac_energy(Z, 2, 1.5, c) - dirac_energy(Z, 2, 0.5, c)
    # one-electron ion: no screening; X2C is exact for the Dirac one-electron problem within the basis
    assert abs(split / exact - 1.0) < 0.01, (split, exact)
    # absolute 2p1/2 energy of the projected block vs Dirac: the projection onto one scalar 2p radial
    # function neglects the j-dependent radial relaxation (1e-4 relative for Ar17+)
    assert abs(e[0] / dirac_energy(Z, 2, 0.5, c) - 1.0) < 1e-4
    # scalar part consistency: spin-averaged 2c one-body block ~ sfx2c one-body block
    hsf_mo = C.T @ hsf @ C
    avg = 0.5 * (h_so[0::2, 0::2] + h_so[1::2, 1::2])
    assert np.max(np.abs(avg - hsf_mo)) < 1e-3
