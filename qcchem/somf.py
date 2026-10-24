"""Spin-orbit mean-field (SOMF / AMFI) two-electron screening in the Breit-Pauli form, reusing the
validated integral drivers of pyscf.prop.gtensor (Neese, JCP 122, 034107 (2005), eq. 19).

Operators (real antisymmetric AO matrices h_x, x = 1..3, "imaginary part of i sigma . (pV x p)"):
    one-electron:  h1e_x   = (alpha^2/4) <mu| (p x V_nuc p)_x |nu>                 (int1e_pnucxp)
    two-electron:  h2e_x   = (alpha^2/4) [ J_x[rho] - 1.5 K_x[rho] ]                (int2e_p1vxp1, SOMF)
                   AMFI    = same, restricted to one-centre (atomic) blocks
The effective one-body SOC operator in a spin-orbital basis is
    H_SO = sum_x (i h_x) (x) sigma_x
(the pyscf matrices already carry the 1/2 of s = sigma/2).  The spin factor was fixed by comparing the
one-electron BP 2p1/2-2p3/2 splitting of hydrogenic ions with the analytic Dirac value: with sigma_x/2 the
ratio was exactly 0.50 (tests/unit/test_somf.py).

Use: the Breit-Pauli *two-electron* mean-field term is added to the X2C one-electron SOC block of
qcchem.soc (X2C-1e + BP-2e-mean-field approximation); the Breit-Pauli one-electron term is used only for
consistency checks against X2C.
"""
from __future__ import annotations

import numpy as np
from pyscf import gto
from pyscf.data import nist
from pyscf.prop.gtensor.uhf import get_jk, get_jk_amfi

ALPHA2 = nist.ALPHA ** 2
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
         np.array([[1, 0], [0, -1]], complex)]


def bp_soc_1e_ao(mol: gto.Mole) -> np.ndarray:
    """(3, nao, nao) real antisymmetric one-electron Breit-Pauli SOC matrices (alpha^2/4 included)."""
    return (ALPHA2 / 4) * mol.intor_asymmetric("int1e_pnucxp", 3)


def bp_soc_2e_ao(mol: gto.Mole, dm_total: np.ndarray, amfi: bool = False) -> np.ndarray:
    """(3, nao, nao) SOMF two-electron mean-field SOC matrices for a (closed-shell-like) total density
    rho (alpha^2/4 included).  rho is split evenly between the two spin channels as required by get_jk."""
    dm0 = (0.5 * dm_total, 0.5 * dm_total)
    vj, vk = (get_jk_amfi if amfi else get_jk)(mol, dm0)
    return (ALPHA2 / 4) * ((vj[0] + vj[1]) - 1.5 * (vk[0] + vk[1]))


def soc_block_spinorb(h_x: np.ndarray, mo: np.ndarray) -> np.ndarray:
    """Transform 3 real antisymmetric AO SOC matrices into the interleaved spin-orbital block
    H[(p,s),(q,s')] = sum_x i (mo^T h_x mo)[p,q] (sigma_x)[s,s'] (complex Hermitian, 2n x 2n)."""
    n = mo.shape[1]
    out = np.zeros((2 * n, 2 * n), dtype=complex)
    for x in range(3):
        hx = mo.T @ h_x[x] @ mo                      # real antisymmetric
        out += np.kron(1j * hx, PAULI[x])
    return 0.5 * (out + out.conj().T)
