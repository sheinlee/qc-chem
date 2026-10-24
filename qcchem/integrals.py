"""PySCF driver: molecule -> SCF -> active space (CASCI/CASSCF, optional AVAS) -> folded integrals.

The folded active-space Hamiltonian is
    H = ecore + sum_pq h1[p,q] a+_p a_q + 1/2 sum_pqrs eri[p,q,r,s] a+_p a+_r a_s a_q   (chemist (pq|rs))
with spatial orbitals p,q,r,s.  `ecore` includes nuclear repulsion and the frozen-core energy, so
exact diagonalisation of this Hamiltonian in the (nalpha, nbeta) sector equals PySCF CASCI e_tot.
"""
from __future__ import annotations

import dataclasses
import json
import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from pyscf import ao2mo, gto, mcscf, scf
from pyscf.mcscf import avas as avas_mod


@dataclass
class ActiveSpace:
    name: str
    h1: np.ndarray            # (ncas, ncas)
    eri: np.ndarray           # (ncas, ncas, ncas, ncas), chemist notation, 8-fold symmetric
    ecore: float
    ncas: int
    nelecas: tuple[int, int]  # (nalpha, nbeta)
    e_hf: float
    e_casci: float
    e_casscf: Optional[float] = None
    meta: dict = field(default_factory=dict)

    @property
    def n_spin_orbitals(self) -> int:
        return 2 * self.ncas

    @property
    def n_electrons(self) -> int:
        return sum(self.nelecas)

    def save(self, path: str) -> None:
        np.savez_compressed(
            path, h1=self.h1, eri=self.eri, ecore=self.ecore, ncas=self.ncas,
            nelecas=np.array(self.nelecas), e_hf=self.e_hf, e_casci=self.e_casci,
            e_casscf=np.nan if self.e_casscf is None else self.e_casscf,
            name=self.name, meta=json.dumps(self.meta, default=_json_default),
        )

    @staticmethod
    def load(path: str) -> "ActiveSpace":
        d = np.load(path, allow_pickle=False)
        e_casscf = float(d["e_casscf"])
        return ActiveSpace(
            name=str(d["name"]), h1=d["h1"], eri=d["eri"], ecore=float(d["ecore"]),
            ncas=int(d["ncas"]), nelecas=tuple(int(x) for x in d["nelecas"]),
            e_hf=float(d["e_hf"]), e_casci=float(d["e_casci"]),
            e_casscf=None if np.isnan(e_casscf) else e_casscf,
            meta=json.loads(str(d["meta"])),
        )


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def build_mol(atom: str, basis, charge: int = 0, spin: int = 0, verbose: int = 0, **kw) -> gto.Mole:
    """basis may be a PySCF basis name/dict, or 'bse:<name>' to fetch from basis-set-exchange."""
    if isinstance(basis, str) and basis.startswith("bse:"):
        import basis_set_exchange as bse
        bname = basis[4:]
        elements = sorted({a.split()[0] for a in atom.replace(";", "\n").splitlines() if a.strip()})
        basis = {el: gto.basis.parse(bse.get_basis(bname, elements=[el], fmt="nwchem")) for el in elements}
    return gto.M(atom=atom, basis=basis, charge=charge, spin=spin, verbose=verbose, **kw)


def run_scf(mol: gto.Mole, x2c: bool = False, conv_tol: float = 1e-10, max_cycle: int = 200,
            frac_occ: bool = False):
    """RHF/ROHF, optionally scalar-relativistic (sfx2c1e).  `frac_occ=True` smears the open-shell
    electrons evenly over the degenerate frontier orbitals (spherically averaged reference for free
    ions: avoids the symmetry-broken core polarisation that splits the 4f levels)."""
    mf = scf.RHF(mol) if mol.spin == 0 else scf.ROHF(mol)
    if x2c:
        mf = mf.sfx2c1e()          # scalar-relativistic exact two-component (spin-free)
    if frac_occ:
        mf = scf.addons.frac_occ(mf)
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    mf.kernel()
    if not mf.converged and not frac_occ:
        mf = mf.newton().run()
    return mf


def orbital_composition(mol: gto.Mole, mo_coeff: np.ndarray, mo_idx, labels: list[str]) -> dict:
    """Mulliken-type AO population of each listed MO on the given AO label patterns (e.g. 'Ce 4f')."""
    s = mol.intor_symmetric("int1e_ovlp")
    ao_labels = mol.ao_labels()
    out = {}
    for i in mo_idx:
        c = mo_coeff[:, i]
        pop = c * (s @ c)          # Mulliken populations per AO (sum = 1 for normalised MO)
        comp = {}
        for lab in labels:
            mask = np.array([lab in al for al in ao_labels])
            comp[lab] = float(pop[mask].sum())
        comp["other"] = float(1.0 - sum(comp.values()))
        out[int(i)] = comp
    return out


def relax_open_shell_reference(mol, mf, labels: list[str], nelecas: int, ncas: Optional[int] = None,
                               minao: str = "minao", threshold: float = 0.2, fix_spin: Optional[float] = None):
    """Preliminary CASSCF in the AVAS space spanned by `labels` plus all singly-occupied orbitals
    (openshell_option=3), with `nelecas` active electrons.  Returns an SCF-like object whose mo_coeff /
    mo_occ describe the CASSCF determinant-like occupation (core doubly occupied, active natural
    orbitals with their CASSCF occupations), suitable as input to a second AVAS."""
    t0 = time.time()
    ncas_p, nelec_p, mo_p = avas_mod.avas(mf, labels, threshold=threshold, minao=minao, openshell_option=3,
                                          canonicalize=True)
    ncas_p = int(ncas_p) if ncas is None else int(ncas)
    na = (nelecas + mol.spin) // 2
    nb = nelecas - na
    mc = mcscf.CASSCF(mf, ncas_p, (na, nb))
    if fix_spin is not None:
        mc.fix_spin_(ss=fix_spin)
    mc.kernel(mo_p)
    # natural orbitals + occupations in the active space
    mo_nat, ci_nat, occ_full = mc.cas_natorb(sort=True)      # occ_full has length nmo
    ncore = mc.ncore
    occ_nat = np.asarray(occ_full)[ncore:ncore + ncas_p]
    mf2 = mf.copy()
    mf2.mo_coeff = mo_nat
    # AVAS needs integer-like occupations to split occupied/virtual: round the natural occupations
    mf2.mo_occ = np.where(occ_full > 1.5, 2.0, np.where(occ_full > 0.5, 1.0, 0.0))
    info = {"labels": labels, "ncas": ncas_p, "nelecas": nelecas, "e_casscf": float(mc.e_tot),
            "converged": bool(mc.converged), "natural_occupations": [round(float(x), 4) for x in occ_nat],
            "n_singly_occupied_after": int(np.sum(mf2.mo_occ == 1.0)),
            "e_scf_before": float(mf.e_tot), "t_seconds": round(time.time() - t0, 1)}
    return mf2, info


def make_active_space(
    name: str,
    mol: gto.Mole,
    ncas: Optional[int] = None,
    nelecas=None,
    x2c: bool = False,
    avas_labels: Optional[list[str]] = None,
    avas_threshold: float = 0.2,
    avas_minao: str = "minao",
    avas_openshell_option: int = 2,
    relax_reference: Optional[dict] = None,
    casscf: bool = False,
    fix_spin: Optional[float] = None,
    composition_labels: Optional[list[str]] = None,
    solve_ci: bool = True,
    ci_max_memory_mb: int = 16000,
    frac_occ: bool = False,
    state_average: int = 0,
    mf=None,
) -> ActiveSpace:
    """Build an active space and fold integrals.

    Either give (ncas, nelecas) for a canonical-orbital CAS around the Fermi level, or `avas_labels`
    (list of AO label patterns such as ['Ce 4f'] or ['Ce 4f', 'Ce 5d', 'Cl 3p']) to let AVAS pick the
    orbitals and their count.  `fix_spin` (S(S+1) value) penalises other spin states in the CI solver.
    """
    t0 = time.time()
    if mf is None:
        mf = run_scf(mol, x2c=x2c, frac_occ=frac_occ)
    meta: dict = {"basis": str(mol.basis) if isinstance(mol.basis, str) else "dict", "x2c": x2c, "frac_occ": frac_occ,
                  "charge": mol.charge, "spin": mol.spin, "nao": mol.nao, "scf_converged": bool(mf.converged),
                  "atom": mol.atom if isinstance(mol.atom, str) else str(mol.atom)}

    if relax_reference:
        # HF often puts the open-shell electron(s) of a lanthanide in 5d/6s instead of 4f.  Relax the
        # reference by a preliminary CASSCF in (target labels + SOMOs), then re-partition occupations.
        mf, relax_info = relax_open_shell_reference(mol, mf, **relax_reference, minao=avas_minao)
        meta["relax_reference"] = relax_info

    mo = mf.mo_coeff
    if avas_labels:
        ncas_avas, nelecas_avas, mo = avas_mod.avas(mf, avas_labels, threshold=avas_threshold, minao=avas_minao,
                                                    openshell_option=avas_openshell_option, canonicalize=True)
        if ncas is None:
            ncas = int(ncas_avas)
        nelecas = int(nelecas_avas) if nelecas is None else nelecas
        meta["avas"] = {"labels": avas_labels, "threshold": avas_threshold, "minao": avas_minao,
                        "ncas": int(ncas_avas), "nelecas": int(nelecas_avas)}
    assert ncas is not None and nelecas is not None, "need ncas/nelecas or avas_labels"

    if isinstance(nelecas, (int, np.integer)):
        na = (int(nelecas) + mol.spin) // 2
        nb = int(nelecas) - na
    else:
        na, nb = int(nelecas[0]), int(nelecas[1])

    mc = (mcscf.CASSCF if casscf else mcscf.CASCI)(mf, ncas, (na, nb))
    mc.max_memory = ci_max_memory_mb
    if fix_spin is not None:
        mc.fix_spin_(ss=fix_spin)
    if casscf and state_average and state_average > 1:
        # equal-weight state averaging over the lowest `state_average` CI roots: for a free ion averaged
        # over a complete multiplet this yields spherically symmetric orbitals by construction
        mc = mc.state_average_([1.0 / state_average] * state_average)
        meta["state_average"] = int(state_average)
    if not solve_ci:
        # integrals only (CI space too large): e_casci is NaN, orbitals are the AVAS/canonical ones
        mc.mo_coeff = mo
        h1, ecore = mc.get_h1eff(mo)
        eri = ao2mo.restore(1, mc.get_h2eff(mo), ncas)
        meta.update(ncore=int(mc.ncore), ci_solved=False, t_seconds=round(time.time() - t0, 2))
        if composition_labels:
            act = range(mc.ncore, mc.ncore + ncas)
            meta["active_orbital_composition"] = orbital_composition(mol, mo, act, composition_labels)
        return ActiveSpace(name=name, h1=np.ascontiguousarray(h1), eri=np.ascontiguousarray(eri),
                           ecore=float(ecore), ncas=int(ncas), nelecas=(na, nb),
                           e_hf=float(mf.e_tot), e_casci=float("nan"), e_casscf=None, meta=meta)
    mc.kernel(mo)
    e_casscf = float(np.min(np.atleast_1d(mc.e_tot))) if casscf else None
    if casscf and state_average and state_average > 1:
        meta["sa_casscf_state_energies"] = [float(e) for e in np.atleast_1d(mc.e_states)]
        meta["sa_casscf_converged"] = bool(mc.converged)
    if casscf:
        # final CASCI in the optimised *natural* orbitals sorted by occupation (same energy), so that the
        # reference determinant used by UCCSD (first na alpha / nb beta active orbitals) is the leading one
        mo, _, occ_full = mc.cas_natorb(sort=True)
        meta["natural_occupations"] = [round(float(x), 4) for x in np.asarray(occ_full)[mc.ncore:mc.ncore + ncas]]
        mc = mcscf.CASCI(mf, ncas, (na, nb))
        if fix_spin is not None:
            mc.fix_spin_(ss=fix_spin)
        mc.kernel(mo)

    h1, ecore = mc.get_h1eff()
    eri = ao2mo.restore(1, mc.get_h2eff(), ncas)
    meta["ncore"] = int(mc.ncore)
    meta["casci_converged"] = bool(getattr(mc, "converged", True))
    meta["t_seconds"] = round(time.time() - t0, 2)
    try:
        meta["spin_square"] = [float(x) for x in mc.fcisolver.spin_square(mc.ci, ncas, (na, nb))]
    except Exception:
        pass
    if composition_labels:
        act = range(mc.ncore, mc.ncore + ncas)
        meta["active_orbital_composition"] = orbital_composition(mol, mc.mo_coeff, act, composition_labels)
        meta["active_orbital_energies"] = [float(x) for x in np.asarray(mf.mo_energy)[list(act)]] \
            if mf.mo_energy is not None and not avas_labels else None

    return ActiveSpace(name=name, h1=np.ascontiguousarray(h1), eri=np.ascontiguousarray(eri),
                       ecore=float(ecore), ncas=int(ncas), nelecas=(na, nb),
                       e_hf=float(mf.e_tot), e_casci=float(mc.e_tot), e_casscf=e_casscf, meta=meta)
