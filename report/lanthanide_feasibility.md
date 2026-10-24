# Quantum algorithms for lanthanide electronic structure: a feasibility study with our own numbers

2026-10-24, second revision: after the independent referee-style review, the two scientific gaps it exposed
(an orbital-consistent CeCl3 ladder with a 4f + σ control, and two-electron screening of the spin-orbit
operator) were addressed. Companion to `technical_report.md`; the chronology, failed attempts and
diagnostics are in `progress.md`. Scalar-relativistic sfx2c1e, cc-pVDZ-DK, PySCF 2.14 / OpenFermion 1.8.1 /
PennyLane 0.45.1 / pyscf-properties 0.1.

## 1. The claim being tested

`PLAN.md` set out the claim that *4f-only lanthanide active spaces are classically trivial; the problem only
becomes hard once ligand orbitals and relativistic/spin-orbit effects enter, and at that point NISQ-style
VQE is already out of reach, leaving only fault-tolerant phase estimation to discuss.* The calculations were
set up so that the numbers, not the claim, decide. Section 6 states what survived.

## 2. Free ions: the reference problem comes before the quantum problem

Ce3+ 4f1, Pr3+ 4f2, Eu3+ 4f6, Gd3+ 4f7 in the 4f-only CAS(n,7) space, 14 qubits.

**What went wrong first.** Hartree-Fock is an unreliable starting point for these ions: ROHF put two of
Gd3+'s electrons into 5d (CASSCF in the proper 4f7 space is 415 mHa lower), converged Pr3+ to a wrong
configuration (78 mHa above CASSCF), and for CeCl3 landed on a 5d1 solution that is 118 mHa *above* the 4f1
ROHF solution (which exists and is found once the density is started from a 4f1 configuration). A
single-determinant reference also breaks spherical symmetry: with ROHF/CASSCF orbitals the seven 4f levels
of Ce3+ were split by 5–14 mHa, Eu3+'s by 10–18 mHa, and two CASSCF runs for Pr3+ landed on solutions 1.5 mHa
apart. A fractional-occupation reference did not help (CASSCF(1,7) leaves the six unoccupied active orbitals
undetermined). **What worked:** state-averaged (SA) CASSCF over the complete high-spin multiplet manifold
(7 states for Ce3+ and Eu3+, 21 for Pr3+). Restoring spherical symmetry costs 3.5 mHa (Ce3+) and 4.5 mHa
(Eu3+) relative to the symmetry-broken ROHF energy, but the sparse ED of the resulting qubit Hamiltonian
then reproduces the exact term structure: Ce3+ 2F seven-fold degenerate; Pr3+ 3H (11) at 0, 3F (7) at
27.19 mHa, 3P (3) at 126.88 mHa (SA-CASSCF state energies and qubit ED agree state by state); Eu3+ 7F (7)
degenerate.

**How hard is 4f-only?** Not hard. The fixed-S_z CI spaces are 14 (Ce3+), 21 (Pr3+), 7 (Eu3+) and 1 (Gd3+)
determinants. UCCSD-VQE is exact where it applies: Ce3+ 7.8e-9 Ha with 6 parameters and 168 CNOTs; Pr3+
−2.3e-8 Ha (below the CASCI reference, whose convergence on the degenerate root is 8e-8; the ED is the
better reference) with 20 parameters and 1400 CNOTs. For Ce3+ every two-electron term vanishes in the N = 1
sector: the physical problem is a 7×7 matrix whatever the Fock-space Hamiltonian looks like.

**The Hamiltonian looks dense, but that is a basis statement.** In SA natural orbitals (degenerate, hence
arbitrarily rotated within the 4f shell) the Ce3+ Hamiltonian has ~3300 JW terms, λ ≈ 55 Ha and two Z2
symmetries (14 → 12 qubits); in the ROHF-based natural orbitals of the same ion it has 845 terms, λ = 16.4 Ha
and four Z2 symmetries (14 → 10 qubits), with the same spectrum. Term counts, λ and the tapering yield are
orbital-choice artefacts for a degenerate shell; λ-minimising orbital rotations (Koridon et al. 2021) are a
first-order lever for Ln Hamiltonians, before any algorithmic choice.

## 3. CeCl3: an active-space ladder on one reference, with the 4f + σ control

Planar D3h CeCl3, Ce–Cl 2.57 Å, 138 basis functions. **Reference protocol (identical for every level):**
ROHF re-converged in the 4f1 configuration (E_ROHF = −10239.411130 Ha), AVAS with ANO-RCC atomic projectors,
CASCI in the AVAS orbitals, CI natural occupations and natural-orbital compositions recorded. The 4f + σ
control is built from the 4f+5d+σ AVAS space by removing the five unoccupied 5d-dominated orbitals, so C3 ⊂
L3 and the two are variationally comparable.

**Why CASCI and not CASSCF.** A fully variational CASSCF does not preserve the chemically defined spaces
(`progress.md`): without a frozen core, CASSCF(11,12) rotates the five Cl σ orbitals into Ce 4d semicore
orbitals (228 mHa of 4d–4f core-valence correlation that has nothing to do with the ligand); with a frozen
core it swaps the six empty 4f orbitals for ligand σ* orbitals (22.8 mHa of σ→σ* pair correlation; active-space
4f weight drops from 7 to 1.04); and for one active electron the unoccupied active orbitals have zero gradient
and drift arbitrarily. A frozen core taken from the state-specific 4f1 ROHF also polarises the 4f manifold
(SA-CASSCF splittings inflate from 0.7–6 to 6–19 mHa). Constraining CASSCF enough to prevent all of this
reintroduces the arbitrariness it was meant to remove, so the comparison below keeps the orbitals fixed and
inspected and lets the CI occupations speak.

| Level | Active space | qubits | E_CASCI − E_ROHF / mHa | CI natural occupations (beyond 2/1/0) | JW terms | λ / Ha |
|---|---|---|---|---|---|---|
| L1 | 4f only, CAS(1,7) | 14 (12) | +3.53 | — (one electron) | 1409 | 56.1 |
| L1-SA | same, SA-CASSCF(7) orbitals, unfrozen core | 14 (12) | +2.69 | — | 3365 | 59.7 |
| L2 | 4f + 5d, CAS(1,12) | 24 (22) | +3.53 | — (one electron) | 13596 | 116.1 |
| **C3 (control)** | **4f + Cl σ donors, no 5d, CAS(11,12)** | 24 (22) | **+3.44** | σ 2.0000 ×5, 4f 1.0000, 4f virtuals 0.0000 ×6 | 13040 | 65.2 |
| L3 | 4f + 5d + Cl σ donors, CAS(11,17) | 34 | **−3.01** | σ 1.997–1.999, 4f 0.9994, **5d 0.0025/0.0025/0.0022/0.0022/0.0010**, 4f virtuals 0.0001 | — | 144.6 |
| L4 | 4f + all nine Cl 3p, no 5d, CAS(19,16) | 32 | +0.02 | Cl 3p 1.9999–2.0000 ×9, 4f 0.9999, 4f virtuals ≤ 0.0001 | — | 92.6 |
| L5 | 4f + 5d + all Cl 3p, CAS(19,21) | 42 | not solved (1.1e11) | — | — | 178.9 |

(Positive entries for L1/L2/C3 are the cost of projecting the ROHF SOMO onto the ANO 4f function; L4's
CASCI, 9.2e7 determinants, ran in 12 minutes on the workstation once the cluster became unreachable.)

What the ladder shows, within the stated protocol.

1. **The ligand electrons do not correlate with the 4f shell at all.** With the five σ donors and the six
   empty 4f orbitals in the same CASCI space (C3), the natural occupations are 2.0000/1.0000/0.0000 to four
   decimals and the energy gain over 4f-only is 0.09 mHa; with all nine Cl 3p orbitals (L4) the 4f virtual
   occupations stay ≤ 1e-4 and the gain is 0.02 mHa. The 4f virtuals are not a correlation degree of freedom
   for the ligand electrons.
2. **Whatever metal–ligand correlation exists in this description goes through 5d.** Adding the five 5d
   orbitals to C3 (giving L3, a nested space) lowers the energy by 6.45 mHa and is the only step on the ladder
   that gets below the ROHF energy; the 5d natural occupations sum to 0.010 electrons, a hundred times the 4f
   virtual occupations but still small in absolute terms. So 5d is necessary in the active-space description
   of CeCl3 (it carries all of the ligand-to-metal correlation/donation), while its magnitude (0.01 e, 6 mHa)
   marks it as a polarisation/donation channel rather than a strongly correlated degree of freedom. This is a
   statement about the CASCI description in this basis, not a general covalency mechanism.
3. **The 4f ligand-field splitting is small and reference-sensitive.** With state-averaged orbitals the 4f
   manifold splits as 0 / 0.68 (×2) / 1.46 (×2) / 5.72 / 5.96 mHa (150–1300 cm⁻¹), the 1-2-2-1-1 pattern of
   f orbitals in D3h; with single-determinant 4f orbitals the same CASCI gives 0 / 7.4 (×2) / 10.7 (×2) /
   20.3 / 21.2 mHa, i.e. residual core polarisation exceeds the ligand field five-fold.
4. **λ and term counts along the ladder.** λ per qubit is roughly constant (4.0 at L1, 4.3 at L3 and L5); L4
   (32 qubits, 92.6 Ha) sits below L2 (24 qubits, 116 Ha) because λ tracks which orbitals are present (5d
   contributes large exchange integrals) more than the qubit count. The JW term count over all systems
   follows N^4.35 in the number of active orbitals (fig. 3).

## 4. Scalar versus two-component: spin-orbit coupling, one-electron and screened

Model: the scalar active space plus a spin-orbit one-body block in the same orbitals. Two versions:
(a) one-electron X2C spin-orbit block (`qcchem/soc.py`), validated against the analytic Dirac 2p fine
structure of hydrogenic Ne9+/Ar17+ (< 1 %); (b) the same plus the Breit-Pauli spin-orbit mean field (SOMF,
Neese 2005) of the state density, built with the integral drivers of `pyscf.prop.gtensor` (`qcchem/somf.py`).
The SOMF convention was fixed by the hydrogenic test (a factor-2 error in the spin operator showed up as a
ratio of exactly 0.50 to Dirac), and the Breit-Pauli one-electron block agrees with the X2C one to 0.5 % on
the Ce 4f shell, cross-validating both implementations. Full N-sector ED on the SA 14-qubit spaces:

| System | one-electron SOC: levels (cm⁻¹) | + SOMF screening: levels (cm⁻¹) | experiment (cm⁻¹) | screening factor |
|---|---|---|---|---|
| Ce3+ 4f1 | 2F5/2 (6) 0; 2F7/2 (8) 5266 | 0; **2397** | 2253 | 0.455 |
| Pr3+ 4f2 | 3H4 (9) 0; 3H5 (11) 5169; 3F2 (5) 7582; 3H6 (13) 10163; 3F3 (7) 10873 | 0; **2259**; 4643 (3H6) | 3H5 2152; 3H6 4389 | 0.466 |
| CeCl3 L1-SA | Kramers doublets 0 / 189 / 1033 (J = 5/2); 5116 / 5250 / 5536 / 6331 (J = 7/2) | 0 / 191 / 996; **2320 / 2454 / 2778 / 3535** | — | 0.455 |

Findings. (i) The one-electron operator overestimates the 2F7/2–2F5/2 interval 2.3×; the two-electron mean
field screens the 4f spin-orbit coupling by a factor 0.46 and brings the intervals to within +5–6 % of
experiment (Ce3+ 2397 vs 2253; Pr3+ 3H5 2259 vs 2152, 3H6 4643 vs 4389). The remaining 5 % is consistent
with a mean-field treatment without dynamic correlation and a double-zeta basis. The reviewer's doubt
("may be a code error") is thus resolved: the transformation was correct and the physics was the screening.
(ii) Structurally, SOC costs ~7.5 % more Pauli strings in the SA orbital basis (182 one-body SOC elements →
~250 strings, 168 of them S_z-non-conserving, real coefficients), < 1 % in λ, and one Z2 symmetry (S_z
parity). Screening does not change this: the SOMF block has the same sparsity pattern, so term counts and λ
are unchanged to within 0.3 %. These percentages are representation-dependent (the SA basis maximises the
scalar term count, §2) and should not be quoted as constants. (iii) The screened SOC stabilisation of the
ground state relative to the scalar barycentre is 6.2 mHa (Ce3+), 12.1 mHa (Pr3+) and 5.7 mHa (CeCl3) —
halved from the one-electron values, but still of the same size as the entire scalar ligand-field spread
(≤ 6 mHa) and as the whole 5d-mediated correlation gain of the ladder (6.5 mHa). The CeCl3 two-component
levels show the ligand field and SOC on equal footing: the three J = 5/2-parentage Kramers doublets span
1000 cm⁻¹ and the J = 7/2 set is pushed up by ~2300 cm⁻¹.

## 5. Resource proxies, and what they do not mean

For each Hamiltonian the records contain λ, the qubitised-walk query count π λ / (2ε) for ε = 1 mHa, a
sparse-LCU Toffoli proxy (queries × Pauli terms) and a logical-qubit proxy. Along the ladder: L1 8.8e4
queries, 1.3e8 Toffolis, 42 logical qubits; L3 2.3e5 queries, ~2e10 Toffolis, 69 qubits; L5 2.8e5 queries,
~6e10 Toffolis, 79 qubits. Caveats: Fock-space 1-norms (for N = 1 spaces most terms vanish in the physical
sector); sparse-LCU numbers, whereas double factorisation / THC reduce both λ and the per-query cost (not
evaluated here; `openfermion.resource_estimates` is the next step); no error-correction overhead. They
compare Hamiltonians to each other; they are not hardware estimates.

## 6. What the numbers say about the original claim

- *4f-only is classically trivial*: confirmed, and sharper than expected — the hard part is the orbital
  reference, not the many-body problem; the Fock-space Hamiltonian's apparent density is an orbital artefact.
- *Ligand orbitals make it hard*: the chemically motivated space is 4f + 5d + σ (34 qubits), still classically
  solvable in minutes, and the ladder now says *why* the 5d orbitals have to be there: without them the ligand
  electrons do not interact with the metal shell at all in the CI (natural occupations 2/1/0), with them a
  0.01-electron donation channel carries 6.5 mHa. The classical wall is reached only at L5 (42 qubits, 1e11
  determinants) or with larger ligands and basis sets.
- *Spin-orbit makes it hard*: not for the Hamiltonian (+7.5 % terms, +0.5 % λ, −1 symmetry, unchanged by
  screening); yes for the physics. A one-electron SOC operator is 2.3× off; a mean-field two-electron
  screening (cheap, validated here) brings Ln multiplet intervals to ~5 %. The difficulty is in building the
  right reference and operator, which is classical pre-processing.
- *NISQ VQE is out of reach there*: supported by the molecular numbers (9024 CNOTs for 12 qubits in a naive
  compilation; UCCSD errors of 2–13 mHa in the multireference regime) but not proven: compact ansätze were
  not run and the CNOT count is unoptimised.
- *Only fault-tolerant QPE remains*: "fault-tolerant QPE is the only route whose cost can even be written down
  for L5, and in the sparse-LCU proxy it is ~6e10 Toffolis before factorisation"; whether it beats DMRG at
  21–30 orbitals is the open question.

The sentence the data support best: **for lanthanides the difficult part comes before the quantum solver —
a physically correct orbital reference (right shell occupation, spherical or ligand-field-adapted 4f
orbitals, 5d in the active space) and a screened relativistic Hamiltonian.** Once those are in hand, the
Hamiltonians of chemically meaningful Ln(III) halide spaces are 24–42 qubits with λ of 60–180 Ha.

## 7. Three research projects, ranked

1. **Two-component lanthanide Hamiltonians with SOMF screening, a spinor (Kramers) basis and
   symmetry-adapted tapering** (feasibility high — the screening step now exists and is validated; novelty
   high; cost low–medium). Map complex spinor Hamiltonians to qubits, derive which Z2/time-reversal
   symmetries survive and how to target a J-manifold; validate against the exact 2F/3H intervals.
2. **ADAPT-VQE and pool design for 4f exchange-dominated Hamiltonians, with λ-minimising orbital
   optimisation** (feasibility high, novelty medium, cost medium): fermionic vs qubit pools on the SA
   Ce3+/Pr3+/Eu3+ and CeCl3 L1 Hamiltonians against ED; §2 shows a 3× spread in λ from orbital choice alone.
3. **Classical-vs-fault-tolerant crossover on the CeCl3 ladder** (feasibility medium, novelty high, cost
   medium–high): DMRG (block2) references for L5 and larger spaces (bigger basis, EuCl3/GdCl3), double-
   factorised and THC λ for the same integrals, and a crossover plot of DMRG cost vs QPE Toffoli count.

## 8. Limitations

cc-pVDZ-DK basis; scalar X2C; one-electron X2C SOC plus Breit-Pauli SOMF screening (no spin-dependent
two-electron terms beyond the mean field; core frozen at the scalar level); AVAS selection thresholds as
stated; CASCI at every CeCl3 level (CASSCF does not preserve the chemically defined spaces, §3); no dynamic
correlation outside the active space; no noise; QPE proxies are 1-norm/sparse-LCU based; L5 has no classical
reference. Term counts, λ and tapering yields depend on the orbital basis for degenerate shells (§2). The
four-component Dirac-Coulomb(-Gaunt) HF benchmark of the 4f screening (`scripts/dhf_benchmark.py`) was
still running when this revision was written; its result will be added to `results/soc_dhf_benchmark.json`.
