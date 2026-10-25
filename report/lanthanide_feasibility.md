# Quantum algorithms for lanthanide electronic structure: a feasibility study with our own numbers

2026-10-24, third revision. Changes since the second: the CeCl3 ladder is now built with the ROHF SOMO kept
unprojected (every level shares one exact baseline), nestedness of the active subspaces is checked explicitly,
the physical conclusions are stated within the protocol that supports them, λ is defined precisely for the
spin-orbit Hamiltonians, and factorised fault-tolerant resource estimates (sparse / SF / DF, OpenFermion
`resource_estimates`) and λ-minimising orbital rotations replace the earlier 1-norm proxies. Companion to
`technical_report.md`; chronology, failed attempts and diagnostics in `progress.md`. Scalar-relativistic
sfx2c1e, cc-pVDZ-DK, PySCF 2.14 / OpenFermion 1.8.1 / PennyLane 0.45.1 / pyscf-properties 0.1.

## 1. The claim being tested

`PLAN.md` set out the claim that *4f-only lanthanide active spaces are classically trivial; the problem only
becomes hard once ligand orbitals and relativistic/spin-orbit effects enter, and at that point NISQ-style
VQE is already out of reach, leaving only fault-tolerant phase estimation to discuss.* The calculations were
set up so that the numbers, not the claim, decide. Section 8 states what survived.

## 2. Free ions: the reference problem comes before the quantum problem

Ce3+ 4f1, Pr3+ 4f2, Eu3+ 4f6, Gd3+ 4f7 in the 4f-only CAS(n,7) space, 14 qubits.

**What went wrong first.** Hartree-Fock is an unreliable starting point for these ions: ROHF put two of
Gd3+'s electrons into 5d (CASSCF in the proper 4f7 space is 415 mHa lower), converged Pr3+ to a wrong
configuration (78 mHa above CASSCF), and for CeCl3 landed on a 5d1 solution that is 118 mHa *above* the 4f1
ROHF solution (which exists and is found once the density is started from a 4f1 configuration). A
single-determinant reference also breaks spherical symmetry: with ROHF/CASSCF orbitals the seven 4f levels
of Ce3+ were split by 5–14 mHa, Eu3+'s by 10–18 mHa, and two CASSCF runs for Pr3+ landed on solutions 1.5 mHa
apart. **What worked:** state-averaged (SA) CASSCF over the complete high-spin multiplet manifold (7 states
for Ce3+ and Eu3+, 21 for Pr3+), at a variational cost of 3.5–4.5 mHa relative to the symmetry-broken ROHF.
The sparse ED of the resulting qubit Hamiltonian then reproduces the exact term structure: Ce3+ 2F seven-fold
degenerate; Pr3+ 3H (11) at 0, 3F (7) at 27.19 mHa, 3P (3) at 126.88 mHa; Eu3+ 7F (7) degenerate.

**How hard is 4f-only?** Not hard. The fixed-S_z CI spaces are 14 (Ce3+), 21 (Pr3+), 7 (Eu3+) and 1 (Gd3+)
determinants. UCCSD-VQE is exact where it applies (Ce3+ 7.8e-9 Ha, 6 parameters, 168 CNOTs; Pr3+ −2.3e-8 Ha
vs CASCI, 20 parameters, 1400 CNOTs). For Ce3+ every two-electron term vanishes in the N = 1 sector, and
correspondingly every factorisation of the two-electron integrals below reproduces its CASCI energy exactly.

**The Hamiltonian looks dense, but that is a basis statement.** In SA natural orbitals the Ce3+ Hamiltonian
has ~3300 JW terms, λ ≈ 55 Ha and two Z2 symmetries (14 → 12 qubits); in the ROHF-based natural orbitals of
the same ion 845 terms, λ = 16.4 Ha and four Z2 symmetries (14 → 10 qubits), with the same spectrum.

## 3. CeCl3: an active-space ladder on one reference, with the 4f + σ control

Planar D3h CeCl3, Ce–Cl 2.57 Å, 138 basis functions. **Protocol (identical at every level):** ROHF
re-converged in the 4f1 configuration (E_ROHF = −10239.411130 Ha); AVAS with ANO-RCC projectors with the
ROHF SOMO kept *unprojected* in the active space (AVAS `openshell_option = 3`), so that a one-electron CAS
reproduces the ROHF energy exactly; CASCI with CI convergence 1e-12; fixed-orbital occupations (diagonal of
the active 1-RDM) and natural occupations recorded. The 4f + σ control (C3) is built from the 4f+5d+σ AVAS
space by removing the five unoccupied 5d-dominated orbitals, so C3 ⊂ L3 exactly.

**Why CASCI and not CASSCF.** A fully variational CASSCF does not preserve the chemically defined spaces
(`progress.md`): without a frozen core, CASSCF(11,12) rotates the five Cl σ orbitals into Ce 4d semicore
orbitals (228 mHa of 4d–4f core-valence correlation); with a frozen core it swaps the six empty 4f orbitals
for ligand σ* orbitals (22.8 mHa of σ→σ* pair correlation, active-space 4f weight 7 → 1.04); for one active
electron the unoccupied active orbitals have zero gradient; and a frozen core from the state-specific 4f1
ROHF polarises the 4f manifold (SA splittings inflate from 0.7–6 to 6–19 mHa). The comparison below therefore
keeps the orbitals fixed and inspected and lets the CI occupations speak.

| Level | Active space | qubits | E_CASCI − E_ROHF / mHa | fixed-orbital occupations of the non-trivial orbitals | JW terms | λ / Ha |
|---|---|---|---|---|---|---|
| L1 | 4f only, CAS(1,7) | 14 (12) | 0.000 | one electron | 1643 | 56.1 |
| L2 | 4f + 5d, CAS(1,12) | 24 (22) | 0.000 | one electron | 13596 | 116.1 |
| **C3 (control)** | **4f + Cl σ donors, no 5d, CAS(11,12)** | 24 (22) | **−0.093** | σ 1.99998 ×5; 4f virtuals 1–3e-5 | 13776 | 65.6 |
| L3 | 4f + 5d + Cl σ donors, CAS(11,17) | 34 | **−6.537** | σ 1.9972–1.9990; **5d 0.0010 / 0.0026 / 0.0026 / 0.0022 / 0.0022 (Σ 0.0105)**; 4f virtuals 0.6–2.7e-4 | — | 146.9 |
| L4 | 4f + all nine Cl 3p, no 5d, CAS(19,16) | 32 | −0.267 | Cl 3p 1.99991–1.99999; 4f virtuals ≤ 1e-4 | — | 92.6 |
| L5 | 4f + 5d + all Cl 3p, CAS(19,21) | 42 | not solved (1.1e11 det.) | — | — | 182.4 |

Nestedness (S-metric projection of one active subspace onto another): C3 ⊂ L3 and L2 ⊂ L3 exactly (1.0000);
L1's 4f orbitals lie 99.5–99.99 % inside C3, L3 and L4. C3 is *not* inside L4 (three of its σ donors have only
56–69 % of their weight in L4's "Cl 3p" projection space, because the 5d-projected donors carry Cl 3s/4p
components), and L3 is not inside L5 (94.9 %), for the same reason; energy differences are therefore only
compared between nested pairs. L4's CASCI (9.2e7 determinants) ran in 12 minutes on the workstation.

What the ladder supports, within this protocol:

1. **Without 5d the ligand electrons barely interact with the 4f shell in the CI.** With five σ donors and
   six empty 4f orbitals active (C3) the occupations are 2/1/0 to 1e-5 and the correlation energy is 0.09 mHa;
   with all nine Cl 3p orbitals (L4) the 4f virtual occupations stay ≤ 1e-4 and the correlation energy is
   0.27 mHa.
2. **In the nested pair C3 → L3, the five 5d orbitals provide the dominant additional correlation:** 6.44 of
   the 6.54 mHa, with a total 5d occupation of 0.0105 e versus ≤ 0.0003 e per 4f virtual. 5d is therefore
   necessary in the active-space description of CeCl3 in this basis; its magnitude (0.01 e, 6 mHa) marks it as
   a donation/polarisation channel rather than a strongly correlated degree of freedom. This is a statement
   about the fixed-orbital CASCI description and the spaces tested, not a general covalency mechanism.
3. **The 4f ligand-field splitting is small and reference-sensitive:** 0 / 0.68 (×2) / 1.46 (×2) / 5.72 /
   5.96 mHa (150–1300 cm⁻¹, the 1-2-2-1-1 pattern of f orbitals in D3h) with state-averaged orbitals, versus
   0 / 7.4 / 7.4 / 10.7 / 10.7 / 20.3 / 21.2 mHa with single-determinant 4f orbitals (core polarisation).
4. **λ tracks which orbitals are present more than the qubit count:** L4 (32 qubits) has λ = 93 Ha, below L2
   (24 qubits, 116 Ha), because 5d brings large exchange integrals.

## 4. Scalar versus two-component: spin-orbit coupling, one-electron and screened

Model: scalar active space plus a spin-orbit one-body block in the same orbitals, (a) one-electron X2C
(`qcchem/soc.py`), (b) plus the Breit-Pauli spin-orbit mean field (SOMF, Neese 2005) of the state density
built with the integral drivers of `pyscf.prop.gtensor` (`qcchem/somf.py`). Validation: the X2C block
reproduces the analytic Dirac 2p fine structure of hydrogenic Ne9+/Ar17+ to < 1 %; the SOMF spin convention
was fixed by the same test (a factor-2 error showed up as a ratio of exactly 0.50); the Breit-Pauli
one-electron block agrees with the X2C block to 0.5 % on Ce 4f. We describe this as a strong cross-validation
of the implementation, not as a proof of its correctness. Full N-sector ED on the SA 14-qubit spaces:

| System | one-electron SOC levels (cm⁻¹) | + SOMF screening (cm⁻¹) | experiment (cm⁻¹) | screening factor |
|---|---|---|---|---|
| Ce3+ 4f1 | 2F5/2 (6) 0; 2F7/2 (8) 5266 | 0; **2397** | 2253 | 0.455 |
| Pr3+ 4f2 | 3H4 (9) 0; 3H5 (11) 5169; 3F2 (5) 7582; 3H6 (13) 10163; 3F3 (7) 10873 | 0; **2259**; 4643 (3H6) | 3H5 2152; 3H6 4389 | 0.466 |
| CeCl3 L1-SA | Kramers doublets 0 / 189 / 1033 (J = 5/2 parentage); 5116 / 5250 / 5536 / 6331 (J = 7/2) | 0 / 191 / 996; **2320 / 2454 / 2778 / 3535** | — | 0.455 |

**λ and term counts, defined.** λ is the sum of |coefficient| over all Pauli strings of the qubit Hamiltonian
except the identity, *including* the one-body spin-orbit terms. Ce3+: scalar 55.21 → one-electron SOC 55.49
(+0.50 %) → SOMF 55.34 (+0.23 %); Pr3+ 50.25 → 50.53 → 50.38; CeCl3 L1-SA 59.00 → 59.23 → 59.11. Pauli
strings: 3297 → 3549 (Ce3+, Pr3+), 3369 → 3617 (CeCl3), unchanged between one-electron and screened SOC
because the SOMF term has the same sparsity pattern (CeCl3: 3609); 168 of the added strings do not conserve
S_z and one Z2 symmetry (S_z parity) is lost. These percentages refer to the SA orbital basis and are
representation-dependent.

Findings. (i) The one-electron operator overestimates the 2F7/2–2F5/2 interval 2.3×; the two-electron mean
field screens the 4f spin-orbit coupling by 0.46 and brings the intervals to within +5–6 % of experiment,
consistent with a mean-field treatment without dynamic correlation in a double-zeta basis. (ii) SOC
changes the size of the Hamiltonian little (+7.5 % strings, +0.2–0.5 % in λ, −1 symmetry, in this basis) but
changes the physics a lot: the screened SOC stabilisation of the ground state relative to the scalar
barycentre is 6.2 (Ce3+), 12.1 (Pr3+) and 5.7 mHa (CeCl3), the same size as the entire scalar ligand-field
spread (≤ 6 mHa) and as the 5d-mediated correlation gain of the ladder (6.5 mHa).

## 5. Fault-tolerant resource estimates with factorised Hamiltonians

OpenFermion `resource_estimates` costings (sparse: Lee et al. 2021; single factorisation: Berry et al.
2019; double factorisation: von Burg et al. 2021), ε = 1 mHa, χ = 10 bits, β = 20 bits; Toffoli counts and
logical qubits exclude surface-code overhead. Error metric: the exact CASCI energy change when the
factorised integrals replace the exact ones (where the CI space allows). `scripts/resource_estimates_ft.py`,
`results/ft_*.json`.

| System (spin-orbitals) | sparse: λ / Toffoli / qubits | SF exact rank: λ / Toffoli / qubits | DF threshold 1e-4: λ / L, Ξ / Toffoli / qubits / ΔE_CASCI |
|---|---|---|---|
| N2 CAS(6,6) (12) | 14.9 / 3.2e7 / 1287 | 10.8 (rank 21) / 1.5e7 / 187 | 9.8 / 20, 92 / not costable at this size / — / 0.000 mHa |
| Ce3+ 4f (14) | 55.2 / 1.2e8 / 1293 | 51.5 (28) / 7.9e7 / 263 | 43.5 / 28, 183 / 1.4e8 / 292 / 0 (N = 1) |
| CeCl3 L1 4f (14) | 56.1 / 1.2e8 / 1293 | 48.6 (28) / 7.5e7 / 263 | 47.0 / 28, 165 / 1.5e8 / 292 / 0 (N = 1) |
| CeCl3 C3 4f+σ (24) | 65.6 / 1.9e8 / 1565 | 55.1 (78) / 2.0e8 / 461 | 48.9 / 75, 693 / 3.2e8 / 651 / 0.007 mHa |
| CeCl3 L3 4f+5d+σ (34) | 146.9 / 5.9e8 / 1839 | 134.7 (153) / 9.0e8 / 861 | 100.3 / 123, 1597 / 1.0e9 / 867 / −0.031 mHa |
| CeCl3 L5 full (42) | 182.4 / 9.1e8 / 1849 | 169.1 (231) / 1.6e9 / 1574 | 104.4 / 140, 2134 / 1.3e9 / 1038 / not solvable |

Threshold dependence (L3): DF at 1e-3 gives λ = 100.2, 8.9e8 Toffolis but a 1.86 mHa CASCI error; 1e-4 is
the first threshold within chemical accuracy; 1e-5 costs 1.08e9 and errors 0.001 mHa. SF needs the full rank
(153) for < 1 mHa; rank 68 errors 2.0 mHa. Sparse truncation at 1e-4 keeps the error at 0.006 mHa with 4047
unique integrals. For C3, SF rank 48 (of 78) already errors 1.3 mHa and DF 1e-3 0.8 mHa.

Reading. (i) Double factorisation reduces λ by 30–45 % relative to the sparse 1-norm for the ligand-containing
spaces (L3: 147 → 100; L5: 182 → 104), as expected, but in this size range the sparse costing is still
competitive or cheaper in Toffolis (its advantage is the small d: 4e3–1e4 unique integrals), at the price of
~1800 logical qubits of QROM workspace versus ~900–1000 for DF. (ii) The chemically motivated CeCl3 space (L3)
costs ~1e9 Toffolis and ~870 logical qubits for 1 mHa with DF — three orders of magnitude below the earlier
Pauli-1-norm proxy (which was a per-query SELECT-by-unary-iteration estimate and is now withdrawn). (iii) The
full 42-qubit space (L5) costs 1.3e9 Toffolis and ~1000 logical qubits at DF 1e-4 — the first point on the
ladder with no classical CASCI reference, and a size at which DMRG is still routine. THC was not evaluated
(requires jax, not installed).

## 6. λ-minimising orbital rotations within the fixed active subspaces

`scripts/lambda_orbital_opt.py` minimises the exact JW 1-norm over orthogonal rotations of the active
orbitals (smoothed |·|, L-BFGS, three starts); the active subspace, CI space and CASCI energy are invariant,
checked by exact CASCI before and after wherever the CI is solvable.

| System | λ before → after | ratio | CASCI invariance | JW terms before → after |
|---|---|---|---|---|
| N2 CAS(6,6) | 14.86 → 13.97 | 0.940 | 0.0 | 382 → 246 |
| Ce3+ 4f (SA) | 55.21 → 52.87 | 0.958 | 0.0 | 3297 → 937 |
| CeCl3 L1 4f | 56.14 → 55.95 | 0.997 | 0.0 | 1643 → 1857 |
| CeCl3 C3 4f+σ | 65.56 → 62.32 | 0.951 | 0.0 | 13776 → 15830 |
| CeCl3 L3 4f+5d+σ | 146.94 → 125.22 | 0.852 | −1.8e-12 Ha | — |
| CeCl3 L5 full | 182.39 → 135.56 | 0.743 | not checkable (unitarity) | — |

The gain grows with the space: nothing for a pure 4f shell (the AVAS orbitals are already near-optimal),
15 % for the 34-qubit and 26 % for the 42-qubit CeCl3 spaces, where the σ/5d/4f mixing leaves room for
rotation. Term counts can go up while λ goes down (C3), so the two metrics must be reported separately.
Combined with double factorisation (which itself profits from more localised orbitals) this is the first
lever to pull before any hardware-side optimisation.

## 7. ADAPT-VQE versus UCCSD: compactness and the multireference regime

`qcchem/adapt.py`: fermionic-pool ADAPT (the same spin-conserving single/double generators that UCCSD uses,
operators may be selected repeatedly), exact statevector, gradient screening by <[H, A]>, full BFGS
re-optimisation with analytic adjoint gradients after every addition; stop when all pool gradients are below
the tolerance or |E - E_exact| < 1e-6 Ha. Circuit cost in the same two gate sets as before: PennyLane's
fermionic (JW-string) excitation gates ("CNOT") and Givens-type particle-conserving gates (2 CNOT per single,
14 per double, "Givens"). `results/adapt_*.json`.

| System | exact-reference error: UCCSD | ADAPT (tolerance 1e-3) | ADAPT (tolerance 1e-4, <= 150 ops) | parameters UCCSD / ADAPT-tight (distinct) | CNOT UCCSD / ADAPT-tight | Givens UCCSD / ADAPT-tight |
|---|---|---|---|---|---|---|
| N2 1.2 A (minimum) | 0.74 mHa | 0.69 mHa (47 ops) | - | 117 / 47 | 9024 / 4112 | 1422 / 658 |
| N2 1.5 A | 1.75 mHa | 1.63 mHa (45) | **4.6e-6 Ha** (93 ops, 51 distinct) | 117 / 93 (51) | 9024 / 8464 | 1422 / 1302 |
| N2 2.0 A | 5.46 mHa | 5.15 mHa (60, cap) | **2.2e-6 Ha** (97, 51 distinct) | 117 / 97 (51) | 9024 / 8816 | 1422 / 1358 |
| N2 2.4 A | 13.4 mHa | 1.18 mHa (58) | **9.2e-7 Ha** (86, 47 distinct) | 117 / 86 (47) | 9024 / 7776 | 1422 / 1204 |
| H2O 0.96 A | 0.098 mHa | 0.098 mHa (24) | - | 92 / 24 | 6976 / 1792 | 1096 / 336 |
| H2O 2.2 A | 3.13 mHa | 3.10 mHa (26) | 3.10 mHa (28, all gradients < 1e-4) | 92 / 28 | 6976 / 1744 | 1096 / 344 |
| LiH (4,4) | 8e-9 | 3.4e-7 (8 ops) | - | 26 / 8 | 1312 / 336 | 268 / 88 |
| Pr3+ 4f2 (SA) | 1e-8 | 7.6e-6 (6 ops) | - | 20 / 6 | 1400 / 488 | 160 / 60 |
| Ce3+, CeCl3 4f1 (SA) | 1e-8 | reference determinant already exact (0 ops) | - | 6 / 0 | 168 / 0 | 12 / 0 |

Reading. (i) With a loose gradient criterion ADAPT lands at the *same* accuracy as UCCSD at every point
(0.69 vs 0.74, 1.63 vs 1.75, 3.10 vs 3.13 mHa ...) with 2-4x fewer parameters and 2-4x fewer CNOTs: at that
level the two ansatze share the same limitation, a single application of each excitation generator. (ii) With
a tight criterion ADAPT keeps going by *re-selecting* operators (51 distinct generators among 93-97) and
reaches 1e-6 Ha at 1.5, 2.0 and 2.4 A - where first-order UCCSD is stuck at 1.7, 5.5 and 13.4 mHa - with
still fewer parameters (86-97 vs 117) and marginally fewer CNOTs (7.8-8.8k vs 9.0k). The accuracy gain in
the multireference regime therefore comes from operator ordering and repetition, not from a smaller circuit.
(iii) H2O at 2.2 A is a genuine failure of this pool: all pool gradients fall below 1e-4 while the energy is
3.1 mHa above the exact (N, S_z) ground state, which is reachable (36 % overlap with the HF determinant,
excited states 6.6 mHa higher); UCCSD stalls at the same 3.13 mHa. A spin-adapted or qubit pool, or
restarts, would be the next thing to try. (iv) For the 4f spaces ADAPT confirms the earlier picture: the
problem is the reference, not the ansatz - once the SA reference is right, 0-6 operators suffice.

## 8. What the numbers say about the original claim

- *4f-only is classically trivial*: confirmed — the hard part is the orbital reference, not the many-body
  problem; the Fock-space Hamiltonian's apparent density is an orbital artefact.
- *Ligand orbitals make it hard*: the chemically motivated space is 4f + 5d + σ (34 qubits), still classically
  solvable in minutes, and within the fixed-orbital protocol the ladder shows why 5d has to be there: without
  it the ligand electrons hardly interact with the 4f shell in the CI, with it a 0.01-electron donation
  channel carries 6.4 mHa. The classical wall is reached only at L5 (42 qubits, 1e11 determinants) or with
  larger ligands and basis sets.
- *Spin-orbit makes it hard*: not for the Hamiltonian (+7.5 % strings, +0.2–0.5 % λ, −1 symmetry, basis-
  dependent); yes for the physics (one-electron SOC 2.3× off; mean-field screening brings Ln multiplet
  intervals to ~5 %).
- *NISQ VQE is out of reach there*: refined rather than confirmed. ADAPT-VQE reaches 1e-6 Ha on the stretched
  N2 points where UCCSD fails by 2-13 mHa, so the multireference regime is not an ansatz wall; but it does so
  with 86-97 operators and ~8000 CNOTs (fermionic compilation) on 12 qubits, and it stalls on one of the two
  H2O points. The cost, not the accuracy, is the obstacle, and it was measured without noise.
- *Only fault-tolerant QPE remains*: for L5 the factorised estimate is ~1.3e9 Toffolis and ~1000 logical
  qubits (DF, 1e-4, before orbital optimisation), not the 6e10 of the earlier proxy; whether that beats DMRG at
  21 orbitals is still the open question, but the number is now one that can be compared.

The sentence the data support best: **for lanthanides the difficult part comes before the quantum solver —
a physically correct orbital reference (right shell occupation, spherical or ligand-field-adapted 4f
orbitals, 5d in the active space) and a screened relativistic Hamiltonian.** Once those are in hand, the
Hamiltonians of chemically meaningful Ln(III) halide spaces are 24–42 qubits, with factorised λ of 50–105 Ha
and 1e8–1e9 Toffolis per 1-mHa phase estimation.

## 9. Three research projects, ranked

1. **Two-component lanthanide Hamiltonians with SOMF screening, a spinor (Kramers) basis and symmetry-adapted
   tapering** (feasibility high, novelty high, cost low–medium).
2. **Pool design for ADAPT (spin-adapted, qubit, λ-optimised orbitals) on the ligand-containing Ln spaces**
   (feasibility high, novelty medium, cost medium): the 4f-only spaces are trivial for ADAPT; the 24-34-qubit
   4f+σ(+5d) spaces are the first non-trivial targets, and the H2O 2.2 A stall shows the spin-orbital pool is
   not enough.
3. **Classical-vs-fault-tolerant crossover on the CeCl3 ladder** (feasibility medium, novelty high, cost
   medium–high): DMRG references for L5 and larger spaces, THC λ, orbital-optimised DF, and a crossover plot of
   DMRG cost vs Toffoli count.

## 10. Limitations

cc-pVDZ-DK basis; scalar X2C; one-electron X2C SOC plus Breit-Pauli SOMF screening (no spin-dependent
two-electron terms beyond the mean field; core frozen at the scalar level); AVAS selection thresholds as
stated; CASCI at every CeCl3 level (CASSCF does not preserve the chemically defined spaces); only nested
pairs are compared energetically; no dynamic correlation outside the active space; no noise; ADAPT run only up to 14 qubits on an exact statevector; resource
estimates are logical-level costings without error-correction overhead and without THC; L5 has no classical
reference. Term counts, λ and tapering yields depend on the orbital basis for degenerate shells. The
four-component Dirac-Coulomb(-Gaunt) HF screening benchmark was still running when this revision was written.
