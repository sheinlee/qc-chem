# Quantum algorithms for lanthanide electronic structure: a feasibility study with our own numbers

2026-10-24. Companion to `technical_report.md`. Every number below comes from `results/*.json` produced by
this repository (scalar-relativistic sfx2c1e, cc-pVDZ-DK, PySCF 2.14 / OpenFermion 1.8.1 / PennyLane 0.45.1).

## 1. The claim being tested

The claim set out in `PLAN.md` was: *4f-only lanthanide active spaces are classically trivial and do not
constitute a quantum-advantage scenario; the problem only becomes hard once ligand orbitals and
relativistic/spin-orbit effects enter, and at that point NISQ-style VQE is already out of reach, leaving
only fault-tolerant phase estimation to discuss.* The calculations were set up so that the numbers, not the
claim, decide. Section 6 states what survived.

## 2. Free ions: the reference problem comes before the quantum problem

Four trivalent ions (Ce3+ 4f1, Pr3+ 4f2, Eu3+ 4f6, Gd3+ 4f7) in the 4f-only CAS(n,7) space, 14 qubits.

**What went wrong first.** Hartree-Fock is unreliable for these ions: ROHF put two of Gd3+'s electrons into 5d
(CASSCF in the proper 4f7 space is 415 mHa lower), converged Pr3+ to a wrong configuration (78 mHa above
CASSCF), and for CeCl3 put the single open-shell electron into Ce 5d rather than 4f (the CASSCF 4f1 state is
118 mHa lower). A single-determinant reference also breaks spherical symmetry: with ROHF/CASSCF orbitals the
seven 4f levels of Ce3+ were split by 5–14 mHa and Eu3+'s by 10–18 mHa, and two CASSCF runs for Pr3+ landed
on solutions 1.5 mHa apart. **Fix:** state-averaged CASSCF over the complete high-spin multiplet manifold
(7 states for Ce3+/Eu3+, 21 for Pr3+). The sparse ED of the resulting qubit Hamiltonian then reproduces the
exact term structure: Ce3+ 2F, seven-fold degenerate to 1e-7 mHa; Pr3+ 3H (11) at 0, 3F (7) at 27.19 mHa,
3P (3) at 126.88 mHa; Eu3+ 7F (7) degenerate. These degeneracies are a sharper validation of the whole
folding/mapping chain than any single energy.

**How hard is 4f-only?** Not hard. The fixed-S_z CI spaces are 14 (Ce3+), 21 (Pr3+), 7 (Eu3+) and 1 (Gd3+)
determinants. UCCSD-VQE is exact where it is applicable: Ce3+ 7.6e-9 Ha with 6 parameters and 168 CNOTs
(3.8 s); Pr3+ 2e-8 Ha with 20 parameters and 1400 CNOTs (149 s, including the restart needed because the
state-averaged orbitals have no leading determinant). The Hamiltonians are nevertheless *dense*: 3285–3365
JW terms and λ = 50–54 Ha for Ce3+/Pr3+ (compare N2 CAS(6,6): 382 terms, λ = 13.8 Ha) because the seven 4f
orbitals are all strongly coupled by the same Coulomb/exchange integrals. λ depends strongly on the orbital
choice — Ce3+ gives 50.7 Ha (ROHF-based orbitals), 54.4 Ha (state-averaged), 16.6 Ha in one symmetry-broken
run — which is the orbital-optimisation lever that Koridon et al. exploit and is worth a project of its own.

## 3. CeCl3: an active-space ladder chosen by orbital composition

Planar D3h CeCl3, Ce–Cl 2.57 Å, 138 basis functions. Reference: ROHF (open-shell electron relaxed into 4f by
a preliminary CASSCF), then AVAS with ANO-RCC atomic projectors. Every level's active orbitals were inspected
(Mulliken weights are stored in the JSON records).

| Level | Active space | Orbitals (composition) | qubits | JW terms | λ / Ha | CASCI E / Ha | CI dimension | Where solved |
|---|---|---|---|---|---|---|---|---|
| L1 | 4f only, CAS(1,7) | 7 × Ce 4f (0.98–1.00) | 14 (12 tapered) | 1549 | 56.0 | −10239.40760 | 7 | local, ED+VQE |
| L1-SA | same, SA-CASSCF(7) | 7 × 4f (0.99) | 14 | 3377 | 59.1 | −10239.40844 | 7 | local, ED+VQE |
| L2 | 4f + 5d, CAS(1,12) | 7 × 4f, 5 × 5d (0.81–0.93) | 24 (22) | 14540 | 118.4 | −10239.40760 | 12 | local, explicit JW statistics |
| L3 | 4f + 5d + Cl σ donors, CAS(11,17) | + 5 occupied Cl 3p/4p orbitals with 5–15 % Ce 5d | 34 | — (λ from integrals) | 147.7 | −10239.41414 | 7.7e7 | local CASCI, 10 min |
| L4 | 4f + all Cl 3p, CAS(19,16) | 7 × 4f + 9 Cl 3p (4–10 % Ce 5d) | 32 | — | 92.1 | −10239.41111 | 9.2e7 | cluster, 8 cores, 6 min |
| L5 | 4f + 5d + all Cl 3p, CAS(19,21) | all of the above | 42 | — | 182.4 | not solved (~1e11) | 1.1e11 | integrals + λ only |

Three things the ladder shows.

1. **Covalency goes through 5d, not 4f.** Adding the five 5d orbitals plus the five σ-donor combinations that
   mix with them (L3) lowers the energy by 6.5 mHa relative to 4f-only, whereas adding all nine Cl 3p orbitals
   without 5d (L4, three more electrons pairs and the same qubit count) gains only 3.5 mHa — L4 is 3.0 mHa
   *above* L3. The chemically meaningful active space for a Ce(III) halide is therefore 4f + 5d + ligand σ, and
   it already costs 34 qubits.
2. **The ligand-field splitting is small and reference-sensitive.** With state-averaged orbitals the 4f
   manifold of CeCl3 splits as 0 / 0.68 (×2) / 1.46 (×2) / 5.72 / 5.96 mHa (150–1300 cm⁻¹), the 1-2-2-1-1
   pattern expected for f orbitals in D3h. With single-determinant-relaxed orbitals the same CASCI gives
   0 / 7.4 (×2) / 10.7 (×2) / 20.3 / 21.2 mHa: residual core polarisation is five times larger than the
   ligand field. Any quantum (or classical) treatment of Ln multiplet structure has to solve this orbital
   problem first; it is not a property of the algorithm.
3. **λ grows faster than the qubit count.** 14 → 34 → 42 qubits takes λ from 56 to 148 to 182 Ha, and the
   JW term count from 1.5e3 (14 qubits) to 1.5e4 (24 qubits); the fitted power law over all systems is
   N^4.2 in the number of active orbitals (fig. 3), i.e. ~2.5e5 terms at 42 qubits. For the 1-norm-based
   QPE proxy (π λ / 2ε, ε = 1 mHa) this is 8.8e4 walk queries at L1, 2.3e5 at L3 and 2.9e5 at L5; the
   sparse-LCU Toffoli proxy is 2.5e8 at L5 (1.6e8 at L3), and the logical-qubit proxy 82 (72). These proxies ignore double
   factorisation/THC, which typically reduce λ by 3–10×, and ignore error correction; they are only for
   comparing Hamiltonians.

## 4. Scalar versus two-component: what spin-orbit coupling does to the qubit Hamiltonian

Model: the scalar active space plus the one-electron X2C spin-orbit block in the same orbitals
(`qcchem/soc.py`; no two-electron/mean-field SOC screening, no spin-dependent two-electron terms). On the
state-averaged 14-qubit spaces:

| System | JW terms scalar → 2c | λ scalar → 2c | Z2 generators | ED levels, 2c (mHa) | ground-state SOC shift |
|---|---|---|---|---|---|
| Ce3+ 4f1 | 3285 → 3537 (×1.077) | 54.37 → 54.64 (×1.005) | 2 → 1 | 2F5/2 (6) at 0, 2F7/2 at 23.99 (5266 cm⁻¹) | −13.7 mHa |
| Pr3+ 4f2 | 3365 → 3617 (×1.075) | 50.56 → 50.89 (×1.007) | 2 → 1 | 3H4 (9) at 0, 3H5 (11) at 23.55, first 3H6 at 34.55 | −27.5 mHa |
| CeCl3 L1-SA | 3377 → 3621 (×1.072) | 59.14 → 59.39 (×1.004) | 2 → 1 | Kramers doublets at 0, 0.86, 4.71; 23.3 | −13.0 mHa |

Findings. (i) The two-component Hamiltonian reproduces the right multiplet multiplicities (6/8 for 2F, 9/11/13
for 3H_J) but overestimates the 2F7/2–2F5/2 interval 2.3× (5266 vs 2253 cm⁻¹ experimental): the missing
two-electron screening is not a small correction for 4f, and a spin-orbit mean-field (AMFI/SOMF) or a full
two-component CASSCF is required before any quantitative Ln spin-orbit work on qubits. (ii) Structurally the
cost of SOC is modest: ~7 % more Pauli terms (the 182 nonzero one-body SOC matrix elements add ~250 Pauli
strings, 168 of them S_z-non-conserving; all coefficients stay real because the purely imaginary spin-orbit
matrix elements map onto XY−YX-type strings), <1 % in λ, and the loss of one Z2 symmetry (S_z parity; particle-number parity survives), i.e. one
fewer tapered qubit. The hard part of relativistic Ln problems is therefore not the Hamiltonian size but the
reference (spinor orbitals, Kramers symmetry, state targeting within J-manifolds), exactly as in classical
two-component CASSCF. (iii) The ground-state SOC stabilisation (13–27 mHa) is an order of magnitude larger
than the ligand-field splitting (≤ 6 mHa) and comparable to the correlation gain from the whole 4f→5d/σ
ladder (6.5 mHa): a scalar-only Ln calculation misses the leading effect.

## 5. Three tiers

**Demonstrable today (this repository).** Everything up to 14 qubits end to end (H2, LiH, N2, H2O, Ln 4f-only,
CeCl3 4f with ligand field, two-component 4f), with VQE at 1e-8 Ha or better where the reference determinant
is right, plus exact Hamiltonian statistics to 24 qubits and exact λ to any size. Cost of the largest VQE run
here: 12 qubits, 117 parameters, 9024 CNOTs, 65 BFGS iterations — and it misses chemical accuracy beyond
1.4 Å on N2.

**Plausible on near-term hardware (not done here).** The 14-qubit Ln 4f Hamiltonians are small enough in qubits
but have 3300 Pauli terms and λ ≈ 50 Ha, i.e. an estimated ≥ 1e6–1e7 shots per energy for mHa precision
under naive grouping; the 1400-CNOT Pr3+ UCCSD circuit is well beyond current two-qubit-gate fidelities.
Compact ansätze (ADAPT, with pools adapted to the 4f exchange structure), measurement reduction (grouping,
classical shadows, or double-factorised measurement) and symmetry verification are the enabling steps, and
all of them can be benchmarked exactly against the ED references produced here.

**Fault-tolerant only.** The chemically meaningful CeCl3 space (L3, 34 qubits, λ = 148 Ha; L5, 42 qubits,
λ = 182 Ha) and any Ln complex with more than one metal or a realistic ligand set. L3 is solvable classically
in ten minutes; L5 is at the edge of DMRG (21 orbitals) and the first level where a quantum solver could in
principle add something — and there the question is purely one of fault-tolerant cost (1-norm 182 Ha before
factorisation).

## 6. What the numbers say about the original claim

- *4f-only is classically trivial*: confirmed, and sharper than expected — the hard part there is the orbital
  reference, not the many-body problem. Quantum solvers add nothing on 4f-only spaces.
- *Ligand orbitals make it hard*: partly. The right ligand channel is 5d-mediated σ donation; with it the
  space is 34 qubits (L3) and still classically solvable in minutes. The classical wall is reached only at L5
  (42 qubits, 1e11 determinants) or with larger ligands/basis sets.
- *Spin-orbit makes it hard*: not for the Hamiltonian (+7 % terms, +0.5 % λ, −1 symmetry); yes for the
  physics (SOC shift ≫ ligand-field splitting, and a one-electron treatment is 2.3× off). The difficulty is
  in building the right two-component reference, which is a classical pre-processing problem.
- *NISQ VQE is out of reach at that point*: supported by the N2 and CeCl3 numbers (9024 CNOTs for 12 qubits;
  ansatz error 2–16 mHa in the multireference regime), but not proven for compact ansätze, which were not run.
- *Only fault-tolerant QPE remains*: the honest version is "fault-tolerant QPE is the only route whose cost
  can even be written down for L5, and that cost (λ = 182 Ha before factorisation) is not small." Whether it
  beats DMRG at 21–30 orbitals is exactly the open question.

## 7. Three research projects, ranked

1. **Two-component lanthanide Hamiltonians with proper spin-orbit mean field and Kramers-adapted tapering**
   (feasibility high, novelty high, cost low–medium). Replace the one-electron SOC block by an AMFI/SOMF or
   two-component CASSCF spinor basis; map complex spinor Hamiltonians to qubits; derive which Z2 (and
   time-reversal) symmetries survive and how to target a J-manifold; validate against the exact 2F/3H
   intervals. Output: a general Ln(III) 4f-shell Hamiltonian generator plus resource statistics.
2. **ADAPT-VQE and pool design for 4f exchange-dominated Hamiltonians** (feasibility high, novelty medium,
   cost medium). Fermionic vs qubit pools on the state-averaged Ce3+/Pr3+/Eu3+ and CeCl3 L1 Hamiltonians,
   measured against ED; question: does the dense 4f exchange structure defeat operator sparsity, and does
   orbital choice (λ-minimising rotations) change the answer?
3. **Classical-vs-fault-tolerant crossover on the CeCl3 ladder** (feasibility medium, novelty high, cost
   medium–high). DMRG (block2) references for L5 and larger spaces (bigger basis, EuCl3/GdCl3), double-
   factorised and THC λ for the same integrals, and an honest crossover plot of DMRG cost vs QPE Toffoli
   count as a function of active orbitals. This is the figure that would answer the title question.

## 8. Limitations

cc-pVDZ-DK basis; scalar X2C; one-electron SOC without screening; AVAS orbital selection with threshold
choices stated in the configs; no dynamic correlation outside the active space; no noise; QPE proxies are
1-norm based. L5 has no classical reference. All free-ion and CeCl3 numbers depend on the orbital reference
at the mHa level, as documented in §2–3.
