# Quantum algorithms for lanthanide electronic structure: a feasibility study with our own numbers

2026-10-24, revised after an independent referee-style review (see `progress.md`). Companion to
`technical_report.md`. Every number comes from `results/*.json` regenerated from one commit (67e551f;
the single cluster run records no SHA because the mirror is not a git checkout). Scalar-relativistic
sfx2c1e, cc-pVDZ-DK, PySCF 2.14 / OpenFermion 1.8.1 / PennyLane 0.45.1.

## 1. The claim being tested

`PLAN.md` set out the claim that *4f-only lanthanide active spaces are classically trivial; the problem only
becomes hard once ligand orbitals and relativistic/spin-orbit effects enter, and at that point NISQ-style
VQE is already out of reach, leaving only fault-tolerant phase estimation to discuss.* The calculations were
set up so that the numbers, not the claim, decide. Section 6 states what survived and what did not.

## 2. Free ions: the reference problem comes before the quantum problem

Ce3+ 4f1, Pr3+ 4f2, Eu3+ 4f6, Gd3+ 4f7 in the 4f-only CAS(n,7) space, 14 qubits.

**What went wrong first.** Hartree-Fock is an unreliable starting point for these ions: ROHF put two of
Gd3+'s electrons into 5d (CASSCF in the proper 4f7 space is 415 mHa lower), converged Pr3+ to a wrong
configuration (78 mHa above CASSCF), and for CeCl3 put the single open-shell electron into Ce 5d rather than
4f (a CASSCF 4f1 state is 118 mHa lower). A single-determinant reference also breaks spherical symmetry:
with ROHF/CASSCF orbitals the seven 4f levels of Ce3+ were split by 5–14 mHa, Eu3+'s by 10–18 mHa, and two
CASSCF runs for Pr3+ landed on solutions 1.5 mHa apart. A fractional-occupation reference did not help
(CASSCF(1,7) leaves the six unoccupied active orbitals undetermined). **What worked:** state-averaged (SA)
CASSCF over the complete high-spin multiplet manifold (7 states for Ce3+ and Eu3+, 21 for Pr3+). Restoring
spherical symmetry has a variational price — the SA CASCI lies 3.5 mHa (Ce3+) and 4.5 mHa (Eu3+) above the
symmetry-broken ROHF energy — but the sparse ED of the resulting qubit Hamiltonian then reproduces the exact
term structure: Ce3+ 2F seven-fold degenerate to 1e-7 mHa; Pr3+ 3H (11) at 0, 3F (7) at 27.19 mHa, 3P (3) at
126.88 mHa (SA-CASSCF state energies and qubit ED agree state by state); Eu3+ 7F (7) degenerate. These
degeneracies validate the folding/mapping chain more sharply than any single energy.

**How hard is 4f-only?** Not hard. The fixed-S_z CI spaces are 14 (Ce3+), 21 (Pr3+), 7 (Eu3+) and 1 (Gd3+)
determinants. UCCSD-VQE is exact where it applies: Ce3+ 7.8e-9 Ha with 6 parameters and 168 CNOTs (6 s);
Pr3+ −2.3e-8 Ha (below the CASCI reference, whose own convergence is 8e-8 on this degenerate state; the ED
reference is the better comparison) with 20 parameters and 1400 CNOTs (108 s). Note that for Ce3+ every
two-electron term vanishes in the N = 1 sector: the physical problem is a 7×7 matrix, whatever the
Fock-space Hamiltonian looks like.

**The Hamiltonian looks dense, but that is a basis statement.** In SA natural orbitals (degenerate, hence
arbitrarily rotated within the 4f shell) the Ce3+ Hamiltonian has 3297 JW terms, λ = 54.6 Ha and two Z2
symmetries (14 → 12 qubits). In the ROHF-based natural orbitals of the same ion it has 845 terms, λ = 16.4 Ha
and four Z2 symmetries (14 → 10 qubits), with the same spectrum. Term counts, λ and the tapering yield are
therefore orbital-choice artefacts for a degenerate shell, and λ-minimising orbital rotations (Koridon et
al. 2021) are a first-order lever for Ln Hamiltonians, before any algorithmic choice.

## 3. CeCl3: an active-space ladder chosen by orbital composition

Planar D3h CeCl3, Ce–Cl 2.57 Å, 138 basis functions. Reference: ROHF (open-shell electron moved into 4f by a
preliminary CASSCF(1,8) in the {4f + SOMO} space, E = −10239.41113 Ha), then AVAS with ANO-RCC atomic
projectors; every level's orbitals were inspected (Mulliken weights are in the JSON records).

| Level | Active space | Orbitals (composition) | qubits | JW terms | λ / Ha | CASCI E / Ha | CI dim. | Solved |
|---|---|---|---|---|---|---|---|---|
| L1 | 4f only, CAS(1,7), CASCI | 7 × Ce 4f (0.98–1.00) | 14 (12 tapered) | 1549 | 56.0 | −10239.40760 | 7 | local, ED+VQE |
| L1-SA | same, SA-CASSCF over 7 states | 7 × 4f (0.99) | 14 (12) | 3365 | 59.6 | −10239.40844 | 7 | local, ED+VQE |
| L2 | 4f + 5d, CAS(1,12), CASCI | 7 × 4f, 5 × 5d (0.81–0.93) | 24 (22) | 14252 | 118.5 | −10239.40760 | 12 | local, JW statistics |
| L3 | 4f + 5d + Cl σ donors, CAS(11,17) | + 5 occupied Cl 3p/4p with 5–15 % Ce 5d | 34 | — | 148.0 | −10239.41414 | 7.7e7 | local CASCI, 10 min |
| L4 | 4f + all Cl 3p, CAS(19,16) | 7 × 4f + 9 Cl 3p (4–10 % Ce 5d) | 32 | — | 92.1 | −10239.41111 | 9.2e7 | cluster, 8 cores, 6 min |
| L5 | 4f + 5d + all Cl 3p, CAS(19,21) | all of the above | 42 | — | 182.7 | not solved (1.1e11) | 1.1e11 | integrals + λ only |

Above 24 qubits λ comes from the integrals (Koridon identity; agrees with the explicit expansion to 1e-9
relative for molecules and free ions and 2e-7 relative for CeCl3 L1/L2 — the residual is not understood and
is immaterial here).

What the ladder shows, and what it cannot show.

1. **Energies across levels are not variationally consistent, so they do not rank "covalency channels".**
   The L3 and L4 CASCI energies (−10239.41414 and −10239.41111 Ha) differ by 3.0 mHa, but the one-electron
   CASSCF(1,8) used to relax the reference is already at −10239.41113 Ha, i.e. 3.5 mHa *below* L1 and equal
   to L4 within 0.02 mHa: the AVAS re-projection of the 4f orbitals after relaxation gives up most of the
   orbital relaxation, and each CASCI level inherits that loss differently. L1-SA (−10239.40844) is 0.8 mHa
   below L1 for the same reason. The honest statement is that the orbital-reference sensitivity (3–4 mHa,
   and 7–21 mHa for the 4f level pattern, see point 2) is as large as the differences between levels. The
   original inference that "covalency goes through 5d rather than 4f" is therefore **withdrawn**; testing it
   needs orbital-consistent (CASSCF) energies at each level and a 4f + σ-without-5d control, which were not
   run. What *can* be said is structural: AVAS with a 10 % Ce-5d projection threshold selects exactly five
   occupied Cl orbitals (the σ-donor combinations with 5–15 % 5d weight), while the same threshold on 4f
   selects none — the ligand orbitals that mix with the metal mix with 5d.
2. **The 4f ligand-field splitting is small and reference-sensitive.** With state-averaged orbitals the 4f
   manifold of CeCl3 splits as 0 / 0.68 (×2) / 1.46 (×2) / 5.72 / 5.96 mHa (150–1300 cm⁻¹), the 1-2-2-1-1
   pattern of f orbitals in D3h. With the single-determinant-relaxed orbitals the same CASCI gives
   0 / 7.4 (×2) / 10.7 (×2) / 20.3 / 21.2 mHa: residual core polarisation exceeds the ligand field five-fold.
   Any treatment of Ln multiplet structure, quantum or classical, must solve this reference problem first.
3. **λ and term counts along the ladder.** λ per qubit is roughly constant (4.0 at L1, 4.3 at L3 and L5) and
   L4 (32 qubits, λ = 92 Ha) sits below L2 (24 qubits, 118 Ha): λ tracks which orbitals are present (5d
   contributes large exchange integrals) more than the qubit count. The JW term count over all systems
   follows N^4.35 in the number of active orbitals (fig. 3, fit over all systems), i.e. of order 2–3e5 explicit terms at 42 qubits.

## 4. Scalar versus two-component: what spin-orbit coupling does to the qubit Hamiltonian

Model: the scalar active space plus the one-electron X2C spin-orbit block in the same orbitals
(`qcchem/soc.py`). The transformation was validated on hydrogenic Ne9+ and Ar17+ (`tests/unit/
test_soc_hydrogenic.py`): in a large even-tempered basis the 2p1/2–2p3/2 splitting from the projected
one-body block agrees with the analytic Dirac value to < 1 %, so the remaining disagreement with experiment
for Ce3+ below is physics (missing two-electron screening), not a transformation error. On the SA 14-qubit
spaces, full N-sector ED (all 14 one-electron states; 45 of 91 two-electron states):

| System | JW terms scalar → 2c | λ scalar → 2c | Z2 gens | two-component levels (mHa) and multiplicities | SOC ground-state shift |
|---|---|---|---|---|---|
| Ce3+ 4f1 | 3297 → 3549 (×1.076) | 54.64 → 54.88 (×1.004) | 2 → 1 | 2F5/2 (6) at 0; 2F7/2 (8) at 23.99 (5266 cm⁻¹; expt 2253) | −13.7 mHa |
| Pr3+ 4f2 | 3333 → 3585 (×1.076) | 51.37 → 51.72 (×1.007) | 2 → 1 | 3H4 (9) at 0; 3H5 (11) at 23.55; 3F2 (5) at 34.55; 3H6 (13) at 46.31; 3F3 (7) at 49.54 | −27.5 mHa |
| CeCl3 L1-SA | 3365 → 3613 (×1.074) | 59.57 → 59.83 (×1.004) | 2 → 1 | three Kramers doublets at 0 / 0.86 / 4.71 (J = 5/2 parentage), four at 23.3 / 23.9 / 25.2 / 28.8 (J = 7/2) | −13.0 mHa |

Findings. (i) The two-component Hamiltonian reproduces the exact J-multiplicities (6/8; 9/11/5/13/7; 3+4
Kramers doublets) and the ordering of the Pr3+ levels including the 3F2 intrusion below 3H6, but the
one-electron spin-orbit operator overestimates the 2F7/2–2F5/2 interval 2.3× (effective ζ4f ≈ 1500 vs ≈ 640
cm⁻¹). The missing two-electron screening is not a small correction for 4f; a spin-orbit mean-field (AMFI/
SOMF) or a two-component CASSCF reference is required before quantitative Ln spin-orbit work on qubits.
(ii) Structurally the cost of SOC is modest in this basis: ~7.5 % more Pauli strings (the 182 nonzero
one-body SOC matrix elements add ~250 strings, 168 of them S_z-non-conserving; all coefficients remain real
because the purely imaginary one-body SOC elements map onto XY−YX-type strings), < 1 % in λ, and the loss of
one Z2 symmetry (S_z parity; particle-number parity survives), i.e. one fewer tapered qubit. (iii) The SOC
stabilisation of the ground state relative to the scalar barycentre (13–27 mHa) is larger than the entire
ligand-field spread of the scalar 4f manifold (≤ 6 mHa) — the familiar SOC ≫ crystal-field hierarchy for
lanthanides, reproduced here on the qubit Hamiltonian — so a scalar-only treatment of a Ln(III) complex
targets the wrong ground state.

## 5. Resource proxies, and what they do not mean

For each Hamiltonian the records contain λ, the qubitised-walk query count π λ / (2ε) for ε = 1 mHa, a
sparse-LCU Toffoli proxy (queries × number of Pauli terms, SELECT by unary iteration) and a logical-qubit
proxy (system + LCU + phase registers). Along the ladder: L1 8.8e4 queries, 1.4e8 Toffolis, 42 logical
qubits; L3 2.3e5 queries, ~2e10 Toffolis (term count from the nonzero-integral proxy), 69 qubits; L5 2.9e5
queries, ~6e10 Toffolis, 79 qubits. Three caveats a reader should not skip: these are Fock-space 1-norms
(for N = 1 spaces most terms vanish in the physical sector); they are sparse-LCU numbers, whereas double
factorisation or THC would reduce both λ and the per-query cost substantially (not evaluated here; OpenFermion's
`resource_estimates` module is the obvious next step); and no error-correction overhead is included. They
compare Hamiltonians to each other; they are not hardware estimates.

## 6. What the numbers say about the original claim

- *4f-only is classically trivial*: confirmed, and sharper than expected — the hard part is the orbital
  reference, not the many-body problem, and the Fock-space Hamiltonian's apparent density is an orbital
  artefact. Quantum solvers add nothing on 4f-only spaces.
- *Ligand orbitals make it hard*: the chemically motivated space (4f + 5d + σ donors) is 34 qubits and still
  classically solvable in minutes; the classical wall is reached only at L5 (42 qubits, 1e11 determinants)
  or with larger ligands and basis sets. Which ligand channel matters energetically could not be decided
  with the present (non-variationally-consistent) ladder.
- *Spin-orbit makes it hard*: not for the Hamiltonian (+7.5 % terms, +0.5 % λ, −1 symmetry); yes for the
  physics (SOC shift ≫ ligand field; one-electron SOC 2.3× off). The difficulty is in the reference.
- *NISQ VQE is out of reach there*: supported by the molecular numbers (9024 CNOTs for 12 qubits in the
  naive fermionic-excitation compilation used here; UCCSD errors of 2–16 mHa in the multireference regime)
  but not proven: compact ansätze were not run and the CNOT count is unoptimised.
- *Only fault-tolerant QPE remains*: the honest version is "fault-tolerant QPE is the only route whose cost
  can even be written down for L5, and in the sparse-LCU proxy it is ~6e10 Toffolis before factorisation";
  whether it beats DMRG at 21–30 orbitals is the open question.

## 7. Three research projects, ranked

1. **Two-component lanthanide Hamiltonians with a proper spin-orbit mean field and Kramers-adapted tapering**
   (feasibility high, novelty high, cost low–medium). Replace the one-electron SOC block by AMFI/SOMF or a
   two-component CASSCF spinor basis; map complex spinor Hamiltonians to qubits; derive which Z2 and
   time-reversal symmetries survive and how to target a J-manifold; validate against the exact 2F/3H
   intervals. Output: a general Ln(III) 4f-shell Hamiltonian generator with resource statistics.
2. **ADAPT-VQE and pool design for 4f exchange-dominated Hamiltonians, with orbital optimisation for λ**
   (feasibility high, novelty medium, cost medium). Fermionic vs qubit pools on the SA Ce3+/Pr3+/Eu3+ and
   CeCl3 L1 Hamiltonians against ED; include λ-minimising orbital rotations, since §2 shows a 3× spread.
3. **Classical-vs-fault-tolerant crossover on the CeCl3 ladder with orbital-consistent references**
   (feasibility medium, novelty high, cost medium–high). CASSCF/DMRG (block2) references at every level and
   beyond (bigger basis, EuCl3/GdCl3), double-factorised and THC λ for the same integrals, and a crossover
   plot of DMRG cost vs QPE Toffoli count as a function of active orbitals.

## 8. Limitations

cc-pVDZ-DK basis; scalar X2C; one-electron SOC without screening; AVAS selection thresholds as stated in
the configs; CASCI (not CASSCF) at L1–L5, hence no variational consistency across levels; no dynamic
correlation outside the active space; no noise; QPE proxies are 1-norm/sparse-LCU based; L5 has no classical
reference. Term counts, λ and tapering yields depend on the orbital basis for degenerate shells (§2).
