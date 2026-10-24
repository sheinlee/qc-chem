# From Active-Space Electronic Structure to Quantum Algorithms: Benchmarks and Resource Estimates for Strongly Correlated and Lanthanide Systems

Technical report, 2026-10-24. All numbers in this document were produced by the code in this repository
(`results/*.json`, `results/summary_table.csv`, figures in `results/figures/`). Software: PySCF 2.14.0,
OpenFermion 1.8.1, PennyLane 0.45.1 (lightning.qubit), SciPy 1.17.1. Hardware: Ryzen 9 5900 workstation
(all VQE/ED/statistics) and one 8-core slurm job on a Xeon Gold 6326 node (the CAS(19,16) CASCI reference).

## 1. Question and scope

Can a variational quantum eigensolver reproduce CASCI/FCI-quality energies when given exactly the same
active-space Hamiltonian, and what does the Hamiltonian cost in qubits, Pauli terms, 1-norm, variational
parameters and circuit depth, from textbook molecules up to a lanthanide complex? The lanthanide part
asks where, along a chemically motivated active-space ladder for CeCl3, the problem leaves the regime
that can be simulated or solved classically today, and what a fault-tolerant phase-estimation
(qubitization-style) algorithm would need on a 1-norm proxy basis.

Everything is a statevector or exact calculation. No hardware noise is modelled; this is a study of the
algorithmic objects (Hamiltonians, ansätze, resource counts), not of NISQ hardware.

## 2. Pipeline and validation

**Active-space Hamiltonian.** PySCF RHF/ROHF (scalar-relativistic sfx2c1e for Ln), CASCI or CASSCF in
canonical or AVAS-selected orbitals; the folded Hamiltonian `H = ecore + Σ h_pq a†a + ½ Σ (pq|rs) a†a†aa`
is obtained from `CASCI.get_h1eff()` / `get_h2eff()`; the chemist two-electron tensor is converted to the
OpenFermion physicist convention with `0.5·eri.transpose(0,2,3,1)` and expanded to interleaved spin orbitals
(2p = α, 2p+1 = β). **Acceptance test:** sparse exact diagonalisation (ED) of the Jordan-Wigner (JW) qubit
Hamiltonian in the (N, S_z) sector reproduces the PySCF CASCI energy: H2 7e-16 Ha, LiH(4,4) 4e-15, N2 CAS(6,6)
7e-14, Ce3+ CAS(1,7) 1e-11, Pr3+ CAS(2,7) 6e-8 (degenerate ground state, eigsh tolerance) — see
`tests/unit/test_hamiltonian.py`.

**Mappings and symmetry.** JW and Bravyi-Kitaev give identical spectra and identical term counts for all
systems. Z2 tapering (PennyLane) removes 2–4 qubits: H2 4→1, LiH(4,4) 8→5, N2/H2O 12→8, Ln 4f spaces 14→12,
CeCl3 4f+5d 24→22; the tapered ground-state energy agrees with CASCI to 1e-13.

**1-norm identity.** The JW Pauli 1-norm λ = Σ|h_j| (identity excluded) is computed in two independent
ways: from the explicit Pauli expansion and directly from the spatial-orbital integrals with the Koridon et
al. (PRR 3, 033127, 2021) expression. For real orbitals the two agree to 1e-13 (molecules) and ~1e-7
relative (lanthanides, limited by coefficient compression at 1e-12). The integral formula is therefore an
*exact* λ for the 32- and 42-qubit CeCl3 Hamiltonians where no Pauli expansion is performed (fig. 3). It is
only a proxy for complex (two-component) integrals, and is labelled as such in the data.

**UCCSD-VQE.** Spin-conserving singles and doubles from an (n_α, n_β) reference determinant (open-shell
capable), disentangled first-order Trotter UCCSD (`qml.UCCSD`), `lightning.qubit` with adjoint gradients,
BFGS (gtol 1e-6) with an L-BFGS-B polish if BFGS stops on precision loss. Unit tests check θ = 0 → E_HF
and a stored reference angle set; integration tests check |E_VQE − E_FCI| < 1e-6 (H2) and 1e-5 (LiH).
Circuit depth and CNOT counts are reported after decomposition to {CNOT, single-qubit rotations}.
Where the first optimisation ends > 1e-5 Ha above the known CASCI energy, up to two restarts from random
angles (σ = 0.2) are made and the best result kept; the number of restarts is recorded. This is a
benchmark luxury (the exact answer is known) and is stated wherever it was used.

## 3. Molecules: where first-order UCCSD stops being chemically accurate

| System | qubits | JW terms | λ / Ha | UCCSD θ | CNOT | depth | UCCSD error at equilibrium |
|---|---|---|---|---|---|---|---|
| H2 / STO-3G (2,2) | 4 (1 tapered) | 14 | 1.89 | 3 | 64 | 100 | 1e-14 (exact) |
| LiH / STO-3G (4,4) | 8 (5) | 192 | 6.82 | 26 | 1312 | 1786 | 8e-9 |
| H2O / STO-3G (8,6) | 12 (8) | 550 | 27.7 | 92 | 6976 | 8878 | 1.0e-4 |
| N2 / STO-3G (6,6) | 12 (8) | 382 | 13.7 | 117 | 9024 | 11496 | 4.2e-4 |

The N2 dissociation curve (fig. 2) is the central molecular result. The optimiser converges at every
point (final gradient norm < 1e-6), so the error is the ansatz's, not the optimiser's: the UCCSD error
grows monotonically from 0.19 mHa (R = 1.0 Å) through the chemical-accuracy threshold at R ≈ 1.4 Å
(1.58 mHa at 1.35 Å, 2.16 mHa at 1.5 Å) to 4.6 mHa at 2.0 Å and 16.5 mHa at 3.0 Å, while the number of
energy evaluations needed rises from 13 to 65 and the wall time per point from 39 s to 377 s. The
Hamiltonian itself does not get harder along the curve (382 terms, λ between 12.6 and 15.9 Ha); what
changes is the multireference character that a single-reference, singly-Trotterised UCCSD cannot
represent. The H2O symmetric stretch shows the same qualitative behaviour with a later onset: errors stay
at 0.06–0.25 mHa up to R = 1.8 Å and reach 3.1 mHa at 2.2 Å (at 2.6 Å the ROHF reference changes
character, λ drops to 18 Ha, and the optimiser needs 260 evaluations to reach 0.7 mHa).

Observation on cost: for the same 12 qubits, the UCCSD circuit for N2 CAS(6,6) is 9024 CNOT deep. This is
the size of a first-order UCCSD circuit for a problem whose classical FCI space is 400 determinants.

## 4. Lanthanide free ions: the reference problem comes first

Ce3+, Pr3+, Eu3+, Gd3+ in 4f-only CAS(n,7), 14 qubits, scalar X2C, cc-pVDZ-DK, AVAS with ANO-RCC projectors.
ROHF is an unreliable starting point (wrong shell occupation for Gd3+ and CeCl3, wrong configuration for
Pr3+, symmetry-broken core polarisation splitting the 4f levels by 5–18 mHa); state-averaged CASSCF over the
full high-spin multiplet manifold gives spherical orbitals, after which sparse ED of the qubit Hamiltonian
reproduces exact term degeneracies (Ce3+ 2F ×7; Pr3+ 3H ×11 at 0, 3F ×7 at 27.19 mHa, 3P ×3 at 126.88 mHa;
Eu3+ 7F ×7). UCCSD-VQE is exact on these spaces (Ce3+ 7.6e-9 Ha, 6 parameters, 168 CNOTs; Pr3+ 2e-8 Ha, 20
parameters, 1400 CNOTs, one restart). The Hamiltonians are dense: 3300 JW terms and λ ≈ 50 Ha at 14 qubits,
and λ varies by a factor 3 with the orbital choice. Details, the two-component comparison and the
interpretation are in `lanthanide_feasibility.md` §2 and §4.

## 5. CeCl3 active-space ladder

Five levels selected by orbital composition (4f; 4f+5d; 4f+5d+σ donors; 4f+all Cl 3p; everything), 14 to 42
qubits, with CASCI references up to 9.2e7 determinants (the CAS(19,16) level ran as the single cluster job).
Key numbers: L1 λ = 56 Ha (1549 terms); L3 (34 qubits) λ = 148 Ha, E = −10239.41414 Ha; L4 (32 qubits)
λ = 92 Ha, E = −10239.41111 Ha, 3.0 mHa above L3; L5 (42 qubits) λ = 182 Ha, no classical reference. The
ligand-field splitting of the 4f manifold with state-averaged orbitals is 0.68–5.96 mHa in the D3h 1-2-2-1-1
pattern. Full table and discussion: `lanthanide_feasibility.md` §3.

## 6. Resource proxies

For every Hamiltonian we report λ and the number of qubitised-walk queries π·λ/(2ε) for ε = 1 mHa, a
Toffoli proxy (queries × ⌈log2 terms⌉ × qubits for a sparse LCU) and a logical-qubit proxy (system +
LCU register + phase register). These are 1-norm proxies; they are not optimised fault-tolerant resource
counts (no double factorisation or THC, no surface-code or physical-qubit assumptions). They are useful
for *comparing* Hamiltonians along the ladder (fig. 1), not for absolute hardware estimates.

## 7. Reproducibility

`scripts/run_system.py configs/<name>.yaml` regenerates any JSON; `scripts/make_plots.py` regenerates the
figures and `results/summary_table.csv`; `scripts/run_soc.py` the two-component comparison;
`pytest tests` the validation suite. `progress.md` documents every failed attempt and the reason.
