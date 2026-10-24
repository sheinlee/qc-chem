# From Active-Space Electronic Structure to Quantum Algorithms: Benchmarks and Resource Estimates for Strongly Correlated and Lanthanide Systems

Technical report, 2026-10-24, revised after an independent referee-style review (the review and the
response are logged in `progress.md`). All numbers were produced by the code in this repository from one
commit (67e551f; `results/*.json` carry the SHA, `results/summary_table.csv`, figures in `results/figures/`).
Software: PySCF 2.14.0, OpenFermion 1.8.1, PennyLane 0.45.1 (lightning.qubit), SciPy 1.17.1. Hardware: a
Ryzen 9 5900 workstation for everything except one 6-minute 8-core job (the CAS(19,16) CASCI reference).

## 1. Question and scope

Can a variational quantum eigensolver reproduce CASCI/FCI-quality energies when given exactly the same
active-space Hamiltonian, and what does that Hamiltonian cost in qubits, Pauli terms, 1-norm, variational
parameters and circuit depth, from textbook molecules to a lanthanide complex? The lanthanide part asks
where, along a chemically motivated active-space ladder for CeCl3, the problem leaves the regime that can be
simulated or solved classically today, and what a qubitised phase-estimation algorithm would need on a
1-norm proxy basis. Everything is a statevector or exact calculation; no hardware noise is modelled.

## 2. Pipeline and validation

**Active-space Hamiltonian.** PySCF RHF/ROHF (scalar-relativistic sfx2c1e for Ln), CASCI or CASSCF in
canonical or AVAS-selected orbitals; `H = ecore + Σ h_pq a†a + ½ Σ (pq|rs) a†a†aa` from `CASCI.get_h1eff()` /
`get_h2eff()`, with the two-electron tensor averaged over its eight permutational symmetries, converted to
the OpenFermion physicist convention with `0.5·eri.transpose(0,2,3,1)` and expanded to interleaved spin
orbitals (2p = α, 2p+1 = β). **Acceptance test:** sparse exact diagonalisation (ED) of the Jordan-Wigner
(JW) Hamiltonian in the (N, S_z) sector reproduces the PySCF CASCI energy: H2 ≤ 2e-15 Ha, LiH(4,4) 4e-15,
N2 CAS(6,6) ≤ 7e-11 over the whole curve, Ce3+ 2e-9, Eu3+ 1e-11; for the degenerate state-averaged Pr3+
space the ED lies 8e-8 Ha *below* CASCI, which is the CI solver's convergence on a degenerate root, and ED is
then the better reference. Molecules are covered by `tests/unit/test_hamiltonian.py`; the lanthanide values
are in the JSON records.

**Mappings and symmetry.** JW and Bravyi-Kitaev give identical spectra and term counts for the four molecules
(BK was not run for the lanthanides). Z2 tapering (PennyLane generators; the symmetry sector is taken from
the actual open- or closed-shell reference determinant) removes 2–4 qubits: H2 4→1, LiH(4,4) 8→5, N2/H2O
12→8, Ln 4f spaces 14→12 (state-averaged orbitals) or 14→10 (ROHF-based orbitals of the same ion), CeCl3
4f+5d 24→22; the tapered ground-state energy agrees with CASCI to 1e-13.

**1-norm identity.** The JW Pauli 1-norm λ (identity excluded) is computed from the explicit expansion and,
independently, from the spatial-orbital integrals with the expression of Koridon et al. (PRR 3, 033127,
2021). With 8-fold-symmetrised integrals the two agree to 1e-14 (H2, LiH), 1e-11 (N2, H2O), 1e-7 (free Ln
ions) and 2e-7 relative (CeCl3 L1/L2); the residual at the CeCl3 level is not understood but immaterial for
resource comparison. The integral formula is used for the 32–42-qubit CeCl3 Hamiltonians, where no
expansion is performed. It is only a proxy for complex (two-component) integrals and is labelled so.

**UCCSD-VQE.** Spin-conserving singles and doubles from an (n_α, n_β) reference determinant, disentangled
first-order Trotter UCCSD (`qml.UCCSD`), `lightning.qubit` with adjoint gradients, BFGS (gtol 1e-6, 500
iterations) with an L-BFGS-B polish on precision loss. Unit tests: θ = 0 reproduces E_HF; a stored reference
angle set reproduces its energy; open-shell excitation counts (117, 20, 6 for N2, Pr3+, Ce3+). Integration
tests: |E_VQE − E_FCI| < 1e-6 (H2) and 1e-5 (LiH). Depth and CNOT counts are reported after PennyLane's
decomposition of each fermionic excitation to {CNOT, single-qubit rotations}; this is a naive compilation and
the counts are upper bounds. Molecular scans use a single start from the HF determinant. For the
state-averaged lanthanide spaces, whose natural orbitals carry no leading determinant, up to two restarts
from random angles (σ = 0.2) are allowed when the first run ends > 1e-5 Ha above the known reference; the
number used is recorded (0 for all free ions, 1 for CeCl3 L1-SA, which had converged to the first excited
state 0.68 mHa above the ground state).

**Spin-orbit block.** The one-electron X2C spin-orbit extension (`qcchem/soc.py`) is validated against the
analytic Dirac 2p fine structure of hydrogenic Ne9+ and Ar17+ in a large even-tempered basis: splitting
within 1 %, absolute 2p1/2 energy within 1e-4 relative (`tests/unit/test_soc_hydrogenic.py`).

## 3. Molecules: where first-order UCCSD stops being chemically accurate

Molecular curves use CASSCF orbitals at every geometry, so "CASCI" below is FCI within the CASSCF-optimised
active space; HF is the RHF energy.

| System | qubits (tapered) | JW terms | λ / Ha | UCCSD θ | CNOT | depth | UCCSD error at the CASCI minimum |
|---|---|---|---|---|---|---|---|
| H2 / STO-3G (2,2) | 4 (1) | 14 | 1.89 | 3 | 64 | 100 | ≤ 2e-13 (exact) |
| LiH / STO-3G (4,4), R = 1.6 Å | 8 (5) | 192 | 6.82 | 26 | 1312 | 1786 | 8e-9 |
| H2O / STO-3G (8,6), R = 0.96 Å | 12 (8) | 550 | 27.7 | 92 | 6976 | 8878 | 9.8e-5 |
| N2 / STO-3G (6,6), R = 1.2 Å | 12 (8) | 382 | 14.9 | 117 | 9024 | 11496 | 7.4e-4 |

The N2 curve (fig. 2) is the central molecular result. From 0.9 to 2.4 Å the optimiser reaches a stationary
point at every geometry (final gradient norm < 1e-6), so the residual error is a property of the ansatz at
that stationary point (not a proof that it is the global minimum; Trotter ordering and other starts were not
explored). The UCCSD error rises from 0.08 mHa at 0.9 Å to 0.74 mHa at the minimum (1.2 Å), crosses the
chemical-accuracy line between 1.35 Å (1.43 mHa) and 1.5 Å (1.75 mHa), dips slightly at 1.7 Å (1.42 mHa),
and then grows to 5.5 mHa at 2.0 Å and 13.4 mHa at 2.4 Å; iterations needed rise from 12 to 55. The 3.0 Å
point is excluded from this statement: there the CASSCF converged to a symmetry-broken orbital set (1778
instead of 382 Pauli terms), and the optimiser hit the 500-iteration cap at 0.22 mHa with a gradient of
5e-5. The Hamiltonian is not what gets harder along the curve (382 terms, λ between 12.2 and 16.3 Ha); the
multireference character is. The H2O symmetric stretch shows the same behaviour with a later onset: 0.06–0.25
mHa up to 1.8 Å, 3.1 mHa at 2.2 Å, 5.8 mHa at 2.6 Å.

Cost observation: for 12 qubits the UCCSD circuit for N2 CAS(6,6) is 9024 CNOTs deep in this compilation,
for a problem whose classical FCI space is 400 determinants. Efficient excitation compilations reduce this
substantially; the number is an upper bound.

## 4. Lanthanide free ions: the reference problem comes first

Ce3+, Pr3+, Eu3+, Gd3+ in 4f-only CAS(n,7), 14 qubits, scalar X2C, cc-pVDZ-DK, AVAS with ANO-RCC projectors.
ROHF is an unreliable start (wrong shell occupation for Gd3+ and CeCl3, wrong configuration for Pr3+,
symmetry-broken core polarisation splitting the 4f levels by 5–18 mHa); state-averaged CASSCF over the full
high-spin multiplet manifold restores spherical orbitals at a variational cost of 3.5–4.5 mHa, after which
sparse ED of the qubit Hamiltonian reproduces exact term degeneracies (Ce3+ 2F ×7; Pr3+ 3H ×11 at 0, 3F ×7 at
27.19 mHa, 3P ×3 at 126.88 mHa; Eu3+ 7F ×7). UCCSD-VQE is exact on these spaces (Ce3+ 7.8e-9 Ha, 6
parameters, 168 CNOTs; Pr3+ −2.3e-8 Ha vs CASCI, 20 parameters, 1400 CNOTs). The Hamiltonian's term count
(845 vs 3297), λ (16.4 vs 54.6 Ha) and tapering yield (4 vs 2 symmetries) for the *same* Ce3+ ion depend on
the orbital rotation chosen within the degenerate 4f shell. Discussion and the two-component comparison:
`lanthanide_feasibility.md` §2 and §4.

## 5. CeCl3 active-space ladder

Five levels selected by orbital composition (4f; 4f+5d; 4f+5d+σ donors; 4f+all Cl 3p; everything), 14 to 42
qubits, CASCI references up to 9.2e7 determinants (the CAS(19,16) level ran as the single cluster job).
Key numbers: L1 λ = 56 Ha (1549 terms); L3 (34 qubits) λ = 148 Ha, E = −10239.41414 Ha; L4 (32 qubits)
λ = 92 Ha, E = −10239.41111 Ha; L5 (42 qubits) λ = 183 Ha, no classical reference. The level energies are
not variationally consistent (CASCI in re-projected orbitals; the relaxation CASSCF(1,8) is already at
−10239.41113 Ha), so they are not used to rank ligand channels. With state-averaged orbitals the 4f
ligand-field splitting is 0.68–5.96 mHa in the D3h 1-2-2-1-1 pattern. Full table and discussion:
`lanthanide_feasibility.md` §3.

## 6. Resource proxies

For every Hamiltonian the records give λ, the qubitised-walk query count π λ / (2ε) for ε = 1 mHa, a
sparse-LCU Toffoli proxy (queries × Pauli terms; SELECT by unary iteration) and a logical-qubit proxy. They
are Fock-space 1-norm numbers without double factorisation / THC and without error-correction overhead,
intended for comparing Hamiltonians along the ladder (fig. 1), not for absolute hardware estimates.

## 7. Reproducibility

`scripts/run_system.py configs/<name>.yaml` regenerates any JSON (records package versions and the git
SHA); `scripts/run_soc.py configs/<name>.yaml [n_levels]` the two-component comparison;
`scripts/make_plots.py` the figures and `results/summary_table.csv`; `pytest tests` the 19-test validation
suite (~40 s). `progress.md` documents every failed attempt, the external review and the changes it caused.
Hard-coded local paths exist only in the convenience wrappers (`scripts/qpy.sh`, `scripts/sync.sh`,
`slurm/`); the package and scripts themselves are path-independent.
