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
within 1 %, absolute 2p1/2 energy within 1e-4 relative (`tests/unit/test_soc_hydrogenic.py`). The
two-electron spin-orbit mean field (SOMF, Breit-Pauli form, `qcchem/somf.py`) reuses the integral drivers of
`pyscf.prop.gtensor`; its spin convention is fixed by the same hydrogenic test (`tests/unit/test_somf.py`,
BP vs Dirac 0.997/0.989 for Ne9+/Ar17+), and the Breit-Pauli one-electron block agrees with the X2C block to
0.5 % on Ce 4f.

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

**ADAPT-VQE** (fermionic pool, exact statevector, `qcchem/adapt.py`; table and discussion in
`lanthanide_feasibility.md` §7): with a loose gradient criterion it reproduces the UCCSD accuracy at every
point with 2-4x fewer parameters and CNOTs; with a tight criterion it reaches 1e-6 Ha at N2 1.5/2.0/2.4 A
(86-97 operators, 47-51 distinct, 7.8-8.8k CNOTs) where UCCSD stays at 1.7/5.5/13.4 mHa, and it stalls at
H2O 2.2 A (3.1 mHa, all pool gradients < 1e-4) exactly where UCCSD does; a spin-adapted pool halves the
operator count on N2 (43 ops, 2e-8 Ha) and a qubit pool escapes the H2O trough (2.7e-5 Ha, 75 Pauli
rotations). A bounded noise study (`scripts/run_noise.py`) is in `lanthanide_feasibility.md` §8.

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

Six levels on one reference (4f1 ROHF, E = −10239.411130 Ha, SOMO kept unprojected so that one-electron
spaces reproduce the ROHF energy exactly; AVAS orbitals; CASCI with CI convergence 1e-12; fixed-orbital and
natural occupations recorded): 4f (0.000 mHa); 4f+5d (0.000); 4f+σ donors without 5d (−0.093; 4f virtual
occupations 1–3e-5); 4f+5d+σ (−6.537; 5d occupations Σ 0.0105 e); 4f+all Cl 3p without 5d (−0.267; 4f
virtuals ≤ 1e-4); everything (42 qubits, no CI reference). Subspace nestedness was checked explicitly and
only nested pairs are compared. A fully variational CASSCF does not preserve these spaces (σ → Ce 4d without
a frozen core, 4f → σ* with one). λ: 56 / 116 / 66 / 147 / 93 / 182 Ha. Full table, nesting data and
discussion: `lanthanide_feasibility.md` §3; spin-orbit (one-electron and SOMF-screened) in §4.

## 6. Resource estimates

Factorised fault-tolerant costings (OpenFermion `resource_estimates`: sparse, single and double
factorisation; ε = 1 mHa) with the exact CASCI energy error of each factorisation where the CI is solvable,
and λ-minimising orbital rotations within the fixed active subspaces (CASCI invariance checked):
`lanthanide_feasibility.md` §5–6, `results/ft_*.json`, `results/lamopt_*.json`. The earlier Pauli-1-norm
proxies (queries × Pauli terms) are superseded; they overestimated the sparse-LCU Toffoli count by about two
orders of magnitude for the 42-qubit space. The JSON records still contain the proxy fields for traceability.

## 7. Reproducibility

`scripts/run_system.py configs/<name>.yaml` regenerates any JSON (records package versions and the git
SHA); `scripts/run_soc.py configs/<name>.yaml [n_levels]` the two-component comparison;
`scripts/make_plots.py` the figures and `results/summary_table.csv`; `scripts/dhf_benchmark.py` the
four-component screening benchmark; `pytest tests` the 21-test validation suite (~45 s). `progress.md` documents every failed attempt, the external review and the changes it caused.
Hard-coded local paths exist only in the convenience wrappers (`scripts/qpy.sh`, `scripts/sync.sh`,
`slurm/`); the package and scripts themselves are path-independent.
