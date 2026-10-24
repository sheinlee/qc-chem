# From Active-Space Electronic Structure to Quantum Algorithms

Benchmarks and resource estimates for strongly correlated and lanthanide systems: PySCF active spaces →
qubit Hamiltonians → UCCSD-VQE / exact diagonalisation / 1-norm resource proxies, from H2 to a 42-qubit
CeCl3 active space, with a one-electron two-component (spin-orbit) extension.

**Read first:** `report/technical_report.md` (pipeline, validation, molecular benchmarks) and
`report/lanthanide_feasibility.md` (the lanthanide ladder, spin-orbit comparison, three-tier assessment,
proposed projects). `progress.md` is the honest log of what failed and why.

## Main results in one paragraph

The folded active-space Hamiltonian reproduces PySCF CASCI to 1e-11 Ha or better after Jordan-Wigner or
Bravyi-Kitaev mapping; the Pauli 1-norm computed from integrals equals the explicit expansion to 1e-13 and
is therefore exact for Hamiltonians too large to expand (32–42 qubits). First-order UCCSD-VQE converges
(gradient < 1e-6) at every point of the N2 CAS(6,6) curve but its intrinsic error rises from 0.2 mHa at
equilibrium to 16.5 mHa at 3.0 Å, crossing chemical accuracy at ~1.4 Å, with a 9024-CNOT circuit for 12
qubits. Lanthanide 4f-only spaces (14 qubits) are classically trivial (7–21 determinants) but the reference
problem is not: HF puts Ln electrons in the wrong shell and breaks spherical symmetry; state-averaged
CASSCF fixes it and the qubit ED then reproduces exact multiplet degeneracies (Pr3+ 3H/3F/3P 11/7/3).
Along a composition-selected CeCl3 ladder the chemically meaningful space is 4f + 5d + ligand σ (34 qubits,
λ = 148 Ha, still solvable classically in 10 minutes); 4f + all Cl 3p (32 qubits) gains less energy than
5d + σ. Spin-orbit coupling adds only ~7 % Pauli terms and < 1 % to λ and removes one Z2 symmetry, but its
energetic effect (13–27 mHa) dwarfs the ligand field (≤ 6 mHa).

## Layout

```
qcchem/        integrals (PySCF, AVAS, SA-CASSCF, reference relaxation), hamiltonian (JW/BK, stats, λ identity,
               tapering), benchmark (sparse ED in (N,Sz) sector), vqe (open-shell UCCSD, lightning+adjoint),
               resources (QPE proxies), soc (one-electron X2C spin-orbit block)
configs/       one YAML per system / active-space level
scripts/       run_system.py (YAML -> JSON), run_soc.py, make_plots.py, qpy.sh (env wrapper), sync.sh (cluster)
results/       JSON per system, integrals/*.npz (h1, eri, ecore, mo_act), summary_table.csv, figures/
tests/         unit (deterministic) and integration (optimizer runs)
report/        technical_report.md, lanthanide_feasibility.md, skill_gap_analysis.md
slurm/         job template used for the one cluster run (CAS(19,16) CASCI)
```

## Reproduce

```
conda env create -f environment.yml          # pyscf 2.14.0, openfermion 1.8.1, pennylane 0.45.1 (+lightning)
conda activate qchem                          # or use scripts/qpy.sh on Windows
pytest tests                                  # 17 tests, ~30 s
python scripts/run_system.py configs/n2_cas66.yaml
python scripts/run_soc.py configs/ce3_4f_sa.yaml
python scripts/make_plots.py
```

Hard limits enforced in code: VQE ≤ 12 qubits (14 for ≤ 2-electron spaces), sparse ED ≤ 16, explicit Pauli
expansion ≤ 24; above that only integral-based λ. Everything ran on a desktop except one 6-minute 8-core
CASCI job.

## Figures

- `results/figures/fig1_ladder.png` qubits / Pauli terms / λ / QPE-proxy queries along the ladder with the
  simulation-regime bands.
- `results/figures/fig2_curves.png` N2 and H2O: energies and E − E_CASCI for HF and UCCSD-VQE.
- `results/figures/fig3_scaling.png` Pauli terms and λ vs active orbitals (N^4.2 fit).

## Known limitations

cc-pVDZ-DK basis, scalar X2C, one-electron SOC without screening (overestimates Ce3+ 2F splitting 2.3×),
AVAS thresholds as stated, no noise model, no ADAPT-VQE yet, 1-norm proxies instead of factorised resource
estimates. See report §8 and `progress.md`.
