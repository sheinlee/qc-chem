# From Active-Space Electronic Structure to Quantum Algorithms

Benchmarks and resource estimates for strongly correlated and lanthanide systems: PySCF active spaces →
qubit Hamiltonians → UCCSD-VQE / exact diagonalisation / 1-norm resource proxies, from H2 to a 42-qubit
CeCl3 active space, with a one-electron two-component (spin-orbit) extension.

**Read first:** `report/technical_report.md` (pipeline, validation, molecular benchmarks) and
`report/lanthanide_feasibility.md` (the lanthanide ladder, spin-orbit comparison, three-tier assessment,
proposed projects). `progress.md` is the honest log of what failed and why.

## Main results in one paragraph

The folded active-space Hamiltonian reproduces PySCF CASCI after Jordan-Wigner or Bravyi-Kitaev mapping
(1e-11 Ha or better for molecules and non-degenerate ions; 8e-8 on a degenerate state-averaged Pr3+ root,
where exact diagonalisation is the better reference); the Pauli 1-norm computed from integrals equals the
explicit expansion to 1e-11 (molecules) and 1e-7 (lanthanides) and is therefore used for Hamiltonians too
large to expand (32–42 qubits). First-order UCCSD-VQE reaches a stationary point (gradient < 1e-6) at every
point of the N2 CAS(6,6) curve but its error rises from 0.08 mHa at 0.9 Å to 0.74 mHa at the minimum and
13.4 mHa at 2.4 Å, crossing chemical accuracy near 1.4 Å, with a 9024-CNOT circuit for 12 qubits in a naive
compilation. Lanthanide 4f-only spaces (14 qubits) are classically trivial (7–21 determinants) but the
reference is not: HF puts Ln electrons in the wrong shell and breaks spherical symmetry; state-averaged
CASSCF fixes it and the qubit ED then reproduces exact multiplet degeneracies (Pr3+ 3H/3F/3P = 11/7/3; with
spin-orbit coupling 3H4/3H5/3F2/3H6/3F3 = 9/11/5/13/7). Along a CeCl3 active-space ladder built on one 4f1 ROHF reference (SOMO unprojected, so the one-electron
levels reproduce ROHF exactly), the ligand electrons hardly interact with the 4f shell in the CI unless Ce 5d
is in the active space (occupations 2/1/0 to 1e-5 without 5d; in the nested 4f+σ → 4f+5d+σ pair the five 5d
orbitals carry 6.4 of 6.5 mHa with a total occupation of 0.01 e), and a fully variational CASSCF does not
preserve the chemically defined spaces (σ drifts into Ce 4d, or empty 4f into σ*). One-electron spin-orbit
coupling adds ~7 % Pauli terms and +0.5 % to λ and removes one Z2 symmetry but overestimates the Ce3+ 2F
splitting 2.3×; adding the Breit-Pauli spin-orbit mean field (convention fixed against Dirac on hydrogenic
ions) screens it by 0.46 and gives 2397 cm⁻¹ (experiment 2253) and Pr3+ 3H5/3H6 at 2259/4643 cm⁻¹
(experiment 2152/4389). Factorised fault-tolerant estimates put the chemically motivated 34-qubit CeCl3 space
at ~1e9 Toffolis and ~870 logical qubits per 1-mHa phase estimation (DF, 1e-4, CASCI error 0.03 mHa) and
the full 42-qubit space at 1.3e9 / ~1000; λ-minimising orbital rotations within the same subspaces lower λ
by 15 % (34 qubits) and 26 % (42 qubits) with the CASCI energy unchanged. ADAPT-VQE with the same operator
pool reaches 1e-6 Ha on the stretched N2 points where UCCSD is 2-13 mHa off, using 86-97 operators (47-51
distinct) and fewer CNOTs than UCCSD, and stalls together with UCCSD at one H2O point.

## Layout

```
qcchem/        integrals (PySCF, AVAS, SA-CASSCF, reference relaxation), hamiltonian (JW/BK, stats, λ identity,
               tapering), benchmark (sparse ED in (N,Sz) sector), vqe (open-shell UCCSD, lightning+adjoint),
               resources (QPE proxies), soc (one-electron X2C spin-orbit block), somf (Breit-Pauli SOMF screening),
               adapt (fermionic-pool ADAPT-VQE on an exact statevector)
configs/       one YAML per system / active-space level
scripts/       run_system.py (YAML -> JSON), run_soc.py, run_adapt.py, resource_estimates_ft.py, lambda_orbital_opt.py, check_nesting.py,
               dhf_benchmark.py, make_plots.py, qpy.sh (env wrapper), sync.sh (cluster)
results/       JSON per system, ft_*.json (resource estimates), lamopt_*.json, soc_*.json, integrals/*.npz, summary_table.csv, figures/
tests/         unit (deterministic) and integration (optimizer runs)
report/        technical_report.md, lanthanide_feasibility.md
slurm/         job template used for the one cluster run (CAS(19,16) CASCI)
```

## Reproduce

```
conda env create -f environment.yml          # pyscf 2.14.0, openfermion 1.8.1, pennylane 0.45.1 (+lightning)
conda activate qchem                          # or use scripts/qpy.sh on Windows
pytest tests                                  # 21 tests, ~45 s
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
- `results/figures/fig3_scaling.png` Pauli terms and λ vs active orbitals (N^4.35 fit).

## Known limitations

cc-pVDZ-DK basis, scalar X2C, one-electron X2C SOC plus mean-field (SOMF) screening only, AVAS thresholds as
stated, CASCI rather than CASSCF across the CeCl3 ladder (CASSCF does not preserve the chemically defined spaces), no noise model, ADAPT-VQE only on <= 14 qubits with a spin-orbital pool, no THC, logical-level resource estimates without error-correction overhead, naive
CNOT compilation. Term counts, λ and tapering yields depend on the orbital basis for
degenerate shells. See `report/lanthanide_feasibility.md` §8 and `progress.md` (developer log, in Chinese).
