# Task brief for Claude Code

Project: **From Active-Space Electronic Structure to Quantum Algorithms: Benchmarks and Resource Estimates for Strongly Correlated and Lanthanide Systems**

Central claim to test with our own numbers: 4f-only lanthanide active spaces are classically trivial; the problem only becomes hard once ligand orbitals and relativistic/spin-orbit effects enter, at which point NISQ-VQE is already out of reach and only fault-tolerant QPE remains worth discussing. A negative result ("no quantum advantage for realistic Ln chemistry today") is a valid outcome if quantified.

Read `PLAN.md` first. Keep `progress.md` updated at the end of every stage: what ran, where, versions, failures and causes, next step.

## Hard engineering and performance constraints (MUST FOLLOW)

1. **Environment.** Local Windows conda env `qchem` (python 3.11, pyscf 2.14.0, openfermion 1.8.1, pennylane 0.45.1, pennylane-lightning, OpenBLAS). Always run through `conda run -n qchem python ...` or an activated shell; calling `python.exe` directly crashes with `0xC06D007F` because `Library\bin` is not on PATH. Pin versions in `environment.yml`.
2. **Simulator and gradients.** PennyLane `lightning.qubit` with `diff_method="adjoint"`; optimizer `scipy.optimize.minimize(method="BFGS", jac=True)`. Never use parameter-shift above 6 qubits.
3. **Qubit ceilings.**
   - Statevector VQE / ADAPT-VQE optimization: max 12 qubits; 14 qubits only for tiny excitation spaces (Ce3+ 4f1, Pr3+ 4f2).
   - Exact qubit diagonalization: sparse `scipy.sparse.linalg.eigsh` only, max 16 qubits. Never build a dense Hamiltonian matrix.
   - Explicit Jordan-Wigner construction and Pauli statistics: max 24 qubits. Above that compute the 1-norm lambda and term counts directly from the integrals without expanding Pauli strings, and extrapolate the O(N^4) scaling.
   - Noise: shot noise up to 12 qubits; depolarizing noise with `default.mixed` max 6 qubits, evaluated at fixed optimal angles theta* over p in [1e-4, 1e-1]. Never re-optimize under noise. Noise is optional and is the first thing to cut.
4. **Active-space integral folding.** Use `mcscf.CASCI.get_h1eff()` and `ao2mo.restore(1, mc.get_h2eff(), ncas)`. Convert chemist (pq|rs) to the OpenFermion physicist two-body tensor via `0.5 * eri.transpose(0, 2, 3, 1)`, then `spinorb_from_spatial`, and include `ecore` explicitly. Acceptance test: OpenFermion sparse ED equals PySCF CASCI `e_tot` within 1e-8 Ha (verified locally on 2026-10-24 for N2 CAS(6,6): 3e-13 Ha, 383 Pauli terms).
5. **Tests.** Separate `tests/unit` (fixed geometry, fixed reference parameters, deterministic statevector, |E - E_ref| < 1e-8) from `tests/integration` (optimizer actually runs, |E_VQE - E_FCI| < 1e-5 Ha, also report whether chemical accuracy 1.6e-3 Ha is reached). Optimizer flakiness must not break unit tests.
6. **cluster cluster (ssh alias `cluster`, env `~/anaconda3/envs/qchem`).** Use it only for what the local machine cannot do: X2C-CASSCF on real Ln complexes, large CASCI, DMRG. Login node: smoke tests only, <= 4 threads, < 4 minutes. Everything else via `sbatch`, 1 node, 4-8 CPUs, `--time` <= 60 min, at most one queued job. Never write large files to `/home`; scratch is `/data/$USER/qc-chem`. After each job dump integrals (`.npz`) and results (`.json/.csv`), pull them to the local machine with `scripts/sync.ps1 pull`, then delete scratch. Keep the local copy authoritative; cluster access is not guaranteed long-term.
7. **Wording.** Write "evaluate the degradation of UCCSD-VQE in the stretched-bond regime"; do not presuppose failure. Report what the numbers show.

## Systems and layers

| Layer | System | Active space | Qubits | Quantum side | Classical reference | Where |
|---|---|---|---|---|---|---|
| Validation | H2 STO-3G | (2,2) | 4 | UCCSD-VQE, JW vs BK, tapering | FCI | local |
| Validation | LiH STO-3G | (2,2), (4,4) frozen core | 4, 8 | same | FCI within CAS | local |
| Core | N2 dissociation, 8-10 points | CAS(6,6) | 12 | UCCSD vs fermionic ADAPT-VQE | HF, CASCI = FCI within CAS(6,6) | local |
| Core | H2O symmetric stretch | CAS(6,5) or (8,6) | 10-12 | same | same | local |
| Ln baseline | Ce3+ 4f1, Pr3+ 4f2 | CAS(1,7), (2,7) | 14 | VQE runs | CASCI; fixed N and Sz sector | local |
| Ln baseline | Eu3+ 4f6, Gd3+ 4f7 | CAS(6,7), (7,7) | 14 | statistics + sparse ED only, no VQE | CASCI, fixed Sz, selected low states | local |
| Core selling point | LnCl3 (fallback LnF3 / [LnCl6]3-), Ln = Ce, Eu; scalar X2C (sfx2c1e), x2c-SVPall | 4f-only (n,7) -> +5d (n,12) -> +ligand (n,15+) | 14 -> 24 -> 30+ | JW statistics <= 24 qubits; lambda from integrals above | CASSCF/CASCI, DMRG if needed | cluster sbatch |
| Extension (S4) | same complex, two-component X2C + SOC (complex Kramers spinors) | same | same | statistics only: nonzero integrals, broken Z2 spin symmetry (only particle parity remains), Pauli term inflation | two-component CASCI if feasible | cluster sbatch |
| Optional | one 3d transition-metal CAS <= (12,12) | | <= 24 | statistics | CASCI | if time |

State in the report: free-ion 4f-only spaces in a fixed Sz sector contain very few determinants (Ce3+: 14; Eu3+ high-spin Sz=3: 7); Hund's rule gives most of the answer. Their value is validating mapping and tapering on degenerate 4f orbitals and providing baseline Hamiltonian statistics, not demonstrating strong correlation. The spin-free model does not reproduce the 7F_J multiplet structure of Eu(III).

## Per-system JSON record

active electrons, active orbitals, spin orbitals, qubits before/after tapering, Pauli terms, lambda = sum |h_j|, VQE parameter count, circuit depth, two-qubit gates, optimizer energy evaluations, energy error vs CASCI/FCI, wall time, machine, package versions.

## Figures

1. Core figure: x-axis H2 -> N2 -> Ce3+ -> Eu3+ -> LnCl3(4f) -> LnCl3(4f+5d) -> LnCl3(4f+ligand) -> two-component LnCl3; log-scale panels for qubits, Pauli terms, lambda, ADAPT operators or UCCSD parameters, estimated QPE cost; shaded regions "classically trivial" / "NISQ-simulatable benchmark" / "fault-tolerant target".
2. E - E_FCI vs bond length for N2 and H2O (HF, UCCSD-VQE, ADAPT-VQE) with the chemical-accuracy line.
3. Scaling: Pauli terms and lambda vs number of active orbitals, fitted exponent, extrapolation to 40-60 qubits.
4. ADAPT vs UCCSD: parameters, CNOTs, depth, energy evaluations.
5. Optional: energy error vs depolarizing rate at fixed theta* (H2, LiH 6 qubits).

## Resource estimation

For each system report N_qubits, N_Pauli, lambda, and for epsilon = 1 mHa a qubitization-style proxy Toffoli count and logical-qubit count. State explicitly: these are Pauli-1-norm proxy estimates, not optimized fault-tolerant resource counts (no double factorization / THC, no surface-code or physical-qubit assumptions); list those as future work.

## Stages

- **S1 (must complete):** H2 -> LiH -> N2. Full pipeline with acceptance test, JW/BK, tapering, sparse ED, UCCSD-VQE, unit/integration tests, N2 curve locally. Deliver results JSON, first resource table, progress.md.
- **S2 (core):** N2 + H2O, UCCSD vs fermionic ADAPT (qubit pool optional). Figures 2 and 4. Noise (figure 5) last, cut first.
- **S3 (chemistry specialization):** Ce3+/Pr3+ VQE and Eu3+/Gd3+ statistics locally; LnCl3 scalar-X2C CASSCF on cluster via sbatch with three active-space sizes, integrals pulled back and analyzed locally. Figures 1 and 3, draft of `report/lanthanide_feasibility.md` with the three-tier classification (demonstrable today / near-term hardware / fault-tolerant only). If CASSCF does not converge, record orbital selection and attempts (AVAS, 4f projection initial guess), fall back to CASCI, and say so.
- **S4 (research-level extension):** two-component X2C + SOC spinor Hamiltonian statistics (report plus one standalone script; time-boxed), QPE proxy analysis, three proposed research projects ranked by feasibility, novelty and cost. Then review the repo as a quantum-algorithms interviewer and a referee; fix README and reports. Final pull; confirm the local copy is complete; clean cluster scratch.

Run the calculations, inspect failures, debug, record them. Do not stop at example scripts.
