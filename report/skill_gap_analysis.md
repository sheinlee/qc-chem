# Skill-gap analysis and 12-week transition plan: quantum computational chemistry

Date: 2026-10-24. Target roles: quantum algorithms scientist / quantum computational chemistry (quantum software companies, hardware vendors' applications teams, pharma/chemicals quantum groups, postdocs in quantum algorithms for chemistry).

Baseline (from the person's own account only): computational chemistry PhD; DFT, ab initio MD, AMOEBA force fields, free-energy calculations; lanthanide coordination chemistry (4f electronic structure, relativistic effects); molecular generative models; Python. A repository built in one day: PySCF active spaces to OpenFermion/PennyLane qubit Hamiltonians (JW, BK, Z2 tapering), UCCSD-VQE on lightning.qubit with adjoint gradients, sparse exact diagonalisation, exact Pauli 1-norm from integrals (Koridon identity), N2/H2O dissociation benchmarks, a CeCl3 active-space ladder to 42 qubits with CASCI references, a one-electron X2C spin-orbit two-component extension, and 1-norm-based QPE resource proxies. Not done: ADAPT-VQE, noise models, double-factorised/THC resource estimation, hardware runs, qubitization circuits, error mitigation.

## 1. Gap analysis

| Area | What such roles typically expect | What the person has | Gap |
|---|---|---|---|
| Mathematics | Linear algebra on operators (spectral theory, Hermitian block structure, polynomial approximation of functions of H); perturbation theory and Trotter error bounds via commutators; Chebyshev/Jacobi-Anger expansions; basic probability (concentration, variance of estimators); low-rank tensor factorisation (eigendecomposition, Cholesky, least-squares tensor fits). | Strong operator/matrix background from quantum chemistry; exact diagonalisation; ERI tensors handled daily. No evidence of formal complexity or approximation-theory training. | Medium |
| Quantum information | Circuit model, Clifford+T, Toffoli/T as cost units, stabiliser formalism, phase estimation, amplitude amplification, density matrices and channels (depolarising, amplitude damping), basic surface-code overhead reasoning. | Z2 tapering (stabiliser-adjacent), state-vector simulation. No channels, no QEC reasoning, no Toffoli accounting. | Medium to large |
| Quantum algorithms | Derive the qubitization walk operator and its Toffoli cost; explain LCU, block encoding, QSVT at the level of a whiteboard interview; implement ADAPT-VQE with a fermionic pool; know why the VQE measurement and optimisation costs are the real bottleneck; explain initial-state-overlap assumptions in QPE; read and reproduce a double-factorisation resource estimate. | UCCSD-VQE with adjoint gradients; 1-norm computed exactly; QPE cost proxies (so the central quantity lambda is understood). No qubitization circuit, no DF/THC, no ADAPT. | Large |
| Electronic structure theory | Active-space selection, CASCI/CASSCF/DMRG as baselines, orbital localisation/rotation effects on lambda, second quantization fluency, knowing where classical methods stop. | Core strength: PySCF active spaces, CASCI references, DFT and AIMD history, lanthanide multi-reference problems. Gap is only in DMRG/selected-CI baselines at scale. | Small |
| Programming / software engineering | Tested, packaged, reproducible code; CI; typed interfaces; profiling of Pauli-sum code; sparse and GPU simulation; benchmark harnesses with fixed seeds and logged environments. | Working repository in a day, adjoint-gradient simulation, sparse methods. Unknown: tests, CI, packaging, benchmark reproducibility. | Small to medium |
| Quantum SDKs | Fluent in Qiskit (the most common in job ads), plus one of PennyLane/Cirq; Qualtran or OpenFermion resource estimation; transpilation and hardware-aware compilation; running jobs on a cloud backend with noise models. | OpenFermion, PennyLane (lightning.qubit). No Qiskit, no hardware, no noise models, no transpilation. | Medium |
| Research skills | Reproduce a published result and state discrepancies; write a clear preprint-style note; ablation discipline; talk to both chemists and computer scientists; arXiv-level output in the quantum-chemistry-on-QC area. | PhD-level research habits; the benchmarks (N2/H2O, CASCI references) show reference-checking discipline. No visible output in this field yet. | Medium |

The central gap is not chemistry or coding; it is fault-tolerant algorithm literacy (qubitization, factorisations, cost accounting) and the lack of a public track record in this subfield.

## 2. Twelve-week plan

Assumes about 15 to 20 hours per week. Each week ends with one artefact committed to the repository (or a `notes/` folder). When a week slips, drop the optional parts, not the output.

| Week | Topic | Concrete output |
|---|---|---|
| 1 | Second quantization and fermion-to-qubit mappings at derivation level: anticommutation, JW, parity, BK (Seeley-Richard-Love), Pauli weight scaling, symmetry-based tapering (Bravyi et al.). | A short derivation note (hand-derived JW and BK images of a_p, a_p^dagger, and of a one-body excitation) plus a test that checks anticommutation relations and Pauli-weight scaling for 4 to 20 modes against the existing code. |
| 2 | Hamiltonian simulation basics: Trotter-Suzuki error, commutator bounds, QPE with Trotterised evolution, gate-count units. | Script: Trotterised QPE on a small tapered Hamiltonian; plot of measured Trotter error against the first/second-order commutator bound; one page on why the bound is loose. |
| 3 | ADAPT-VQE with a fermionic (singles/doubles) pool: commutator-based gradient screening, operator selection loop, reuse of adjoint differentiation. | `adapt.py` module with tests; N2 at stretched geometry (e.g. a 3 to 4 Angstrom point on the existing dissociation curve) showing energy error versus iteration against UCCSD and exact diagonalisation. |
| 4 | Qubit-ADAPT (Pauli-string pool, Tang et al.) and benchmarking. Count CNOTs, parameters, and gradient-measurement cost honestly. | Benchmark figure and table: error versus CNOT count and versus number of parameters for UCCSD, fermionic ADAPT, qubit-ADAPT on N2 and on CeCl3 ladder rungs up to about 16 qubits. |
| 5 | LCU, block encoding, qubitization (Low-Chuang): walk operator W, eigenphase relation, Toffoli cost of SELECT and PREPARE for a generic Pauli sum. | Derivation note with a hand resource count for H2/STO-3G; numerical check on a 4 to 6 qubit Hamiltonian that the dense walk-operator eigenphases equal +/- arccos(E_k / lambda). |
| 6 | Chemistry-specific qubitization: sparse (Babbush 2018), single factorisation, double factorisation (Berry 2019, von Burg 2021). Implement DF and truncation-error control. | `double_factorization.py`: DF of the ERI tensor, truncation threshold versus energy error, and lambda_DF versus the exact Pauli 1-norm from the repo for N2, H2O, CeCl3 rungs. |
| 7 | Reproduce a published tensor-hypercontraction estimate (Lee et al. 2021) on one repository Hamiltonian: THC fit, lambda_THC, logical qubits, Toffoli count per walk step and total. Check OpenFermion's resource-estimate module (if installed) as an independent cross-check, not a substitute for the hand derivation. | Reproduction report: the paper's cost formulas implemented from the text, applied to N2 (or H2O) active space, with a side-by-side with the paper's small-system numbers where available and a list of unexplained differences. |
| 8 | Error-mitigation basics (ZNE, symmetry verification, readout mitigation, Temme et al.) and noise models; one small run. | Noise-model study on the tapered H2 or a 4 to 6 qubit N2 Hamiltonian in a noisy simulator, then one run on a cloud device (e.g. IBM open plan, subject to current access) or an emulator; table of raw, symmetry-verified, ZNE-corrected energies with error bars. |
| 9 | Relativistic and spin-orbit Hamiltonians on qubits: two-component spinor basis, complex integrals, loss of Sz and spin symmetry, remaining symmetries (particle number, time reversal/Kramers), effect on tapering and term count. | Note plus code: for the CeCl3 X2C ladder, term count, qubit count after tapering, and Pauli 1-norm for scalar versus spin-orbit Hamiltonians; Kramers-pair check on the spectrum. |
| 10 | Fault-tolerant cost of spin-orbit Hamiltonians: DF/THC with complex two-component integrals (or the real-form embedding), cost comparison versus scalar-relativistic. | Table and figure: lambda and Toffoli estimates for scalar versus X2C two-component CeCl3 at fixed active-space size. This is the core of portfolio project 1. |
| 11 | Engineering and SDK breadth: tests, CI, packaging, pinned environments, reproducible benchmark scripts; cross-check one small circuit (e.g. UCCSD for H2 or a qubitization step) in Qiskit and one resource count in Qualtran or an equivalent tool. | Repository with CI badge, `pytest` suite, README with one-command reproduction of all headline figures, and a Qiskit cross-check script. |
| 12 | Write-up and talk. | A 12 to 15 slide talk and a 6 to 10 page arXiv-style write-up of the main result (spin-orbit lambda/resource scaling, with the ADAPT benchmark as a secondary section); practise once with a colleague. |

### Reading list (16 papers, 3 textbooks)

Papers:
1. McClean et al., "The theory of variational hybrid quantum-classical algorithms" (2016): VQE formalism and the measurement/optimisation problem.
2. Seeley, Richard, Love, "The Bravyi-Kitaev transformation for quantum computation of electronic structure" (2012): the mapping derivation you need for week 1.
3. Bravyi, Gambetta, Mezzacapo, Temme, "Tapering off qubits to simulate fermionic Hamiltonians" (2017): the theory behind your Z2 tapering.
4. Grimsley et al., "An adaptive variational algorithm for exact molecular simulations on a quantum computer" (2019): ADAPT-VQE.
5. Tang et al., "Qubit-ADAPT-VQE" (2021): hardware-efficient pool and its CNOT trade-offs.
6. Reiher et al., "Elucidating reaction mechanisms on quantum computers" (2017): the FeMoco benchmark that framed the field and set the cost vocabulary.
7. Low and Chuang, "Hamiltonian simulation by qubitization" (2019): the walk operator and its eigenphase relation.
8. Babbush et al., "Encoding electronic spectra in quantum circuits with linear T complexity" (2018): first chemistry qubitization with explicit Toffoli counts.
9. Berry et al., "Qubitization of arbitrary basis quantum chemistry leveraging sparsity and low rank factorization" (2019): introduces single/double factorisation costs.
10. von Burg et al., "Quantum computing enhanced computational catalysis" (2021): DF with careful resource accounting on Ru catalysts.
11. Lee et al., "Even more efficient quantum computations of chemistry through tensor hypercontraction" (2021): the THC estimate you will reproduce.
12. Koridon et al., "Orbital transformations to reduce the 1-norm of the electronic structure Hamiltonian" (2021): extends what the repository already does.
13. Lee et al., "Evaluating the evidence for exponential quantum advantage in ground-state quantum chemistry" (2023): the sceptical counterweight; essential for any honest crossover study.
14. Temme, Bravyi, Gambetta, "Error mitigation for short-depth quantum circuits" (2017): ZNE and probabilistic error cancellation.
15. Huang, Kueng, Preskill, "Predicting many properties of a quantum system from very few measurements" (2020): classical shadows, for measurement reduction.
16. Saue, "Relativistic Hamiltonians for chemistry: a primer" (2011): connects your relativistic knowledge to the Hamiltonian you are encoding.

Textbooks:
- Nielsen and Chuang, *Quantum Computation and Quantum Information*: circuits, QPE, channels; the shared vocabulary of the field.
- Helgaker, Jorgensen, Olsen, *Molecular Electronic-Structure Theory*: second quantization and orbital-rotation machinery, which every mapping and factorisation paper assumes.
- Dyall and Faegri, *Introduction to Relativistic Quantum Chemistry*: for the spin-orbit/X2C angle that is the person's differentiator.

## 3. Portfolio projects (ranked)

### Project 1: Spin-orbit lanthanide Hamiltonians and fault-tolerant resource estimates (extension of the repository)
- **Question:** How much do spin-orbit coupling and two-component structure change the 1-norm, factorisation rank, qubit count after tapering and Toffoli cost of simulating Ln complexes, compared with a scalar-relativistic Hamiltonian at the same active-space size?
- **Deliverables:** (i) a validated two-component pipeline (X2C spinors to qubit Hamiltonian, checked against CASCI/two-component references on small cases); (ii) DF and THC cost estimates for CeCl3 rungs, scalar versus spin-orbit; (iii) symmetry analysis (Kramers, number parity) and its effect on tapering; (iv) extension to a second 4f count if feasible (e.g. a 4f2 or 4f7 ion, only if the SCF converges credibly); (v) a reproducible benchmark and a preprint.
- **Why not a toy:** real 4f active spaces up to 42 qubits, an independent classical reference at every rung, and costs compared on equal footing. The existing literature on resource estimates is dominated by light-element systems and scalar Hamiltonians.
- **Effort:** about 6 to 8 weeks part-time beyond the 12-week plan's weeks 7 to 10.
- **Lanthanide link:** direct; the 4f spin-orbit splitting is much larger than chemical accuracy, so a scalar treatment cannot be the right model for the question.

### Project 2: ADAPT-VQE and ansatz compactness on strongly correlated active spaces
- **Question:** For stretched N2 and small CeCl3 active spaces, how many CNOTs and measurements do fermionic ADAPT and qubit-ADAPT need to reach 1.6 mHa, compared with UCCSD, and where does ADAPT's gradient-measurement overhead cancel its compactness?
- **Deliverables:** pool implementations (fermionic, qubit, a symmetry-adapted variant); convergence and plateau analysis; CNOT, parameter, and shot-cost tables; a comparison of ground-state overlap reached (relevant for QPE state preparation); a noisy-simulator robustness check.
- **Why not a toy:** most published ADAPT results stop at LiH/H6-class systems and report iterations rather than total measurement cost; running adjoint-gradient benchmarks up to 16 to 20 qubits with the full cost accounting is a defensible contribution, and negative results (ADAPT stalls in near-degenerate spaces) are informative.
- **Effort:** about 4 weeks part-time (largely weeks 3 and 4 plus polishing).
- **Lanthanide link:** near-degenerate 4f manifolds are exactly where single-reference-like ansatzes struggle; this tests which excitations ADAPT actually selects.

### Project 3: DMRG versus QPE crossover study
- **Question:** For the CeCl3 active-space ladder and the N2/H2O sets, at what active-space size and target accuracy does the estimated fault-tolerant QPE cost (from lambda and DF/THC Toffoli counts, including initial-state overlap) fall below the cost of a DMRG calculation reaching the same accuracy?
- **Deliverables:** DMRG energies versus bond dimension (e.g. via PySCF-interfaced DMRG software such as block2) at fixed accuracy targets; a cost model on both sides, with assumptions stated; a crossover plot and sensitivity analysis (orbital basis, localisation, overlap).
- **Why not a toy:** it takes the sceptical argument of Lee et al. (2023) and tests it on a class of systems (open-shell 4f, relativistic) they do not cover. The likely honest result may be that DMRG remains competitive throughout the ladder; that is still a publishable and employer-relevant answer.
- **Effort:** about 5 weeks part-time; depends on Project 1's cost model.
- **Lanthanide link:** whether the multiplet and ligand-field structure of Ln complexes is low-entanglement (DMRG-friendly) or not is an open, chemically meaningful question.

(Alternative to Project 3 if DMRG tooling proves painful: classical shadows or grouped measurement for 4f Hamiltonians, using the repository's Pauli-sum code to quantify variance reduction.)

## 4. What will not transfer, and what is an unfair advantage

**Will mostly not transfer:**
- AMOEBA/force-field parametrisation, AIMD sampling, and free-energy workflows: these are classical-simulation skills with no counterpart in quantum-algorithm roles. At most they show scientific range.
- Molecular generative models: mildly relevant for pharma groups interested in hybrid workflows; irrelevant to algorithm-scientist hiring.
- Broad DFT workflow experience: useful for sanity checks, not a deciding factor.
- A one-day repository is evidence of speed, not of depth: interviewers will probe qubitization and resource counts, and that is exactly what the plan targets. Competitors from physics and computer science will often be stronger on complexity theory; do not try to win there.

**Unfair advantage:**
- Relativistic electronic structure and 4f chemistry: very few quantum-algorithms researchers can say what X2C does, why spin-orbit coupling breaks the usual spin symmetries, or whether a given active space is chemically sensible. Almost all resource estimates in the literature use light-element, scalar Hamiltonians.
- Active-space judgement: lambda, factorisation rank, and state overlap all depend on orbital choice; someone who can make and defend that choice produces more credible resource numbers.
- Reference discipline: CASCI baselines at every rung of the ladder is the habit that quantum-algorithms papers often lack.
- Domain credibility with chemists at pharma/chemicals employers, who need a translator more than another circuit designer.

Caveat: this advantage only counts if it is demonstrated in public code and a write-up; the portfolio, not the CV line, is what changes the outcome.
