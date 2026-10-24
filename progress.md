# progress

- 2026-10-24 上午：计划 v1（cluster 为主）。cluster 连通、`~/anaconda3/envs/qchem` 建好并跑通 H2/N2 FCI。
- 2026-10-24 下午：两轮外部评审后改为计划 v2（PLAN.md、TASK_BRIEF.md）。本机 conda 环境 `qchem` 建好：pyscf 2.14.0（conda-forge win-64）、openfermion 1.8.1、pennylane 0.45.1、pennylane-lightning。
  - 踩坑：conda 的 numpy/MKL 直接调 python.exe 时所有 BLAS 调用崩（0xC06D007F），根因是 `Library\bin` 不在 PATH；必须 `conda activate` 或 `conda run -n qchem`。期间已把 BLAS 换成 OpenBLAS（保留）。
  - 本机验证：N2/STO-3G CASCI(6,6) 积分折叠进 OpenFermion（`0.5*eri.transpose(0,2,3,1)` + `spinorb_from_spatial` + ecore），12 qubit 383 个 Pauli 项，稀疏 ED 与 CASCI 差 3e-13 Ha，1.3 s；lightning.qubit 可用。
  - 验证脚本：scratchpad 的 local_test.py（S1 时搬进 tests/）。
- 下一步：S0（environment.yml、sync.ps1、git init），然后 S1。
