# progress

- 2026-10-24 上午：计划 v1（cluster 为主）。cluster 连通、`~/anaconda3/envs/qchem` 建好并跑通 H2/N2 FCI。
- 2026-10-24 下午：两轮外部评审后改为计划 v2（PLAN.md、TASK_BRIEF.md）。本机 conda 环境 `qchem` 建好：pyscf 2.14.0（conda-forge win-64）、openfermion 1.8.1、pennylane 0.45.1、pennylane-lightning。
  - 踩坑：conda 的 numpy/MKL 直接调 python.exe 时所有 BLAS 调用崩（0xC06D007F），根因是 `Library\bin` 不在 PATH；必须 `conda activate` 或 `conda run -n qchem`。期间已把 BLAS 换成 OpenBLAS（保留）。
  - 本机验证：N2/STO-3G CASCI(6,6) 积分折叠进 OpenFermion（`0.5*eri.transpose(0,2,3,1)` + `spinorb_from_spatial` + ecore），12 qubit 383 个 Pauli 项，稀疏 ED 与 CASCI 差 3e-13 Ha，1.3 s；lightning.qubit 可用。
  - 验证脚本：scratchpad 的 local_test.py（S1 时搬进 tests/）。
- 下一步：S0（environment.yml、sync.ps1、git init），然后 S1。

## 2026-10-24 S0 / S1 / S3 执行日志

- **S0 完成**：environment.yml、scripts/sync.sh(+ps1) tar-over-ssh 推拉、scripts/qpy.sh（本机 conda 环境包装，解决 PATH 问题）、git 初始化。cluster 往返验证通过。踩坑：sync 里 `~`/`$HOME` 被本地展开；cluster 的 ~/.bashrc 在 cd 时打印目录列表污染 tar 流（改用 `tar -C`）。
- **S1 完成**（验收全部满足）：
  - 积分折叠 ED == CASCI：H2 7e-16、LiH(4,4) 4e-15、N2 CAS(6,6) 7e-14、拉伸 N2 通过。
  - Koridon 1-范数恒等式与显式 JW 一致到 1e-13（小分子）；镧系离子 ~3e-6 绝对偏差（相对 1e-7，compress 截断引起），仍标为 exact。
  - JW/BK 同谱；Z2 tapering 保能量（H2 4→1 qubit，LiH 8→5，N2 12→8，Ln 14→12）。
  - UCCSD-VQE（lightning+adjoint+BFGS）：H2 1e-14、LiH(4,4) 1e-8 到 4e-6、N2 平衡点 4.2e-4 Ha（梯度 3e-7 已收敛 → 这是一阶 Trotter UCCSD 的本征误差，不是优化器）。N2 12 qubit UCCSD 分解后 117 参数、9024 CNOT、深度 11496、单点 42 s。
  - 测试：tests/unit 15 个（折叠、恒等式、谱、tapering、θ=0=HF、参考角回归、开壳层激发）+ tests/integration 2 个，全部通过。
  - 扫描：h2_sto3g（8 点）、lih_cas44（6 点）完成；n2_cas66（10 点）、h2o_cas86（9 点）后台进行中。
- **S3 进行中**：
  - 基组：PySCF 自带库无 x2c-SVPall（basis-set-exchange 可取），改用自带的 cc-pVDZ-DK（Ce–Gd 有 DK3 版本）+ sfx2c1e。
  - AVAS 默认 minao 无 Ce → 用 ANO-RCC 作投影参考（avas_minao: ano）。
  - **HF 对镧系的典型毛病**：Gd3+ ROHF 把 2 个电子放进 5d；CeCl3 ROHF 把 SOMO 放在 Ce 5d 而非 4f（AVAS 报 nelec=-1）；Pr3+ ROHF 组态错误（CASSCF 低 78 mHa）。对策：自由离子一律 CASSCF（高自旋下近乎单行列式，便宜）；CeCl3 加 relax_reference（先在 4f+SOMO 空间 CASSCF，再 AVAS）。
  - Ce3+ 4f-only：CAS(1,7) 14 qubit，1585 Pauli 项，λ=50.68，VQE 误差 7e-11（6 参数、168 CNOT、1.6 s）。
  - Eu3+ 4f-only：Sz=3 扇区 7 个行列式，ED 与 CASCI 差 2e-12；7 个能级分裂 10–18 mHa —— 来自 ROHF 单组态核心极化破坏球对称。已加 frac_occ（球平均参考）选项，待重跑对比。
  - Pr3+ 4f-only：VQE 6 mHa 误差、梯度已收敛 → 参考行列式不是主导组态。修复：CASSCF 后按自然占据排序活性轨道。重跑中。
  - relax_reference 首版 bug（cas_natorb 返回全长占据数）已修，CeCl3 L1–L5 重跑中。
