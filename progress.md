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
- **cluster 作业**（唯一一次）：cecl3_L4_4f_cl3p_ci，CAS(19,16) CASCI，9.2e7 行列式，node4 8 核 6 分钟，E = −10239.411113 Ha。第一次提交因 cluster 环境缺 pyyaml 失败（已 pip 装上）。结果已拉回本机，/data/$USER/qc-chem 已清空，队列为空。
  - 能量梯子（CASCI，Ha）：L1 4f(1,7) −10239.40760；L3 4f+5d+σ(11,17) −10239.41414；L4 4f+Cl3p(19,16) −10239.41111。**L4 比 L3 高 3.0 mHa**：加 5 个 5d 比加 9 个 Cl 3p 更重要，Cl→Ce 给电子的通道是 5d，不是 4f。
- **S3 完成（本机）**：
  - 自由离子改用**态平均 CASSCF**（Ce/Eu 7 态、Pr 21 态）后轨道球对称，ED 复现精确多重态简并：Ce3+ 2F×7；Pr3+ 3H×11 (0)、3F×7 (27.19 mHa)、3P×3 (126.88 mHa)；Eu3+ 7F×7。frac_occ 路线无效（CASSCF(1,7) 中未占据活性轨道能量无关、被任意旋转），已删除。
  - VQE 重启策略：态平均轨道无主导行列式时 BFGS 会落到激发态（CeCl3 L1-SA 卡在 0.68 mHa = 第一激发态），从随机角度重启后 1.9e-10。记录在 JSON 的 n_restarts_used。
  - CeCl3 梯子 L1–L5 全部完成（L4 的 CASCI 在 cluster）。L2 阈值 0.5 才得到干净的 4f×7+5d×5。
  - 二分量扩展 `qcchem/soc.py`：标量活性空间 + 一电子 X2C SOC 块。Ce3+ 2F7/2–2F5/2 = 5266 cm⁻¹（实验 2253，无二电子屏蔽高估 2.3×）；Pr3+ 3H4(9)/3H5(11) 多重度正确；项数 +7%、λ +0.5%、Z2 生成元 2→1。PySCF 踩坑：自旋无关 X2C 类在 `pyscf.x2c.sfx2c1e.SpinFreeX2C`。
  - 报告：report/technical_report.md、report/lanthanide_feasibility.md、report/skill_gap_analysis.md（后台 agent 写，引文年份未核）、README.md；图 results/figures/fig1–3。
- **未做**（按用户指示延后）：S2 的 ADAPT-VQE 与噪声；双因子化/THC 资源估算；二电子 SOC 屏蔽。
- **偏离 PLAN 的地方及原因**：(1) CeCl3 L3 的 CASCI（7.7e7 行列式）在本机跑了 10 分钟而不是 cluster，因为本机能跑且更符合"少用 cluster"；(2) 主图中自由离子用态平均版本替代原 ROHF 参考版本，后者保留为对照；(3) 噪声图（fig5）和 ADAPT 图（fig4）未做。
