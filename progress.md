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

## 2026-10-24 S4：独立评审与修订

以审稿人 + 面试官身份的独立评审（后台 agent，只读）给出"major revision / 进入技术面"的结论。采纳并修复的项目：
- **代码错误**：开壳层体系的 tapering 扇区用了闭壳层假设（改为按 (na,nb) 参考行列式计算生成元本征值；Ce3+ ROHF 轨道下由 14→12 变为 14→10）；`ao2mo.restore` 得到的积分只有 ~1e-7 的 8-fold 对称，Koridon 恒等式偏差由 3e-5 降到 1e-7（CeCl3 L2 仍 2e-7 相对，原因未明，对资源估算无影响）；Toffoli proxy 改为 queries × 项数（稀疏 LCU 的 SELECT 一元迭代 ~L）；>24 qubit 的项数 proxy 去掉 ×8；SA-CASSCF 的 e_casscf 改记最低态而非平均；结果 JSON 记录 git SHA。
- **科学结论**：撤回"共价走 5d 而非 4f"的能量排序结论（各层 CASCI 参考态不一致，relax CASSCF(1,8) 已低于 L1 3.5 mHa、与 L4 持平）；Pr3+ 二分量第 21 个能级原误标为 3H6，补算 45 个态后指认为 3F2（5 重），3H6（13 重）在 46.3 mHa；"λ 比 qubit 数增长快"撤回；N2 "单调"措辞改为实际数据；"CASCI 到 1e-11"按体系区分。
- **SOC 验证**：新增类氢 Ne9+/Ar17+ 测试，偶调和大基组下 2p1/2–2p3/2 分裂与 Dirac 解析值差 <1%，证明变换正确，Ce3+ 2.3× 高估是缺二电子屏蔽的物理。
- **一致性**：所有结果从同一 commit（67e551f）重新生成；分子曲线改用 CASSCF 轨道（RHF 轨道在拉伸区变性，N2 CASCI 曲线此前在 2.4 Å 非单调）、分子不用重启；技能差距分析移出仓库（../qc-chem-notes）。
- **新结果中需要说明的点**：N2 3.0 Å 的 CASSCF 收敛到对称破缺轨道（1778 项），VQE 500 步未收敛（0.22 mHa，梯度 5e-5），该点从"ansatz 误差"讨论中剔除。
- **未采纳/留作后续**：DF/THC λ（openfermion.resource_estimates）、各层 CASSCF 与 4f+σ 对照、二电子 SOC 屏蔽、ADAPT-VQE、误差棒与噪声、CNOT 编译优化。

## 2026-10-24 晚：参考态一致的 CeCl3 梯子（优先级 1）——CASSCF 路线的结果

- 诊断 1：CeCl3 的 4f1 ROHF 解存在，E = −10239.411130 Ha，与 relax CASSCF(1,8) 完全一致；此前记录的 E_HF（−10239.2928）是 SCF 先落入的 5d1 局部解，比 4f1 解高 118 mHa。现在 relax 步之后把 ROHF 重新收敛到 4f1 组态并以它为参考（JSON 里 relax_reference.reference 字段）。
- 诊断 2（活性空间漂移，两种）：
  - 不冻芯的 CASSCF(11,12)（4f+σ）把 5 个 Cl σ 占据轨道旋转成 Ce 4d 半芯轨道（成分 4d 0.99），能量比 ROHF 低 228 mHa，全是 4d–4f 芯价关联，不是配体效应。
  - 冻芯的 CASSCF(11,12) 守住了 σ，但把 6 个空 4f 轨道换成配体 σ*（活性空间 Ce 4f 权重 7→1.04），降低 22.8 mHa，是 σ→σ* 配体对关联，与 4f 无关；自然占据 σ 1.994–1.998、虚轨道 0.003–0.0045。
  - 一个活性电子的 CASSCF(1,n)（C1、C2）对未占据活性轨道无梯度，它们漂成 5d/6s 杂轨道；能量恒等于 ROHF。
  - 冻芯的 SA-CASSCF(1,7)：芯取自 4f1 ROHF，对其他 6 个 4f 是极化的，4f 分裂假性放大到 6–19 mHa；不冻芯的 SA-CASSCF（L1_sa，0.7–6 mHa）才是配体场分裂的正确做法。
- 结论：对 CeCl3，完全变分的 CASSCF 守不住按化学定义的活性空间（不冻芯漂向 4d，冻芯漂向 σ*），"orbital-consistent CASSCF ladder" 不可行，强行加约束只是换一种任意性。可辩护的比较是：**同一个 4f1 ROHF 参考 + AVAS 轨道 + CASCI**，每层报告 CI 自然占据数和自然轨道成分；对照 4f+σ（去掉 5d 虚轨道）同样做 CASCI。已据此重做（L1、L1_sa、L2、L3、L5 本机，L4 CASCI cluster），CASSCF 版本的结果移到 results/superseded/。
- cluster 上的 C4 CASSCF(11,17) 作业（112914）只用来取宏迭代计时，其轨道结果不采用。

## 2026-10-24 晚：二电子 SOC 屏蔽（优先级 2）
- PySCF 主包无 SOMF/AMFI；`pyscf-properties`（pyscf.prop.gtensor.uhf）有基于 int2e_p1vxp1 的 SOMF（Neese 2005 式 19）与 AMFI 驱动，直接复用（qcchem/somf.py）。
- 约定校准：按 σ/2 组装时类氢离子 2p 分裂恰为 Dirac 的 0.50，说明 pyscf 矩阵已含 1/2；改为 σ 后 Ne9+ 0.997、Ar17+ 0.989（BP 为 O(α²)）。tests/unit/test_somf.py。
- 交叉验证：Ce3+ 4f 的 Breit–Pauli 一电子分裂 5292 cm⁻¹ 与 X2C 一电子 5266 cm⁻¹ 差 0.5%（Frobenius 比 1.005），两个独立实现一致。
- 屏蔽：SOMF 二电子平均场因子 0.458（AMFI 单原子时相同）。Ce3+ 2F7/2–2F5/2：X2C-1e 5266 → X2C-1e+SOMF 2397 cm⁻¹（实验 2253，+6%）。Pr3+：3H5 2258–2264（实验 2152）、3H6 ≈4640（实验 4389）cm⁻¹。态平均密度版本（去掉 5 cm⁻¹ 的非球对称残差）重跑中。
- 模型表述：X2C 一电子 SOC + Breit–Pauli 二电子平均场（SOMF），密度取 CASCI/SA 态密度。
- 四分量 DHF/DCG 基准仍在本机跑（Gaunt 项慢），作为独立的屏蔽估计。

## 2026-10-24 晚：一致参考上的 CASCI 梯子（优先级 1 的最终做法）——结果
同一个 4f1 ROHF 参考（E = −10239.411130）、AVAS 轨道、CASCI；每层报告 CI 自然占据数与自然轨道成分：
- L1 4f-only CAS(1,7)：E = −10239.407597（比 ROHF 高 3.53 mHa，是 AVAS 投影丢掉 SOMO 弛豫的代价；一个电子无关联）。
- C3 对照 4f+σ（无 5d）CAS(11,12)：E = −10239.407691；自然占据 2.000×5 / 1.000 / 0.000×6（到 1e-4），即 Cl σ 电子完全不向 4f 虚轨道激发；相对 L1 只低 0.09 mHa。
- L3 4f+5d+σ CAS(11,17)：E = −10239.414137，比 C3 低 6.45 mHa（嵌套空间，变分可比）；5d 自然占据 0.0025/0.0025/0.0022/0.0022/0.0010（合计 0.010 e），4f 虚轨道 0.0001 each。
- 结论（限定在当前活性空间描述内）：Cl σ 对的关联/给予通道是 Ce 5d，不是 4f；去掉 5d 后配体电子与金属之间没有任何可观的关联；5d 占据 ~0.01 e，量级上是极化/给予通道而非强关联自由度。"5d 是否必须进入活性空间"——是（6.4 mHa，全部的金属–配体关联能），"4f 虚轨道是否必须"——否（0.09 mHa）。这不是一般性的共价机理结论。
- λ：L1 56.1、C3 65.2、L2 116.1、L3 144.6、L4 92.6、L5 178.9 Ha（L4/L5 由积分恒等式）。
- L4 CASCI（CAS(19,16)）因 cluster 不通改在本机跑（内存上限 22 GB）。

## 2026-10-24 晚：三项闭合检查
1. **嵌套一致性**：原因找到——AVAS 把 ROHF 的 SOMO 重新投影到 ANO 4f（损失 3.53 mHa），各层级 CI 以不同程度"赎回"这个损失，所以"降能 0.09 / 0.02 mHa"是 SOMO 弛豫假象且基准不同。修正：AVAS openshell_option=3（SOMO 原样进入活性空间），L1 与 L2 的 CASCI 能量严格等于 ROHF（0.0000 mHa）；CI 收敛阈值 1e-12；所有层级共用一个 ROHF 基准（漂移 0 nHa）。新数字（E_CASCI − E_ROHF，mHa）：L1 0、L2 0、C3(4f+σ) −0.093、L4(4f+Cl3p) −0.267、L3(4f+5d+σ) −6.537。子空间包含检查（S 度量投影）：C3 ⊂ L3 与 L2 ⊂ L3 严格嵌套（1.0000）；C3 的三个 5d 投影 σ 轨道只有 56–69% 落在 L4 的"Cl 3p"AVAS 占据空间内（含 Cl 3s/4p 成分），所以 C3 与 L4 不嵌套，不能直接比；L3 对 L5 同理（94.9%）。有效的能量比较只在嵌套对上做。
2. **物理结论收紧**：固定轨道占据（CI 一阶约化密度对角元）：C3 4f 虚轨道 1e-5–3e-5；L3 5d 0.00096–0.00259（合计 0.0105）、4f 虚轨道 6e-5–2.7e-4；L4 4f 虚轨道 ≤1e-4。措辞改为"在当前固定轨道与所测空间内，5d 提供主要的新增关联贡献（6.4 mHa，0.01 e），4f 虚轨道的贡献 ≤0.3 mHa"，不写"全部"或"完全惰性"。
3. **λ 定义**：λ 为排除恒等项后全部 Pauli 系数绝对值之和，含 SOC 一体项。Ce3+：标量 55.209 → 一电子 SOC 55.485（+0.50%）→ SOMF 55.335（+0.23%）；Pr3+ 50.245 → 50.531 → 50.378；CeCl3 L1-SA 58.999 → 59.231 → 59.105。项数一电子与 SOMF 相同（3549/3549；CeCl3 3617/3609）。SOC 实现表述为"经类氢 Dirac 校准与 BP/X2C 交叉验证（0.5%）"，不写"彻底排除代码错误"。
- 下一步已启动：OpenFermion resource_estimates 的稀疏/SF/DF 估算（N2、Ce3+、L1、C3、L3、L5；C3/L3 带精确 CASCI 误差），λ 轨道优化（同一活性子空间内，CASCI 不变检查）。THC 需 jax，未装。

## 2026-10-24 夜：分解资源估算与 λ 轨道优化（完成）
- `scripts/resource_estimates_ft.py`：OpenFermion resource_estimates 的稀疏（Lee 2021）/SF（Berry 2019）/DF（von Burg 2021）成本，ε=1 mHa、χ=10、β=20；误差用分解积分下的精确 CASCI 能量变化（C3、L3 可解）。THC 需 jax，未装。SF/DF 的 QROM 成本函数对太小的空间或秩会 SystemExit，已捕获记为 costing_error。
- 关键数字：L3（34 自旋轨道）DF 1e-4：λ 100.3、L=123、Ξ=1597、1.0e9 Toffoli、867 逻辑 qubit、CASCI 误差 −0.031 mHa（1e-3 时误差 1.86 mHa 不够）；L5（42）DF 1e-4：λ 104.4、1.3e9 Toffoli、1038 qubit；稀疏法在这个尺寸下 Toffoli 反而相当或更低（d 只有 4e3–1e4），但需 ~1800 个 QROM 工作 qubit。此前的"查询数×项数"proxy（L5 5.5e10）撤回，高估约两个量级。
- `scripts/lambda_orbital_opt.py`（Koridon 思路，固定活性子空间内正交旋转，平滑 |·| + L-BFGS，3 次起点）：N2 0.940、Ce3+ 0.958、L1 0.997、C3 0.951、L3 0.852、L5 0.743；CASCI 能量不变（≤2e-12，L5 因 CI 不可解只能靠幺正性）。项数可能上升（C3 13776→15830），λ 与项数要分开报。L3 的优化 56 分钟主要花在两次 CASCI 检查上。
- 报告第三版：lanthanide_feasibility.md §3（嵌套、固定轨道占据、措辞）、§4（λ 定义含 SOC）、§5（分解估算）、§6（λ 优化）；technical_report.md §5–6；README。
- 仍在跑：四分量 DHF/DCG 基准（本机，已 >6 h）。集群连不通后不再探测。
