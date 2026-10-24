# From Active-Space Electronic Structure to Quantum Algorithms: Benchmarks and Resource Estimates for Strongly Correlated and Lanthanide Systems

计划 v2，2026-10-24。v1 经两轮外部评审后修订；修订理由见文末"v1 → v2 变更"。

## 0. 核心主张（一句话）

**4f-only 的镧系活性空间对经典方法太容易，不构成量子优势场景；只有加入配体轨道并考虑相对论/自旋轨道耦合后问题才真正变难，而此时 NISQ-VQE 已不可行，只剩容错 QPE 值得讨论。本项目用自己的计算给出这条边界上的定量数字。**

结论即使是"目前量子计算对真实镧系化学没有优势"也算成功，只要数字是自己算出来的。

## 1. 算力与原则

| 资源 | 现状（10/24 已验证） | 用途 |
|---|---|---|
| 本机 Windows，Ryzen 9 5900 / 32 GB / 无 CUDA | conda 环境 `qchem`：pyscf 2.14.0、openfermion 1.8.1、pennylane 0.45.1、pennylane-lightning、BLAS 已换成 OpenBLAS。**必须 `conda activate qchem` 或 `conda run -n qchem` 运行**，直接调 python.exe 会因 Library\bin 不在 PATH 而崩溃（0xC06D007F） | **主开发环境**。所有 ≤ 14 qubit 的计算、全部 pytest、作图、报告 |
| cluster（ssh 别名 cluster，Xeon Gold 6326 登录节点 32 核 251 GB；slurm 节点 32 核 120 GB；无 GPU） | 环境 `~/anaconda3/envs/qchem` 同版本 | **只做本机做不动的经典参考**：真实镧系配合物 X2C-CASSCF、大活性空间 CASCI、必要时 block2 DMRG |

cluster 使用原则（硬约束）：
1. 集群访问不保证长期可用。代码与结果以本机为主副本，集群只有镜像与临时产物。
2. 登录节点只做 smoke test：≤ 4 线程、< 4 分钟。任何 benchmark 一律 `sbatch`，1 节点、4 到 8 核、`--time` ≤ 60 分钟，同时最多排 1 个作业，不用 tmux 挂通宵。
3. 不写 `/home`，scratch 用 `/data/$USER/qc-chem/`，用完即删。
4. 每个 sbatch 作业结束立刻把积分 `.npz`、结果 `.json/.csv` 拉回本机，写 `progress.md`。
5. 额度紧时先停昂贵的验证步骤。

## 2. 目录布局

```
本机  C:\Users\shich\project\qc-chem\        # git 主副本
      ├── PLAN.md  TASK_BRIEF.md  progress.md  environment.yml
      ├── qcchem/        # integrals, active_space, mapping, ansatz, vqe, adapt, benchmark, resources
      ├── configs/       # 每个体系一个 yaml
      ├── scripts/       # run_system.py, make_plots.py, sync.ps1
      ├── slurm/         # sbatch 模板
      ├── results/       # 每体系一个 json；积分 npz（小于 5 MB 的进 git）；csv；png
      ├── tests/         # unit（固定参数、确定性）与 integration（真跑优化器）分开
      └── report/        # technical_report.md, lanthanide_feasibility.md
cluster  ~/qc-chem/            # 代码镜像
         /data/$USER/qc-chem/   # scratch
```

## 3. 工程硬约束（写进 TASK_BRIEF，Claude Code 必须遵守）

| 项 | 上限 / 规定 | 理由 |
|---|---|---|
| 模拟器 | PennyLane `lightning.qubit` + `diff_method="adjoint"`，优化器 `scipy.optimize.minimize(method="BFGS", jac=True)`。> 6 qubit 禁用 parameter-shift | 比 parameter-shift 快一到两个量级 |
| VQE / ADAPT-VQE 优化 | ≤ 12 qubit；14 qubit 仅限小激发空间（Ce³⁺ 4f¹、Pr³⁺ 4f²） | CPU 上 16 到 20 qubit 的 UCCSD 一次优化要小时到天 |
| qubit 空间精确对角化 | 稀疏 `eigsh`，≤ 16 qubit。**从不构造稠密矩阵** | 24 qubit 稠密矩阵 4 PiB；经典参考本来就有 PySCF CASCI |
| 哈密顿量构造与统计 | Jordan-Wigner 显式构造 ≤ 24 qubit；更大的活性空间**直接从积分算 λ 和项数**（不展开 Pauli 串），外推 O(N⁴) | 40 qubit 的 JW 展开要数十 GB 内存与数小时 |
| 噪声 | 只做 shot noise（≤ 12 qubit）和去极化（`default.mixed`，≤ 6 qubit）。固定无噪声最优参数 θ*，扫 p ∈ [1e-4, 1e-1]，**不在噪声下重跑优化** | 密度矩阵维度 4^N；这不是本项目主线 |
| 积分折叠 | `mcscf.CASCI.get_h1eff()`、`ao2mo.restore(1, get_h2eff(), ncas)`，二体张量 `0.5 * eri.transpose(0,2,3,1)` 再 `spinorb_from_spatial`，`ecore` 显式带上。**验收：OpenFermion 稀疏 ED 与 PySCF CASCI `e_tot` 差 < 1e-8 Ha**（10/24 本机实测 N₂ CAS(6,6)：3e-13） | 化学家记号 (pq\|rs) 转物理学家记号是最常见的能量对不上的 bug |
| 测试 | **unit**：固定几何 + 固定参考参数 + 确定性 statevector，\|E − E_ref\| < 1e-8。**integration**：真跑优化器，\|E_VQE − E_FCI\| < 1e-5 Ha，另报告是否达到化学精度 1.6e-3 Ha | 优化器偶然失败不应让 CI 变红 |
| 措辞 | 写 "evaluate the degradation of UCCSD-VQE in the stretched-bond regime"，不预设"失效" | 不为叙事强行制造结果 |

## 4. 体系与分层

| 层 | 体系 | 活性空间 | qubit | 量子侧 | 经典参考 | 跑在哪 |
|---|---|---|---|---|---|---|
| 验证 | H₂ STO-3G | (2,2) | 4 | UCCSD-VQE、JW/BK 对拍、tapering | FCI | 本机 |
| 验证 | LiH STO-3G | (2,2)、(4,4) 冻芯 | 4、8 | 同上 | FCI within CAS | 本机 |
| **核心** | N₂ 解离曲线 | CAS(6,6) | 12 | UCCSD vs 费米子池 ADAPT-VQE | HF、CASCI = FCI within CAS(6,6) | 本机 |
| 核心 | H₂O 对称拉伸 | CAS(6,5) 或 (8,6) | 10、12 | 同上 | 同上 | 本机 |
| 镧系基线 | Ce³⁺ 4f¹、Pr³⁺ 4f² | CAS(1,7)、(2,7) | 14 | VQE 实跑（激发空间很小） | CASCI；明确 N 与 S_z 扇区 | 本机 |
| 镧系基线 | Eu³⁺ 4f⁶、Gd³⁺ 4f⁷ | CAS(6,7)、(7,7) | 14 | 只做 Hamiltonian 统计与稀疏 ED，**不做 VQE** | CASCI，固定 S_z，挑选低能态 | 本机 |
| **核心卖点** | 真实小配合物（首选 LnCl₃，备选 LnF₃ / [LnCl₆]³⁻；Ln = Ce、Eu），标量 X2C（sfx2c1e），x2c-SVPall 基组 | 4f-only (n,7) → +5d (n,12) → +配体 σ/π (n,15+) | 14 → 24 → 30+ | ≤ 24 构造 JW 并统计；> 24 从积分算 λ | CASSCF / CASCI，必要时 DMRG | **cluster sbatch** |
| 扩展 | 同一配合物，二分量 X2C + SOC（复数 Kramers 旋量） | 同上 | 同上 | 只统计：非零积分数、Z₂ 对称性破缺（只剩粒子数奇偶）、Pauli 项数膨胀倍数 | 二分量 CASCI 若可行 | cluster sbatch |
| 可选 | 一个 3d 过渡金属小 CAS | ≤ (12,12) | ≤ 24 | 统计 | CASCI | 若时间有余 |

诚实声明（写进报告）：自由离子 4f-only 活性空间在固定 S_z 扇区里行列式数极少（Ce³⁺ 14 个，Eu³⁺ 高自旋 S_z = 3 扇区 7 个），Hund 规则基本给出答案；这一层的价值是验证映射与 tapering 在简并 4f 轨道下的行为，并给出哈密顿量统计基线，而非展示强关联。非相对论、无 SOC 的模型不重现 Eu(III) 的 ⁷F_J 多重态结构。

每个体系统一输出 JSON：活性电子数、活性轨道数、自旋轨道数、tapering 前后 qubit 数、Pauli 项数、λ = Σ|h_j|、VQE 参数量、电路深度、双比特门数、优化器能量评估次数、相对 CASCI/FCI 的能量误差、耗时、运行机器与版本。

## 5. 主图

**图 1（核心）**：横轴 H₂ → N₂ → Ce³⁺ → Eu³⁺ → LnCl₃ (4f) → LnCl₃ (4f+5d) → LnCl₃ (4f+配体) → 二分量 LnCl₃。纵轴（多面板，对数）：qubit 数、Pauli 项数、λ、ADAPT 算符数或 UCCSD 参数量、估算 QPE 代价。用底色标三个区域：classically trivial / NISQ-simulatable benchmark / fault-tolerant target。

**图 2**：N₂ 与 H₂O 的 E − E_FCI vs 键长，方法 HF / UCCSD-VQE / ADAPT-VQE，附化学精度线。

**图 3**：标度图。Pauli 项数与 λ 随活性轨道数 N_orb 的曲线，拟合指数，外推到 40 到 60 qubit。

**图 4**：ADAPT vs UCCSD 的参数量、CNOT 数、深度、能量评估次数对比。

**图 5（若做）**：固定 θ* 下能量误差随去极化率 p 的退化曲线（H₂、LiH 6 qubit）。

## 6. 资源估算的口径

对每个体系给：N_qubits、N_Pauli、λ，目标精度 ε = 1 mHa，按 qubitization 的 proxy 公式给 Toffoli 数量级（∝ λ/ε × 每步代价）与逻辑 qubit 数。报告明确写：**这是基于 Pauli 1-范数的 proxy 估算，不是经过 double factorization / THC 优化的容错资源计数；逻辑与物理 qubit、surface code 假设留作未来工作。**

## 7. 会话安排

**S0（手动，≈20 分钟，不烧额度）**
- 本机环境已建好并验证。补 `environment.yml`（导出 pin 版本）、`scripts/sync.ps1`（push 代码 / pull results），测一次往返。
- `git init`，首次 commit。

**S1 — 必须完成：pipeline 成立**
H₂ → LiH → N₂。PySCF → active space → 积分折叠（验收 < 1e-8）→ JW/BK → Z₂ tapering → 稀疏 ED → UCCSD-VQE（lightning + adjoint）。unit 与 integration 测试分开。N₂ 曲线 8 到 10 点全在本机跑。产出 results/*.json、第一版资源表、progress.md。

**S2 — 项目核心：UCCSD vs ADAPT-VQE**
N₂ + H₂O。费米子池 ADAPT（qubit 池可选）。输出能量误差、参数量、能量评估次数、Pauli 项数、深度、双比特门、耗时。图 2、图 4。噪声（图 5）放最后，额度紧第一个砍。

**S3 — 化学特色：镧系**
- 本机：Ce³⁺、Pr³⁺ 14 qubit VQE；Eu³⁺、Gd³⁺ 统计与 ED。
- cluster sbatch：LnCl₃ 标量 X2C-CASSCF，活性空间三档递增，积分 npz 拉回本机后在本机做 JW 统计与 λ。
- 图 1、图 3 初版。可行性报告初稿：今天可演示 / 近期硬件 / 仅容错 三档。
- CASSCF 不收敛时如实记录轨道选择与尝试，允许报告"未收敛"。

**S4 — 研究级扩展与审查**
- 二分量 X2C + SOC 旋量哈密顿量统计（只要报告与一个独立小脚本，不要求完整代码）。
- QPE proxy 资源分析与三个候选课题（按可行性、新颖性、代价排序）。
- 让 Fable 以量子算法岗面试官 + 审稿人身份审 repo；后台 agent 并行做技能差距分析。
- 最后一次 pull，确认本机为完整副本，清 cluster scratch。

## 8. 风险

- cluster 中断：只影响 S3 的真实配合物参考；S1、S2 完全在本机。
- 4f 轨道在 CASSCF 迭代中旋出活性空间：用 AVAS 或按 4f 投影选初猜，失败则退到 CASCI。
- 二分量 PySCF API（spinor 积分、复数 InteractionOperator）调试可能吞额度：限定为 S4，设时间盒，超时只写报告。
- 项目的 H₂/N₂/ADAPT 部分是成熟结果，科研新意全在镧系层；若 S3 没有干净数字，作品集价值仍在，发表价值不在。

## v1 → v2 变更

1. 主开发环境从 cluster 改为本机（conda-forge 已有 Windows 版 PySCF，10/24 实测可用）；cluster 只做真实镧系配合物的经典参考。cluster 占用因此大幅下降。
2. 删除"24 qubit 精确对角化"；稀疏 ED ≤ 16 qubit，更大的用 PySCF CASCI。
3. VQE 上限从 20 降到 12（14 仅限小激发空间）；强制 lightning + adjoint。
4. JW 显式构造 ≤ 24 qubit；更大活性空间从积分直接算 λ。
5. 噪声降为可选、≤ 6 qubit、固定 θ*。
6. Eu³⁺ 从"VQE 复现多重态"改为"统计与 ED，明确 S_z 扇区"，VQE 改用 Ce³⁺/Pr³⁺；加诚实声明。
7. 镧系主线改为同一真实配合物上活性空间逐级增大，形成标度图；标量 X2C 与二分量 SOC 分开，后者为 S4 扩展。
8. 测试拆 unit / integration；N₂ 措辞改为"评估退化"。
9. 过渡金属降为可选。
10. 登录节点只做 smoke test，benchmark 全走 sbatch。
11. 标题升级，核心主张收窄为一句。
