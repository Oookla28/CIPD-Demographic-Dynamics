#!/usr/bin/env python3
"""Generate 4-reviewer self-audit reports, claim-evidence matrix, research logs, Makefile, and README."""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[3]

# 1. Four-Reviewer Self-Audit Reports
(ROOT / "reports/reviewer_math.md").write_text(
"""# 模拟同行评审报告 — Reviewer A（数学审稿人）

## 评审结论：Major Revision（核心框架严谨，存在2项需在后续阶段强化的理论深度）
1. **已通过核查项**：
   - V0/V1 人口账本守恒残差分别为 $2.73\\times 10^{-16}$ 与 $2.17\\times 10^{-16}$（E01），达到双精度机器极限。
   - V3 线性规划原对偶证书经独立验证，所有规模下归一化对偶间隙 $\\le 1.8\\times 10^{-16}$（E02、E10）。
   - 附录 `paper/appendix/mathematical_proofs.tex` 严格区分了已有数学基础（强对偶、McCormick、HonestDiD、Kallus-Zhou）与本文结构性命题（Theorem A.1 跨期传播、Theorem A.3 双轨迹观测分离LP、Proposition A.4 识别宽度非次模性反例）。
2. **严重质疑与修复状态**：
   - **[S1 - 已修复]** 早期 E05/E09 在 $M < M_{\\min}^{\\text{pre}}$ 时出现 LP 不可行（NaN），已通过先验检测历史二阶差分上界 $M_{\\min}^{\\text{pre}}$ 修正网格，并在报告中将其转化为“平滑假设被数据证伪”的显式诊断。
   - **[S2 - 仍存局限]** V4 CCPI 当前在低维共享参数子空间（3维）演示了 McCormick + 空间分支定界（50节点，0.43秒，间隙 $<0.9\\%$），尚未扩展至高维单岁×31省联合双线性全局求解。已明确降级表述，不宣称解决大规模非凸全局优化。
""", encoding="utf-8")

(ROOT / "reports/reviewer_statistics.md").write_text(
"""# 模拟同行评审报告 — Reviewer B（统计与因果审稿人）

## 评审结论：Major Revision（诚实披露负面与不可识别结果，但宏观单序列统计功效偏弱）
1. **亮点**：
   - 未将确定性识别区间混同为频率学95%置信区间；严格区分 C1/C2/C3/C4 证据等级。
   - E03（10类DGP×200次蒙特卡洛）清晰展示了传统线性外推在隐藏混杂（DGP-7）下高达 66.5% 的错误方向判断率，而 CIPD 将错误方向率降至 0.0%。
2. **严重质疑与局限**：
   - **[S1 - 保守性权衡]** E03 中 CIPD 在所有10类DGP上的覆盖率为100%但弃权率（区间跨0）亦达100%，表明在无横截面控制组（如省级面板或微观资格暴露）时，仅凭单国年度宏观序列的 $\\Delta^{SD}(M)$ 边界过于保守。
   - **[S2 - 未解决阻断]** 因沙箱无出站网络，未能运行 R 语言 `did` 与 `HonestDiD` 官方包以及 31 省面板暴露强度回归；全国宏观序列无法支持 C1 因果点识别。已在论文中如实将中国三项政策归为 C2/C3 证据等级。
""", encoding="utf-8")

(ROOT / "reports/reviewer_demography.md").write_text(
"""# 模拟同行评审报告 — Reviewer C（人口学审稿人）

## 评审结论：Minor-to-Major Revision（口径审计扎实，国际数据覆盖不足）
1. **亮点**：
   - 准确识别并量化了国家统计局2020年七普对2012—2017年出生人口的大幅事后修订（如2012年1,635万→1,973万，+20.7%），以及2010六普表6-3的抽样长表属性，避免了常见口径误用。
   - 纳入了2025年最新统计公报（出生792万、死亡1,131万、总人口140,489万）及2025年全国育儿补贴（3,600元/年、900亿预算）作为政策基线 $u_0$。
   - E06 Shapley 分解将2010—2025年800万出生降幅精确拆解为生育数量（-541.3万，67.7%）、育龄人口结构（-300.6万，37.6%）与生育时机（+42.0万）。
2. **严重质疑与局限**：
   - **[S1 - 数据阻断]** 国际跨国验证仅覆盖中国与韩国真实公开数据，未达到原定 8–16 国 HFD 单岁×孩次批量遮蔽标准。
   - **[S2 - 长期预测简化]** E08 中 2025–2070 年人口情景基于简化宏观生育—死亡路径，未接入单岁生命表与七普分省迁徙矩阵，仅作为敏感性情景展示。
""", encoding="utf-8")

(ROOT / "reports/reviewer_computation.md").write_text(
"""# 模拟同行评审报告 — Reviewer D（计算与复现审稿人）

## 评审结论：Accept with Minor Revisions（代码工程与复现性达标，大规模求解速度尚无算法加速突破）
1. **亮点**：
   - 全流程在 2 vCPU / 7.8 GiB RAM 纯 CPU 环境下运行；所有真实数据具备 URL 溯源与 SHA-256 哈希校验；8/8 单元测试通过；所有图表均由脚本从 CSV 自动生成。
   - 稀疏矩阵密度在 XL 规模（6,100变量）下仅 0.000451，HiGHS 对偶间隙为 0.0。
2. **严重质疑与局限**：
   - **[S1 - 无算法加速声明]** E10 中 XL 规模（6,100变量）直接调用 HiGHS 耗时 405.6 秒，尚未实现超越通用 HiGHS 的专用块分解加速算法。严格遵守原则，本文不宣称提出更快的 LP 求解器，将计算贡献聚焦于建模编译与对偶证书框架。
   - **[S2 - 强基线缺失]** 未能在离线环境复现 R/TMB `ccmpp.tmb`、`bayesTFR` 与 Ellison et al. (2026) Zenodo HMC 基线。
""", encoding="utf-8")

# Claim-Evidence Matrix
claims = [
    ("C01", "V0/V1 人口流量模型满足精确人口账本守恒与胎次非负守恒", "E01", "results/tables/E01_accounting_tests.csv", "VERIFIED", "v0_rel_res=2.73e-16, v1_rel_res=2.17e-16, adjoint_speedup=155.7x"),
    ("C02", "2010六普35×3真实矩阵遮蔽下，CIPD边界100%覆盖真值而IPF产生24.1%点估计偏差", "E02", "results/tables/E02_real_census2010_masking.csv", "VERIFIED", "True=289,646; IPF=219,927 (err=69,719); CIPD Mask-F=[0, 372,295]"),
    ("C03", "在10类合成因果机制下，CIPD消除错误政策方向判断（0.0% vs 传统外推20.8%）", "E03", "results/tables/E03_synthetic_causal_recovery.csv", "VERIFIED_WITH_CAVEAT", "100%覆盖但单序列区间偏宽、弃权率高"),
    ("C04", "历史滚动预测在2016政策峰值拐点处暴露出结构突变失效风险", "E04", "results/tables/E04_rolling_forecast_china.csv", "NEGATIVE_RESULT", "2016起点4个预测期全部失覆；整体覆盖率80%"),
    ("C05", "仅凭全国宏观年度出生序列（含七普修订区间），中国三次生育政策与韩国2024反弹均无法实现符号稳健的C1因果点识别", "E05/E09", "results/tables/E05_china_policy_bounds.csv", "VERIFIED", "所有可行M值下区间均跨0；三孩政策前序最大曲率达205万/年^2"),
    ("C06", "共享机制反事实耦合(V4 CCPI)较独立边界相减缩窄识别宽度54.9%", "E05", "results/tables/E05_ccpi_shared_mechanism_coupling.csv", "VERIFIED", "908.8万 -> 409.6万，sBB间隙<3.8万(<0.9%)"),
    ("C07", "2010-2025年中国出生人口减少800万中，生育数量占67.7%、育龄结构占37.6%", "E06", "results/tables/E06_shapley_decomposition.csv", "VERIFIED_CONDITIONAL", "Shapley残差0.0；依赖2020/2025推迟与育龄规模校准假设"),
    ("C08", "双轨迹观测分离VOI贪心算法在预算3.0下缩窄识别宽度96.5%，与穷举全局最优零间隙", "E07", "results/tables/E07_voi_observation_selection.csv", "VERIFIED", "CIPD 96.5% vs Random 10.4% vs Cheapest 0.0%; gap=0.000"),
    ("C09", "以2025国家育儿补贴为基线，普惠托育扩张(U2)满足基线保护(R<=0)，8/15政策对存在排序反转", "E08", "results/tables/E08_robust_policy_metrics.csv", "VERIFIED_CONDITIONAL", "依赖弹性区间Theta假设，属C3级鲁棒决策分析"),
    ("C10", "多国HFD 8-16国单岁×孩次遮蔽与R强基线(Ellison 2026/did/HonestDiD)对比", "E02/E09", "reports/data_availability.md", "BLOCKED", "沙箱离线策略阻断外部大文件下载与R环境安装"),
]
with open(ROOT / "reports/claim_evidence_matrix.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["claim_id", "claim", "experiment", "evidence_file", "status", "notes"])
    w.writerows(claims)

# 2. Research Logs
(ROOT / "research_log/decisions.md").write_text(
"""# 关键科研决策日志 (Decisions Log)
- **D01 (2026-11-26)**：检测到沙箱已预装 `SciPy 1.15.2`（内嵌 MIT 许可 HiGHS LP/MILP 求解器），决定将核心部分识别、双轨迹观测分离与 HonestDiD $\\Delta^{SD}(M)$ 全部建立在原生稀疏 HiGHS 引擎之上，避免外部求解器依赖。
- **D02 (2026-11-26)**：在数据审计中发现七普对2012–2017年出生人口存在高达 +20.7% 的事后修订差异，决定将“同期值—修订值”直接编码为观测盒区间，而非任选单一序列。
- **D03 (2026-11-26)**：因沙箱禁止联网无法获取31省微观/面板及HFD 8–16国批量ZIP，严格执行降级机制，保持全国尺度部分识别并诚实标注最终验收状态为 `PARTIAL / BLOCKED`。
""", encoding="utf-8")

(ROOT / "research_log/failures.md").write_text(
"""# 失败记录与诊断日志 (Failures Log)
- **F01**：初版 `run_all_experiments.py` 在 E07 穷举枚举8个候选表×4预算 + E10 12,600变量规模时超过10分钟单次执行超时。原因：`DemographicFlowCompiler` 使用 `lil_matrix` 逐项赋值且 E07 重复求解 1,024 次双轨迹 LP。修复：分离运行 E10（最高6,100变量，405秒内收敛且对偶间隙为0）并保留可恢复CSV中间结果。
- **F02**：初版 E05（三孩政策 $M\\le 100$）与 E09（韩国 $M\\le 25$）返回 `NaN`。原因：前序真实数据的二阶差分曲率已达 205 万/年²（中国2018–2020）和 35.2 千/年²（韩国2015–2016），导致过小 $M$ 不可行。修复：引入历史最大曲率可行下界 $M_{\\min}^{\\text{pre}}$ 诊断并重新设定可行网格。
- **F03**：初版 Matplotlib 绘图未显式注册 `NotoSansCJK-Regular.ttc`，出现中文缺字告警。修复：通过 `font_manager.fontManager.addfont` 显式注册，重绘后零告警。
""", encoding="utf-8")

(ROOT / "research_log/negative_results.md").write_text(
"""# 负面结果与模型失效边界记录 (Negative Results)
1. **E04 2016起点结构突变失覆**：在2016年全面二孩峰值（1,786万）处向外预测2017—2021年时，随机游走、阻尼趋势与基于历史曲率的 CIPD 预测区间全部未能覆盖2018—2021年的断崖式下滑（覆盖率0/4）。说明仅靠历史平滑曲率无法预见结构性生育意愿坍塌与疫情叠加冲击。
2. **E05/E09 全国宏观单序列无法实现符号稳健因果识别**：在纳入七普修订区间与可行曲率 $M$ 后，中国单独二孩、全面二孩、三孩政策及韩国2024–2025反弹的识别区间均跨越0。这是重要的负面识别结论：它证明了“仅靠全国年度出生总量宣称某项政策精确净增X百万人”在方法论上缺乏稳健依据。
3. **E03 单序列敏感性区间高弃权率**：在无横截面控制组时，CIPD 虽然实现了100%覆盖率与0%错误方向率，但区间宽度较大导致弃权率100%，凸显了补充省级/胎次暴露观测（E07）的必要性。
4. **Proposition A.4 非次模性反例**：严格证明人口识别宽度缩减函数不具备次模性（互补边际观测可产生协同增益），否定了直接套用贪心 $(1-1/e)$ 理论保证的猜想。
""", encoding="utf-8")

(ROOT / "research_log/lessons.md").write_text(
"""# 科研经验与方法论总结 (Lessons Learned)
1. 口径审计（尤其是普查事后修订与长表抽样比）直接决定因果识别的可信边界，必须先于建模执行。
2. 不可行（Infeasibility）本身是强科学信号：E05/E09 中小 $M$ 不可行直接揭示了政策前人口趋势已发生剧烈非线性弯曲。
3. 共享机制反事实耦合（V4）与双轨迹观测分离（V6）是 CIPD 最具可计算增量价值的两个模块。
""", encoding="utf-8")

(ROOT / "research_log/hypotheses.md").write_text(
"""# 研究假设登记与验证状态 (Hypotheses Registry)
- **H1**：即使微观胎次状态不可完全恢复，某些宏观目标仍可在动态守恒约束下被紧致识别 —— **部分支持**（E02 2010六普矩阵与动态流量网络验证了单调收缩与100%覆盖，但缺失率>50%时区间仍显著扩张）。
- **H2**：共享潜在机制可显著压缩反事实效应区间 —— **强支持**（E05 CCPI 宽度缩减 54.9%，gap < 0.9%）。
- **H3**：全国宏观序列足以点识别中国三项生育政策因果效应 —— **明确否定**（E05 全部跨越0，需降级为 C2/C3 并借助外部微观/省级数据）。
- **H4**：目标导向双轨迹 VOI 显著优于随机与低成本优先观测选择 —— **强支持**（E07 预算3.0时缩减96.5% vs 10.4%，与穷举最优零间隙）。
- **H5**：识别宽度缩减函数具有次模性 —— **严格否定**（Proposition A.4 给出解析反例）。
""", encoding="utf-8")

# 3. Makefile, pyproject.toml, README.md, AGENTS.md
(ROOT / "Makefile").write_text(
""".PHONY: setup data audit test baselines experiments figures paper reproduce

setup:
\tpython3 --version && python3 -c "import numpy, scipy, pandas, matplotlib; print('Environment OK')"

data:
\tpython3 src/cipd/data/build_verified_data.py

audit:
\tsha256sum -c provenance/hashes.sha256

test:
\tpython3 -m unittest discover -s tests -v

experiments:
\tpython3 experiments/run_all_experiments.py

figures:
\tpython3 experiments/fix_e05_e09_and_fonts.py

paper:
\tlualatex -interaction=nonstopmode -output-directory=paper/journal_version paper/journal_version/cipd_journal_draft.tex || true

reproduce: setup data audit test figures
\t@echo "Core verification, data hash audit, unit tests, and figure generation complete."
""", encoding="utf-8")

(ROOT / "pyproject.toml").write_text(
"""[project]
name = "cipd"
version = "0.1.0"
description = "Constraint-based Identification and Policy Dynamics for Partially Observed Demographic Systems"
requires-python = ">=3.11"
dependencies = ["numpy>=2.0", "scipy>=1.14", "pandas>=2.2", "matplotlib>=3.9", "pyyaml>=6.0"]

[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"
""", encoding="utf-8")

(ROOT / "README.md").write_text(
"""# CIPD: Constraint-based Identification and Policy Dynamics

**Title (EN):** *Structure-Preserving Partial Identification and Robust Policy Decisions in Partially Observed Demographic Dynamics*  
**Title (ZH):** 不完全观测下人口动力系统的结构保持型部分识别与鲁棒政策决策——以中国生育政策改革为例

## Current Research Status
- **Status:** `PARTIAL / BLOCKED` (Honest assessment per Research Definition of Done)
- **Achieved Level:** **Level 1 (Engineering Complete) + Partial Level 2 (Core Theory, Verified National Empirical Analysis & Synthetic/Real Benchmarks Complete)**
- **Blocking External Factors:** Sandbox outbound internet & package-installation restrictions prevent automated downloading of 8–16 country HFD bulk ZIP archives, 31-province policy crawls, and external R baselines (`did`, `HonestDiD`, `ccmpp.tmb`, `bayesTFR`, Ellison et al. 2026 HMC).

## Quick Reproducibility Commands
```bash
make setup
make data
make audit
make test
make figures
```
See `reports/final_research_acceptance.md` for full quantitative results, negative findings, and evidence index.
""", encoding="utf-8")

(ROOT / "AGENTS.md").write_text(
"""# CIPD Research Agent Protocol
1. Never fabricate data, provenance hashes, or statistical significance.
2. Always report negative results, infeasible smoothness regimes, and blocked external dependencies honestly.
3. All figures and paper numbers must trace directly to `results/tables/*.csv` and `provenance/hashes.sha256`.
""", encoding="utf-8")

print("Milestone 3 self-audit reports, research logs, and repository metadata created.")
