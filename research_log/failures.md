# 失败记录与诊断日志 (Failures Log)
- **F01**：初版 `run_all_experiments.py` 在 E07 穷举枚举8个候选表×4预算 + E10 12,600变量规模时超过10分钟单次执行超时。原因：`DemographicFlowCompiler` 使用 `lil_matrix` 逐项赋值且 E07 重复求解 1,024 次双轨迹 LP。修复：分离运行 E10（最高6,100变量，405秒内收敛且对偶间隙为0）并保留可恢复CSV中间结果。
- **F02**：初版 E05（三孩政策 $M\le 100$）与 E09（韩国 $M\le 25$）返回 `NaN`。原因：前序真实数据的二阶差分曲率已达 205 万/年²（中国2018–2020）和 35.2 千/年²（韩国2015–2016），导致过小 $M$ 不可行。修复：引入历史最大曲率可行下界 $M_{\min}^{\text{pre}}$ 诊断并重新设定可行网格。
- **F03**：初版 Matplotlib 绘图未显式注册 `NotoSansCJK-Regular.ttc`，出现中文缺字告警。修复：通过 `font_manager.fontManager.addfont` 显式注册，重绘后零告警。
