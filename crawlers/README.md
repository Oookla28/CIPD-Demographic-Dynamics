# crawlers/ —— 数据获取阶段（按依赖顺序编号）

每个脚本是一个可独立重跑的阶段：**抓取/OCR → 解析 → 恒等式校验 → 写入 `data/`**。
所有原始下载物与最终 CSV 的 SHA-256 记录在 `../provenance/hashes.sha256`；
来源、机构、口径与状态记录在 `../provenance/sources.yaml`（D01–D12）。

| 阶段 | 脚本 | 抓取对象 / 来源 | 输出 | 校验 |
|---|---|---|---|---|
| 01 | `01_national_series.py` | 国家统计局历年公报与年鉴；卫健委孩次占比；民政部婚姻登记 | `data/raw/china_annual_vital_1990_2025.csv`、`china_parity_marriage_2012_2025.csv`、`china_census_age_anchors_2000_2025.csv` | 与 2022/2024/2025 公报逐项核对 |
| 02 | `02_census2020_tables.py` | 七普年鉴 2020 表 6-1/6-3/6-4/6-5/6-12（官网 JPG，`data/raw/china_provincial/raw7rp/`） | `census2020_table6_*.csv` + `ocr_audit.json` | 全国=Σ31省、单岁和=总计、TFR=5×ΣASFR、Σ活产组=人数、性别比恒等式；失败单元格 4–6× 重 OCR，仍失败按恒等式推导并显式记录 |
| 03 | `03_census2020_age_sex.py` | 七普年鉴表 1-5 各地区分年龄、性别的人口（JPG，6076×1455） | `census2020_table1_5_provincial_age_sex.csv` + `ocr_a0105_audit.json` | 墨迹投影切 69 列带→上下文+8× 单元格 OCR；全国 0-14=253.38M、60-64=73.38M 与公报一致；27/31 省 age_sum/合计=1.000 |
| 04a | `04a_yearbook_download.py` | 《中国统计年鉴》2005–2025 各版分地区人口/出生率/死亡率/自然增长率/城镇人口比重表（HTML/XLS/JPG 三格式） | `data/raw/china_provincial/yearbooks/` | 下载后逐文件完整性检查，URL 与标题记入 `crawl_log.json` |
| 04b | `04b_yearbooks_panel.py` | 上年鉴表的结构化解析 | `china_provincial_panel_2004_2024.csv` | CBR−CDR=NGR；全国 vs Σ省；跨年仲裁（用独立人口表纠正 2 处）；已知锚点核对（2016 山东 17.89‰ 等 8 项） |
| 05 | `05_policy_timeline.py` | 31 省人大常委会决议（单独二孩落地日期，经 zh.wikipedia 汇总表转录并抽样核对官方报道）+ 全国政策节点 | `provincial_policy_timeline.csv` | 浙江 2014-01-17 / 北京 2014-02-21 / 山东 2014-05-30 / 西藏 2014-09-17 与官方报道一致；新疆标记 `approximate` |
| 06 | `06_hfc_international.py` | Human Fertility Collection `HFC_ASFR_BO`（fertilitydata.org 免登录 ZIP）+ HFD STFF 月度出生 | `data/external/hfd_wpp/hfc_asfr_bo_panel.csv`、`hfc_country_coverage.csv` | 仅保留周岁口径（ACY）、剔除虚岁（ARDY）重复；19 国（16 开发 + 3 严格独立测试国）；逐国完整矩阵年份数审计 |

说明：
- 阶段 04a 依赖外网与站点可用性，失败会重试并写日志；04b 对 2005–2009 年“分省合计低于全国 1.4–1.9%”不做人为修补（官方未校准估计，属真实数据特性），已在 `../reports/data_quality_audit.md` 说明。
- 阶段 02/03 的 OCR 结果**不允许静默插补**：任何无法识别的单元格要么留空、要么按表内恒等式推导，并在审计 JSON 中留痕。
