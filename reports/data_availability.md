# 数据可获取性审计报告（Data Availability，Level-3 更新）

**审计日期**：2026-11-26
**运行环境**：本地 macOS arm64（可联网；Python 3.12 venv + R 4.5.3 便携安装）

## 1. 已核实并入库的真实公开数据（全部 SHA-256 锁定，共 104 文件）

### 中国
1. **中国 1990-2025 年度人口与自然变动序列**（`data/raw/china_annual_vital_1990_2025.csv`）：与 NBS 2022/2024/2025 公报逐项核对。
2. **中国 2012-2025 分孩次出生占比与结婚登记**（`data/raw/china_parity_marriage_2012_2025.csv`）。
3. **2010 六普长表 6-3**（单岁 15-49 × 3 孩次）与 **6-5**（31 省活产子女分布，本轮新增）：官方 HTML 表直接解析。
4. **2020 七普年鉴 6-1/6-3/6-4/6-5/6-12**（`data/raw/china_provincial/`）：官网 JPG 经 OCR + 结构恒等式校验（详见 `data_quality_audit.md`）。
5. **31 省 2004-2024 面板**（人口/城镇化/CBR/CDR/NGR）：《中国统计年鉴》2005-2025 版 HTML/Excel/JPG 三格式解析。
6. **省级政策时间线**（`data/raw/china_provincial/provincial_policy_timeline.csv`）：31 省单独二孩精确到日 + 全国三节点。
7. **2025 国家育儿补贴**（policy_codebook.csv）。

### 国际（19 国）
8. **HFC ASFR by birth order**（fertilitydata.org，免登录）：16 开发国 + 3 独立测试国（GBRTENW/CHE/TWN），单岁 × {1,2,3+}，1975-2025 为主。
9. **HFD STFF 月度出生**（约 50 国，2000-2025）。
10. **UN WPP2024**（PPgp 官方 R 镜像 `wpp2024`/`wpp2024extra`）。
11. **Ellison et al. (2026) 复现代码**（Zenodo 15179793）。

## 2. 受阻/降级项（如实记录）

| 项目 | 原因 | 降级处理 |
|---|---|---|
| HFD 批量 ZIP（8-16 国全长序列） | 官网改版后需登录账号 | 使用 HFC（同一 HFD/MPIDR 体系、免登录开放数据）；HFC 覆盖 58 国，本项目的 19 国均在其中 |
| UN WPP 官方 API / bulk CSV | API 返回 401（需密钥）、bulk 链接改版 404 | 使用 PPgp 官方镜像 R 包（数据同源） |
| data.stats.gov.cn easyquery API | 403（WAF 拦截） | 改走统计年鉴网页表（OCR/HTML/Excel） |
| 七普 A0105（1-5 分省年龄×性别） | 源 JPG 数字仅 ~6px 高 | ✅ **本轮补齐**：墨迹投影切列 + 8× 单元格 OCR；全国行锚点精确、27/31 省完整（`census2020_table1_5_provincial_age_sex.csv`；4 省部分缺口已标注） |
| 七普 A0110（1-10 分月出生） | 同上 | 保持 PARTIAL（原图哈希存档）；2020 出生总量已由表 6-3 承载 |
| 省级 2020 年 CBR | 普查年无分省抽样表 | E05 将 2020 作为缺失年；2020 出生总量由七普 6-3 承载 |
| 省级 2025 年数据 | 省级公报爬取未执行 | 全国 2025 序列完整；省级面板止于 2024 |
| ccmpp.tmb 完整多国基线 | 仅中国完成（WPP2024 单岁人口 + 冻结 ASFR） | 中国 E04 已完成；多国 CCM 未跑（记录于报告） |
| 省级交错政策（育儿假/补贴）落地日期库 | 无权威汇总表可爬 | E05 使用连续强度/时点设计，并把“全国同年推行→无真正交错采纳”作为识别限制明确写入报告与论文 |

## 3. 复现与溯源

- `provenance/sources.yaml`：12 个来源（D01-D12）含 URL、机构、口径、状态。
- `provenance/hashes.sha256`：104 个文件，`make audit` 通过。
- `provenance/code_sources.yaml`：外部软件栈全部登记（精确版本/commit）。
- OCR 审计：`data/raw/china_provincial/ocr_audit.json`、`yearbooks/panel_audit.json`。
