# CIPD — Constraint-based Identification and Policy Dynamics

**Structure-Preserving Partial Identification and Robust Policy Decisions in Partially Observed Demographic Dynamics**
**不完全观测下人口动力系统的结构保持型部分识别与鲁棒政策决策——以中国生育政策改革为例**

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Data hashes: 85 verified](https://img.shields.io/badge/data%20hashes-85%20SHA--256%20verified-brightgreen)](provenance/hashes.sha256)
[![Unit tests: 8/8](https://img.shields.io/badge/unit%20tests-8%2F8%20passing-brightgreen)](tests/test_cipd_core.py)
[![R baselines: optional](https://img.shields.io/badge/R%20baselines-optional-lightgrey)](provenance/code_sources.yaml)

> **TL;DR** — CIPD does not try to produce another point forecast of births. It compiles age–parity demographic accounting into a sparse **flow polytope** and answers three questions with auditable linear programs: *What range can the target be bounded in?* (primal–dual certificates), *Which missing table is worth collecting first?* (two-trajectory value of information), and *Which policy cannot be beaten in the worst case?* (baseline-safe minimax regret). Validated on 19 countries of real age×parity fertility data, 31 Chinese provinces, and a 2,160-row reconstruction benchmark; every number below is regenerated from `results/tables/`.

---

## 1. Project & Positioning / 项目定位

CIPD is the research codebase behind a Chinese course paper and an English manuscript in preparation; **the manuscripts themselves are not part of this repository**. It is built around one idea: when micro mechanisms are unobservable and official statistics are revised, the honest object is not a point estimate but a **certified interval plus a decision rule**.

The engine models the female population as a time-expanded network of age × parity flows. With conservation equations `Az = b`, rate boxes and observation intervals, the identified set is a polytope

$$
\mathcal{P} = \{ z \ge 0 : Az = b,\ \underline{y} \le Hz \le \overline{y} \},
$$

and every reported bound is the solution of an LP whose primal and dual residuals are checked independently:

$$
\underline{Y} = \min_{z \in \mathcal{P}} c^\top z, \qquad
\overline{Y} = \max_{z \in \mathcal{P}} c^\top z .
$$

The framework has seven layers (V0–V7): cohort-component accounting (V0), the age–parity flow network (V1), the flow-polytope compiler (V2), goal-oriented partial identification with primal–dual certificates (V3), shared-mechanism counterfactual coupling (V4), causal sensitivity with HonestDiD and its embedding into the flow polytope (V5), two-trajectory value of information (V6), and baseline-safe robust decisions (V7).

| Status | Value |
|---|---|
| Code | Python 3.12 (uv-managed venv); no compiled extension required |
| External baselines | R 4.5.3 *optional*, invoked as isolated subprocesses (`experiments/rmods/`) |
| Licence | MIT for code; external data keep their publishers' terms |
| Reproducibility | 85 data/provenance files SHA-256 locked; 8/8 unit tests; deterministic seeds |
| Results shipped | 41 CSV tables, 11 publication figures, run logs with runtime/memory |

---

## 2. Key Findings / 核心发现

All values below come from the shipped CSVs; see §6 for the file mapping.

| Exp. | Question | Headline result (real data) | Evidence |
|---|---|---|---|
| **E01** | Is the accounting exact? | Conservation residual `2.7e-16` (V0) / `2.2e-16` (V1); adjoint gradient matches finite differences to `5e-6` at **155.7×** speed-up | `E01_accounting_tests.csv` |
| **E02** | Reconstruct masked age×parity matrices? | 432 country×mask×seed instances, 16 countries, 5 methods. CIPD MAE **0.00062–0.00078** for 10–90% random masking (MARE 4.4–4.9%), **2.7–16.9× better than every baseline at ≥50% masking**; certified intervals cover truth 88.0–96.2% | `E02_multicountry_summary.csv`, `E02_multicountry_cipd_summary.csv` |
| **E03** | Recover causal signs on synthetic DGPs? | 10 DGPs × 200 Monte-Carlo runs: CIPD 100% coverage, **0% wrong-sign**; naive trend extrapolation wrong in 20.8% (66.5% on the hidden-confounder DGP) | `E03_synthetic_causal_recovery.csv` |
| **E04** | Strict rolling forecasts? | 587 forecasts × methods, origins 2005/2010/2013/2016/2020. Point accuracy: ETS 0.1179, ARIMA 0.1198, **CIPD 0.1221**, persistence 0.1235. CIPD worst-case bands cover **98.5%** | `E04_rolling_multicountry_*.csv` |
| **E05** | Identify China's three reforms? | 31-province panels, continuous-intensity event studies; **HonestDiD robust CIs cross 0 for every M∈{0,…,3}**; embedding the causal band into the polytope shrinks 2016–2020 births by **56.1%** and 2017 second-parity births by **88.0%** (dual gap ≤5.8e-16) | `E05_honestdid_sensitivity.csv`, `E05_ccpi_did_embedding.csv` |
| **E06** | Decompose the birth decline? | 2010–2025 fall of 800.0 万: quantum **67.7%**, age structure **37.6%**, tempo interaction `+42.0 万`; identity residual 0.0 | `E06_shapley_decomposition.csv` |
| **E07** | Which table to buy first? | Two-trajectory VOI at budget 3.0 shrinks width **96.5%** with **zero gap** to exhaustive optimum (random selection: 10.4%) | `E07_voi_observation_selection.csv` |
| **E08** | Robust policy under the 2025 subsidy baseline? | Universal childcare (`U2`) is the only baseline-safe option: `R(u) = −0.046 ≤ 0` across 64 scenarios; **8/15** pairwise rankings reverse; 2050 population 12.05–12.54 亿 | `E08_robust_policy_metrics.csv`, `E08_pairwise_dominance.csv` |
| **E09** | Does it generalise to unseen countries? | 3 hold-out countries never used for tuning: CHE `0.00061` / 96.9%, GBRTENW `0.00068` / 89.2%, TWN `0.00139` / 100% — all inside the 16-country development range (`0.00038–0.00188`) | `E09_holdout_vs_dev.csv`, `E09_loco_summary.csv` |
| **E10** | Can the solver be accelerated? | Same accuracy (Δobjective `1e-12`): **5.10–5.18×** on the 9,555-variable 26-target sequence via basis reuse; province-block dual decomposition **1.19–1.96×** on the 100,905-variable coupled instance (relative gap `5.5e-9`) | `E10_acceleration.csv` |

### Negative results & scope boundaries / 负面发现与边界

These are kept deliberately, because they bound what the numbers above mean:

1. **Light masking favours interpolation.** At 10% random missingness linear interpolation beats CIPD (`0.00056` vs `0.00062`); when an entire parity channel is missing, temporal interpolation wins (`0.00040` vs `0.00079`). CIPD's advantage appears at ≥50% masking or structural (age-block) gaps.
2. **Bayesian "strong" baselines under-cover.** Ellison et al. (2026) GAM+HMC adapted to period matrices reaches MAE 0.0030–0.0091 but its nominal 90% interval covers only **47.5–81.9%**. `bayesTFR` (origins 2005/2010 only, 5-year resolution) scores MAE 0.1773 on matched cells with **67.2%** coverage; its 2016 origin fails on an internal `column 'year'` defect.
3. **Frozen-rate CCM fails in a fertility transition.** `ccmpp.tmb` (single-age, annual steps, WPP2024 inputs) yields MAE 392.3 万 / MAPE 34.8% on Chinese births, worse than persistence (257.0 万 / 24.7%): the echo-boom cohort structure predicts rising births while births collapsed after 2017.
4. **No point forecaster beats persistence.** Across 587 forecasts the spread between ETS/ARIMA/CIPD/persistence is <5%, and every MASE > 1. Certified worst-case bands are *not* predictive distributions: 98.5% coverage with the worst WIS is a feature of the object, not a bug.
5. **Policy effects are not point-identified.** National-series and 31-province evidence both fail HonestDiD robustness (every interval crosses 0); the 2021 three-child break is indistinguishable from the 2018 placebo break, hence evidence grade **C3**.
6. **The acceleration claim is partial.** Target-cone presolve is useless on the fully time-coupled network (0.84×); the 100,905-variable decomposition is **1.19–1.96×**, so the 2× claim for that instance is **withdrawn** (single-threaded runs: 0.49×).
7. **Data extraction has known holes.** 七普 Table 1-5: national anchors exact (0–14 = 253.38 M, 60–64 = 73.38 M), 27/31 provinces close exactly, but 吉林/黑龙江/浙江/广西 have 4–6% local OCR gaps kept as-is (no imputation); Table 1-10 is not structured; provincial CBR for 2020 (census year) is missing.

---

## 3. Repository Structure / 仓库结构

```text
CIPD/
├── src/cipd/                     # Core engine V0–V7 (data-independent, unit-tested)
│   ├── demography/ccm.py         #   V0/V1 cohort-component, age–parity flows, adjoint, Shapley
│   ├── observation/flow_compiler.py   # V2 time-expanded flow polytope (sparse constraints)
│   ├── identification/partial_id.py   # V3 goal-oriented bounds + primal–dual certificates
│   ├── counterfactual/coupler.py      # V4 shared-mechanism coupling (McCormick + sBB)
│   ├── causality/causal_bounds.py     # V5 Fréchet eligibility + Delta^SD(M) LP
│   ├── information_value/voi.py       # V6 two-trajectory observation separation
│   ├── optimization/robust_policy.py  # V7 baseline-safe minimax regret
│   └── data/                     # Shared data tooling (OCR engine, parsers, dictionaries)
├── crawlers/                     # Data acquisition stages 01–06 (see crawlers/README.md)
│   ├── 01_national_series.py          # NBS annual series + parity shares + marriages
│   ├── 02_census2020_tables.py        # 7th-census long-form tables 6-1/6-3/6-4/6-5/6-12
│   ├── 03_census2020_age_sex.py       # 7th-census table 1-5 provincial age x sex
│   ├── 04a_yearbook_download.py       # China Statistical Yearbook 2005–2025 tables
│   ├── 04b_yearbooks_panel.py         # 31-province annual panel with cross-year arbitration
│   ├── 05_policy_timeline.py          # 31-province Dandu dates + national policy nodes
│   └── 06_hfc_international.py        # HFC age x parity panels, 19 countries (auto-unzip)
├── experiments/
│   ├── level2/run_all.py         # E01–E10 baseline suite + Fig01–Fig06
│   ├── level3/                   # Cross-country, provincial causality, acceleration, ablations
│   │   ├── e02_multicountry_masking.py   e02_run_ellison.py
│   │   ├── e04_rolling_multicountry.py   e04_merge_bayestfr.py
│   │   ├── e05_provincial_did.py         e05_ccpi_did_embedding.py   e05_robustness.py
│   │   ├── e09_loco.py                   e10_acceleration.py          ablations_full.py
│   │   └── make_level3_figures.py        # Fig07–Fig11 (pure table -> figure)
│   └── rmods/                    # R baselines, isolated subprocesses (GPL lives here)
│       ├── e02_ellison_gam.R     #   Ellison et al. (2026) P-spline GAM + HMC (Stan)
│       ├── e04_bayestfr.R        #   bayesTFR native WPP2019 projections
│       ├── e04_ccm_china.R       #   ccmpp.tmb cohort-component baseline
│       ├── e05_did_honest.R      #   did::att_gt + HonestDiD sensitivity grid
│       └── e05_honestdid_single.R
├── data/
│   ├── raw/                      # China: NBS series, 2010/2020 census tables, yearbook tables
│   └── external/                 # HFC 19-country panels, STFF, WPP2024 mirror, R baseline code
├── provenance/                   # sources.yaml, hashes.sha256, code_sources.yaml, policy_codebook
├── results/
│   ├── tables/                   # 41 CSV result tables (every paper number comes from here)
│   ├── figures/                  # Fig01–Fig11 (PNG, 180 dpi)
│   └── logs/                     # runtime, peak memory, solver status, seeds
├── reports/                      # acceptance, data audit, claim-evidence matrix, mock reviews
├── research_log/                 # decisions, failures, negative results
├── tests/test_cipd_core.py       # 8 unit tests (V0–V7 mathematical properties)
├── Makefile                      # setup / data / crawl / audit / test / level2 / level3 / figures
└── requirements.txt              # pinned Python stack
```

---

## 4. Quick Start & Reproduction / 安装与复现

### 4.1 Prerequisites

- Python **3.12** (uv recommended), ~2 GB free disk for the venv.
- Optional: **R 4.5.3**, only to re-run the external baselines (`bayesTFR`, `did`, `HonestDiD`, `ccmpp.tmb`, Ellison GAM). Without R, `make level3` still runs E02/E04/E05/E09/E10 and all ablations; R-dependent steps print a notice.
- Linux or macOS; on Windows use WSL2 (paths assume POSIX).

### 4.2 From zero to verified results

```bash
# 1) Python environment
uv venv --python 3.12 .venv
uv pip install -r requirements.txt

# 2) Verify the frozen data artifacts (85 SHA-256 hashes)
make audit          # expected: every line "OK"

# 3) Core mathematical tests (no data needed)
make test           # expected: Ran 8 tests ... OK

# 4) Level-3 experiments (behind the paper's tables)
make level3         # ~25–40 min; writes results/tables/* and results/logs/*

# 5) All figures, pure table -> figure (no hidden numbers)
make figures        # regenerates Fig01–Fig11
```

Each experiment is standalone and prints its own summary, for example:

```bash
python experiments/level3/e05_ccpi_did_embedding.py
# total_births_2016_2020: P0 width 124,084,886 -> PDiD 54,458,205 (-56.1%)
# p2_births_2017      : P0 width 16,065,195 -> PDiD  1,926,978 (-88.0%)
```

### 4.3 Optional: R baselines

```bash
# R 4.5.3 + CRAN binaries (Posit Package Manager):
#   did 2.1.2, HonestDiD 0.2.8 (needs Rglpk+GLPK), bayesTFR, ccmpp.tmb, rstan, TMB
Rscript experiments/rmods/e04_bayestfr.R results/tables/E04_bayestfr_raw.csv "2005,2010,2016" 600
Rscript experiments/rmods/e05_did_honest.R results/logs/e05_did_input.csv \
        results/logs/e05_did_output.json results/logs/e05_honestevent.json
```

GPL components (`did`, `HonestDiD`, `bayesTFR`, `Rglpk`, `GLPK`) are **never linked** into this repository; they run as separate processes through `experiments/rmods/`. Exact versions/commits: `provenance/code_sources.yaml` and `THIRD_PARTY_NOTICES.md`.

### 4.4 Reproducibility contract

| Item | How it is fixed |
|---|---|
| Randomness | E02 mask seeds `101–505` (five seeds for masks A–E; Ellison uses 101 only, for HMC cost); E03 Monte-Carlo seed `20261126`; all other experiments are deterministic LPs |
| Data version | 85 files SHA-256 locked in `provenance/hashes.sha256` (data + provenance only; result tables are outputs); `make audit` fails loudly on drift |
| Solver accuracy | every bound reports primal/dual residuals and relative duality gap; independent re-check in `partial_id._verify_lp_certificate` |
| Figures | generated only from `results/tables/*.csv` — no hand-typed numbers |
| Logs | `results/logs/experiment_run_log.json` records runtime, peak memory, solver status per experiment |

---

## 5. Data Provenance & Licensing / 数据来源与授权

All raw data come from public official channels; nothing is simulated or imputed under an official label.

| ID | Source | Content used |
|---|---|---|
| D01–D03 | 国家统计局 (NBS) communiqués & yearbooks | 1990–2025 national vital series; 2012–2025 parity shares & marriages; census age anchors |
| D04 | 国家卫生健康委员会 (NHC) | 2025 national childcare subsidy (3,600 CNY/child/year) |
| D05 | Statistics Korea / KOSIS | 2014–2025 vital statistics (external validation) |
| D07 | 中国人口普查年鉴-2020 (7th census, JPG tables OCR'd) | Tables 6-1/6-3/6-4/6-5/6-12, 1-5 |
| D08 | 中国统计年鉴 2005–2025 | 31-province panel: population, urbanisation, CBR/CDR/NGR |
| D09 | Human Fertility Collection | 19 countries, single-age × parity ASFR (16 development + 3 hold-out) |
| D10 | HFD STFF | Monthly births, ~50 countries |
| D11 | UN WPP2024 (via PPgp R mirror) | External reference only — never used for fitting/tuning in E04 |
| D12 | Provincial legislature records (transcribed, spot-checked) | 31-province Dandu-2nd-child effective dates |

- Full URLs, retrieval dates, units and audit notes: [`provenance/sources.yaml`](provenance/sources.yaml).
- Extraction quality (OCR identities, rescued cells, derived values): [`data/raw/china_provincial/ocr_audit.json`](data/raw/china_provincial/ocr_audit.json), [`reports/data_quality_audit.md`](reports/data_quality_audit.md).
- **Size policy:** raw HFC text dumps (43–88 MB each) are not committed; the ZIP archives (28 MB total) are. `crawlers/06_hfc_international.py` auto-unzips them, or run `make crawl` to rebuild every raw artifact. The committed repo stays ≈55 MB.
- **Licensing:** code is MIT (see §8). External data remain © their publishers (NBS, NHC, KOSIS, HFC/MPIDR, UN) and appear here only as derived, citable research artifacts; check each publisher's terms before further redistribution.

---

## 6. Artifact Mapping / 实验与图表索引

### 6.1 Figures (paper numbering -> file -> generator)

| Paper Fig. | File | Produced by | Reads |
|---|---|---|---|
| Fig. 1 | `Fig01_china_population_policy_timeline.png` | `experiments/level2/run_all.py` | `data/raw/*`, `policy_codebook.csv` |
| Fig. 2 | `Fig02_E02_masking_identification.png` | `level2/run_all.py` | `E02_real_census2010_masking.csv` |
| Fig. 3 | `Fig07_E02_multicountry_headtohead.png` | `level3/make_level3_figures.py` | `E02_multicountry_summary.csv` |
| Fig. 4 | `Fig09_E09_loco_generalization.png` | `level3/make_level3_figures.py` | `E09_holdout_vs_dev.csv` |
| Fig. 5 | `Fig03_E03_synthetic_causal_recovery.png` | `level2/run_all.py` | `E03_synthetic_causal_recovery.csv` |
| Fig. 6 | `Fig08_E04_rolling_forecast_headtohead.png` | `level3/make_level3_figures.py` | `E04_rolling_multicountry_summary.csv` |
| Fig. 7 | `Fig04_E05_china_policy_sensitivity.png` | `level2/run_all.py` | `E05_china_policy_bounds.csv` |
| Fig. 8 | `Fig10_E05_provincial_event_studies_honestdid.png` | `level3/make_level3_figures.py` | `E05_event_study_*.csv`, `E05_honestdid_sensitivity.csv` |
| Fig. 9 | `Fig11_E10_E05CCPI_embedding.png` | `level3/make_level3_figures.py` | `E10_acceleration.csv`, `E05_ccpi_did_embedding.csv` |
| Fig. 10 | `Fig06_E06_decomposition_and_E08_projections.png` | `level2/run_all.py` | `E06_shapley_decomposition.csv`, `E08_china_longterm_scenarios_2025_2070.csv` |
| Fig. 11 | `Fig05_E07_VOI_and_E08_robust_policy.png` | `level2/run_all.py` | `E07_voi_observation_selection.csv`, `E08_robust_policy_metrics.csv` |

### 6.2 Tables (paper -> result files)

- Main-text tables 1–13 read directly from `results/tables/E*.csv` (e.g. Table 1 ← `data/raw/china_census_age_anchors_2000_2025.csv`; E02 tables ← `E02_multicountry_summary.csv`; E05 tables ← `E05_honestdid_sensitivity.csv`, `E05_ccpi_did_embedding.csv`; E10 table ← `E10_acceleration.csv`).
- The paper's appendix tables (18–27) are generated from these same CSVs by a local LaTeX-package script, so no number is typed twice; the manuscript source is kept out of this repository.
- Full experiment → script → table → metric matrix: `reports/claim_evidence_matrix.csv`.

---

## 7. Citation / 引用

```bibtex
@software{zhang2026cipd,
  author       = {Zhang, Jinghan},
  title        = {{CIPD}: Constraint-based Identification and Policy Dynamics ---
                  Structure-Preserving Partial Identification and Robust Policy
                  Decisions in Partially Observed Demographic Dynamics},
  year         = {2026},
  version      = {1.0.0},
  license      = {MIT},
  url          = {https://github.com/Oookla28/CIPD-Demographic-Dynamics}
}
```

This repository ships [`CITATION.cff`](CITATION.cff), so GitHub shows a **"Cite this repository"** button in the right sidebar. If you use only the benchmark tables, please cite the version tag you used.

---

## 8. License / 开源协议

- **Code:** MIT — see [`LICENSE`](LICENSE).
- **Third-party R/GPL dependencies:** used as separate executables (`Rscript ...`); no GPL code is copied, linked, or vendored into this repository. Details: [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md), [`provenance/code_sources.yaml`](provenance/code_sources.yaml).
- **Data:** © respective official publishers (see §5); derived tables are provided for reproducibility of the published analysis.

> Note: if the GitHub repository page currently shows a GPL-3.0 licence tab, set the repository licence to **MIT** (Settings → License) after replacing `LICENSE` with the MIT text — GPL obligations attach only to the separately invoked R packages, not to this codebase.

---

### Topics

`demography` · `partial-identification` · `causal-inference` · `robust-optimization` · `linear-programming` · `china-fertility-policy` · `reproducible-research`
