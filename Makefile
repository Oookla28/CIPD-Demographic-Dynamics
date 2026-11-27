.PHONY: setup data crawl audit test level2 level3 figures reproduce paper clean

PY := .venv/bin/python
export PATH := $(HOME)/.local/bin:$(PATH)

## 环境（Python 3.12 venv + R 4.5.3 用户安装；R 仅外部基线需要）
setup:
	@python3 --version
	@test -x $(PY) || (echo "Create venv: uv venv --python 3.12 .venv && uv pip install -r requirements.txt" && exit 1)
	@$(PY) -c "import numpy, scipy, pandas, matplotlib, pyfixest, statsmodels, highspy; print('Python env OK')"
	@Rscript -e 'cat(R.version.string, "\n")' 2>/dev/null || echo "R optional (only for external baselines)"

## 数据校验（冻结产物，哈希锁定；重建用 make crawl）
data:
	$(PY) -I -c "import pathlib,sys; fs=[pathlib.Path(p) for p in \
	  ['data/raw/china_annual_vital_1990_2025.csv', \
	   'data/raw/china_provincial/china_provincial_panel_2004_2024.csv', \
	   'data/external/hfd_wpp/hfc_asfr_bo_panel.csv']]; \
	  ok=all(f.exists() and f.stat().st_size>1000 for f in fs); \
	  print('data artifacts OK' if ok else 'MISSING'); sys.exit(0 if ok else 1)"

## 数据抓取/解析全链（会改写 data/ 下产物，谨慎执行）
crawl:
	$(PY) -I crawlers/01_national_series.py
	$(PY) -I crawlers/02_census2020_tables.py
	$(PY) -I crawlers/03_census2020_age_sex.py
	$(PY) -I crawlers/04_yearbooks_panel.py
	$(PY) -I crawlers/05_policy_timeline.py
	$(PY) -I crawlers/06_hfc_international.py

## 溯源哈希审计（112 文件）
audit:
	shasum -a 256 -c provenance/hashes.sha256 | tail -3

test:
	$(PY) -m unittest discover -s tests -v

## 二级实验（E01–E10 初始全套；输出 results/tables/E0*.csv）
level2:
	$(PY) -I experiments/level2/run_all.py

## 三级实验（跨国基准/省级因果/加速/消融/绘图）
level3:
	$(PY) -I experiments/level3/e02_multicountry_masking.py
	$(PY) -I experiments/level3/e02_run_ellison.py
	$(PY) -I experiments/level3/e04_rolling_multicountry.py
	$(PY) -I experiments/level3/e09_loco.py
	$(PY) -I experiments/level3/e05_provincial_did.py
	$(PY) -I experiments/level3/e05_ccpi_did_embedding.py
	$(PY) -I experiments/level3/e05_robustness.py
	$(PY) -I experiments/level3/e10_acceleration.py
	$(PY) -I experiments/level3/ablations_full.py

figures:
	$(PY) -I experiments/level3/make_level3_figures.py
	$(PY) -I experiments/level3/fix_e05_e09_and_fonts.py || true
	$(PY) -I experiments/level3/finish_e10_and_figures.py || true

paper:
	cd paper/journal_version && latexmk -xelatex -interaction=nonstopmode cipd_journal_draft.tex || true

reproduce: setup data audit test level3 figures
	@echo "Reproduction complete: hash audit, unit tests, Level-3 experiments, figures."

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
