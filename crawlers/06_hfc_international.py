#!/usr/bin/env python3
"""Parse Human Fertility Collection (HFC) ASFR-by-birth-order data into tidy
per-country single-age x parity panels for E02/E04/E09.

Input:  data/external/hfd_wpp/raw_hfc/HFC_ASFR_BO/HFC_ASFR_BO.txt
Output: data/external/hfd_wpp/hfc_asfr_bo_panel.csv   (long: country,year,age,asfr_p1,p2,p3p)
        data/external/hfd_wpp/hfc_country_coverage.csv (coverage audit)

Source note: HFC (fertilitydata.org) redistributes official national vital
statistics; SourceType column is retained for provenance.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/external/hfd_wpp/raw_hfc"
SRC = RAW / "HFC_ASFR_BO/HFC_ASFR_BO.txt"


def _ensure_extracted():
    """原始 ZIP 随仓库分发，解压后的 .txt 不入库；缺文件时按需解压到约定子目录。"""
    import zipfile
    target = SRC
    if target.exists():
        return
    zp = RAW / "HFC_ASFR_BO.zip"
    if not zp.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zp) as z:
        members = [m for m in z.namelist() if m.endswith(".txt") and "__MACOSX" not in m]
        if members:
            with z.open(members[0]) as src:
                target.write_bytes(src.read())
            print(f"[stage 06] extracted HFC_ASFR_BO.zip -> {target.relative_to(ROOT)}")


_ensure_extracted()
OUT = ROOT / "data/external/hfd_wpp"

# Study design: development (tuning allowed) vs strictly independent hold-out
# test countries (never used in any tuning; E09 final evaluation only).
DEV_COUNTRIES = ["JPN", "KOR", "DEUTNP", "SWE", "NOR", "FIN", "DNK", "AUT",
                 "CZE", "NLD", "HUN", "ITA", "ESP", "USA", "CAN", "POL"]
HOLDOUT_COUNTRIES = ["GBRTENW", "CHE", "TWN"]
ALL_COUNTRIES = DEV_COUNTRIES + HOLDOUT_COUNTRIES


def main():
    rows = list(csv.reader(open(SRC), delimiter=","))
    hdr = [h.strip() for h in rows[0]]
    idx = {h: i for i, h in enumerate(hdr)}
    data = []
    for r in rows[1:]:
        if len(r) != len(hdr):
            continue
        d = {h: r[i].strip() for h, i in idx.items()}
        data.append(d)

    out_rows = []
    cov = defaultdict(lambda: {"years": set(), "ages": set(), "sources": set(),
                               "n_cells": 0, "n_missing": 0, "agedef": set()})
    per_year_ages = defaultdict(set)
    for d in data:
        c = d["Country"]
        if c not in ALL_COUNTRIES:
            continue
        if d["Region"] != "." or d["Urban"] != "." or d["Origin"] != ".":
            continue  # national totals only
        # keep completed-years age definition (ACY); if a country only publishes
        # ARDY (age reached in year), keep ARDY and flag it
        agedef = d["AgeDef"]
        if d["AgeInt"] != "1":
            # allow age-15 open-interval marker rows (-99) as the age-15 value
            if not (d["Age"] == "15" and d["AgeInt"] == "-99"):
                continue
        try:
            age = int(d["Age"])
            year = int(d["Year1"])
        except ValueError:
            continue
        if not (15 <= age <= 49):
            continue
        cov[c]["agedef"].add(agedef)
        vals = {}
        for k, col in [("asfr_total", "ASFR"), ("asfr_p1", "ASFR1"), ("asfr_p2", "ASFR2"),
                       ("asfr_p3", "ASFR3"), ("asfr_p3p", "ASFR3P")]:
            v = d[col].strip()
            vals[k] = float(v) if v not in (".", "", "-99") else None
        # Aggregate ASFR3 and ASFR3P -> asfr_p3p (prefer explicit 3P, else 3)
        if vals["asfr_p3p"] is None and vals["asfr_p3"] is not None:
            vals["asfr_p3p"] = vals["asfr_p3"]
        cov[c]["years"].add(year)
        cov[c]["ages"].add(age)
        cov[c]["sources"].add(d["SourceType"])
        per_year_ages[(c, year)].add(age)
        if vals["asfr_total"] is None:
            cov[c]["n_missing"] += 1
        cov[c]["n_cells"] += 1
        out_rows.append({
            "country": c, "role": "dev" if c in DEV_COUNTRIES else "holdout",
            "year": year, "age": age, "age_def": agedef,
            "asfr_total": vals["asfr_total"],
            "asfr_p1": vals["asfr_p1"],
            "asfr_p2": vals["asfr_p2"],
            "asfr_p3p": vals["asfr_p3p"],
            "source_type": d["SourceType"],
        })

    # coverage audit: how many of ages 15-49 present per year
    cov_rows = []
    for c in sorted(cov):
        years = sorted(cov[c]["years"])
        full_years = [y for y in years if set(range(15, 50)) <= per_year_ages[(c, y)]]
        cov_rows.append({
            "country": c, "role": "dev" if c in DEV_COUNTRIES else "holdout",
            "n_years": len(years), "first_year": years[0], "last_year": years[-1],
            "n_full_years_15_49": len(full_years),
            "n_cells": cov[c]["n_cells"], "n_missing_total": cov[c]["n_missing"],
            "sources": ";".join(sorted(cov[c]["sources"])),
            "age_defs": ";".join(sorted(cov[c]["agedef"])),
        })

    # dedupe on (country, year, age): prefer ACY, then strict AgeInt=1 rows
    best = {}
    for r in out_rows:
        key = (r["country"], r["year"], r["age"])
        rank = (0 if r["age_def"] == "ACY" else 1)
        if key not in best or rank < best[key][0]:
            best[key] = (rank, r)
    out_rows = [v[1] for v in best.values()]
    out_rows.sort(key=lambda r: (r["country"], r["year"], r["age"]))
    with open(OUT / "hfc_asfr_bo_panel.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    with open(OUT / "hfc_country_coverage.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cov_rows[0].keys()))
        w.writeheader()
        w.writerows(cov_rows)

    print(f"panel rows: {len(out_rows)}  countries: {len(cov)}")
    for cr in cov_rows:
        print(f"  {cr['country']:8} {cr['role']:8} years {cr['first_year']}-{cr['last_year']}"
              f" n_years={cr['n_years']} full15-49={cr['n_full_years_15_49']} src={cr['sources'][:20]}")


if __name__ == "__main__":
    main()
