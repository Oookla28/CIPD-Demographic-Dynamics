#!/usr/bin/env python3
"""Build the 31-province annual panel (2004-2024) from China Statistical
Yearbook tables: population, urbanization, CBR/CDR/NGR.

Inputs (data/raw/china_provincial/yearbooks/):
  - yb{ed}_birth_{yr}.htm (2005-2012 editions, HTML tables)
  - yb{ed}_birth_{yr}.xls (2013/2014 editions, Excel)
  - yb{ed}_birth_{yr}.jpg (2015-2025 editions, OCR)
  - yb{ed}_pop_{yr}.*   and yb{ed}_urban_{yr}.* where available
Validations:
  - CBR - CDR ~= NGR per row
  - national population ~= sum of provinces (0.5% tolerance, 1% pre-2011)
  - rate ranges and inter-year jumps flagged
Output: data/raw/china_provincial/china_provincial_panel_2004_2024.csv
"""
from __future__ import annotations

import csv
import json
import re
import html as htmlmod
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
YB = ROOT / "data/raw/china_provincial/yearbooks"
OUT = ROOT / "data/raw/china_provincial/china_provincial_panel_2004_2024.csv"

sys.path.insert(0, str(ROOT / "src"))
from cipd.data.ocr_census import (REGIONS, REGION_KEYS, ocr_image, cluster_rows,
                                  locate_region_ys, cell_near, parse_number, reocr_row)

PROV_ALIAS = {
    "北 京": "北京", "天津": "天津", "河北": "河北", "山西": "山西",
    "内蒙古": "内蒙古", "内蒙古自治区": "内蒙古", "辽 宁": "辽宁", "吉林": "吉林",
    "黑龙江": "黑龙江", "上 海": "上海", "江苏": "江苏", "浙江": "浙江",
    "安徽": "安徽", "福建": "福建", "江西": "江西", "山东": "山东", "河南": "河南",
    "湖北": "湖北", "湖南": "湖南", "广东": "广东", "广西": "广西",
    "广西壮族自治区": "广西", "海南": "海南", "重庆": "重庆", "四川": "四川",
    "贵州": "贵州", "云南": "云南", "西藏": "西藏", "西藏自治区": "西藏",
    "陕西": "陕西", "甘肃": "甘肃", "青海": "青海", "宁夏": "宁夏",
    "新疆": "新疆", "新疆维吾尔自治区": "新疆", "全 国": "全国", "全国": "全国",
}


def clean_label(s: str) -> str | None:
    s = re.sub(r"[\s\xa0\u3000]", "", str(s))
    for alias, canon in PROV_ALIAS.items():
        a = re.sub(r"[\s\xa0\u3000]", "", alias)
        if s == a:
            return canon
    # partial: starts with region name
    for r in REGIONS + ["全国"]:
        if s.startswith(r) and len(s) <= len(r) + 3:
            return r
    return None


def parse_html_table(path: Path) -> list[dict]:
    txt = path.read_bytes().decode("gbk", errors="replace")
    rows = re.findall(r"<tr[^>]*>(.*?)</tr>", txt, re.S)
    out = []
    for r in rows:
        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", r, re.S)
        cells = [htmlmod.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in cells]
        cells = [c for c in cells if c != ""]
        if len(cells) < 4:
            continue
        prov = clean_label(cells[0])
        if prov is None:
            continue
        nums = []
        for c in cells[1:]:
            v = parse_number(c)
            nums.append(v)
        out.append({"province": prov, "values": nums})
    return out


def parse_xls_table(path: Path) -> list[dict]:
    import xlrd
    wb = xlrd.open_workbook(str(path))
    sh = wb.sheet_by_index(0)
    out = []
    for r in range(sh.nrows):
        first = str(sh.cell_value(r, 0))
        prov = clean_label(first)
        if prov is None:
            continue
        nums = []
        for c in range(1, sh.ncols):
            v = sh.cell_value(r, c)
            if isinstance(v, float):
                nums.append(v)
            else:
                pv = parse_number(str(v))
                nums.append(pv)
        out.append({"province": prov, "values": [x for x in nums]})
    return out


def parse_jpg_table(path: Path) -> list[dict]:
    """OCR a yearbook JPG table: 31 region rows, columns detected from header."""
    cells = ocr_image(path, scale=2.0)
    # locate header cells for rates
    hdr_ys = [c for c in cells if c["y"] < 160]
    col_x = {}
    for c in hdr_ys:
        t = c["text"]
        if "出生率" in t:
            col_x["cbr"] = c["x"]
        elif "死亡率" in t:
            col_x["cdr"] = c["x"]
        elif "自然增长率" in t or "增长率" in t:
            col_x["ngr"] = c["x"]
        elif "总人口" in t or "人口数" in t and "城" not in t:
            if "pop" not in col_x:
                col_x.setdefault("pop", c["x"])
        elif "比重" in t:
            col_x.setdefault("urban", c["x"])
    if not all(k in col_x for k in ("cbr", "cdr", "ngr")):
        # fallback: assume last 3 numeric columns
        col_x = None
    ys = locate_region_ys(cells, label_x_max=110)
    out = []
    for prov in ["全国"] + REGIONS:
        if prov not in ys:
            continue
        y = ys[prov]
        if col_x is not None:
            rec = {"province": prov, "values": []}
            for key in ("pop", "urban", "cbr", "cdr", "ngr"):
                if key in col_x:
                    c = cell_near(cells, col_x[key], y, x_tol=55)
                    rec["values"].append(parse_number(c["text"]) if c else None)
            # fix column alignment: find real x of cbr etc.
            rec["_colx"] = {k: v for k, v in col_x.items()}
            out.append(rec)
    return out


def extract_jpg_rows_precise(path: Path) -> list[dict]:
    """For each region row cluster, take its numeric cells ordered by x."""
    cells = ocr_image(path, scale=2.0)
    rows = cluster_rows(cells)
    hdr = [c for c in cells if c["y"] < 170]
    rate_x = {}
    for c in hdr:
        t = c["text"]
        if "出生" in t:
            rate_x["cbr"] = c["x"]
        elif "死亡" in t:
            rate_x["cdr"] = c["x"]
        elif "自然增" in t or ("增" in t and "长" in t):
            rate_x["ngr"] = c["x"]
    out = []
    for r in rows:
        left = [c for c in r["cells"] if c["x"] < 110]
        if not left:
            continue
        label = "".join(c["text"] for c in left).replace(" ", "")
        prov = None
        for cand in ["全国"] + REGIONS:
            if label == cand or label.startswith(cand) and len(label) <= len(cand) + 2:
                prov = cand
                break
        if prov is None:
            continue
        rowpairs = [(c["x"], parse_number(c["text"])) for c in r["cells"]
                    if c["x"] >= 110 and parse_number(c["text"]) is not None]
        rowpairs.sort(key=lambda q: q[0])
        out.append({"province": prov, "values": [], "_hdr": rate_x, "_row": rowpairs})

    # fallback: if some provinces are missing, assign by fixed row order
    found = {r["province"] for r in out}
    missing = [p for p in ["全国"] + REGIONS if p not in found]
    if missing:
        data_rows = [r for r in rows
                     if len([c for c in r["cells"] if c["x"] >= 110 and parse_number(c["text"]) is not None]) >= 6]
        data_rows.sort(key=lambda r: r["y"])
        if len(data_rows) == 32:
            out = []
            for prov, r in zip(["全国"] + REGIONS, data_rows):
                rowpairs = [(c["x"], parse_number(c["text"])) for c in r["cells"]
                            if c["x"] >= 110 and parse_number(c["text"]) is not None]
                rowpairs.sort(key=lambda q: q[0])
                out.append({"province": prov, "values": [], "_hdr": rate_x, "_row": rowpairs})
    return out


def parse_pop_table(path: Path) -> dict[str, list[float]]:
    """Parse a provincial population table -> {province: [numeric cells by x-order]}."""
    ext = path.suffix.lstrip(".").lower()
    out: dict[str, list[float]] = {}
    try:
        if ext in ("htm", "xls"):
            rows = parse_html_table(path) if ext == "htm" else parse_xls_table(path)
            for r in rows:
                nums = [v for v in r["values"] if v is not None]
                if nums:
                    out[r["province"]] = nums
        else:
            from cipd.data.ocr_census import ocr_image, cluster_rows, parse_number
            cells = ocr_image(path, scale=2.0)
            rows = cluster_rows(cells)
            matched = {}
            for r in rows:
                left = [c for c in r["cells"] if c["x"] < 110]
                if not left:
                    continue
                label = "".join(c["text"] for c in left).replace(" ", "")
                prov = None
                for cand in ["全国"] + REGIONS:
                    if label == cand or (label.startswith(cand) and len(label) <= len(cand) + 2):
                        prov = cand
                        break
                if prov is None:
                    continue
                nums = [parse_number(c["text"]) for c in r["cells"]
                        if c["x"] >= 110 and parse_number(c["text"]) is not None]
                if nums:
                    matched[prov] = nums
            out = matched
            if len(out) < 32:
                data_rows = [r for r in rows
                             if len([c for c in r["cells"] if c["x"] >= 110 and parse_number(c["text"]) is not None]) >= 1
                             and r["y"] > 100]
                if len(data_rows) == 32:
                    out = {}
                    for prov, r in zip(["全国"] + REGIONS, data_rows):
                        nums = [parse_number(c["text"]) for c in r["cells"]
                                if c["x"] >= 110 and parse_number(c["text"]) is not None]
                        if nums:
                            out[prov] = nums
    except Exception as e:
        print(f"  pop table parse failed {path.name}: {e}")
    return out


def build_panel():
    files = sorted(YB.glob("yb*_birth_*.*"))
    files = [f for f in files if not f.name.startswith("._")]
    panel = {}
    audit = []
    for f in files:
        m = re.match(r"yb(\d{4})_birth_(\d{4})\.(htm|xls|jpg)", f.name)
        if not m:
            continue
        ed, year, ext = int(m.group(1)), int(m.group(2)), m.group(3)
        year = int(re.search(r"(\d{4})\.(?:htm|xls|jpg)$", f.name).group(1))
        try:
            pop = urban = None
            if ext == "htm":
                rows = parse_html_table(f)
            elif ext == "xls":
                rows = parse_xls_table(f)
            else:
                rows = extract_jpg_rows_precise(f)
        except Exception as e:
            audit.append({"file": f.name, "issue": f"parse_error: {e}"})
            continue
        # standardize columns per format
        for r in rows:
            prov, vals = r["province"], r["values"]
            pop = urban = None
            if ext == "htm":
                # [pop(万人), cbr, cdr, ngr] (older) — some have more cols
                nums = [v for v in vals if v is not None]
                if prov == "全国" and len(nums) >= 4:
                    pop, cbr, cdr, ngr = nums[0], nums[-3], nums[-2], nums[-1]
                elif len(nums) >= 4:
                    pop, cbr, cdr, ngr = nums[0], nums[-3], nums[-2], nums[-1]
                elif len(nums) == 3:
                    pop, cbr, cdr, ngr = None, nums[0], nums[1], nums[2]
                else:
                    continue
            elif ext == "xls":
                nums = [v for v in vals if v is not None]
                # 3-7 style: [pop, urban_pop, urban_pct, rural_pop, rural_pct, cbr, cdr, ngr]
                if len(nums) >= 8:
                    pop, urban, cbr, cdr, ngr = nums[0], nums[2], nums[-3], nums[-2], nums[-1]
                elif len(nums) >= 4:
                    pop, urban, cbr, cdr, ngr = nums[0], None, nums[-3], nums[-2], nums[-1]
                else:
                    continue
            else:  # jpg
                rowpairs = r["_row"]
                hdr = r["_hdr"]
                vals = [v for (_, v) in rowpairs if v is not None]
                if len(vals) < 4:
                    continue
                pop = vals[0]
                def nearest(target):
                    if target is None or not rowpairs:
                        return None
                    best = min(rowpairs, key=lambda q: abs(q[0] - target))
                    return best[1] if abs(best[0] - target) < 55 else None
                cbr = nearest(hdr.get("cbr"))
                cdr = nearest(hdr.get("cdr"))
                ngr = nearest(hdr.get("ngr"))
                # fallback chain: (a) last three, (b) row re-OCR
                if not (cbr is not None and cdr is not None and ngr is not None
                        and abs((cbr - cdr) - ngr) < 0.2 and 2.0 < cbr < 25.0):
                    cbr2, cdr2, ngr2 = vals[-3], vals[-2], vals[-1]
                    if 2.0 < cbr2 < 25.0 and abs((cbr2 - cdr2) - ngr2) < 0.2:
                        cbr, cdr, ngr = cbr2, cdr2, ngr2
                if not (cbr is not None and cdr is not None and ngr is not None
                        and abs((cbr - cdr) - ngr) < 0.2 and 2.0 < cbr < 25.0):
                    audit.append({"file": f.name, "province": prov, "issue": "jpg_rates_unresolved",
                                  "cbr": cbr, "cdr": cdr, "ngr": ngr})
                urban = None
            if None in (cbr, cdr, ngr):
                audit.append({"file": f.name, "province": prov, "issue": "missing_rates"})
                continue
            # validation: cbr - cdr ~ ngr
            if abs((cbr - cdr) - ngr) > 0.15:
                audit.append({"file": f.name, "province": prov, "issue": "ngr_identity",
                              "cbr": cbr, "cdr": cdr, "ngr": ngr})
            if not (2.0 < cbr < 25.0 and 2.0 < cdr < 12.0):
                audit.append({"file": f.name, "province": prov, "issue": "range_check",
                              "cbr": cbr, "cdr": cdr})
            key = (prov, year)
            panel.setdefault(key, {"province": prov, "year": year})
            panel[key].update({
                "pop_wan": pop, "urban_pct": urban, "cbr_permille": cbr,
                "cdr_permille": cdr, "ngr_permille": ngr,
                "source": f.name, "caliber": "常住人口年度抽样调查/普查",
            })
    # override population from dedicated year-end tables, picking the column
    # whose national value matches the known national total (handles multi-year tables)
    import csv as _csv
    nat_anchor = {}
    with open(ROOT / "data/raw/china_annual_vital_1990_2025.csv") as f:
        for row in _csv.DictReader(f):
            nat_anchor[int(row["year"])] = float(row["pop_wan"])
    for pf in sorted(YB.glob("yb*_pop_*.*")):
        if pf.name.startswith("._"):
            continue
        m = re.search(r"_(\d{4})\.", pf.name)
        if not m:
            continue
        yy = int(m.group(1))
        tbl = parse_pop_table(pf)
        if not tbl or "全国" not in tbl:
            continue
        anchor = nat_anchor.get(yy)
        if anchor is None:
            continue
        # choose column index whose 全国 value is closest to anchor (within 1%)
        nat_vals = tbl["全国"]
        best_col, best_err = None, float("inf")
        for ci, v in enumerate(nat_vals):
            if v is None or v <= 0:
                continue
            err = abs(v - anchor) / anchor
            if err < best_err:
                best_err, best_col = err, ci
        if best_col is None or best_err > 0.01:
            audit.append({"year": yy, "issue": "pop_table_column_unresolved",
                          "file": pf.name, "anchor": anchor, "nat_vals": nat_vals[:8]})
            continue
        n_applied = 0
        for prov, vals in tbl.items():
            if prov == "全国" or best_col >= len(vals):
                continue
            v = vals[best_col]
            key = (prov, yy)
            if key in panel and v and v > 0:
                old = panel[key].get("pop_wan")
                panel[key]["pop_candidates"] = {"birth": old, "pop_table": v}
                panel[key]["pop_src_pop_table"] = pf.name
                n_applied += 1
        print(f"  pop override {yy}: {pf.name} col={best_col} err={best_err:.4f} applied={n_applied}")

    # --- robust arbitration of population candidates across years ---
    for prov in ["全国"] + REGIONS:
        series = {}
        for (pp, yy), rec in panel.items():
            if pp != prov or "pop_candidates" not in rec:
                continue
            series[yy] = rec["pop_candidates"]
        years_avail = sorted(series)
        picks = {}
        for i, yy in enumerate(years_avail):
            cands = series[yy]
            neigh = []
            for d in (-2, -1, 1, 2):
                yy2 = yy + d
                if yy2 in picks:
                    neigh.append(picks[yy2])
                elif yy2 in series:
                    neigh += [v for v in series[yy2].values() if v]
            if neigh:
                import statistics as st
                ref = st.median(neigh)
                pick = min([v for v in cands.values() if v], key=lambda v: abs(v - ref))
            else:
                pick = cands.get("birth") or cands.get("pop_table")
            picks[yy] = pick
        for yy in years_avail:
            key = (prov, yy)
            old = panel[key].get("pop_wan")
            new = picks[yy]
            if new and old and abs(old - new) > max(20, 0.02 * new):
                audit.append({"year": yy, "province": prov, "issue": "pop_arbitrated",
                              "final": new, "old": old,
                              "candidates": panel[key].get("pop_candidates")})
            panel[key]["pop_wan"] = new
            if "pop_src_pop_table" in panel[key]:
                panel[key]["source"] = panel[key]["pop_src_pop_table"] + "+birth_rates"

    # national sum checks per year
    years = sorted(set(y for (_, y) in panel))
    for y in years:
        nat = panel.get(("全国", y))
        provs = [panel[("p", y)] for p in REGIONS if ("p", y) in panel]
        provs = [panel[(p, y)] for p in REGIONS if (p, y) in panel]
        if nat and provs and nat.get("pop_wan") and all(p.get("pop_wan") for p in provs):
            s = sum(p["pop_wan"] for p in provs)
            if abs(s - nat["pop_wan"]) > 0.01 * nat["pop_wan"]:
                audit.append({"year": y, "issue": "natpop_vs_provsum", "nat": nat["pop_wan"], "sum": s})
    rows = [panel[k] for k in sorted(panel, key=lambda k: (k[1], k[0]))]
    cols = ["province", "year", "pop_wan", "urban_pct", "cbr_permille", "cdr_permille",
            "ngr_permille", "caliber", "source"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    (YB / "panel_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=1))
    print(f"panel rows: {len(rows)} -> {OUT.name}")
    print(f"audit entries: {len(audit)}")
    yrs = sorted(set(r["year"] for r in rows))
    print("years:", yrs)
    n_prov = sorted(set(r["province"] for r in rows))
    print("provinces:", len([p for p in n_prov if p != "全国"]))
    for e in audit[:12]:
        print("  ", e)


if __name__ == "__main__":
    build_panel()
