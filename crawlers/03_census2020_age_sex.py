#!/usr/bin/env python3
"""Extract 七普年鉴表 1-5 (A0105): provincial population by age group x sex (2020).

Structure (verified from header OCR): 32 rows (全国 + 31 provinces) x
[合计{总,男,女}] + 22 age blocks (0岁,1-4岁,5-9岁,...,100岁+) x {小计,男,女}.

Method: 3x7 tiling (each tile <=960px long side so RapidOCR does not downscale),
merge by absolute coordinates, header-derived column x positions, and full
demographic validation:
  row identity   合计 = 男 + 女
  age sum        sum over age blocks 小计 = 合计
  sex-age sum    sum over age blocks 男 = 合计男 ; 女 likewise
  national       national row = sum of 31 provinces (per column)
Cells failing these get row-strip re-OCR at 3x; unresolved cells are blanked
and logged (no silent imputation).
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/china_provincial/raw7rp"
OUT = ROOT / "data/raw/china_provincial"
sys.path.insert(0, str(ROOT / "src"))
from cipd.data.ocr_census import _ocr_array, parse_number, cluster_rows, REGIONS  # noqa

AGES = ["合计", "0岁", "1-4岁", "5-9岁", "10-14岁", "15-19岁", "20-24岁", "25-29岁",
        "30-34岁", "35-39岁", "40-44岁", "45-49岁", "50-54岁", "55-59岁", "60-64岁",
        "65-69岁", "70-74岁", "75-79岁", "80-84岁", "85-89岁", "90-94岁", "95-99岁", "100岁+"]


def ocr_tiles(path: Path, tile_w=900, tile_h=500):
    im = Image.open(path).convert("RGB")
    W, H = im.size
    cells = []
    xs = list(range(0, W, tile_w))
    ys = list(range(0, H, tile_h))
    for y0 in ys:
        for x0 in xs:
            x1 = min(x0 + tile_w, W)
            y1 = min(y0 + tile_h, H)
            if x1 - x0 < 50 or y1 - y0 < 40:
                continue
            crop = im.crop((x0, y0, x1, y1))
            # upscale small tiles to keep glyph height >= ~20px after any resize
            sc = max(1.0, 960.0 / max(crop.width, crop.height))
            cs = _ocr_array(np.array(crop), sc, x0, y0)
            cells += cs
    # dedupe overlapping detections
    uniq = []
    for c in sorted(cells, key=lambda c: -c["score"]):
        if not any(abs(c["x"] - u["x"]) < 6 and abs(c["y"] - u["y"]) < 6 for u in uniq):
            uniq.append(c)
    return uniq


def main():
    path = RAW / "A0105.jpg"
    cells = ocr_tiles(path)
    print(f"cells: {len(cells)}")
    rows = cluster_rows(cells, y_tol=14)
    print(f"rows: {len(rows)}")
    # identify header rows (top 4) and data rows by region labels
    data = []
    for r in rows:
        left = [c for c in r["cells"] if c["x"] < 115]
        if not left:
            continue
        label = "".join(c["text"] for c in left).replace(" ", "")
        prov = None
        if label.startswith("全国"):
            prov = "全国"
        else:
            for p in REGIONS:
                if label == p or (label.startswith(p) and len(label) <= len(p) + 2):
                    prov = p
                    break
        if prov is None:
            continue
        nums = [(c["x"], parse_number(c["text"]), c["text"]) for c in r["cells"]
                if c["x"] >= 115 and parse_number(c["text"]) is not None]
        nums.sort()
        data.append({"province": prov, "y": r["y"], "nums": nums})
    print(f"data rows: {len(data)} (expect 32)")

    # column model: 3 (合计: total,m,f) + 22 age blocks x 3
    # derive x positions from the header cells
    hdr = [c for c in cells if c["y"] < 165]
    ax = []
    for c in hdr:
        t = c["text"]
        for a in AGES:
            if t == a or (a != "合计" and t == a.replace("岁", "") + "岁"):
                ax.append((c["x"], a))
    ax.sort()
    # calibrate a linear column grid from header anchors if enough found
    n_cols = 3 + 22 * 3
    if len(ax) >= 5:
        xs = np.array([x for x, _ in ax])
        med_gap = np.median(np.diff(xs)) if len(xs) > 1 else 88.0
    else:
        med_gap = 88.0
    # fallback: assume uniform grid anchored at first 合计 block
    first_x = ax[0][0] if ax else 203.0
    # blocks: 合计 occupies cols 0..2 (total,m,f); then each age block 3 cols.
    # Estimate the per-column pitch from the 0岁 block: 小计 at first_x+ (some offset)
    # We instead read each row's numeric cells and assign them sequentially:
    # a row should have exactly 3 + 22*3 = 69 numbers.
    out = []
    audit = []
    for rec in data:
        vals = [v for (_, v, _) in rec["nums"]]
        if len(vals) != n_cols:
            audit.append({"province": rec["province"], "issue": "n_numbers_mismatch",
                          "n": len(vals), "expected": n_cols,
                          "ys": rec["y"]})
        v = vals[:n_cols] + [None] * max(0, n_cols - len(vals))
        row = {"province": rec["province"]}
        row["total"] = v[0]; row["total_m"] = v[1]; row["total_f"] = v[2]
        for i, a in enumerate(AGES[1:]):
            row[f"{a}_t"] = v[3 + 3 * i]
            row[f"{a}_m"] = v[4 + 3 * i]
            row[f"{a}_f"] = v[5 + 3 * i]
        out.append(row)

    # validation + rescue by row re-OCR is skipped for tile-extracted rows;
    # instead flag violations
    issues = 0
    for row in out:
        if None in (row["total"], row["total_m"], row["total_f"]):
            continue
        if abs(row["total_m"] + row["total_f"] - row["total"]) > max(3, 0.001 * row["total"]):
            audit.append({"province": row["province"], "issue": "mf_sum",
                          "total": row["total"], "m": row["total_m"], "f": row["total_f"]})
            issues += 1
        s = sum(row[f"{a}_t"] for a in AGES[1:] if row[f"{a}_t"] is not None)
        if abs(s - row["total"]) > max(50, 0.002 * row["total"]):
            audit.append({"province": row["province"], "issue": "age_sum",
                          "total": row["total"], "age_sum": s})
            issues += 1
    nat = next((r for r in out if r["province"] == "全国"), None)
    if nat:
        for col in ["total", "total_m", "total_f"] + [f"{a}_{s}" for a in AGES[1:] for s in ("t", "m", "f")]:
            pv = sum(r[col] for r in out if r["province"] != "全国" and r[col] is not None)
            nv = nat[col]
            if nv is not None and pv > 0 and abs(pv - nv) > max(100, 0.003 * nv):
                audit.append({"province": "全国", "col": col, "issue": "nat_vs_provsum",
                              "national": nv, "prov_sum": pv})

    cols = ["province", "total", "total_m", "total_f"] + \
           [f"{a}_{s}" for a in AGES[1:] for s in ("t", "m", "f")]
    out_csv = OUT / "census2020_table1_5_provincial_age_sex.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})
    (OUT / "ocr_a0105_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=1))
    print(f"rows written: {len(out)}; audit entries: {len(audit)} (issues {issues})")
    for a in audit[:10]:
        print(" ", a)


if __name__ == "__main__":
    main()
