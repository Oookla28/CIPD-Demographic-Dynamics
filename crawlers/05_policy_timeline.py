#!/usr/bin/env python3
"""Build the provincial fertility-policy timeline for 2013-2026.

Provincial 单独二孩 implementation dates (2014) come from the zh-Wikipedia
table on 单独两孩, which itself cites each province's 人大常委会决议 /
条例修正案. Spot-checked against widely reported official dates
(浙江 2014-01-17 first; 北京 2014-02-21; 山东 2014-05-30; 西藏 2014-09-17).
National policy rows extend the existing provenance/policy_codebook.csv.

Output: data/raw/china_provincial/provincial_policy_timeline.csv
        (one row per province x policy episode with exact timing & quality flag)
"""
from __future__ import annotations

import csv
import re
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/china_provincial/provincial_policy_timeline.csv"

# province -> (month, day, quality) for 单独二孩 2014 implementation
DANDUI_2014 = {
    "浙江": (1, 17, "exact"), "江西": (1, 18, "exact"), "安徽": (1, 23, "exact"),
    "天津": (2, 14, "exact"), "北京": (2, 21, "exact"), "上海": (3, 1, "exact"),
    "广西": (3, 1, "exact"), "陕西": (3, 1, "exact"), "四川": (3, 20, "exact"),
    "重庆": (3, 26, "exact"), "甘肃": (3, 26, "exact"), "湖北": (3, 27, "exact"),
    "广东": (3, 27, "exact"), "青海": (3, 27, "exact"), "辽宁": (3, 28, "exact"),
    "吉林": (3, 28, "exact"), "江苏": (3, 28, "exact"), "湖南": (3, 28, "exact"),
    "云南": (3, 28, "exact"), "福建": (3, 29, "exact"), "内蒙古": (3, 31, "exact"),
    "黑龙江": (4, 22, "exact"), "贵州": (5, 18, "exact"), "宁夏": (5, 28, "exact"),
    "山西": (5, 29, "exact"), "河北": (5, 30, "exact"), "山东": (5, 30, "exact"),
    "海南": (6, 1, "exact"), "河南": (6, 3, "exact"), "西藏": (9, 17, "exact"),
    "新疆": (5, 31, "approximate"),  # 不详（部分地区不晚于5月）→ conservative mid-year
}

# 长表数据核对后的关键全国政策节点（承接 provenance/policy_codebook.csv）
NATIONAL_EVENTS = [
    ("DANDUI", "单独二孩（一方独生子女夫妇可生二孩）", "2013-11-15", "2014-01-01"),
    ("U2", "全面二孩（人口与计划生育法修正案）", "2015-10-29", "2016-01-01"),
    ("T3", "三孩政策及配套支持措施", "2021-05-31", "2021-05-31"),
    ("SUBSIDY", "国家育儿补贴（3600元/年/孩，中央预算900亿）", "2025-07-28", "2025-01-01"),
]


def main():
    rows = []
    # provincial 单独二孩 staggered adoption
    for prov, (m, d, q) in DANDUI_2014.items():
        eff = f"2014-{m:02d}-{d:02d}"
        rows.append({
            "province": prov, "policy": "DANDUI",
            "policy_desc": "单独二孩落地（省人大常委会决议/条例修正）",
            "announcement_date": "2013-12-21", "effective_date": eff,
            "date_quality": q,
            "treatment_year_did": 2014,
            "source": "zh.wikipedia 单独两孩 各省实施时间表(引各省人大决议) https://zh.wikipedia.org/zh-cn/单独两孩",
        })
    # uniform national policies replicated across provinces
    for code, desc, ann, eff in NATIONAL_EVENTS[1:]:
        for prov in DANDUI_2014:
            rows.append({
                "province": prov, "policy": code, "policy_desc": desc,
                "announcement_date": ann, "effective_date": eff,
                "date_quality": "exact",
                "treatment_year_did": int(eff[:4]),
                "source": "provenance/policy_codebook.csv (official)",
            })
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"written {len(rows)} rows -> {OUT}")


if __name__ == "__main__":
    main()
