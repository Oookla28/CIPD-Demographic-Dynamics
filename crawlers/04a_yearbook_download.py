#!/usr/bin/env python3
"""Crawl China Statistical Yearbook provincial population tables (2004-2024 data years).

For each yearbook edition, downloads:
  - provincial birth/death/natural-growth rates table (title contains 出生率)
  - provincial year-end population table (年末人口数)
  - provincial urbanization table (城镇人口比重) when available
Handles both HTML (older editions) and JPG (2019+) table formats.
Retries + polite delays; logs every fetch to provenance/download-style log.
"""
from __future__ import annotations

import re
import sys
import urllib.parse
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/china_provincial/yearbooks"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://www.stats.gov.cn/sj/ndsj/{ed}/"
EDITIONS = {
    2004: "2005", 2005: "2006", 2006: "2007", 2007: "2008", 2008: "2009",
    2009: "2010", 2010: "2011", 2011: "2012", 2012: "2013", 2013: "2014",
    2014: "2015", 2015: "2016", 2016: "2017", 2017: "2018", 2018: "2019",
    2019: "2020", 2020: "2021", 2021: "2022", 2022: "2023", 2023: "2024",
    2024: "2025",
}

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0 Safari/537.36")


def fetch(url: str, tries: int = 4, sleep: float = 3.0) -> bytes | None:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except Exception as e:
            print(f"  retry {i+1} for {url}: {e}", flush=True)
            time.sleep(sleep * (i + 1))
    return None


def get_toc(edition: str) -> list[tuple[str, str]]:
    for left in ("left.htm", "left.html"):
        data = fetch(BASE.format(ed=edition) + left)
        if data:
            try:
                txt = data.decode("gbk", errors="replace")
            except Exception:
                txt = data.decode("utf-8", errors="replace")
            pairs = re.findall(r'href=["\']([^"\']*)["\'][^>]*>([^<]+)<', txt)
            return [(u, t.strip()) for u, t in pairs if t.strip()]
    return []


def main():
    log = []
    targets = {
        "birth": ("出生率", ("死亡率",)),
        "pop": ("年末人口数", ()),
        "urban": ("城镇人口比重", ()),
    }
    for data_year, ed in EDITIONS.items():
        print(f"=== edition {ed} (data year {data_year})", flush=True)
        toc = get_toc(ed)
        if not toc:
            log.append({"edition": ed, "data_year": data_year, "status": "toc_failed"})
            continue
        # pick best-matching table for each target
        for key, (must, extras) in targets.items():
            cands = []
            for u, t in toc:
                if must in t and "地区" in t:
                    # skip tables that are clearly subpopulations
                    if key != "urban" and any(x in t for x in ("民族", "受教育", "户", "性别", "婚姻", "文盲", "城乡构成")):
                        continue
                    if key == "urban" and "城乡构成" in t:
                        continue
                    score = sum(1 for x in extras if x in t)
                    cands.append((score, u, t))
            if not cands and key == "birth":
                # fall back: 城乡构成和出生率 table
                for u, t in toc:
                    if "出生率" in t and "地区" in t:
                        cands.append((0, u, t))
            if not cands:
                log.append({"edition": ed, "data_year": data_year, "key": key, "status": "no_table"})
                continue
            cands.sort(key=lambda c: -c[0])
            _, url, title = cands[0]
            url = url.replace("\\", "/").lstrip("/")
            full = BASE.format(ed=ed) + urllib.parse.quote(url)
            ext = url.split(".")[-1].lower()
            dest = OUT / f"yb{ed}_{key}_{data_year}.{ext}"
            if dest.exists() and dest.stat().st_size > 1000:
                print(f"  [{key}] cached: {dest.name}", flush=True)
                log.append({"edition": ed, "data_year": data_year, "key": key,
                            "url": full, "title": title, "file": str(dest), "status": "cached"})
                continue
            data = fetch(full)
            if data is None or len(data) < 1000:
                log.append({"edition": ed, "data_year": data_year, "key": key,
                            "url": full, "title": title, "status": "download_failed"})
                print(f"  [{key}] FAILED: {full}", flush=True)
                continue
            dest.write_bytes(data)
            print(f"  [{key}] {title[:40]} -> {dest.name} ({len(data)}B)", flush=True)
            log.append({"edition": ed, "data_year": data_year, "key": key,
                        "url": full, "title": title, "file": str(dest), "status": "ok"})
            time.sleep(1.5)
        time.sleep(1.0)
    import json
    (OUT / "crawl_log.json").write_text(json.dumps(log, ensure_ascii=False, indent=1))
    n_ok = sum(1 for l in log if l.get("status") in ("ok", "cached"))
    print(f"done: {n_ok}/{len(log)} ok; log -> {OUT/'crawl_log.json'}")


if __name__ == "__main__":
    main()
