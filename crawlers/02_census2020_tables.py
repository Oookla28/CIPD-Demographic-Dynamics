#!/usr/bin/env python3
"""阶段 02：七普年鉴（2020）长表生育表格 OCR 提取。

输入（已下载原图，SHA-256 锁定于 provenance/hashes.sha256）:
  data/raw/china_provincial/raw7rp/B0601.jpg  —— 表 6-1 各地区分性别、孩次的出生人口
  data/raw/china_provincial/raw7rp/B0603.jpg  —— 表 6-3 全国育龄妇女分年龄、孩次的生育状况
  data/raw/china_provincial/raw7rp/B0604.jpg  —— 表 6-4 各地区育龄妇女年龄别生育率
  data/raw/china_provincial/raw7rp/B0605.jpg  —— 表 6-5 各地区按活产子女数分的 15-64 岁妇女人数
  data/raw/china_provincial/raw7rp/B0612.jpg  —— 表 6-12 全国按年龄分的平均活产/存活子女数

输出:
  data/raw/china_provincial/census2020_table6_1_provincial_parity_births.csv
  data/raw/china_provincial/census2020_table6_3_age_parity_national.csv
  data/raw/china_provincial/census2020_table6_4_provincial_asfr.csv
  data/raw/china_provincial/census2020_table6_5_provincial_ceb.csv
  data/raw/china_provincial/census2020_table6_12_age_mean_ceb.csv
  data/raw/china_provincial/ocr_audit.json         （逐表校验/救援/推导记录）

校验规则（失败即记录，不做静默插补）:
  6-1: 每省每组 小计=男+女；合计=Σ孩次；全国=Σ31省；性别比=男/女×100
  6-3: 单岁求和=5岁组=总计；妇女人数×生育率=出生数；孩次分项求和=出生数
  6-4: TFR=5×ΣASFR
  6-5: Σ活产组=妇女人数；全国=Σ省
  6-12: 男+女=合计；平均活产=活产总数/妇女人数
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cipd.data.ocr_census import (extract_b0601, extract_b0603, extract_b0604,
                                  extract_b0605, extract_b0612, AUDIT)

if __name__ == "__main__":
    extract_b0601()
    extract_b0603()
    extract_b0604()
    extract_b0605()
    extract_b0612()
    import json
    (ROOT / "data/raw/china_provincial/ocr_audit.json").write_text(
        json.dumps(AUDIT, ensure_ascii=False, indent=1))
    print("[stage 02] census 2020 tables done")
