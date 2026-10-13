#!/usr/bin/env python3
"""Assemble one workbook (and a folder of CSVs) per Supplementary Table from
docs/SupplementaryTable_Map.csv. Values are never edited; text is stored as text.

usage: 11_assemble_supplement.py --map MAP --run RUN_DIR [--repo REPO] [--upstream UP]
                                 [--handcurated DIR] --out OUT_DIR
"""
import argparse, csv, hashlib, os, shutil, sys
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font

ap = argparse.ArgumentParser()
ap.add_argument("--map", required=True); ap.add_argument("--run", required=True)
ap.add_argument("--repo"); ap.add_argument("--upstream"); ap.add_argument("--handcurated")
ap.add_argument("--out", required=True)
a = ap.parse_args()
roots = {"run": a.run, "repo": a.repo, "upstream": a.upstream, "handcurated": a.handcurated}
m = pd.read_csv(a.map, dtype=str, keep_default_na=False)
m = m[m.include.str.upper() == "Y"]
missing = []
def md5(p):
    h = hashlib.md5(); h.update(open(p, "rb").read()); return h.hexdigest()

def read(path):
    sep = "\t" if path.endswith(".tsv") else ","
    return pd.read_csv(path, sep=sep, na_values=["NA", ""], keep_default_na=False, low_memory=False)

for t, g in m.groupby("table", sort=False):
    tdir = os.path.join(a.out, f"Supplementary_Table_{int(t):02d}")
    os.makedirs(os.path.join(tdir, "csv"), exist_ok=True)
    title = g.table_title.iloc[0]
    wb = Workbook(); ws0 = wb.active; ws0.title = "README"
    ws0.append([f"Supplementary Table {t}. {title}"]); ws0["A1"].font = Font(bold=True)
    ws0.append([]); ws0.append(["sheet", "source file", "rows", "columns", "md5 of source", "note"])
    for c in ws0[3]: c.font = Font(bold=True)
    for _, r in g.iterrows():
        root = roots.get(r.root)
        src = os.path.join(root, r.source_file) if root else None
        if not src or not os.path.exists(src):
            missing.append((t, r.sheet, r.source_file)); continue
        if src.endswith(".xlsx"):          # hand-curated workbook: copy unchanged
            shutil.copy(src, os.path.join(tdir, f"Supplementary_Table_{int(t):02d}.xlsx"))
            ws0.append([r.sheet, r.source_file, "", "", md5(src), r.note + " (workbook is the source itself)"])
            wb = None; continue
        df = read(src)
        df.to_csv(os.path.join(tdir, "csv", f"{r.sheet}.csv"), index=False)
        ws = wb.create_sheet(r.sheet[:31]); ws.append(list(df.columns))
        for row in df.itertuples(index=False):
            ws.append([None if pd.isna(v) else (v.item() if hasattr(v, "item") else v) for v in row])
        ws.freeze_panes = "A2"
        ws0.append([r.sheet, r.source_file, len(df), df.shape[1], md5(src), r.note])
    if wb is not None:
        for col, w in zip("ABCDEF", (30, 70, 8, 9, 34, 60)): ws0.column_dimensions[col].width = w
        wb.save(os.path.join(tdir, f"Supplementary_Table_{int(t):02d}.xlsx"))
for r_, d_, f_ in os.walk(a.out, topdown=False):
    if not os.listdir(r_): os.rmdir(r_)
if missing:
    print("MISSING SOURCES:"); [print("  ", x) for x in missing]; sys.exit(1)
print("assembled", m.table.nunique(), "tables,", len(m), "sheets")
