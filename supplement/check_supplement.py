#!/usr/bin/env python3
"""Consistency checks across the assembled Supplementary Tables.
usage: 12_check_supplement.py --out ASSEMBLED_DIR --run RUN_DIR [--map MAP]
Exit code 1 if any check fails. Expected values are numbers quoted in the Supplementary Tables text."""
import argparse, glob, os, sys
import pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--run", required=True)
ap.add_argument("--map"); a = ap.parse_args()
def S(t, sheet):
    return pd.read_csv(os.path.join(a.out, f"Supplementary_Table_{t:02d}", "csv", f"{sheet}.csv"), low_memory=False)
fails = 0
def chk(name, got, want):
    global fails
    ok = got == want; fails += (not ok)
    print(("PASS " if ok else "FAIL ") + f"{name}: got {got}, expected {want}")

x = S(5, "Status_crossings")
chk("T5 genes in Figure 7", x.symbol.nunique(), 58)
chk("T5 rows == Concordance age-axis rows are a subset", len(S(7, "Age_axis_crossings")) <= len(x), True)
chk("T6 unique panel genes x genotype", len(S(6, "Trajectories").drop_duplicates(["genotype", "symbol"])) <= 2 * 263, True)
i = S(8, "Interaction_tests"); d = i.drop_duplicates(["contrast", "gene_id"])
chk("T8 distinct tests", len(d), 483)
chk("T8 differ between ages", int((d.outcome == "genotype effect differs between ages").sum()), 26)
chk("T8 equivalent", int((d.outcome == "no difference detected (bounded)").sum()), 7)
cg = S(11, "Categorical_genes"); fs = S(9, "Fold_sweep").set_index("min_fold")
chk("T9 fold-2 total == T11 rows", int(fs.loc[2, "categorical_total"]), len(cg))
chk("T9 fold-2 distinct == T11 distinct", int(fs.loc[2, "distinct_calls"]), len(cg.drop_duplicates(["contrast", "symbol"])))
chk("T9 fold-2 unique genes == T11 genes", int(fs.loc[2, "unique_genes"]), cg.symbol.nunique())
chk("T9 gate-cost rows", len(S(9, "Gate_cost")), 8)
pan = sum(len(S(15, s)) for s in ("Axon_guidance_panel", "Myelination_panel", "Vascular_lymphatic_panel"))
chk("T15 panel sizes 92+25+177", pan, 294)
ec = S(10, "Expression_calls")
chk("T10 genes in expression calls", ec.gene_id.nunique(), 263)
chk("T12 trigeminal libraries", len(S(12, "Sample_table")), 14)
chk("T14 panel genes in WT concordance", len(S(14, "WT_status_concordance")), 263)
chk("T14 null draws", len(S(14, "WT_null_distribution")), 2000)
chk("T17 program-level tests", len(S(17, "Program_level")), 15)
def xl(t): return glob.glob(os.path.join(a.out, f"Supplementary_Table_{t:02d}", "*.xlsx"))[0]
t19 = pd.read_excel(xl(18), sheet_name="Libraries", header=3)
t19 = t19[t19["Study"].astype(str).str.contains("GSE")]
chk("T18 libraries (48 + 6)", len(t19), 54)
chk("T18 sheets", sorted(pd.ExcelFile(xl(18)).sheet_names), ["Libraries", "Raw file checksums", "Software and reference", "Strandedness"])
t13 = pd.read_excel(xl(13), header=None)
ids = set(t13.iloc[:, 0].astype(str)); lab = set(t19.Library.astype(str))
chk("T13 IDs cover the 48 T18 study IDs", len(set(t19[t19.Study.str.contains("this study", case=False)].Library) - ids), 0)
chk("T18 workbook says Supplementary, not Supplemental",
    any("Supplemental" in str(v) for ws in pd.read_excel(xl(18), sheet_name=None, header=None).values() for v in ws.values.ravel()), False)
t13h = t13.copy(); cols = [str(c) for c in t13h.iloc[1]]; cols[0] = "ID"; t13h.columns = cols; t13h = t13h.iloc[2:].set_index("ID")
j = t13h.join(t19.set_index("Library"), how="inner"); num = lambda c: pd.to_numeric(j[c])
chk("T13 vs T18 agree: libraries compared", len(j), 48)
for x, sx, y, sy, tol in [("read_pairs_M", 1e6, "Read pairs sequenced", 1, 6000), ("trim_both_surviving_pct", 1, "Pairs retained after trimming (%)", 100, 0.006),
        ("hisat2_overall_alignment_pct", 1, "Overall alignment rate (%)", 100, 0.006), ("assigned_fragments_M", 1e6, "Fragments assigned to genes", 1, 6000),
        ("library_size_fragments", 1, "Fragments in matrix, reprocessed (matrix of record)", 1, 0)]:
    chk(f"T13 {x} == T18 (rounding tol {tol})", int(((num(x) * sx - num(y) * sy).abs() > tol).sum()), 0)
chk("T6 steps rows", len(S(6, "Steps_per_gene")), 1194)
chk("T16 DE rows == T8-style panel family", len(S(16, "DE_results")), 1213)
chk("T20 libraries in depth table", len(S(20, "Depth_per_sample")), 48)
chk("T19 sheets", len(glob.glob(os.path.join(a.out, "Supplementary_Table_19", "csv", "*.csv"))), 7)
# every run CSV should be assigned or listed as deliberately unassigned
if a.map:
    m = pd.read_csv(a.map, dtype=str); used = set(m.source_file[m.root == "run"])
    allf = {os.path.relpath(p, a.run) for p in glob.glob(os.path.join(a.run, "**/*.csv"), recursive=True)}
    un = sorted(allf - set(m.source_file[m.root == "run"]))
    chk("run CSVs not named anywhere in the map (include=N rows count as named)", un, [])
print("\nFAILED:", fails); sys.exit(1 if fails else 0)
