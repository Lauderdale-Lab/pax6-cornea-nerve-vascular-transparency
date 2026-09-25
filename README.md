# Nerve and vascular abnormalities precede loss of corneal transparency in *Pax6*-haploinsufficient mice — analysis pipeline

R code for the curated-panel RNA-seq analysis reported in:

> **Nerve and vascular abnormalities precede loss of corneal transparency in *Pax6*-haploinsufficient mice**
> Sneha K. Mohan, James D. Lauderdale
> Department of Cellular Biology, University of Georgia, Athens, GA 30602, USA

Repository: <https://github.com/Lauderdale-Lab/pax6-cornea-nerve-vascular-transparency> · Archived release: `<ZENODO_DOI>` · Contact: James D. Lauderdale, <jdlauder@uga.edu>

Licensed under the MIT License — see [`LICENSE`](LICENSE).

---

## What this pipeline does

Three curated gene panels — axon guidance (92 genes), myelination (25) and
vascular/lymphatic (177) — are examined in wild-type and *Pax6*<sup>Sey-Neu</sup>/+
mouse cornea at P3/P4, P15 and adult, with adult mutant corneas separated into
transparent and opaque.

The design answers three questions, and is built so that it can decline to
answer where the data do not support one:

1. **Which genes can be called expressed, not expressed, or neither.** Exact
   Poisson limits on counts pooled within a group; a call that changes when any
   single library is removed is reported as indeterminate rather than kept.
2. **Which genes differ, and by how much.** Only genes confidently expressed in
   both groups enter a differential model. Genes confidently expressed in one
   group and silent — or separated by non-overlapping intervals — in the other
   are reported as categorical differences instead, with no model fitted.
3. **What a negative result means.** Every non-significant comparison is
   reported as either *bounded*, with the excluded effect size stated, or
   *undetermined*. Undetermined is not evidence of no effect and is never
   reported as one.

A fourth step checks the first two against a conventional analysis. Every
status call is looked up in the same negative-binomial model fitted to every
panel gene without the expression-status gate, and a call the replicate-aware
model does not confirm is drawn as an **open symbol** in Figures 7 and 8B
rather than removed. See *Concordance* below.

Every library is also screened for tissue other than central cornea carried
at dissection (adherent iris, angle, lens, retina) and for its share of stroma
and epithelium (`01c`, Supplementary Table 6). The screen excludes nothing; it
is reported, and `07` lists any status call whose reads come mostly from a
flagged library.

## Running it

1. Clone the repository.
2. Download the study's count matrix and save it as
   `RNAseq_Quantification_Matrices/Combined_2025_2026/gene_counts_featureCounts_all48.txt`
   (see the README in that folder for the source and the MD5 checksum).
3. In R, from the root of the clone:

```r
setwd("<path to the clone>")
## Optional: the cross-dataset check (09) runs only when both are set.
## Sys.setenv does not survive an R restart.
Sys.setenv(PAX6_DUNCAN_COUNTS = "RNAseq_Quantification_Matrices/Duncan_GSE183742/gene_counts_featureCounts.txt",
           PAX6_DUNCAN_META   = "Metadata_and_Sample_Guide/Duncan_GSE183742/GSE183742_sample_metadata.csv")
source("scripts/run_all.R")
```

The working directory is the project root: inputs are read from it and outputs
written under it. The scripts are found in `scripts/` automatically; set
`PAX6_SCRIPT_DIR` only if they are kept somewhere else, and
`PAX6_PROJECT_ROOT` only if the data are not under the working directory.
`run_all.R` sets `PAX6_SCRIPT_DIR` for the rest of the R session, so restart R
before running a different copy of the pipeline.

### Reproducing the published results

A run from a clone should reproduce the published results exactly. Three
checks confirm it:

- **Inputs.** `RunProvenance.txt` records the MD5 of the count matrix and
  metadata it read. Both should match
  `docs/RunProvenance_v1.0.0_2026-09-24.txt`, the record of the run behind
  the manuscript. That file also lists the R and package versions used.
- **Panels.** The three curated panel master lists are committed in
  `Gene_Panels_and_Reference_Lists/master_lists/`. `00b` rebuilds each panel
  from its curated CSV and stops the run if the result differs from the
  committed list by a single gene, for example because of a
  different `org.Mm.eg.db` version.
- **Derived files.** The run rewrites `scripts/PanelPrimaryAssignment.R` and
  the three master lists. Afterwards, `git diff` should show no change except
  the "Written" date line in `PanelPrimaryAssignment.R`.

Numbers quoted in figure captions are regenerated in `Figure_CaptionNumbers.txt`
on every run and can be compared with the manuscript directly.

Scripts are numbered in execution order and share one R session by design:
`05` needs transcriptome-wide expression status, `06` needs the shared-gene
assignment, `07` needs the fitted models from `04` and `06`, and `09` uses the
marker identifiers and the per-library rule defined in `01c`, none of which is
written to disk. The runner stops at the first failure, reports any other
numbered script present but not run, records a step whose external inputs are
not set as skipped rather than failed, and writes `RunProvenance.txt` recording
the input files, every resolved parameter, the step timings and the full
`sessionInfo()`.

| Script | Purpose |
|---|---|
| `00_config.R` | every analysis parameter; computes nothing |
| `00b_panel_universes.R` | builds the three curated panels |
| `01_load.R` | metadata, counts, pool sex composition, panels |
| `01b_qc_library_depth.R` | library complexity at matched sequencing depth |
| `01c_contamination_indices.R` | off-target tissue and stromal/epithelial share per library (Supplementary Table 6); defines the rule `09` applies to GSE183742 |
| `02_shared_gene_assignment.R` | primary panel for genes on more than one panel; writes `PanelPrimaryAssignment.R` |
| `03_expression_status.R` | expression status and assessability |
| `04_differential.R` | differential expression, equivalence, omnibus, interactions |
| `05_program_level.R` | abundance-matched program-level permutation tests |
| `06_developmental.R` | developmental trajectories within genotype |
| `07_concordance.R` | status calls against a replicate-aware analysis of the same data; calls led by a flagged library; threshold sweep; QC figures including the sample PCA |
| `08_figures.R` | the manuscript's four RNA-seq figures (optional; needs `patchwork`) |
| `09_cross_dataset_replication.R` | the calls against GSE183742 (optional; runs only when `PAX6_DUNCAN_COUNTS` and `PAX6_DUNCAN_META` are set) |
| `run_all.R` | driver and provenance |

`07` precedes the figures because they read its verdicts. `09` is a check that
nothing consumes.

`PanelPrimaryAssignment.R` is generated by `02` and sourced by `00_config.R`; it
is committed so that every shared-gene assignment can be read without running
anything. Do not edit it by hand: change `PANEL_PRECEDENCE` or
`PANEL_PRIMARY_OVERRIDE` in `00_config.R` and re-run. A copy is also written
with each set of outputs.

Every parameter can be overridden by an environment variable of the same name
prefixed `PAX6_`, which is how sensitivity analyses are run without editing
code. `00_config.R` is the single place to read to see every choice the
analysis makes.

### Editing a curated panel

`00b_panel_universes.R` rebuilds each panel from its curated CSV and compares
it with the committed master list in `Gene_Panels_and_Reference_Lists/master_lists/`. Any difference stops the
run, because every denominator in the paper depends on the panel sizes. That
is the intended behaviour for an *unintended* change. For an intended edit,
run once with the guard off, which rewrites the master list, then return to
strict mode:

```r
Sys.setenv(PAX6_PANEL_STRICT = "FALSE")
source("scripts/run_all.R")
Sys.unsetenv("PAX6_PANEL_STRICT")
```

The rewritten master lists are then committed with the curated CSV change.


## Outputs

Written to `results/` in the project root, or to `results/<name>/` when
`PAX6_OUT_SUBDIR` is set to `<name>` (useful for keeping a sensitivity run
apart from the main one):

| File | Contents |
|---|---|
| `PanelUniverses_BuildReport.txt` | how each curated panel was built, and any difference from the saved master list |
| `QC_LibraryDepth_*.csv` / `.png` | detection at matched depth, depth series, and each library's share of its pooled group |
| `Contamination_Indices.csv` | off-target tissue indices, composition and flags, every library |
| `SuppTable6_Contamination.csv` | the published subset of the above |
| `SharedGeneAssignment.csv`, `PanelPrimaryAssignment.R` | the primary panel of each shared gene, and why |
| `ExpressionCalls.csv` | status, pooled counts and Poisson limits per group |
| `Assessability.csv` | what statement each gene supports, per comparison |
| `CategoricalGenes.csv` | genes on in one group and off in the other |
| `DE_Results.csv` | tested genes, bounded vs undetermined negatives |
| `InteractionResults.csv` | whether the genotype effect changes with age |
| `AdultOmnibusLRT.csv` | any difference among the three adult groups |
| `ProgramLevel.csv` | panel-level, abundance-matched |
| `Developmental_*.csv` | trajectories and status crossings |
| `Concordance/` | the concordance check, see below |
| `CrossDataset/` | the calls against GSE183742, the dissection-margin check, and GSE183742's libraries judged by the Supplementary Table 6 rule (only when `09` runs) |
| `RunProvenance.txt` | inputs, parameters, timings, session info |

Figures, from `08_figures.R`:

| File | Contents |
|---|---|
| `Figure7.pdf` / `.png` | genes whose expression status changes with age, within genotype; open symbols are crossings not confirmed by the replicate-aware model |
| `Figure8.pdf` / `.png` | program level (A) and genes crossing between expressed and silent (B); open symbols as in Figure 7 |
| `Figure9.pdf` / `.png` | quantitative change: axon guidance and myelination |
| `Figure10.pdf` / `.png` | quantitative change: vascular and lymphatic |
| `<figure>_SourceData.csv` | the plotted values, one file per figure; Figures 7 and 8 carry `confirmed_by_replicate_model` |
| `Figure_CaptionNumbers.txt` | every count quoted in a caption, including the open-symbol genes by name, regenerated each run |

### Concordance

`07_concordance.R` asks whether the status calls would survive an analysis a
reader already trusts. It takes the negative-binomial models that `04` and `06`
fitted, applies them to every panel gene without the expression-status gate,
recomputes adjusted P values over that larger family, and asks whether each
status crossing (Figure 7) and each categorical call (Figure 8B) is also
significant in the same direction. Where edgeR is installed, the genotype
contrasts are refitted with its quasi-likelihood pipeline on the identical
design as a second engine. It also re-derives the crossings over a grid of
expression thresholds and minimum folds, reports how much of each pooled count
its single largest library supplied, lists every status call whose higher
side comes more than half from one library flagged in `01c`, and writes three
routine diagnostics. One of these is a genome-wide sample PCA, drawn as
measured and with the batch effect removed for display. It is descriptive
quality control only: no principal component is interpreted, selected or
tested, and no result depends on it.

It changes no call and no number upstream, but `08_figures.R` reads its
verdicts: an unconfirmed call is drawn as an open symbol, and the caption file
names those genes.

| File (in `Concordance/`) | Contents |
|---|---|
| `Concordance_AgeAxis.csv` | every Figure 7 crossing with its ungated verdict, log2 fold change, adjusted P and top-library share |
| `Concordance_AgeAxis_StandardSignificant.csv` | ungated calls at the categorical fold that are not crossings, and the route the pipeline sent each down |
| `Concordance_GenotypeAxis.csv` | every gene by comparison, with assessability, ungated verdict and top-library share |
| `Concordance_GenotypeAxis_Summary.csv` | counts by comparison, route and verdict |
| `Concordance_GenotypeAxis_FamilyMoves.csv` | quantitatively tested genes whose `changed` call moves with the correction family |
| `Concordance_edgeR.csv` | the second engine's results (if edgeR is installed) |
| `Concordance_FlaggedLibraryCalls.csv` | status calls whose largest contributor on the higher side is a library flagged in `01c`, with its share |
| `Concordance_ThresholdStability.csv` | crossings and Figure 7 row bands over the threshold grid |
| `Concordance_AgeAxis.png`, `Concordance_GenotypeAxis.png` | bounded fold against ungated log2 fold change, discordant genes labelled |
| `Figure7_DevelopmentalDE_Nerve.*`, `Figure7_DevelopmentalDE_Vascular.*` | supplementary to Figure 7: the ungated within-genotype age contrasts as conventional heatmaps, one colour scale, Figure 7 crossings marked. Named after the figure they support (`PAX6_FIG_DEVELOPMENT`, default `Figure7`); the supplementary number is assigned at assembly |
| `QC_PValueHistograms.png`, `QC_SampleCorrelation.png` | P-value histograms for every fitted contrast; sample–sample Spearman correlation on the panel genes |
| `QC_SamplePCA.png` / `.pdf`, `_SourceData.csv`, `_VarianceExplained.csv` | sample PCA of the corneal libraries (top 500 variable genes, variance-stabilised), as measured and batch-removed; R² of each of PCs 1–5 with age, group, batch and pool sex composition |
| `Concordance_Summary.txt` | every number printed by the script, ending with the sentences for Methods |

### The supplementary combined comparison

Opaque mutant cornea differs from wild type in genotype **and** in transparency
at once, so that comparison cannot separate the two and is absent from the main
figures. To produce it as a supplementary figure:

```r
Sys.setenv(PAX6_FIG7_CONFOUNDED = "TRUE")
source("scripts/08_figures.R")
Sys.unsetenv("PAX6_FIG7_CONFOUNDED")
```

That run writes `Figure8_Confounded.pdf` / `.png`, its source data and its own
caption numbers, and touches none of the main figures. The name is derived
from the main figure it is a variant of (`PAX6_FIG_PROGRAM_MAIN`, default
`Figure8`), so renumbering the main figure renames the supplement with it; its
final supplementary number is assigned when the figures are assembled. It draws
the program-level figure only: adding the combined comparison changes the
heatmap colour limits, so heatmaps from that run would not be comparable with
Figures 9 and 10. (The switch is still named `PAX6_FIG7_CONFOUNDED` for
compatibility with earlier runs.)

### Counting conventions in the outputs

Most tables carry one row per gene per **panel**, so a gene curated onto two
panels appears under each with the same adjusted P value — that is deliberate,
so a reader consulting one panel's table sees every gene on it. Totals taken
**across** panels must therefore be de-duplicated, and the scripts print both
numbers wherever they differ: `03_expression_status.R` reports distinct gene by
contrast calls beside the panel-row count, and `06_developmental.R` does the
same for status crossings and trajectories. `primary_panel` is carried in the
developmental outputs so the de-duplication can be reproduced downstream.

## Change log

**Version 1.0.0 (2026-09-24) is the first public release.** The entries below
record changes made during development, before release, so that any output
produced earlier can be matched to the code that produced it.

**2026-09-25 — Runs from a clone; outputs moved to `results/`.** `run_all.R`
and `00_config.R` now find the scripts in `scripts/` when run from the
repository root, and the project root is the working directory unless
`PAX6_PROJECT_ROOT` is set (the machine-specific fallback paths were removed).
Outputs are written to `results/` (were
`Pax6_dtu_results/Module3_Standalone_Combined_2025_2026/Paper_Exports/ThreeState_Release/`),
and `PAX6_OUT_SUBDIR` now names a folder under `results/`. The three panel
master lists moved to `Gene_Panels_and_Reference_Lists/master_lists/` and are
renamed `AxonGuidance_Master_List.csv`, `Myelination_Master_List.csv` and
`VascularLymphatic_Master_List.csv` (were `GenesOfInterest_Master_List.csv`
in one folder per panel); they are committed, so a fresh clone is checked
against the published panels. In the code, `PANEL_DIRS` is now
`PANEL_MASTER_FILES`.
Two error messages and two comments naming renamed scripts were corrected. No
analysis code or number changed.

**2026-09-24 — Off-target tissue screen added to the pipeline; sample PCA.**
`01c_contamination_indices.R` is new and required. It reproduces
Supplementary Table 6 inside the pipeline: per-library indices for pigmented
tissue, retina, lens, angle and conjunctiva, and stromal and epithelial share,
each judged against the same-age corneal libraries excluding opaque mutants and
the library itself. Stroma-poor is now defined as both stromal controls (Kera,
Angptl7) more than 4-fold below that reference, replacing an absolute rule
(Kera below 100 CPM); on the current data both rules flag the same three
libraries. The marker list and every threshold moved to `00_config.R` and are
shared with the dissection-margin check in `09`, whose own list had drifted: it
no longer counts *Pmel* (expressed in every cornea) or *Chi3l1* (induced with
inflammation in opaque cornea), and its lens and retina sets are now `01c`'s.
`09` also judges each GSE183742 library by the Supplementary Table 6 rule.
`07` gains a list of status calls whose higher side is more than half one
flagged library, and a descriptive genome-wide sample PCA. The confounded
supplement from `08` is now named after its main figure (`Figure8_Confounded`,
was `FigureS3_Confounded`), and so are the developmental heatmaps from `07`
(`Figure7_DevelopmentalDE_*`, were `FigureS4a/b_DevelopmentalDE_*`). The study
title in every file now matches the manuscript. `run_all.R` now reports any stray numbered script,
not only those beginning with `0`. No analysis number changed.

**2026-09-19 — Concordance step added; figures and cross-dataset check renumbered.**
`07_concordance.R` is new and required. The figure script is now
`08_figures.R` (was `07_figures.R`) and the cross-dataset check
`09_cross_dataset_replication.R` (was `08_…`), so that the figures run after
the step they read. In `06_developmental.R` the within-genotype model list was
renamed from `fits` to `dev_fits`; it had been overwriting the object of the
same name left by `04_differential.R`, which nothing read. No analysis number
changed. Figures 7 and 8B now draw calls not confirmed by the replicate-aware
model as open symbols; in the current run that is 18 of 144 crossings and 9 of
45 categorical calls.

**2026-09-18 — Nrg1 Ensembl ID corrected in the myelination panel.**
`Curated_Myelination_Genes.csv` listed *Nrg1* as ENSMUSG00000118541, an
Ensembl entry that receives 11 read pairs across all 34 cornea libraries, so
every earlier run reported *Nrg1* as not expressed in every group. The entry
was replaced by the canonical *Nrg1* gene, ENSMUSG00000062991 (45,381 read
pairs). Panel size is unchanged (25). All 263 panel genes were screened for
the same fault (a same-symbol Ensembl entry carrying more than five times the
reads); *Nrg1* was the only one.

Effect on outputs: *Nrg1* is expressed in all nine group cells and is lower
in mutant cornea (log2FC −2.30 at P3/P4, −0.94 transparent adult vs WT,
−0.72 opaque vs WT, all adjusted P < 0.001; undetermined at P15 and opaque vs
transparent). Six rows were added to `DE_Results.csv`; no other fold change,
adjusted P value, outcome or expression-status call changed. The myelination
Benjamini–Hochberg family grew by one gene per contrast, moving adjusted P
values in that family by at most 0.04 with no call crossing 0.05.

## Requirements

R with `DESeq2`, `SummarizedExperiment`, `AnnotationDbi`, `org.Mm.eg.db`,
`dplyr`, `tidyr`, `tibble`, `readr`, `purrr`, `stringr`, `janitor`,
`data.table`, `rlang` and `ggplot2`; `patchwork` for the figures; `edgeR` for
the second-engine check in `07_concordance.R`, which is skipped with a message
if it is absent; `limma` (installed with `edgeR`) for the batch-removed view of
the sample PCA, which is likewise skipped without it. Exact versions used for the published results are recorded in
`RunProvenance.txt`.

The study's count matrix is not included in this repository; sequencing data
and the count matrix are deposited at GEO `<ACCESSION>`. The count matrices
were produced from the raw reads by the pipeline in
[pax6-cornea-trigeminal-rnaseq-upstream](https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream).
The independent dataset used by `09` is GEO GSE183742; its count matrix,
produced by the same pipeline, is included here.

## Methods references

`docs/REFERENCES.md` lists every statistical method used and what it is cited for.
Each script carries a numbered `STATISTICAL BASIS` block in its header with
inline markers at the implementing line.

## Citing

See `CITATION.cff`, or use the "Cite this repository" button on GitHub.

## Licence

MIT — see [`LICENSE`](LICENSE).
