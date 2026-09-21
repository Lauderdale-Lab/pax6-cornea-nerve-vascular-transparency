# Pax6 Mouse Cornea and Trigeminal RNA-seq

Analysis pipeline and code for bulk RNA-seq studies of Pax6 mutant and control mouse cornea and trigeminal ganglion tissues.

## Overview

This repository contains scripts and workflows used for:

- RNA-seq quality control
- Read alignment
- Count matrix generation
- Differential expression analysis (DESeq2)
- PCA and visualization
- Cornea and trigeminal ganglion comparative analyses

## Repository Structure

metadata/
scripts/
config/
results/
figures/
docs/

## Contact

James D. Lauderdale
Department of Cellular Biology
University of Georgia

## Citation

If you use this repository, please cite the associated manuscript and Zenodo DOI.




# Nerve remodeling in a Pax6 model of keratopathy — analysis pipeline

R code for the curated-panel RNA-seq analysis reported in:

> **Nerve remodeling in a Pax6 model of keratopathy**
> Sneha K. Mohan, James D. Lauderdale
> Department of Cellular Biology, University of Georgia, Athens, GA 30602, USA

Repository: `<REPO_URL>` · Archived release: `<ZENODO_DOI>` · Contact: `<CONTACT_EMAIL>`

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

## Running it

```r
setwd("<project root>")
source("run_all.R")
```

Scripts are numbered in execution order and share one R session by design:
`05` needs transcriptome-wide expression status, `06` needs the shared-gene
assignment, and `07` needs the fitted models from `04` and `06`, none of which
is written to disk. The runner stops at the first failure, reports any other
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
| `02_shared_gene_assignment.R` | primary panel for genes on more than one panel |
| `03_expression_status.R` | expression status and assessability |
| `04_differential.R` | differential expression, equivalence, omnibus, interactions |
| `05_program_level.R` | abundance-matched program-level permutation tests |
| `06_developmental.R` | developmental trajectories within genotype |
| `07_concordance.R` | status calls against a replicate-aware analysis of the same data; threshold sweep; QC figures |
| `08_figures.R` | the manuscript's four RNA-seq figures (optional; needs `patchwork`) |
| `09_cross_dataset_replication.R` | the calls against GSE183742 (optional; runs only when `PAX6_DUNCAN_COUNTS` and `PAX6_DUNCAN_META` are set) |
| `run_all.R` | driver and provenance |

`07` precedes the figures because they read its verdicts. `09` is a check that
nothing consumes.

Every parameter can be overridden by an environment variable of the same name
prefixed `PAX6_`, which is how sensitivity analyses are run without editing
code. `00_config.R` is the single place to read to see every choice the
analysis makes.

### Editing a curated panel

`00b_panel_universes.R` rebuilds each panel from its curated CSV and compares
it with the master list saved by the previous run. Any difference stops the
run, because every denominator in the paper depends on the panel sizes. That
is the intended behaviour for an *unintended* change. For an intended edit,
run once with the guard off, which rewrites the master list, then return to
strict mode:

```r
Sys.setenv(PAX6_PANEL_STRICT = "FALSE")
source("run_all.R")
Sys.unsetenv("PAX6_PANEL_STRICT")
```


## Outputs

Written to `Paper_Exports/ThreeState_Release/`:

| File | Contents |
|---|---|
| `ExpressionCalls.csv` | status, pooled counts and Poisson limits per group |
| `Assessability.csv` | what statement each gene supports, per comparison |
| `CategoricalGenes.csv` | genes on in one group and off in the other |
| `DE_Results.csv` | tested genes, bounded vs undetermined negatives |
| `InteractionResults.csv` | whether the genotype effect changes with age |
| `AdultOmnibusLRT.csv` | any difference among the three adult groups |
| `ProgramLevel.csv` | panel-level, abundance-matched |
| `Developmental_*.csv` | trajectories and status crossings |
| `Concordance/` | the concordance check, see below |
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
its single largest library supplied, and writes two routine diagnostics.

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
| `Concordance_ThresholdStability.csv` | crossings and Figure 7 row bands over the threshold grid |
| `Concordance_AgeAxis.png`, `Concordance_GenotypeAxis.png` | bounded fold against ungated log2 fold change, discordant genes labelled |
| `FigureS4a_DevelopmentalDE_Nerve.*`, `FigureS4b_DevelopmentalDE_Vascular.*` | the ungated within-genotype age contrasts as conventional heatmaps, one colour scale, Figure 7 crossings marked |
| `QC_PValueHistograms.png`, `QC_SampleCorrelation.png` | P-value histograms for every fitted contrast; sample–sample Spearman correlation on the panel genes |
| `Concordance_Summary.txt` | every number printed by the script, ending with the sentences for Methods |

### The supplementary combined comparison

Opaque mutant cornea differs from wild type in genotype **and** in transparency
at once, so that comparison cannot separate the two and is absent from the main
figures. To produce it as a supplementary figure:

```r
Sys.setenv(PAX6_FIG7_CONFOUNDED = "TRUE")
source("08_figures.R")
Sys.unsetenv("PAX6_FIG7_CONFOUNDED")
```

That run writes `FigureS3_Confounded.pdf` / `.png`, its source data and its own
caption numbers, and touches none of the main figures. It draws the
program-level figure only: adding the combined comparison changes the heatmap
colour limits, so heatmaps from that run would not be comparable with Figures 9
and 10.

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
if it is absent. Exact versions used for the published results are recorded in
`RunProvenance.txt`.

Input data are not included in this repository; sequencing data are deposited
at `<ACCESSION>`.

## Methods references

`REFERENCES.md` lists every statistical method used and what it is cited for.
Each script carries a numbered `STATISTICAL BASIS` block in its header with
inline markers at the implementing line.

## Citing

See `CITATION.cff`, or use the "Cite this repository" button on GitHub.

## Licence

MIT — see [`LICENSE`](LICENSE).
