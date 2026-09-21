###############################################################################
## 07_concordance.R -- does a standard analysis agree with the three-state calls?
##
## Analysis pipeline for:
##   Nerve remodeling in a Pax6 model of keratopathy
##   Sneha K. Mohan, James D. Lauderdale
##
## James D. Lauderdale, PhD  (ORCID 0000-0001-7503-0528)
## Department of Cellular Biology, University of Georgia
## Athens, GA 30602, USA
##
## Repository : <REPO_URL>
## Archived   : <ZENODO_DOI>
## Licence    : MIT (see LICENSE)
## Contact    : <CONTACT_EMAIL>
##
## Run the pipeline with run_all.R. Scripts are numbered in execution order and
## share one R session by design; see run_all.R for why.
###############################################################################

###############################################################################
## 07_concordance.R
##
## The three-state framework is an assembly of standard parts, but the assembly
## is this study's own. A reader is entitled to ask whether its calls would
## survive an analysis they already trust. This script answers that with the
## data in hand. It changes no call and no number upstream, but 08_figures.R
## READS ITS OUTPUT: a crossing or categorical call that the replicate-aware
## model does not confirm is drawn as an open symbol in Figures 7 and 8B.
## That is why this step precedes the figures.
##
## THREE CHECKS
##
##   1  AGE AXIS. Every status crossing in Figure 7 is looked up in the SAME
##      within-genotype DESeq2 model that 06_developmental.R fitted, but with
##      the assessability gate removed: every panel gene the model could fit is
##      tested, and adjusted P values are recomputed over that larger family.
##      A crossing is CONFIRMED when the ungated model calls it significant in
##      the same direction. The reverse question is asked too: which genes does
##      the ungated model call significant at the categorical fold that are NOT
##      crossings, and which route the pipeline sent them down instead.
##
##   2  GENOTYPE AXIS. The same for the categorical calls (Figure 8B), against
##      the ungated versions of Models A and B from 04_differential.R. Also:
##      how many `changed` calls among quantitatively tested genes move when
##      the correction family is enlarged by removing the gate, and which
##      `cannot assess` genes the ungated model would have called. Where edgeR
##      is installed, the genotype contrasts are refitted with its
##      quasi-likelihood pipeline on the identical design, so that agreement
##      does not rest on one implementation.
##
##   3  THRESHOLDS. The crossings are re-derived over a grid of ON/OFF
##      thresholds and minimum folds, and the three Figure 7 row bands are
##      checked for how many genes keep their band. The local implementation
##      of the rule must first reproduce 03_expression_status.R exactly, or the
##      script stops -- the same acceptance test 08 uses.
##
## TWO ROUTINE DIAGNOSTICS a reviewer expects and the pipeline did not yet
## write: P-value histograms for every fitted contrast (all genes in the model,
## not just the panels), and a sample-sample correlation heatmap over the
## panel genes.
##
## ONE SUPPLEMENTARY FIGURE: the ungated within-genotype age contrasts as a
## conventional log2-fold-change heatmap of every panel gene the standard
## model calls significant, with the Figure 7 crossings marked on it. This is
## the figure that lets Figure 7 be described as a status-level view of a
## standard result rather than as a method of its own.
##
## WHAT THIS SCRIPT CANNOT DO. It cannot make the crossings and the ungated
## model independent -- they use the same libraries, so agreement is
## expected and disagreement is what carries information. It is a check that
## the framework's restrictions did not manufacture a result, not evidence
## that the result is right. 09_cross_dataset_replication.R is the test on
## independent data.
##
## STATISTICAL BASIS
##   [1] Negative-binomial GLMs with empirical-Bayes dispersion and Wald tests,
##       as fitted in 04 and 06. Love MI, Huber W, Anders S (2014) Genome
##       Biology 15:550.
##   [2] Quasi-likelihood F-tests on negative-binomial GLMs, the second engine.
##       Lund SP, Nelson D, McCarthy DJ, Smyth GK (2012) Stat Appl Genet Mol
##       Biol 11:Article 8; Chen Y, Lun ATL, Smyth GK (2016) F1000Research
##       5:1438.
##   [3] TMM normalisation, edgeR's default. Robinson MD, Oshlack A (2010)
##       Genome Biology 11:R25.
##   [4] Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300, applied
##       within primary panel and contrast as in 04 and 06, but over the
##       UNGATED family.
##   [5] Exact Poisson limits and interval non-overlap, restated locally for
##       the threshold sweep. Garwood F (1936) Biometrika 28:437-442; Schenker
##       N, Gentleman JF (2001) Am Stat 55:182-186.
##
## Inputs : session objects from 01_load.R, 03_expression_status.R,
##          04_differential.R and 06_developmental.R, run in that order in one
##          session; and from OUT_ROOT: Developmental_StatusCrossings.csv,
##          Developmental_Steps.csv, Assessability.csv, DE_Results.csv
## Outputs: Concordance/Concordance_AgeAxis.csv
##          Concordance/Concordance_AgeAxis_StandardSignificant.csv
##          Concordance/Concordance_GenotypeAxis.csv
##          Concordance/Concordance_GenotypeAxis_Summary.csv
##          Concordance/Concordance_edgeR.csv               (if edgeR present)
##          Concordance/Concordance_ThresholdStability.csv
##          Concordance/Concordance_AgeAxis.png
##          Concordance/Concordance_GenotypeAxis.png
##          Concordance/<FIG_DEV_NERVE>.pdf/.png/_SourceData.csv
##          Concordance/<FIG_DEV_VASCULAR>.pdf/.png/_SourceData.csv
##          Concordance/QC_PValueHistograms.png
##          Concordance/QC_SampleCorrelation.png
##          Concordance/Concordance_Summary.txt
##
## Every top-level object here is prefixed c9_ / C9_. The scripts share one
## global environment and a bare name would overwrite an upstream object.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr)
  library(purrr); library(ggplot2)
})

c9_need <- c("OUT_ROOT", "count_mat", "meta", "lib_sizes", "panel_tbl", "panel_genes",
             "cell_stats", "cells_unique",
             "pair_results", "fit_A", "fit_B", "panel_primary", "assessability",
             "raw_steps", "step_assessability", "SERIES", "DEV_INTERVALS")
c9_missing <- c9_need[!vapply(c9_need, exists, logical(1))]
if (length(c9_missing)) {
  stop("07_concordance.R needs objects left in the session by 01, 03, 04 and 06, ",
       "run in that order. Missing: ", paste(c9_missing, collapse = ", "))
}

C9_OUT <- file.path(OUT_ROOT, "Concordance")
dir.create(C9_OUT, showWarnings = FALSE, recursive = TRUE)
c9_path <- function(...) file.path(C9_OUT, ...)

## Supplementary figure names. Parameters, as in 08_figures.R, because the
## supplement is renumbered as often as the main figures. The heatmap is split
## the way Figures 9 and 10 are: nerve panels in one, vascular in the other,
## because together they exceed a printed page. Both share one colour scale.
C9_FIG_DEV_NERVE    <- Sys.getenv("PAX6_FIG_DEV_NERVE",    "FigureS4a_DevelopmentalDE_Nerve")
C9_FIG_DEV_VASCULAR <- Sys.getenv("PAX6_FIG_DEV_VASCULAR", "FigureS4b_DevelopmentalDE_Vascular")

c9_read <- function(f) {
  p <- file.path(OUT_ROOT, f)
  if (!file.exists(p)) stop(f, " not found in ", OUT_ROOT, ". Run the pipeline through 06 first.")
  readr::read_csv(p, show_col_types = FALSE)
}

## Every sentence worth quoting is both printed and kept for the summary file.
c9_log <- character(0)
c9_say <- function(...) { t <- paste0(...); message(t); c9_log <<- c(c9_log, t) }
c9_show <- function(df) print(as.data.frame(df), row.names = FALSE)

## One row per gene, with the panel it is corrected in. Every family below is
## defined on primary panel so that a shared gene is counted once, exactly as
## in 04 and 06.
c9_gene_primary <- panel_tbl %>%
  dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel), by = "symbol") %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::distinct(gene_id, symbol, primary_panel)
if (nrow(c9_gene_primary) != length(panel_genes)) {
  stop("Primary-panel resolution gave ", nrow(c9_gene_primary), " genes for ",
       length(panel_genes), " panel genes. The assignment in 04 is inconsistent.")
}

## The ungated standard analysis: every panel gene the model fitted, corrected
## within primary panel and the given family columns. [4]
c9_standard <- function(df, family_cols) {
  df %>%
    dplyr::inner_join(c9_gene_primary, by = "gene_id") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("primary_panel", family_cols)))) %>%
    dplyr::mutate(padj_std = stats::p.adjust(pvalue, method = "BH"),
                  n_family_std = sum(!is.na(pvalue))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outcome_std = dplyr::case_when(
      is.na(pvalue)  ~ "no P value (outlier)",
      padj_std < FDR ~ "significant",
      TRUE           ~ "not significant"))
}

C9_NOT_IN_MODEL <- "not in model (count pre-filter)"
C9_LFC_FLOOR <- log2(CATEGORICAL_MIN_FOLD)

## Library ids per group, by the cell names 03 used.
c9_ids_all <- stats::setNames(cells_unique$ids, cells_unique$cell)

## How much of a group's pooled count for one gene came from its single
## largest library. Status is called on POOLED counts, so a library carrying
## most of the reads carries the call; the leave-one-out check catches that
## only when dropping the library flips the status, which it need not. A
## negative-binomial model sees the replicate spread and does not agree. With
## five libraries an even share is 0.2; above 0.5 one library IS the call.
c9_top_share <- function(gene_id, cell) {
  k <- count_mat[gene_id, c9_ids_all[[cell]]]
  if (sum(k) == 0) return(NA_real_)
  round(max(k) / sum(k), 3)
}
c9_top_share_vec <- function(gene_ids, cells) {
  unname(mapply(c9_top_share, gene_ids, cells))
}

c9_say(strrep("=", 74))
c9_say("07_concordance -- the three-state calls against a standard analysis")
c9_say(strrep("=", 74))

## ---------------------------------------------------------------------------
## 1. Age axis: the Figure 7 crossings
## ---------------------------------------------------------------------------

c9_say("")
c9_say("=== 1. Age axis: Figure 7 crossings against the ungated within-genotype model ===")

c9_cross <- c9_read("Developmental_StatusCrossings.csv") %>%
  dplyr::filter(panel == primary_panel)

## raw_steps is 06's own DESeq2 output for EVERY gene that passed the count
## pre-filter, before the assessability gate chose which to report. [1]
c9_std_age <- c9_standard(raw_steps, c("genotype", "step"))

c9_age <- c9_cross %>%
  dplyr::left_join(
    c9_std_age %>% dplyr::select(genotype, step, gene_id, log2FoldChange, lfcSE,
                                 padj_std, outcome_std),
    by = c("genotype", "step", "gene_id")) %>%
  dplyr::mutate(
    outcome_std = dplyr::coalesce(outcome_std, C9_NOT_IN_MODEL),
    direction_agrees = dplyr::case_when(
      is.na(log2FoldChange)      ~ NA,
      direction == "on with age" ~ log2FoldChange > 0,
      TRUE                       ~ log2FoldChange < 0),
    verdict = dplyr::case_when(
      outcome_std == "significant" & direction_agrees %in% TRUE ~ "confirmed",
      outcome_std == "significant"                              ~ "significant, opposite direction",
      TRUE                                                      ~ outcome_std),
    cell_earlier = unname(mapply(function(g, s)
      SERIES[[g]][[DEV_INTERVALS$earlier[match(s, DEV_INTERVALS$step)]]], genotype, step)),
    cell_later = unname(mapply(function(g, s)
      SERIES[[g]][[DEV_INTERVALS$later[match(s, DEV_INTERVALS$step)]]], genotype, step)),
    top_library_share_earlier = c9_top_share_vec(gene_id, cell_earlier),
    top_library_share_later   = c9_top_share_vec(gene_id, cell_later)) %>%
  dplyr::select(genotype, step, primary_panel, symbol, gene_id, direction, assessable,
                cpm_earlier, cpm_later, bounded_fold, log2FoldChange, lfcSE,
                padj_std, outcome_std, verdict,
                top_library_share_earlier, top_library_share_later) %>%
  dplyr::arrange(genotype, step, dplyr::desc(bounded_fold))
readr::write_csv(c9_age, c9_path("Concordance_AgeAxis.csv"))

c9_age_tab <- c9_age %>%
  dplyr::count(genotype, verdict) %>%
  tidyr::pivot_wider(names_from = verdict, values_from = n, values_fill = 0)
c9_show(c9_age_tab)

c9_n_confirmed <- sum(c9_age$verdict == "confirmed")
c9_say(sprintf("  %d of %d crossings are significant in the same direction under the ungated model (%.1f%%).",
               c9_n_confirmed, nrow(c9_age), 100 * c9_n_confirmed / nrow(c9_age)))
c9_say(sprintf("  Family sizes: the pipeline corrected %s genes per genotype x interval x panel; ungated, %s.",
               paste(range(c9_read("Developmental_Steps.csv")$n_in_correction_family), collapse = "-"),
               paste(range(c9_std_age$n_family_std), collapse = "-")))

c9_discordant <- c9_age %>% dplyr::filter(verdict != "confirmed")
if (nrow(c9_discordant)) {
  c9_say("  Crossings the ungated model does NOT confirm:")
  c9_show(c9_discordant %>% dplyr::transmute(
    genotype, step, symbol, direction, cpm_earlier, cpm_later, bounded_fold,
    log2FC = round(log2FoldChange, 2), padj_std = signif(padj_std, 2), verdict,
    top_share_earlier = top_library_share_earlier,
    top_share_later = top_library_share_later))
  c9_say("  A crossing that is `not significant` here sits where the pooled Poisson")
  c9_say("  interval and the replicate-aware model disagree. `top_share` is the")
  c9_say("  fraction of the pooled count supplied by the single largest library at")
  c9_say("  that age (even share with five libraries: 0.2). Above 0.5, one library")
  c9_say("  is the call. Quote such genes at status level only, or not at all.")
  c9_n_dominated <- sum(pmax(c9_discordant$top_library_share_earlier,
                             c9_discordant$top_library_share_later, na.rm = TRUE) > 0.5)
  c9_say(sprintf("  %d of the %d unconfirmed crossings have one library supplying more than half the pooled count at one age.",
                 c9_n_dominated, nrow(c9_discordant)))
  c9_n_dominated_ok <- sum(pmax(c9_age$top_library_share_earlier,
                                c9_age$top_library_share_later, na.rm = TRUE) > 0.5 &
                           c9_age$verdict == "confirmed")
  c9_say(sprintf("  For comparison, %d of the %d confirmed crossings do.",
                 c9_n_dominated_ok, c9_n_confirmed))
} else {
  c9_say("  Every crossing is confirmed.")
}

## The reverse: what the standard model calls at the categorical fold that
## Figure 7 does not show, and where the pipeline sent each one.
c9_dev_steps <- c9_read("Developmental_Steps.csv") %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::select(genotype, step, gene_id, outcome_pipeline = outcome)

c9_age_rev <- c9_std_age %>%
  dplyr::filter(outcome_std == "significant", abs(log2FoldChange) >= C9_LFC_FLOOR) %>%
  dplyr::left_join(
    step_assessability %>% dplyr::select(genotype, step, gene_id, assessable,
                                         cpm_earlier, cpm_later, bounded_fold),
    by = c("genotype", "step", "gene_id")) %>%
  dplyr::left_join(c9_dev_steps, by = c("genotype", "step", "gene_id")) %>%
  dplyr::mutate(pipeline_route = dplyr::case_when(
    startsWith(assessable, "status crossing") ~ "status crossing (in Figure 7)",
    assessable == "test quantitatively" ~
      paste0("tested quantitatively: ", dplyr::coalesce(outcome_pipeline, "not in output")),
    TRUE ~ assessable)) %>%
  dplyr::select(genotype, step, primary_panel, symbol, gene_id, log2FoldChange, padj_std,
                cpm_earlier, cpm_later, bounded_fold, assessable, pipeline_route) %>%
  dplyr::arrange(genotype, step, pipeline_route, dplyr::desc(abs(log2FoldChange)))
readr::write_csv(c9_age_rev, c9_path("Concordance_AgeAxis_StandardSignificant.csv"))

c9_say("")
c9_say(sprintf("  Ungated model, significant at >= %.1f-fold: %d gene x genotype x interval calls. Where the pipeline sent them:",
               CATEGORICAL_MIN_FOLD, nrow(c9_age_rev)))
c9_show(c9_age_rev %>% dplyr::count(genotype, pipeline_route) %>%
          tidyr::pivot_wider(names_from = genotype, values_from = n, values_fill = 0))
c9_say("  `tested quantitatively` rows are genes expressed at both ages: a real change")
c9_say("  that is not a change of STATUS, so Figure 7 correctly omits them and the")
c9_say("  pipeline's own outcome for each is shown. `cannot assess` rows are what the")
c9_say("  crossing rule misses: neither age confidently expressed, or the intervals")
c9_say("  overlap at the chosen fold. Those are in the supplementary heatmap.")

## ---------------------------------------------------------------------------
## 2. Genotype axis: the categorical calls
## ---------------------------------------------------------------------------

c9_say("")
c9_say("=== 2. Genotype axis: categorical calls against the ungated Models A and B ===")

## pair_results is 04's own output for every fitted gene, before the gate. [1]
c9_std_geno <- pair_results %>%
  dplyr::filter(contrast != HINGE[["A"]]) %>%
  c9_standard("contrast")

c9_assess <- assessability %>%
  dplyr::filter(contrast != HINGE[["A"]]) %>%
  dplyr::inner_join(dplyr::select(c9_gene_primary, gene_id, primary_panel), by = "gene_id") %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::select(contrast, gene_id, symbol, primary_panel, assessable, call_a, call_b,
                cpm_a, cpm_b, bounded_fold, higher_in, cell_a)

c9_de <- c9_read("DE_Results.csv") %>%
  dplyr::filter(reportable, panel == primary_panel) %>%
  dplyr::select(contrast, gene_id, outcome_pipeline = outcome, padj_pipeline = padj)

c9_geno <- c9_assess %>%
  dplyr::left_join(
    c9_std_geno %>% dplyr::select(contrast, gene_id, log2FoldChange, lfcSE, padj_std,
                                  outcome_std, n_family_std),
    by = c("contrast", "gene_id")) %>%
  dplyr::left_join(c9_de, by = c("contrast", "gene_id")) %>%
  dplyr::mutate(
    outcome_std = dplyr::coalesce(outcome_std, C9_NOT_IN_MODEL),
    direction_agrees = dplyr::case_when(
      is.na(log2FoldChange) | is.na(higher_in) ~ NA,
      higher_in == cell_a                      ~ log2FoldChange > 0,
      TRUE                                     ~ log2FoldChange < 0),
    route = dplyr::case_when(
      startsWith(assessable, "categorical") ~ "categorical",
      assessable == "test quantitatively"   ~ "tested quantitatively",
      TRUE                                  ~ assessable),
    verdict = dplyr::case_when(
      outcome_std == "significant" & direction_agrees %in% TRUE ~ "significant, same direction",
      outcome_std == "significant"                              ~ "significant, opposite direction",
      TRUE                                                      ~ outcome_std)) %>%
  dplyr::arrange(contrast, route, dplyr::desc(bounded_fold))

c9_geno <- c9_geno %>%
  dplyr::mutate(
    cell_b = assessability$cell_b[match(paste(contrast, gene_id),
                                        paste(assessability$contrast, assessability$gene_id))],
    top_library_share_a = c9_top_share_vec(gene_id, cell_a),
    top_library_share_b = c9_top_share_vec(gene_id, cell_b))
readr::write_csv(c9_geno, c9_path("Concordance_GenotypeAxis.csv"))

## 2a. Categorical calls
c9_cat <- c9_geno %>% dplyr::filter(route == "categorical")
c9_say(sprintf("  2a. %d categorical calls (distinct gene x comparison):", nrow(c9_cat)))
c9_show(c9_cat %>% dplyr::count(contrast, verdict) %>%
          tidyr::pivot_wider(names_from = verdict, values_from = n, values_fill = 0))
c9_n_cat_ok <- sum(c9_cat$verdict == "significant, same direction")
c9_say(sprintf("  %d of %d categorical calls are significant in the same direction under the ungated model (%.1f%%).",
               c9_n_cat_ok, nrow(c9_cat), 100 * c9_n_cat_ok / nrow(c9_cat)))
c9_cat_bad <- c9_cat %>% dplyr::filter(verdict != "significant, same direction")
if (nrow(c9_cat_bad)) {
  c9_say("  Categorical calls the ungated model does NOT confirm:")
  c9_show(c9_cat_bad %>% dplyr::transmute(
    contrast, symbol, cpm_a, cpm_b, bounded_fold,
    log2FC = round(log2FoldChange, 2), padj_std = signif(padj_std, 2), verdict,
    top_share_a = top_library_share_a, top_share_b = top_library_share_b))
  c9_say("  `top_share` as in section 1: the largest library's fraction of the pooled")
  c9_say("  count in each group. Above 0.5, that library is the call.")
}

## 2b. Quantitatively tested genes: does enlarging the family move `changed`?
c9_tq <- c9_geno %>%
  dplyr::filter(route == "tested quantitatively", !is.na(outcome_pipeline)) %>%
  dplyr::mutate(changed_pipeline = outcome_pipeline == "changed",
                changed_std = outcome_std == "significant")
## The per-gene flags are combined BEFORE summarise(): inside summarise a
## column that has just been totalled is the total from that line on, so
## `sum(changed_pipeline & ...)` after `changed_pipeline = sum(...)` would be
## evaluated against a scalar. That is how a first version reported more genes
## lost than were ever called.
c9_tq <- c9_tq %>%
  dplyr::mutate(lost = changed_pipeline & !changed_std,
                gained = !changed_pipeline & changed_std)
c9_tq_tab <- c9_tq %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    n = dplyr::n(),
    changed_pipeline = sum(changed_pipeline),
    changed_ungated = sum(changed_std),
    lost_when_ungated = sum(lost),
    gained_when_ungated = sum(gained),
    .groups = "drop")
c9_tq_moved <- c9_tq %>%
  dplyr::filter(lost | gained) %>%
  dplyr::transmute(contrast, symbol, primary_panel, log2FoldChange,
                   padj_pipeline, padj_std,
                   moved = ifelse(lost, "changed only under pipeline family",
                                  "changed only under ungated family"))
readr::write_csv(c9_tq_moved, c9_path("Concordance_GenotypeAxis_FamilyMoves.csv"))
c9_say("")
c9_say("  2b. Quantitatively tested genes: `changed` under the pipeline's family vs the ungated family")
c9_show(c9_tq_tab)
if (nrow(c9_tq_moved)) c9_show(c9_tq_moved)
c9_say("  The only difference between the two columns is the size of the correction")
c9_say("  family. Genes that move are borderline under either frame and should not")
c9_say("  carry a sentence on their own.")

## 2c. What the gate hid that a standard analysis would have called.
c9_hidden <- c9_geno %>%
  dplyr::filter(startsWith(route, "cannot assess"), outcome_std == "significant",
                abs(log2FoldChange) >= C9_LFC_FLOOR)
c9_say("")
c9_say(sprintf("  2c. `cannot assess` genes the ungated model calls significant at >= %.1f-fold: %d",
               CATEGORICAL_MIN_FOLD, nrow(c9_hidden)))
if (nrow(c9_hidden)) {
  c9_show(c9_hidden %>% dplyr::count(contrast, assessable))
  c9_say("  These are genes neither route reports. They are listed in")
  c9_say("  Concordance_GenotypeAxis.csv and belong in a supplementary table if the")
  c9_say("  count is more than a handful; each is a gene the framework declined to")
  c9_say("  call and a conventional analysis would have.")
}

c9_geno_summary <- c9_geno %>%
  dplyr::count(contrast, route, verdict) %>%
  dplyr::arrange(contrast, route, verdict)
readr::write_csv(c9_geno_summary, c9_path("Concordance_GenotypeAxis_Summary.csv"))

## 2d. A second engine. edgeR quasi-likelihood on the identical design, so that
## agreement between the categorical calls and "a standard analysis" does not
## rest on DESeq2 alone. [2] [3]
if (requireNamespace("edgeR", quietly = TRUE)) {
  c9_say("")
  c9_say("=== 2d. Second engine: edgeR quasi-likelihood, same samples, same design ===")

  c9_edger_fit <- function(fit) {
    cd <- fit$col_data
    terms <- c(if (nlevels(cd$batch) > 1) "batch",
               if (USE_SEX_COVARIATE) "male_frac", "cell_group")
    design <- stats::model.matrix(
      stats::as.formula(paste("~", paste(terms, collapse = " + "))), cd)
    ## Same pre-filter as 04, so both engines see the same genes.
    cm <- count_mat[, rownames(cd), drop = FALSE]
    cm <- cm[rowSums(cm >= 10) >= 2, , drop = FALSE]
    ## edgeR 4 renamed calcNormFactors to normLibSizes and warns on the old
    ## name; either computes TMM factors.
    c9_norm <- if ("normLibSizes" %in% getNamespaceExports("edgeR"))
      edgeR::normLibSizes else edgeR::calcNormFactors
    y <- c9_norm(edgeR::DGEList(cm))                          # [3]
    y <- edgeR::estimateDisp(y, design)
    list(fit = edgeR::glmQLFit(y, design), design = design, reference = fit$reference)
  }

  c9_edger_test <- function(ef, label, level_a, level_b) {
    cn <- colnames(ef$design)
    idx_of <- function(level) {
      if (identical(level, ef$reference)) return(NA_integer_)
      i <- match(paste0("cell_group", level), cn)
      if (is.na(i)) stop("No design column for cell_group level ", level)
      i
    }
    v <- numeric(length(cn))
    ia <- idx_of(level_a); ib <- idx_of(level_b)
    if (!is.na(ia)) v[ia] <- 1
    if (!is.na(ib)) v[ib] <- -1
    tt <- edgeR::glmQLFTest(ef$fit, contrast = v)$table
    tibble::tibble(contrast = label, gene_id = rownames(tt),
                   logFC_edger = tt$logFC, pvalue = tt$PValue)
  }

  c9_ef <- list(A = c9_edger_fit(fit_A), B = c9_edger_fit(fit_B))
  c9_contr <- CONTRASTS[CONTRASTS$label != HINGE[["A"]], , drop = FALSE]
  c9_edger <- purrr::pmap_dfr(
    list(c9_contr$model, c9_contr$label, c9_contr$level_a, c9_contr$level_b),
    function(model, label, a, b) c9_edger_test(c9_ef[[model]], label, a, b)) %>%
    c9_standard("contrast") %>%
    dplyr::select(contrast, gene_id, logFC_edger, pvalue_edger = pvalue,
                  padj_edger = padj_std, outcome_edger = outcome_std)
  readr::write_csv(c9_edger, c9_path("Concordance_edgeR.csv"))

  c9_two <- c9_geno %>%
    dplyr::inner_join(c9_edger, by = c("contrast", "gene_id")) %>%
    dplyr::mutate(edger_agrees_direction = dplyr::case_when(
      is.na(higher_in)    ~ NA,
      higher_in == cell_a ~ logFC_edger > 0,
      TRUE                ~ logFC_edger < 0))

  ## Engine agreement over every panel gene both could test.
  c9_eng <- c9_two %>%
    dplyr::filter(outcome_std %in% c("significant", "not significant"),
                  outcome_edger %in% c("significant", "not significant")) %>%
    dplyr::group_by(contrast) %>%
    dplyr::summarise(
      n = dplyr::n(),
      sig_deseq2 = sum(outcome_std == "significant"),
      sig_edger = sum(outcome_edger == "significant"),
      both = sum(outcome_std == "significant" & outcome_edger == "significant"),
      pct_agree = round(100 * mean(outcome_std == outcome_edger), 1),
      r_logFC = round(stats::cor(log2FoldChange, logFC_edger, use = "complete.obs"), 3),
      .groups = "drop")
  c9_say("  Engine agreement over panel genes fitted by both, ungated:")
  c9_show(c9_eng)

  c9_cat2 <- c9_two %>% dplyr::filter(route == "categorical")
  c9_n_cat_edger <- sum(c9_cat2$outcome_edger == "significant" &
                          c9_cat2$edger_agrees_direction %in% TRUE)
  c9_say(sprintf("  Categorical calls significant in the same direction under edgeR QL: %d of %d fitted (%d of %d overall).",
                 c9_n_cat_edger, nrow(c9_cat2), c9_n_cat_edger, nrow(c9_cat)))
} else {
  c9_say("")
  c9_say("  2d. edgeR is not installed; the second-engine check was skipped.")
  c9_say("      install.packages(\"BiocManager\"); BiocManager::install(\"edgeR\")")
}

## ---------------------------------------------------------------------------
## 3. Threshold sweep for the crossings
## ---------------------------------------------------------------------------

c9_say("")
c9_say("=== 3. Do the Figure 7 crossings depend on the thresholds? ===")

## The rule, restated on the panel rows only. Status is per gene and does not
## depend on any other gene, so restricting to the panel changes nothing and
## makes the leave-one-out loop cheap. Library sizes are the FULL-matrix sizes,
## as in 03; deriving them from the panel rows would be wrong. [5]
c9_cm <- count_mat[panel_genes, , drop = FALSE]
c9_limits <- function(ids) {
  k <- rowSums(c9_cm[, ids, drop = FALSE])
  L <- sum(lib_sizes[ids])
  list(lcl = ifelse(k == 0, 0, stats::qgamma(POISSON_ALPHA, shape = k)) / L * 1e6,
       ucl = stats::qgamma(1 - POISSON_ALPHA, shape = k + 1) / L * 1e6)
}
c9_call <- function(ids, t_hi, t_lo) {
  l <- c9_limits(ids)
  s <- rep("indeterminate", length(l$lcl))
  s[l$lcl >= t_hi] <- "expressed"
  s[l$ucl < t_lo & l$lcl < t_hi] <- "not_expressed"
  s
}
c9_status <- function(ids, t_hi, t_lo) {
  full <- c9_call(ids, t_hi, t_lo)
  if (STABILITY_GATE && length(ids) >= 3) {
    loo <- vapply(ids, function(d) c9_call(setdiff(ids, d), t_hi, t_lo),
                  character(length(full)))
    full[rowSums(loo != full) > 0] <- "indeterminate"
  }
  full
}

c9_series_cells <- unique(unlist(SERIES))
c9_ids <- stats::setNames(cells_unique$ids, cells_unique$cell)[c9_series_cells]

## ACCEPTANCE TEST. The local rule must reproduce 03's calls exactly at the
## configured thresholds, or nothing below can be interpreted.
for (cl in c9_series_cells) {
  if (!identical(cell_stats[[cl]]$gene_id, panel_genes)) {
    stop("cell_stats[['", cl, "']] is not in panel_genes order; cannot compare.")
  }
  got <- c9_status(c9_ids[[cl]], T_HI, T_LO)
  ref <- cell_stats[[cl]]$call
  if (any(got != ref)) {
    stop("Local status rule disagrees with 03_expression_status.R in group ", cl,
         " for ", sum(got != ref), " gene(s). The threshold sweep would be uninterpretable.")
  }
}
c9_say(sprintf("  acceptance: local rule reproduces 03's calls in all %d groups of the series.",
               length(c9_series_cells)))

c9_crossings_at <- function(t_hi, t_lo, min_fold) {
  st  <- lapply(c9_ids, c9_status, t_hi = t_hi, t_lo = t_lo)
  lim <- lapply(c9_ids, c9_limits)
  purrr::map_dfr(names(SERIES), function(gt) {
    purrr::pmap_dfr(DEV_INTERVALS[, c("step", "later", "earlier")],
                    function(step, later, earlier) {
      a <- SERIES[[gt]][[later]]; b <- SERIES[[gt]][[earlier]]
      A <- st[[a]]; B <- st[[b]]
      bf <- pmax(lim[[a]]$lcl / lim[[b]]$ucl, lim[[b]]$lcl / lim[[a]]$ucl)
      absent <- (A == "expressed" & B == "not_expressed") |
                (A == "not_expressed" & B == "expressed")
      separated <- xor(A == "expressed", B == "expressed") & !absent & bf >= min_fold
      tibble::tibble(genotype = gt, step = step, gene_id = panel_genes[absent | separated])
    })
  })
}

## Second acceptance test: at the configured settings the sweep must give back
## exactly the crossings 06 wrote.
c9_ref <- c9_crossings_at(T_HI, T_LO, CATEGORICAL_MIN_FOLD)
c9_key <- function(df) paste(df$genotype, df$step, df$gene_id)
if (!setequal(c9_key(c9_ref), c9_key(c9_cross))) {
  stop("The threshold sweep at the configured settings does not reproduce ",
       "Developmental_StatusCrossings.csv (", nrow(c9_ref), " vs ", nrow(c9_cross),
       " crossings). Investigate before reading the sweep.")
}
c9_say("  acceptance: sweep at the configured settings reproduces 06's crossings exactly.")

c9_band <- function(cr) {
  cr %>% dplyr::distinct(gene_id, genotype) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(band = if (dplyr::n() > 1) "both" else genotype[1], .groups = "drop")
}
c9_ref_band <- c9_band(c9_cross)

c9_grid <- tidyr::expand_grid(
  tibble::tibble(t_hi = c(0.5, 1, 2), t_lo = c(0.1, 0.2, 0.4)),
  min_fold = c(1.5, 2, 3))

c9_stab <- purrr::pmap_dfr(c9_grid, function(t_hi, t_lo, min_fold) {
  cr <- c9_crossings_at(t_hi, t_lo, min_fold)
  bd <- c9_band(cr)
  j <- c9_ref_band %>% dplyr::left_join(bd, by = "gene_id", suffix = c("_ref", ""))
  tibble::tibble(
    t_hi = t_hi, t_lo = t_lo, min_fold = min_fold,
    is_configured = t_hi == T_HI & t_lo == T_LO & min_fold == CATEGORICAL_MIN_FOLD,
    n_crossings = nrow(cr), n_genes = nrow(bd),
    ref_genes = nrow(c9_ref_band),
    ref_genes_retained = sum(!is.na(j$band)),
    ref_genes_same_band = sum(!is.na(j$band) & j$band == j$band_ref),
    genes_added = sum(!bd$gene_id %in% c9_ref_band$gene_id))
})
readr::write_csv(c9_stab, c9_path("Concordance_ThresholdStability.csv"))
c9_show(c9_stab)
c9_say("  `ref_genes_same_band` is the number of the configured run's genes that keep")
c9_say("  their Figure 7 row band (wild type only / both / Sey only) at that setting.")
c9_say("  A band that is stable across the grid can be described; one that is not")
c9_say("  should be described only at the configured setting, with that stated.")

## ---------------------------------------------------------------------------
## 4. Figures
## ---------------------------------------------------------------------------

c9_theme <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
                 strip.text = ggplot2::element_text(face = "bold", size = 8),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 plot.title = ggplot2::element_text(face = "bold", size = 10),
                 plot.subtitle = ggplot2::element_text(size = 7, colour = "grey30"))

C9_PANEL_COLOUR <- c("Axon guidance" = "#0072B2", "Myelination" = "#009E73",
                     "Vascular/Lymphatic" = "#D55E00")
C9_VERDICT_COLOUR <- c("confirmed" = "#0072B2",
                       "significant, same direction" = "#0072B2",
                       "significant, opposite direction" = "#D55E00",
                       "not significant" = "#E69F00",
                       "no P value (outlier)" = "grey40")

c9_save <- function(plot, name, width, height) {
  p <- c9_path(paste0(name, ".png"))
  ggplot2::ggsave(p, plot, width = width, height = height, units = "in", dpi = 300)
  if (!file.exists(p) || file.info(p)$size == 0) stop(name, " was not written.")
  invisible(p)
}

## 4a. Age-axis concordance
c9_age_plot <- c9_age %>% dplyr::filter(!is.na(log2FoldChange)) %>%
  dplyr::mutate(genotype = factor(genotype, levels = names(SERIES)))
c9_p_age <- ggplot2::ggplot(c9_age_plot,
    ggplot2::aes(bounded_fold, log2FoldChange, colour = verdict)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = c(-C9_LFC_FLOOR, C9_LFC_FLOOR),
                      linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_point(size = 1.6, alpha = 0.85) +
  ggplot2::geom_text(data = c9_age_plot %>% dplyr::filter(verdict != "confirmed"),
                     ggplot2::aes(label = symbol), size = 2.3, vjust = -0.8,
                     show.legend = FALSE) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = C9_VERDICT_COLOUR, name = NULL) +
  ggplot2::facet_wrap(~ genotype) +
  ggplot2::labs(
    x = "Minimum fold difference from the Poisson intervals (log scale)",
    y = expression(log[2]~fold~change*","~ungated~model),
    title = "Figure 7 crossings under a standard within-genotype model",
    subtitle = sprintf("%d of %d crossings confirmed; %d not in the model (count pre-filter). Dashed: %.1f-fold.",
                       c9_n_confirmed, nrow(c9_age),
                       sum(c9_age$outcome_std == C9_NOT_IN_MODEL), CATEGORICAL_MIN_FOLD)) +
  c9_theme
c9_save(c9_p_age, "Concordance_AgeAxis", 7, 4.2)

## 4b. Genotype-axis concordance, categorical calls only
c9_cat_plot <- c9_cat %>% dplyr::filter(!is.na(log2FoldChange))
c9_p_geno <- ggplot2::ggplot(c9_cat_plot,
    ggplot2::aes(bounded_fold, log2FoldChange, colour = verdict)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_point(size = 1.6, alpha = 0.85) +
  ggplot2::geom_text(data = c9_cat_plot %>%
                       dplyr::filter(verdict != "significant, same direction"),
                     ggplot2::aes(label = symbol), size = 2.3, vjust = -0.8,
                     show.legend = FALSE) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = C9_VERDICT_COLOUR, name = NULL) +
  ggplot2::facet_wrap(~ contrast, ncol = 3) +
  ggplot2::labs(
    x = "Minimum fold difference from the Poisson intervals (log scale)",
    y = expression(log[2]~fold~change*","~ungated~model),
    title = "Categorical calls under the ungated standard model",
    subtitle = sprintf("%d of %d calls significant in the same direction", c9_n_cat_ok, nrow(c9_cat))) +
  c9_theme
c9_save(c9_p_geno, "Concordance_GenotypeAxis", 7, 5)

## 4c. Supplementary: the ungated age contrasts as a conventional heatmap
C9_STEP_LABEL <- c(P15_vs_P3_P4 = "P3/P4\nto P15", Adult_vs_P15 = "P15\nto adult")
C9_GENOTYPE_LABEL <- c(WT = "Wild type", Sey = "Sey")

c9_heat_all <- c9_std_age %>%
  dplyr::filter(step %in% names(C9_STEP_LABEL))
c9_heat_genes <- c9_heat_all %>%
  dplyr::filter(outcome_std == "significant") %>%
  dplyr::group_by(primary_panel, symbol, gene_id) %>%
  dplyr::summarise(m = max(abs(log2FoldChange)), .groups = "drop") %>%
  dplyr::arrange(factor(primary_panel, levels = names(C9_PANEL_COLOUR)), dplyr::desc(m))

c9_cross_keys <- paste(c9_cross$genotype, c9_cross$step, c9_cross$gene_id)
c9_heat <- c9_heat_all %>%
  dplyr::semi_join(c9_heat_genes, by = "gene_id") %>%
  dplyr::mutate(
    shown = ifelse(outcome_std == "significant", log2FoldChange, NA_real_),
    is_crossing = paste(genotype, step, gene_id) %in% c9_cross_keys,
    symbol = factor(symbol, levels = rev(c9_heat_genes$symbol)),
    step = factor(C9_STEP_LABEL[step], levels = unname(C9_STEP_LABEL)),
    genotype = factor(C9_GENOTYPE_LABEL[genotype], levels = unname(C9_GENOTYPE_LABEL)),
    primary_panel = factor(primary_panel, levels = names(C9_PANEL_COLOUR)))
## ONE colour limit for both halves, taken over every gene drawn in either, so
## that a tile of a given colour means the same fold change in S4a and S4b.
c9_lim <- max(abs(c9_heat$shown), na.rm = TRUE)

c9_dev_heat <- function(keep, name, title) {
  h <- c9_heat %>% dplyr::filter(primary_panel %in% keep) %>%
    dplyr::mutate(primary_panel = droplevels(primary_panel))
  genes <- c9_heat_genes %>% dplyr::filter(primary_panel %in% keep)
  if (!nrow(h)) { c9_say("  ", name, ": no genes to draw -- skipped."); return(invisible(0L)) }
  p <- ggplot2::ggplot(h, ggplot2::aes(step, symbol, fill = shown)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_point(data = h %>% dplyr::filter(is_crossing),
                        shape = 16, size = 0.9, colour = "black") +
    ggplot2::facet_grid(primary_panel ~ genotype, scales = "free", space = "free") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-c9_lim, c9_lim),
                                  na.value = "grey92",
                                  name = expression(log[2]~fold~change)) +
    ggplot2::labs(x = NULL, y = NULL, title = title,
                  subtitle = paste(
                    "Every panel gene significant in either interval, ungated DESeq2 within genotype.",
                    "Grey: tested, not significant. Blank: not in the model.",
                    "Dot: a Figure 7 status crossing. Colour scale shared with the companion figure.",
                    sep = "\n")) +
    c9_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.5),
                   panel.grid = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0, size = 7))
  ## Titles are two lines: ggplot2 clips a title at the device edge exactly as
  ## it clips a subtitle, and one line does not fit at 5.0 in. The extra line
  ## is paid for here.
  hgt <- min(48, 2.5 + 0.12 * nrow(genes))
  c9_save(p, name, 5.0, hgt)
  ggplot2::ggsave(c9_path(paste0(name, ".pdf")), p, width = 5.0, height = hgt, units = "in")
  readr::write_csv(
    h %>% dplyr::transmute(genotype, interval = step, panel = primary_panel, symbol,
                           log2FoldChange, padj = padj_std, outcome = outcome_std,
                           figure7_crossing = is_crossing),
    c9_path(paste0(name, "_SourceData.csv")))
  c9_say(sprintf("  %s: %d genes (%s); %d of the %d crossings drawn fall on a significant tile.",
                 name, nrow(genes), paste(keep, collapse = " + "),
                 sum(h$is_crossing & !is.na(h$shown)), sum(h$is_crossing)))
  invisible(nrow(genes))
}

c9_say("")
c9_say(sprintf("  Supplementary heatmaps: %d genes significant in at least one adjacent interval, either genotype, drawn on one colour scale (limit %.1f log2).",
               nrow(c9_heat_genes), c9_lim))
c9_dev_heat(c("Axon guidance", "Myelination"), C9_FIG_DEV_NERVE,
            "Developmental change within genotype, standard model\nAxon guidance and myelination")
c9_dev_heat("Vascular/Lymphatic", C9_FIG_DEV_VASCULAR,
            "Developmental change within genotype, standard model\nVascular and lymphatic")
c9_adj <- c9_cross %>% dplyr::filter(step %in% names(C9_STEP_LABEL))
c9_n_adjacent <- nrow(c9_adj)
c9_say(sprintf("  Across both: of the %d crossings in the two adjacent intervals (the whole-series column is not drawn), %d are drawn and %d fall on a significant tile.",
               c9_n_adjacent, sum(c9_heat$is_crossing), sum(c9_heat$is_crossing & !is.na(c9_heat$shown))))
## The heatmaps show only genes the ungated model calls significant in SOME
## adjacent interval, so a crossing in a gene it never calls has no row to be
## drawn on. Those crossings are named here, because a reader who counts the
## dots against Figure 7 will otherwise find the figures short.
c9_not_drawn <- c9_adj %>%
  dplyr::filter(!gene_id %in% c9_heat_genes$gene_id) %>%
  dplyr::left_join(c9_age %>% dplyr::select(genotype, step, gene_id, verdict),
                   by = c("genotype", "step", "gene_id"))
if (nrow(c9_not_drawn)) {
  c9_say(sprintf("  %d crossing(s) are not drawn because the gene is significant in no adjacent interval under the ungated model:",
                 nrow(c9_not_drawn)))
  c9_show(c9_not_drawn %>% dplyr::transmute(genotype, step, symbol, direction, verdict))
}

## 4d. P-value histograms, every fitted contrast, all genes in the model
c9_pv <- dplyr::bind_rows(
  pair_results %>% dplyr::filter(contrast != HINGE[["A"]]) %>%
    dplyr::transmute(test = contrast, pvalue),
  raw_steps %>% dplyr::transmute(test = paste0(genotype, ": ", step), pvalue)) %>%
  dplyr::filter(!is.na(pvalue))
c9_p_pv <- ggplot2::ggplot(c9_pv, ggplot2::aes(pvalue)) +
  ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "grey55", colour = "white",
                          linewidth = 0.2) +
  ggplot2::facet_wrap(~ test, scales = "free_y", ncol = 4) +
  ggplot2::labs(x = "Unadjusted P value", y = "Genes",
                title = "P-value histograms, all genes in each model",
                subtitle = paste("Expected: flat to the right with a spike at zero.",
                                 "A hump near 1 or a U shape signals a dispersion or filtering problem.")) +
  c9_theme
c9_save(c9_p_pv, "QC_PValueHistograms", 8, 2.2 * ceiling(dplyr::n_distinct(c9_pv$test) / 4))

## 4e. Sample-sample correlation over the panel genes
C9_GROUP_ORDER <- c("P3_P4_WT", "P3_P4_Sey", "P15_WT", "P15_Sey",
                    "Adult_WT", "Adult_Sey_T", "Adult_Sey_O")
c9_samples <- meta %>%
  dplyr::filter(tissue == TISSUE, laboratory %in% LABORATORY, !is.na(age_group)) %>%
  dplyr::mutate(group = dplyr::case_when(
    age_group != "Adult" ~ paste0(age_group, "_", genotype),
    genotype == "WT"     ~ "Adult_WT",
    transparency == "T"  ~ "Adult_Sey_T",
    transparency == "O"  ~ "Adult_Sey_O",
    TRUE                 ~ NA_character_)) %>%
  dplyr::filter(!is.na(group), sample_id %in% colnames(count_mat)) %>%
  dplyr::arrange(factor(group, levels = C9_GROUP_ORDER), sample_id)
c9_lcpm <- log2(t(t(count_mat[panel_genes, c9_samples$sample_id, drop = FALSE]) /
                    lib_sizes[c9_samples$sample_id]) * 1e6 + 1)
c9_cor <- stats::cor(c9_lcpm, method = "spearman")
c9_cor_long <- as.data.frame(as.table(c9_cor), stringsAsFactors = FALSE)
names(c9_cor_long) <- c("a", "b", "rho")
c9_cor_long$a <- factor(c9_cor_long$a, levels = c9_samples$sample_id)
c9_cor_long$b <- factor(c9_cor_long$b, levels = rev(c9_samples$sample_id))
c9_bounds <- cumsum(table(factor(c9_samples$group, levels = C9_GROUP_ORDER)))
c9_bounds <- c9_bounds[c9_bounds < nrow(c9_samples)] + 0.5
c9_p_cor <- ggplot2::ggplot(c9_cor_long, ggplot2::aes(a, b, fill = rho)) +
  ggplot2::geom_tile() +
  ggplot2::geom_vline(xintercept = c9_bounds, colour = "black", linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = nrow(c9_samples) + 1 - c9_bounds,
                      colour = "black", linewidth = 0.3) +
  ggplot2::scale_fill_viridis_c(name = "Spearman rho") +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Sample-sample correlation, panel genes (log2 CPM)",
                subtitle = paste("Samples ordered by group:",
                                 paste(C9_GROUP_ORDER, collapse = ", "))) +
  c9_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,
                                                     size = 5),
                 axis.text.y = ggplot2::element_text(size = 5),
                 legend.position = "right")
c9_save(c9_p_cor, "QC_SampleCorrelation", 7.5, 7)

## ---------------------------------------------------------------------------
## 5. Summary
## ---------------------------------------------------------------------------

c9_say("")
c9_say(strrep("=", 74))
c9_say("For Methods, generated from this run:")
c9_say(sprintf("  Of %d developmental status crossings, %d (%.0f%%) were also significant in",
               nrow(c9_age), c9_n_confirmed, 100 * c9_n_confirmed / nrow(c9_age)))
c9_say("  the same direction when the within-genotype model was applied to every")
c9_say("  panel gene without the expression-status gate.")
c9_say(sprintf("  Of %d categorical genotype calls, %d (%.0f%%) were significant in the same",
               nrow(c9_cat), c9_n_cat_ok, 100 * c9_n_cat_ok / nrow(c9_cat)))
c9_say("  direction under the ungated model.")
c9_say("  State both numbers and name the exceptions; a check that is only quoted")
c9_say("  when it agrees is not a check.")
c9_say(strrep("=", 74))

writeLines(c9_log, c9_path("Concordance_Summary.txt"))
message("\n07_concordance complete. Outputs in: ", C9_OUT)
