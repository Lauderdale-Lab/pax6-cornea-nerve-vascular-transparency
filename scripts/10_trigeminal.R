###############################################################################
## 10_trigeminal.R -- trigeminal ganglion, curated panels
##
## Analysis pipeline for:
##   Nerve and vascular abnormalities precede loss of corneal transparency
##   in Pax6-haploinsufficient mice
##   Sneha K. Mohan, James D. Lauderdale
##
## James D. Lauderdale, PhD  (ORCID 0000-0001-7503-0528)
## Department of Cellular Biology, University of Georgia
## Athens, GA 30602, USA
##
## Repository : https://github.com/Lauderdale-Lab/pax6-cornea-nerve-vascular-transparency
## Archived   : https://doi.org/10.5281/zenodo.22965379
## Licence    : MIT (see LICENSE)
## Contact    : James D. Lauderdale, jdlauder@uga.edu
##
## Run the pipeline with run_all.R. Scripts are numbered in execution order and
## share one R session by design; see run_all.R for why.
###############################################################################

###############################################################################
## 10_trigeminal.R
##
## Asks whether the trigeminal ganglia of the same animals show genotype-
## associated differences in the three curated panels, and whether the
## libraries could have shown them.
##
## WHY NOT THE THREE-STATE FRAMEWORK. Every trigeminal group holds two
## libraries, and 03_expression_status.R declines to call any group smaller
## than MIN_GROUP_N (3), because the leave-one-out stability check is
## meaningless below that. The trigeminal arm is therefore a conventional
## replicate-aware analysis, with the corneal conventions applied wherever
## they carry over:
##   * the same pre-filter (>= 10 counts in >= 2 libraries of the fit);
##   * one factor combining age, genotype and transparency, as in the corneal
##     developmental and opacity models;
##   * unshrunk fold changes with 95% confidence intervals;
##   * an equivalence test at the same bound, so a non-significant gene is
##     reported as bounded or undetermined, never as unchanged;
##   * Benjamini-Hochberg correction within primary panel and contrast, with
##     no independent filtering, each gene corrected once.
##
## TWO DEPARTURES FROM THE CORNEAL MODELS, both forced by the design:
##   * No batch term. All trigeminal libraries were sequenced in one batch.
##   * No male-fraction covariate by default. With two libraries per group it
##     would spend most of the remaining residual degrees of freedom. Each
##     library's male fraction is written to the sample table so the balance
##     can be read directly. TRIGEMINAL_SEX_COVARIATE turns it on.
##
## MODELS
##   D  developmental : WT and transparent-mutant ganglia at P3/P4, P15 and
##                      adult; genotype compared within each age.
##   O  opacity       : adult WT, transparent-mutant and opaque-mutant
##                      ganglia; all three pairs, plus an omnibus LRT. The
##                      adult transparent-vs-WT contrast is reported from this
##                      model, as in the cornea.
##   C  positive control : WT ganglia across the three ages, LRT, all genes
##                      passing the pre-filter. The headline trigeminal result
##                      is a negative; C shows whether the same libraries and
##                      code detect differential expression where it exists.
##
## TISSUE MARKERS. Mbp and Mag are made by both Schwann cells and central
## oligodendrocytes, so a ganglion carrying more nerve root or brainstem would
## show them higher for a reason unrelated to genotype. Every library's CPM for
## the markers in TRIGEMINAL_TISSUE_MARKERS is written out so that reading can
## be checked against the data rather than assumed.
##
## Requires: 00_config.R, 01_load.R (meta, count_mat, lib_sizes, panel_tbl),
##           01c_contamination_indices.R (ot_ids_by_symbol) and
##           04_differential.R (panel_primary).
## Outputs (results/Trigeminal/):
##   Trigeminal_SampleTable.csv      libraries used, with male fraction
##   Trigeminal_TissueMarkers.csv    central and Schwann-cell marker CPM per library
##   Trigeminal_PanelResults.csv     one row per gene per panel per contrast
##   Trigeminal_PanelSummary.csv     outcomes by panel and contrast
##   Trigeminal_GenomeWide.csv       every contrast over all filtered genes
##   Trigeminal_Summary.txt          the numbers quoted in the manuscript
##
## STATISTICAL BASIS
##   [1] DESeq2 negative-binomial GLM, Wald and likelihood-ratio tests.
##       Love MI, Huber W, Anders S (2014) Genome Biol 15:550.
##   [2] Equivalence by DESeq2's lfcThreshold with altHypothesis = "lessAbs",
##       as in 04_differential.R.
##   [3] Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
###############################################################################

if (!exists("count_mat") || !exists("panel_tbl")) {
  stop("Source 00_config.R and 01_load.R before 10_trigeminal.R.")
}
if (!exists("ot_ids_by_symbol")) {
  stop("10_trigeminal.R needs ot_ids_by_symbol from 01c_contamination_indices.R; ",
       "run the pipeline with run_all.R.")
}
if (!exists("panel_primary")) {
  stop("10_trigeminal.R needs panel_primary from 04_differential.R; run the ",
       "pipeline with run_all.R.")
}

suppressPackageStartupMessages({
  library(DESeq2)
})

TG_OUT <- out_path("Trigeminal")
dir.create(TG_OUT, showWarnings = FALSE, recursive = TRUE)
tg_path <- function(...) file.path(TG_OUT, ...)

tg_log <- character(0)
say <- function(...) {
  txt <- paste0(...)
  message(txt)
  tg_log <<- c(tg_log, txt)
}

## ---------------------------------------------------------------------------
## 1. Libraries
## ---------------------------------------------------------------------------

tg <- meta %>%
  dplyr::filter(tissue == TRIGEMINAL_TISSUE,
                laboratory %in% LABORATORY,
                sample_id %in% colnames(count_mat),
                !is.na(age_group)) %>%
  dplyr::mutate(
    tg_group = dplyr::case_when(
      genotype == "WT"                        ~ "WT",
      genotype == "Sey" & transparency == "T" ~ "SeyT",
      genotype == "Sey" & transparency == "O" ~ "SeyO",
      TRUE                                    ~ NA_character_),
    cell = paste(age_group, tg_group, sep = "_"))

if (!nrow(tg)) stop("No trigeminal libraries found in the metadata.")
if (any(is.na(tg$tg_group))) {
  stop("Trigeminal libraries without a genotype/transparency group: ",
       paste(tg$sample_id[is.na(tg$tg_group)], collapse = ", "))
}

say("Trigeminal libraries: ", nrow(tg))
say(paste(utils::capture.output(print(table(tg$age_group, tg$tg_group))),
          collapse = "\n"))
if (nlevels(droplevels(factor(tg$batch))) > 1) {
  say("NOTE: trigeminal libraries span more than one batch; this script does ",
      "not model batch and should be revisited.")
} else {
  say("All trigeminal libraries are from batch ", as.character(tg$batch[1]),
      "; no batch term is fitted.")
}

readr::write_csv(
  tg %>% dplyr::select(sample_id, age_group, genotype, transparency, tg_group,
                       batch, male_frac) %>%
    dplyr::arrange(age_group, tg_group, sample_id),
  tg_path("Trigeminal_SampleTable.csv"))

## ---------------------------------------------------------------------------
## 1b. Tissue markers: central carry-over and Schwann-cell content
## ---------------------------------------------------------------------------

tm <- TRIGEMINAL_TISSUE_MARKERS
tm_missing <- setdiff(tm$symbol, names(ot_ids_by_symbol))
if (length(tm_missing)) {
  stop("Trigeminal tissue marker(s) not found among the count-matrix genes: ",
       paste(tm_missing, collapse = ", "),
       ". Correct TRIGEMINAL_TISSUE_MARKERS in 00_config.R.")
}
tm$gene_id <- vapply(tm$symbol, function(s) ot_ids_by_symbol[[s]][1], character(1))

tg_ordered <- tg %>% dplyr::arrange(age_group, factor(tg_group, c("WT", "SeyT", "SeyO")),
                                    sample_id)
tm_cpm <- t(t(count_mat[tm$gene_id, tg_ordered$sample_id, drop = FALSE]) /
              lib_sizes[tg_ordered$sample_id]) * 1e6
rownames(tm_cpm) <- tm$symbol

tissue_markers <- tg_ordered %>%
  dplyr::select(sample_id, age_group, tg_group) %>%
  dplyr::bind_cols(tibble::as_tibble(round(t(tm_cpm), 2)))
readr::write_csv(tissue_markers, tg_path("Trigeminal_TissueMarkers.csv"))

say("")
say("Tissue markers, CPM (", paste(tm$symbol, tm$gene_id, sep = " = ", collapse = "; "), "):")
say(paste(utils::capture.output(print(as.data.frame(tissue_markers), row.names = FALSE)),
          collapse = "\n"))

## ---------------------------------------------------------------------------
## 2. Fitting
## ---------------------------------------------------------------------------

fit_model <- function(samples, label) {
  cm <- count_mat[, samples$sample_id, drop = FALSE]
  cm <- cm[rowSums(cm >= 10) >= 2, , drop = FALSE]
  cd <- as.data.frame(samples)
  rownames(cd) <- cd$sample_id
  cd$cell <- factor(cd$cell)
  terms <- c(if (isTRUE(TRIGEMINAL_SEX_COVARIATE)) "male_frac", "cell")
  fml <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- stats::model.matrix(fml, cd)
  if (qr(mm)$rank < ncol(mm)) stop("Model ", label, " design is not full rank.")
  say("")
  say("[", label, "] ", ncol(cm), " libraries, ", nrow(cm),
      " genes after pre-filter, design ", deparse(fml))
  dds <- DESeq2::DESeq(
    DESeq2::DESeqDataSetFromMatrix(cm, cd, design = fml), quiet = TRUE)
  list(dds = dds, fml = fml, n = ncol(cm), label = label)
}

## Wald test and matching equivalence test for one contrast. [1] [2]
test_pair <- function(fit, label, a, b) {
  arg <- c("cell", a, b)
  res <- DESeq2::results(fit$dds, contrast = arg, alpha = FDR)
  eq  <- DESeq2::results(fit$dds, contrast = arg, alpha = FDR,
                         lfcThreshold = EQUIV_BOUND, altHypothesis = "lessAbs")
  stopifnot(identical(rownames(res), rownames(eq)))
  tibble::tibble(model = fit$label, contrast = label, gene_id = rownames(res),
                 baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
                 lfcSE = res$lfcSE, pvalue = res$pvalue,
                 pvalue_equiv = eq$pvalue)
}

transparent <- tg %>% dplyr::filter(tg_group %in% c("WT", "SeyT"))
adult       <- tg %>% dplyr::filter(age_group == "Adult")
wt          <- tg %>% dplyr::filter(tg_group == "WT")

fit_D <- fit_model(transparent, "D developmental")
fit_O <- fit_model(adult, "O opacity")

pairs <- dplyr::bind_rows(
  test_pair(fit_D, "SeyT_vs_WT_P3_P4", "P3_P4_SeyT", "P3_P4_WT"),
  test_pair(fit_D, "SeyT_vs_WT_P15",   "P15_SeyT",   "P15_WT"),
  test_pair(fit_O, "SeyT_vs_WT_Adult", "Adult_SeyT", "Adult_WT"),
  test_pair(fit_O, "SeyO_vs_WT_Adult", "Adult_SeyO", "Adult_WT"),
  test_pair(fit_O, "SeyO_vs_SeyT_Adult", "Adult_SeyO", "Adult_SeyT"))

## Omnibus: any difference among the three adult groups.
reduced_O <- if (isTRUE(TRIGEMINAL_SEX_COVARIATE)) ~ male_frac else ~ 1
dds_O_lrt <- DESeq2::DESeq(fit_O$dds, test = "LRT", reduced = reduced_O, quiet = TRUE)
res_O_lrt <- DESeq2::results(dds_O_lrt, alpha = FDR)
omnibus <- tibble::tibble(model = fit_O$label, contrast = "Adult_omnibus_LRT",
                          gene_id = rownames(res_O_lrt),
                          baseMean = res_O_lrt$baseMean,
                          log2FoldChange = NA_real_, lfcSE = NA_real_,
                          pvalue = res_O_lrt$pvalue, pvalue_equiv = NA_real_)

## Positive control: WT age effect, every filtered gene. [1]
fit_C <- fit_model(wt, "C positive control")
reduced_C <- if (isTRUE(TRIGEMINAL_SEX_COVARIATE)) ~ male_frac else ~ 1
dds_C_lrt <- DESeq2::DESeq(fit_C$dds, test = "LRT", reduced = reduced_C, quiet = TRUE)
res_C <- DESeq2::results(dds_C_lrt, alpha = FDR)
control_padj <- stats::p.adjust(res_C$pvalue, method = "BH")
control <- tibble::tibble(
  contrast   = "WT_age_LRT (positive control)",
  n_tested   = sum(!is.na(res_C$pvalue)),
  n_significant = sum(control_padj < FDR, na.rm = TRUE),
  min_padj   = suppressWarnings(min(control_padj, na.rm = TRUE)))

all_tests <- dplyr::bind_rows(pairs, omnibus)

## ---------------------------------------------------------------------------
## 3. Panel results, corrected within primary panel and contrast [3]
## ---------------------------------------------------------------------------

## Every panel gene appears in every contrast, so a gene removed by the
## pre-filter is reported as not tested rather than silently absent.
contrast_models <- dplyr::distinct(all_tests, contrast, model)
panel_rows <- panel_tbl %>%
  dplyr::select(panel, gene_id, symbol) %>%
  dplyr::cross_join(contrast_models) %>%
  dplyr::left_join(dplyr::select(all_tests, -model),
                   by = c("contrast", "gene_id")) %>%
  dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel),
                   by = "symbol") %>%
  dplyr::mutate(primary_panel = dplyr::coalesce(primary_panel, panel))

corrections <- panel_rows %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::group_by(correction_panel = panel, contrast) %>%
  dplyr::mutate(padj = stats::p.adjust(pvalue, method = "BH"),
                padj_equiv = stats::p.adjust(pvalue_equiv, method = "BH"),
                n_in_correction_family = sum(!is.na(pvalue))) %>%
  dplyr::ungroup() %>%
  dplyr::select(contrast, gene_id, correction_panel, padj, padj_equiv,
                n_in_correction_family)

panel_results <- panel_rows %>%
  dplyr::left_join(corrections, by = c("contrast", "gene_id")) %>%
  dplyr::mutate(
    ci_low  = log2FoldChange - stats::qnorm(1 - FDR / 2) * lfcSE,
    ci_high = log2FoldChange + stats::qnorm(1 - FDR / 2) * lfcSE,
    outcome = dplyr::case_when(
      is.na(pvalue)                               ~ "not tested (below pre-filter)",
      !is.na(padj) & padj < FDR                   ~ "changed",
      contrast == "Adult_omnibus_LRT"             ~ "no difference detected",
      !is.na(padj_equiv) & padj_equiv < FDR       ~ "no change detected (bounded)",
      TRUE                                        ~ "undetermined")) %>%
  dplyr::arrange(panel, contrast, padj)
readr::write_csv(panel_results, tg_path("Trigeminal_PanelResults.csv"))

panel_summary <- panel_results %>%
  dplyr::group_by(panel, contrast) %>%
  dplyr::summarise(
    n_panel_genes = dplyr::n(),
    n_tested      = sum(!is.na(pvalue)),
    changed       = sum(outcome == "changed"),
    bounded       = sum(outcome == "no change detected (bounded)"),
    undetermined  = sum(outcome %in% c("undetermined", "no difference detected")),
    min_padj      = suppressWarnings(min(padj, na.rm = TRUE)),
    .groups = "drop") %>%
  dplyr::mutate(min_padj = ifelse(is.finite(min_padj), min_padj, NA_real_))
readr::write_csv(panel_summary, tg_path("Trigeminal_PanelSummary.csv"))

## ---------------------------------------------------------------------------
## 4. Genome-wide context for every contrast
## ---------------------------------------------------------------------------

genome_wide <- all_tests %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(padj_gw = stats::p.adjust(pvalue, method = "BH")) %>%
  dplyr::summarise(n_tested = sum(!is.na(pvalue)),
                   n_significant = sum(padj_gw < FDR, na.rm = TRUE),
                   min_padj = suppressWarnings(min(padj_gw, na.rm = TRUE)),
                   .groups = "drop") %>%
  dplyr::bind_rows(control) %>%
  dplyr::mutate(min_padj = ifelse(is.finite(min_padj), min_padj, NA_real_))
readr::write_csv(genome_wide, tg_path("Trigeminal_GenomeWide.csv"))

## ---------------------------------------------------------------------------
## 5. The numbers the manuscript quotes
## ---------------------------------------------------------------------------

## Distinct genes, not panel rows: a gene on two panels is counted once.
distinct_sig <- function(ctr) {
  panel_results %>%
    dplyr::filter(contrast %in% ctr, outcome == "changed") %>%
    dplyr::distinct(contrast, gene_id, symbol)
}
min_padj_of <- function(ctr) {
  v <- panel_results$padj[panel_results$contrast == ctr]
  if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
}
genotype_ctr <- c("SeyT_vs_WT_P3_P4", "SeyT_vs_WT_P15", "SeyT_vs_WT_Adult")
opacity_ctr  <- c("SeyO_vs_WT_Adult", "SeyO_vs_SeyT_Adult", "Adult_omnibus_LRT")

say("")
say(strrep("=", 74))
say("TRIGEMINAL -- numbers for the manuscript")
say(strrep("=", 74))
say("Libraries per group:")
say(paste(utils::capture.output(print(table(tg$age_group, tg$tg_group))),
          collapse = "\n"))
say("")
say("Genotype comparisons of transparent tissue (panel genes, distinct):")
for (ctr in genotype_ctr) {
  sig <- distinct_sig(ctr)
  say(sprintf("  %-20s %d significant; smallest adjusted P %s%s", ctr, nrow(sig),
              signif(min_padj_of(ctr), 3),
              if (nrow(sig)) paste0(" (", paste(sig$symbol, collapse = ", "), ")") else ""))
}
say("Comparisons involving opaque animals (panel genes, distinct):")
for (ctr in opacity_ctr) {
  sig <- distinct_sig(ctr)
  say(sprintf("  %-20s %d significant%s", ctr, nrow(sig),
              if (nrow(sig)) paste0(" (", paste(sig$symbol, collapse = ", "), ")") else ""))
}
tested_pairs <- panel_results %>%
  dplyr::filter(contrast %in% genotype_ctr, !is.na(pvalue)) %>%
  dplyr::distinct(contrast, gene_id, outcome)
say("")
say("Outcomes over the genotype comparisons of transparent tissue (distinct ",
    "gene x contrast): ",
    paste(names(table(tested_pairs$outcome)), table(tested_pairs$outcome),
          sep = " ", collapse = "; "))
say("")
say("Positive control (WT age effect, all filtered genes): ",
    control$n_significant, " of ", control$n_tested,
    " genes significant at adjusted P < ", FDR, ".")
say("")
say("For Methods: trigeminal models carry no batch term (single batch)",
    if (isTRUE(TRIGEMINAL_SEX_COVARIATE)) "." else
      " and no male-fraction covariate (two libraries per group).")

writeLines(tg_log, tg_path("Trigeminal_Summary.txt"))
message("\nTrigeminal outputs: ", TG_OUT)
