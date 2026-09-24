###############################################################################
## 06_developmental.R -- developmental trajectories within genotype
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
## Repository : https://github.com/Lauderdale-Lab/pax6-mouse-cornea-trigeminal-genesets
## Archived   : <ZENODO_DOI>
## Licence    : MIT (see LICENSE)
## Contact    : James D. Lauderdale, jdlauder@uga.edu
##
## Run the pipeline with run_all.R. Scripts are numbered in execution order and
## share one R session by design; see run_all.R for why.
###############################################################################

###############################################################################
## 06_developmental.R
##
## Describes how each curated-panel gene behaves across developmental time,
## separately within each genotype.
##
## WHAT THIS SCRIPT DELIBERATELY DOES NOT DO
##
## It does not sort genes into a vocabulary of trajectory shapes and then report
## how many fit each one. An earlier version of this analysis classified every
## gene as increasing, decreasing, peaking at P15, troughing at P15, stable, or
## undetermined. Two things were wrong with that. The shapes were chosen in
## advance, so a gene doing something else could only be recorded as a failure
## to fit rather than as what it was. And "stable" required all three intervals
## to be individually bounded, a condition so demanding that it was met by at
## most one gene in any panel -- a category that never fires is not a finding
## about the tissue, it is a property of the rule.
##
## What is reported instead is the outcome of each developmental interval on its
## own terms: changed, bounded, or undetermined, exactly as in the genotype
## contrasts. A shape description is derived from those outcomes mechanically at
## the end, and it asserts a peak or a trough ONLY when both adjacent intervals
## are individually significant and run in opposite directions. Everything else
## is described by what the intervals actually said.
##
## It does not compare trajectories between genotypes. The earlier version
## produced a table of genes whose class differed between wild type and mutant,
## and that table was dominated by genes moving into and out of "undetermined" --
## differences in resolution, not in biology. A gene can change class on a small
## shift in one interval. The question "does the genotype effect change with
## age" has a proper test, it is a difference of differences, and it is fitted
## in 04_differential.R. This script points there rather than offering a
## substitute that reads like an answer.
##
## THE MUTANT TRAJECTORY USES TRANSPARENT ADULTS ONLY. Opaque adult cornea
## differs from transparent adult cornea in the same genotype, so including
## opaque libraries in a developmental series would confound age with disease
## state and the resulting "developmental" change would partly be opacity. The
## opacity axis is a separate comparison and is fitted in 04_differential.R.
##
## DEVELOPMENTAL STATUS CROSSINGS are reported alongside the quantitative
## results. A gene that is confidently silent at one age and confidently
## expressed at another has changed categorically, and no fold change is needed
## to say so. The criterion is the same one used for the genotype contrasts:
## one age confidently expressed, and Poisson intervals separated by at least
## the configured fold difference.
##
## STATISTICAL BASIS
##   [1] Negative-binomial generalised linear models with empirical-Bayes
##       dispersion shrinkage; Wald tests on the fitted coefficients.
##       Love MI, Huber W, Anders S (2014) Genome Biology 15:550.
##   [2] Equivalence testing against the composite null |log2 fold change| >=
##       bound, so that a negative interval is reported as bounded rather than
##       as absence. Schuirmann DJ (1987) J Pharmacokinet Biopharm 15:657-680.
##   [3] Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300, applied
##       within primary panel, genotype and interval.
##   [4] Exact Poisson limits for the status crossings, as in
##       03_expression_status.R. Garwood F (1936) Biometrika 28:437-442.
##
## Inputs : 00_config.R, 01_load.R, 03_expression_status.R, 04_differential.R
##          -- run in that order in one session.
## COUNTING CONVENTION. Every table here is one row per gene per PANEL it
## belongs to, matching 03 and 04, so a gene on two panels appears twice. That
## is right for a per-panel table and wrong for a total: summing a per-panel
## table across panels double-counts the thirty shared genes. Totals printed by
## this script are therefore taken over distinct crossings, and both numbers are
## shown.
##
## Outputs: Developmental_Steps.csv, Developmental_Summary.csv,
##          Developmental_Trajectories.csv, Developmental_StatusCrossings.csv
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(purrr)
  library(DESeq2)
})

if (!exists("cell_stats")) stop("Source 03_expression_status.R first.")
if (!exists("panel_primary")) {
  stop("Source 04_differential.R first -- this script uses the same primary-panel ",
       "assignment so that a gene on two panels is corrected once here too.")
}

## ---------------------------------------------------------------------------
## 1. The developmental series, per genotype
## ---------------------------------------------------------------------------

## Every cornea at P3/P4 and P15 is transparent, and mutant adults are
## restricted to transparent so that age is not confounded with opacity.
SERIES <- list(
  WT  = c(P3_P4 = "P3_P4:WT",  P15 = "P15:WT",  Adult = "Adult:WT"),
  Sey = c(P3_P4 = "P3_P4:Sey", P15 = "P15:Sey", Adult = "Adult:Sey"))

## Ordered intervals. Named DEV_INTERVALS rather than anything shorter because
## this script is sourced into the global environment alongside a runner that
## keeps its own step table; two objects called STEPS in one session is a
## collision waiting to happen.
## The two adjacent ones describe the shape; the third spans
## the whole series and is reported but never used to infer a shape, because a
## gene that rises and falls back can look unchanged across the span.
DEV_INTERVALS <- tibble::tribble(
  ~step,              ~later,   ~earlier, ~adjacent,
  "P15_vs_P3_P4",     "P15",    "P3_P4",  TRUE,
  "Adult_vs_P15",     "Adult",  "P15",    TRUE,
  "Adult_vs_P3_P4",   "Adult",  "P3_P4",  FALSE)

missing_cells <- setdiff(unlist(SERIES), names(cell_stats))
if (length(missing_cells)) {
  stop("Developmental series needs group(s) that were not called: ",
       paste(missing_cells, collapse = ", "))
}

dev_pool <- meta %>%
  dplyr::filter(tissue == TISSUE, laboratory %in% LABORATORY,
                transparency == "T", !is.na(age_group),
                !sample_id %in% DE_EXCLUDE)

fit_genotype <- function(gt) {
  samples <- dev_pool %>% dplyr::filter(genotype == gt)
  col_data <- as.data.frame(samples %>% tibble::column_to_rownames("sample_id"))
  col_data$age_group <- stats::relevel(factor(col_data$age_group, levels = AGE_LEVELS),
                                       ref = "P3_P4")
  col_data$batch <- droplevels(factor(col_data$batch))

  terms <- c(if (nlevels(col_data$batch) > 1) "batch",
             if (USE_SEX_COVARIATE) "male_frac", "age_group")
  fml <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- try(stats::model.matrix(fml, col_data), silent = TRUE)
  if (inherits(mm, "try-error") || qr(mm)$rank < ncol(mm)) {
    stop("Developmental design for ", gt, " is not full rank: ", deparse(fml))
  }

  cm <- count_mat[, rownames(col_data), drop = FALSE]
  cm <- cm[rowSums(cm >= 10) >= 2, , drop = FALSE]

  message("  [", gt, "] ", ncol(cm), " libraries, ", nrow(cm), " genes  design ",
          deparse(fml))
  print(table(col_data$age_group))

  DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(cm, col_data, design = fml),
                quiet = TRUE)
}

message("\n=== Developmental models, fitted within genotype ===")
## Named dev_fits, not fits: 04_differential.R leaves an object called `fits`
## in the shared session, and this script used to overwrite it. Nothing read
## the clobbered copy, but a later script reaching for 04's models by that
## name would have got these instead.
dev_fits <- stats::setNames(lapply(names(SERIES), fit_genotype), names(SERIES))

## ---------------------------------------------------------------------------
## 2. Which intervals a gene can be tested across
## ---------------------------------------------------------------------------

## A quantitative comparison is only interpretable where the gene is confidently
## expressed at BOTH ages. Where it is not, the interval is either a status
## crossing (section 4) or simply unresolvable, and it is reported as such
## rather than tested and reported as null.
step_assessability <- purrr::map_dfr(names(SERIES), function(gt) {
  purrr::pmap_dfr(list(DEV_INTERVALS$step, DEV_INTERVALS$later, DEV_INTERVALS$earlier),
                  function(step, later, earlier) {
    A <- cell_stats[[SERIES[[gt]][[later]]]]
    B <- cell_stats[[SERIES[[gt]][[earlier]]]]
    bounded_fold <- pmax(A$lcl / B$ucl, B$lcl / A$ucl)
    one_expressed <- xor(A$call == "expressed", B$call == "expressed")
    tibble::tibble(
      genotype = gt, step = step, gene_id = A$gene_id,
      call_later = A$call, call_earlier = B$call,
      cpm_later = A$pooled_cpm, cpm_earlier = B$pooled_cpm,
      bounded_fold = round(bounded_fold, 3),
      assessable = dplyr::case_when(
        A$call == "expressed" & B$call == "expressed" ~ "test quantitatively",
        (A$call == "expressed" & B$call == "not_expressed") |
        (A$call == "not_expressed" & B$call == "expressed")
          ~ "status crossing (absent at one age)",
        A$call == "not_expressed" & B$call == "not_expressed" ~ "off at both ages",
        one_expressed & bounded_fold >= CATEGORICAL_MIN_FOLD
          ~ "status crossing (interval separated)",
        TRUE ~ "cannot assess"))
  })
})

## ---------------------------------------------------------------------------
## 3. Quantitative outcome per interval
## ---------------------------------------------------------------------------

test_step <- function(gt, step, later, earlier) {
  dds <- dev_fits[[gt]]
  res <- DESeq2::results(dds, contrast = c("age_group", later, earlier), alpha = FDR)
  eq  <- DESeq2::results(dds, contrast = c("age_group", later, earlier), alpha = FDR,
                         lfcThreshold = EQUIV_BOUND, altHypothesis = "lessAbs")
  stopifnot(identical(rownames(res), rownames(eq)))
  tibble::tibble(genotype = gt, step = step, gene_id = rownames(res),
                 log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                 pvalue = res$pvalue, pvalue_equiv = eq$pvalue)
}

raw_steps <- purrr::pmap_dfr(
  tidyr::expand_grid(genotype = names(SERIES), DEV_INTERVALS) %>%
    dplyr::select(genotype, step, later, earlier),
  function(genotype, step, later, earlier) test_step(genotype, step, later, earlier))

testable <- step_assessability %>%
  dplyr::filter(assessable == "test quantitatively") %>%
  dplyr::select(genotype, step, gene_id) %>%
  dplyr::inner_join(panel_tbl, by = "gene_id", relationship = "many-to-many") %>%
  dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel), by = "symbol")

steps_joined <- testable %>%
  dplyr::inner_join(raw_steps, by = c("genotype", "step", "gene_id"),
                    relationship = "many-to-many")

## Corrected once per gene, in its primary panel, within genotype and interval.
## The value is then shown wherever the gene appears. [3]
corrections <- steps_joined %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::group_by(correction_panel = panel, genotype, step) %>%
  dplyr::mutate(padj = stats::p.adjust(pvalue, method = "BH"),
                padj_equiv = stats::p.adjust(pvalue_equiv, method = "BH"),
                n_in_correction_family = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(genotype, step, gene_id, correction_panel, padj, padj_equiv,
                n_in_correction_family)

dev_steps <- steps_joined %>%
  dplyr::left_join(corrections, by = c("genotype", "step", "gene_id")) %>%
  dplyr::mutate(
    ci_low  = log2FoldChange - stats::qnorm(1 - FDR / 2) * lfcSE,
    ci_high = log2FoldChange + stats::qnorm(1 - FDR / 2) * lfcSE,
    largest_effect_not_excluded = round(pmax(abs(ci_low), abs(ci_high)), 3),
    direction = dplyr::case_when(log2FoldChange > 0 ~ "up", log2FoldChange < 0 ~ "down",
                                 TRUE ~ NA_character_),
    outcome = dplyr::case_when(
      is.na(padj) & is.na(padj_equiv)       ~ "not testable",
      !is.na(padj) & padj < FDR             ~ "changed",
      !is.na(padj_equiv) & padj_equiv < FDR ~ "no change detected (bounded)",
      TRUE                                  ~ "undetermined")) %>%
  dplyr::arrange(panel, genotype, step, padj)
readr::write_csv(dev_steps, out_path("Developmental_Steps.csv"))

dev_summary <- dev_steps %>%
  dplyr::count(panel, genotype, step, outcome) %>%
  tidyr::pivot_wider(names_from = outcome, values_from = n, values_fill = 0)
readr::write_csv(dev_summary, out_path("Developmental_Summary.csv"))
message("\n=== Outcome per developmental interval ===")
print(as.data.frame(dev_summary))

## What an undetermined interval could still be hiding, so that the largest
## category in this analysis is reported with a magnitude rather than as a gap.
undet <- dev_steps %>% dplyr::filter(outcome == "undetermined")
if (nrow(undet)) {
  message("\n=== Size of what `undetermined` leaves open, in log2 units ===")
  print(as.data.frame(undet %>% dplyr::group_by(panel, genotype) %>%
    dplyr::summarise(n = dplyr::n(),
                     median_effect_not_excluded = round(stats::median(
                       largest_effect_not_excluded, na.rm = TRUE), 2),
                     max_effect_not_excluded = round(max(
                       largest_effect_not_excluded, na.rm = TRUE), 2),
                     .groups = "drop")))
  message("  These genes are not evidence of developmental stability. Each row's")
  message("  median is the fold change that could not be ruled out for a typical")
  message("  gene in that group. Compare it with the equivalence bound of ",
          EQUIV_BOUND, ".")
}

## ---------------------------------------------------------------------------
## 4. Status crossings across development
## ---------------------------------------------------------------------------

## One row per panel MEMBERSHIP, as everywhere else in this pipeline: a gene on
## two panels appears under both, so that a reader consulting one panel's table
## sees it. `primary_panel` is carried alongside so that any count taken ACROSS
## panels can be de-duplicated without re-deriving the assignment -- see the
## totals below, and 07_figures.R, which colours one row per gene.
crossings <- step_assessability %>%
  dplyr::filter(startsWith(assessable, "status crossing")) %>%
  dplyr::inner_join(panel_tbl, by = "gene_id", relationship = "many-to-many") %>%
  dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel),
                   by = "symbol") %>%
  dplyr::mutate(direction = ifelse(cpm_later > cpm_earlier,
                                   "on with age", "off with age")) %>%
  dplyr::select(genotype, step, panel, primary_panel, symbol, gene_id, assessable,
                direction, call_earlier, call_later, cpm_earlier, cpm_later,
                bounded_fold) %>%
  dplyr::arrange(panel, genotype, step, dplyr::desc(bounded_fold))
readr::write_csv(crossings, out_path("Developmental_StatusCrossings.csv"))

message("\n=== Status crossings: genes switching on or off with age ===")
if (nrow(crossings)) {
  print(as.data.frame(crossings %>% dplyr::count(genotype, step, panel, direction)))

  ## THE TABLE ABOVE COUNTS PANEL MEMBERSHIP AND MUST NOT BE SUMMED. Adding its
  ## cells across panels counts a shared gene once per panel it sits on. Any
  ## total quoted outside a single panel is taken here instead, over distinct
  ## gene x genotype x interval crossings. The two numbers are printed together
  ## so that a figure in the paper can never be traced to the wrong one.
  distinct_crossings <- crossings %>%
    dplyr::distinct(genotype, step, symbol, direction)
  message("  ", nrow(crossings), " panel-membership row(s); ",
          nrow(distinct_crossings), " distinct crossing(s) over ",
          dplyr::n_distinct(distinct_crossings$symbol), " gene(s).")
  message("  ", sum(distinct_crossings$direction == "on with age"),
          " turn on with age, ",
          sum(distinct_crossings$direction == "off with age"), " turn off.",
          "  (Distinct crossings, not panel rows.)")
} else {
  message("  None at a minimum fold difference of ", CATEGORICAL_MIN_FOLD, ".")
}

## ---------------------------------------------------------------------------
## 5. Shape, derived from the intervals rather than imposed on them
## ---------------------------------------------------------------------------

## The two adjacent intervals are collapsed into a literal pair of outcomes, and
## the description is a direct translation of that pair. A peak or a trough is
## asserted only when both adjacent intervals are individually significant and
## run in opposite directions -- the only circumstance in which the data show a
## reversal rather than merely permit one.
## `primary_panel` is carried through for the same reason the crossings table
## carries it: this is one row per gene per PANEL, so a shared gene appears
## twice with the same shape, and any count taken across panels needs to be
## de-duplicated. Without the column that has to be re-derived downstream.
adjacent <- dev_steps %>%
  dplyr::filter(step %in% DEV_INTERVALS$step[DEV_INTERVALS$adjacent]) %>%
  dplyr::select(panel, primary_panel, genotype, symbol, gene_id, step, outcome,
                direction, largest_effect_not_excluded)

shape_of <- function(o1, d1, o2, d2) {
  code <- function(o, d) dplyr::case_when(
    is.na(o) ~ "absent", o == "changed" ~ d,
    o == "no change detected (bounded)" ~ "flat", TRUE ~ "unresolved")
  a <- code(o1, d1); b <- code(o2, d2)
  dplyr::case_when(
    a == "absent" | b == "absent"      ~ "not assessable across both intervals",
    a == "up"   & b == "up"            ~ "rises across both intervals",
    a == "down" & b == "down"          ~ "falls across both intervals",
    a == "up"   & b == "down"          ~ "rises then falls (peak at P15)",
    a == "down" & b == "up"            ~ "falls then rises (trough at P15)",
    a == "flat" & b == "flat"          ~ "no change detected across either interval",
    a %in% c("up", "down") & b == "flat" ~ paste0(a, " to P15 only"),
    a == "flat" & b %in% c("up", "down") ~ paste0(b, " after P15 only"),
    TRUE                               ~ "not resolved")
}

trajectories <- adjacent %>%
  dplyr::select(-largest_effect_not_excluded) %>%
  tidyr::pivot_wider(names_from = step, values_from = c(outcome, direction)) %>%
  dplyr::mutate(shape = shape_of(
    .data[["outcome_P15_vs_P3_P4"]], .data[["direction_P15_vs_P3_P4"]],
    .data[["outcome_Adult_vs_P15"]], .data[["direction_Adult_vs_P15"]])) %>%
  dplyr::arrange(panel, genotype, shape, symbol)
readr::write_csv(trajectories, out_path("Developmental_Trajectories.csv"))

message("\n=== Shape, translated from the two adjacent intervals ===")
print(as.data.frame(trajectories %>% dplyr::count(panel, genotype, shape) %>%
  tidyr::pivot_wider(names_from = genotype, values_from = n, values_fill = 0)))
message("  Counted by panel membership: ", nrow(trajectories), " rows over ",
        dplyr::n_distinct(paste(trajectories$genotype, trajectories$symbol)),
        " distinct gene x genotype pairs. A shared gene appears under each")
message("  panel it belongs to, with the same shape, so THIS TABLE MUST NOT BE")
message("  SUMMED ACROSS PANELS. Filter on panel == primary_panel for a total.")
message("  These labels carry no model. Each is a restatement of what the two")
message("  adjacent intervals returned, and `not resolved` means at least one of")
message("  them was undetermined -- which is a statement about this experiment,")
message("  not about the gene.")

message("\n", strrep("=", 74))
message("06_developmental complete. Outputs in: ", OUT_ROOT)
message("")
message("DO NOT COMPARE THESE SHAPES BETWEEN GENOTYPES. A shape is a threshold")
message("summary, so a gene can change shape on a small shift in one interval,")
message("and most cross-genotype differences in shape are differences in")
message("resolution. Whether the genotype effect itself changes with age is")
message("tested directly in 04_differential.R, as a difference of differences,")
message("and InteractionResults.csv is the file that answers it.")
message(strrep("=", 74))
