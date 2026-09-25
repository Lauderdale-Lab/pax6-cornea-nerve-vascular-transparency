###############################################################################
## 03_expression_status.R -- expressed / not expressed / indeterminate, and assessability
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
## 03_expression_status.R
##
## Assigns every curated-panel gene, in every group, one of three states --
## expressed, not expressed, or indeterminate -- and then decides, for each
## contrast, what kind of statement the data can support about that gene.
##
## The three-state design exists because a two-state rule is forced to guess on
## genes the data cannot resolve, and the guess is what moves when a single
## library is dropped. Making "cannot tell" a reported outcome converts a hidden
## bias into a stated one.
##
## WHAT THIS SCRIPT DECIDES, per gene and contrast:
##
##   test quantitatively                confidently expressed in both groups;
##                                      passed to the differential model
##   categorical (absent in one group)  confidently expressed in one group and
##                                      confidently silent in the other; needs
##                                      no model
##   categorical (interval separated)   confidently expressed in one group, the
##                                      other indeterminate, and the two Poisson
##                                      intervals do not overlap by at least the
##                                      required fold difference
##   off in both                        confidently silent in both
##   cannot assess (...)                reported with the reason, never as a
##                                      negative result
##
## The third class is the one that distinguishes this from a two-state
## analysis. A gene at 0.3 CPM in wild type and 20 CPM in opaque cornea cannot
## reach "confidently silent" at any sample size this study has, because the
## upper Poisson limit at low counts stays above the OFF threshold. Under a
## rule that requires both groups to be confidently called, such a gene is
## discarded however large the gap. Testing whether the two intervals overlap
## asks the question the data can actually answer.
##
## The class is deliberately conservative in three ways. One group must be
## confidently EXPRESSED -- two silent groups can also have non-overlapping
## intervals and calling those categorical would be meaningless. The fold
## difference is the conservative bound (lower limit of the higher group over
## upper limit of the lower group), not the ratio of point estimates. And the
## threshold's effect is printed in every run, so the choice is never hidden.
##
## This class does NOT change the set passed to the differential model, and so
## changes no fold change, no P value and no adjusted P value anywhere in the
## analysis. It moves genes out of an undifferentiated "cannot assess" bucket
## into a named, evidenced category.
##
## STATISTICAL BASIS
##   [1] Exact (Garwood) confidence limits for a Poisson rate, computed from
##       the gamma quantile function: lower = qgamma(a/2, k), upper =
##       qgamma(1-a/2, k+1), for k pooled counts over a pooled library size.
##       Garwood F (1936) Biometrika 28:437-442; the gamma-quantile form as
##       used here, Ulm K (1990) Am J Epidemiol 131:373-375.
##   [2] Counts are pooled across replicates before the limits are computed, so
##       the criterion depends on total counts and total depth rather than on
##       how many libraries those were split across, and is invariant to group
##       size.
##   [3] Leave-one-out stability, in the jackknife sense: a conclusion that
##       depends on which single observation was retained is not reported as a
##       conclusion. Efron B (1979) Ann Stat 7:1-26.
##   [4] Non-overlap of two confidence intervals is a conservative criterion
##       for a difference -- more conservative than a direct test at the same
##       level -- which suits a class intended to support categorical claims
##       without a model. Schenker N, Gentleman JF (2001) Am Stat 55:182-186.
##
## Inputs : 00_config.R, 01_load.R
## Outputs: ExpressionCalls.csv, ExpressionSummary.csv, GateCost_ByCell.csv,
##          Assessability.csv, Assessability_Summary.csv, CategoricalGenes.csv,
##          Categorical_FoldThresholdSensitivity.csv
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(purrr)
})

if (!exists("count_mat")) stop("Source 00_config.R and 01_load.R first.")

## ---------------------------------------------------------------------------
## 1. The rule
## ---------------------------------------------------------------------------

## Exact Poisson limits on pooled counts, in CPM. [1] [2]
poisson_limits <- function(ids) {
  k <- rowSums(count_mat[, ids, drop = FALSE])
  L <- sum(lib_sizes[ids])
  list(k = k, L = L,
       cpm = k / L * 1e6,
       lcl = ifelse(k == 0, 0, stats::qgamma(POISSON_ALPHA, shape = k)) / L * 1e6,
       ucl = stats::qgamma(1 - POISSON_ALPHA, shape = k + 1) / L * 1e6)
}

status_from_limits <- function(lim) {
  out <- rep("indeterminate", length(lim$lcl))
  out[lim$lcl >= T_HI] <- "expressed"
  out[lim$ucl <  T_LO & lim$lcl < T_HI] <- "not_expressed"
  out
}

status_ungated <- function(ids) {
  stats::setNames(status_from_limits(poisson_limits(ids)), rownames(count_mat))
}

## Any status that changes when one library is removed becomes indeterminate. [3]
status_gated <- function(ids) {
  full <- status_ungated(ids)
  if (!STABILITY_GATE || length(ids) < 3) return(full)
  loo <- vapply(ids, function(drop) status_ungated(setdiff(ids, drop)),
                character(nrow(count_mat)))
  full[rowSums(loo != full) > 0] <- "indeterminate"
  full
}

## ---------------------------------------------------------------------------
## 2. Groups
## ---------------------------------------------------------------------------

call_pool <- meta %>%
  dplyr::filter(tissue == TISSUE, laboratory %in% LABORATORY,
                !sample_id %in% DETECTION_EXCLUDE)
if (!nrow(call_pool)) {
  stop("No libraries remain after filtering on tissue and laboratory. ",
       "Check those two metadata columns.")
}

## Genotype within age, transparent cornea only. Every cornea at P3/P4 and P15
## is transparent, so this axis is genotype at matched developmental stage.
groups_age <- call_pool %>%
  dplyr::filter(transparency == "T") %>%
  dplyr::mutate(cell = paste0(age_group, ":", genotype)) %>%
  dplyr::group_by(cell, age_group, group = genotype) %>%
  dplyr::summarise(ids = list(sample_id), n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(axis = "age")

## Wild type, transparent mutant and opaque mutant in the adult.
groups_opacity <- call_pool %>%
  dplyr::filter(age_group == "Adult", !is.na(group3)) %>%
  dplyr::mutate(cell = paste0("Adult:", group3)) %>%
  dplyr::group_by(cell, age_group, group = group3) %>%
  dplyr::summarise(ids = list(sample_id), n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(axis = "opacity")

## "Adult:Sey" on the age axis and "Adult:Sey_T" on the opacity axis are the
## same libraries under two names, as are the two "Adult:WT" entries. The
## status is therefore computed once per distinct set of libraries and reused.
groups_all <- dplyr::bind_rows(groups_age, groups_opacity)
groups <- groups_all %>% dplyr::filter(n >= MIN_GROUP_N)

dropped <- groups_all %>% dplyr::filter(n < MIN_GROUP_N)
if (nrow(dropped)) {
  message("\nGROUPS BELOW n = ", MIN_GROUP_N, ", not called. Every contrast ",
          "using them is dropped as well:")
  print(as.data.frame(dropped %>% dplyr::select(axis, cell, n)))
}

message("\nGroups:")
print(as.data.frame(groups %>% dplyr::select(axis, cell, n)))

cells_unique <- groups %>% dplyr::distinct(cell, .keep_all = TRUE)

## ---------------------------------------------------------------------------
## 3. Status and limits, per group
## ---------------------------------------------------------------------------

## One tibble per group, rows in panel-gene order, so downstream comparisons
## are element-wise and need no join.
cell_stats <- stats::setNames(lapply(seq_len(nrow(cells_unique)), function(i) {
  ids  <- cells_unique$ids[[i]]
  lim  <- poisson_limits(ids)
  keep <- match(panel_genes, rownames(count_mat))
  tibble::tibble(
    gene_id      = panel_genes,
    n_reps       = length(ids),
    k            = as.integer(lim$k[keep]),
    lib_total    = lim$L,
    pooled_cpm   = round(lim$cpm[keep], 4),
    lcl          = round(lim$lcl[keep], 4),
    ucl          = round(lim$ucl[keep], 4),
    call         = status_gated(ids)[panel_genes],
    call_ungated = status_ungated(ids)[panel_genes])
}), cells_unique$cell)

## Full transcriptome status, needed by the program-level test downstream.
cell_status_genomewide <- stats::setNames(
  lapply(seq_len(nrow(cells_unique)), function(i) status_gated(cells_unique$ids[[i]])),
  cells_unique$cell)

expression_calls <- groups %>%
  dplyr::select(axis, cell, age_group, group) %>%
  dplyr::mutate(stats = cell_stats[cell]) %>%
  tidyr::unnest(stats) %>%
  dplyr::inner_join(panel_tbl, by = "gene_id", relationship = "many-to-many")
readr::write_csv(expression_calls, out_path("ExpressionCalls.csv"))

expression_summary <- expression_calls %>%
  dplyr::count(panel, axis, cell, call) %>%
  tidyr::pivot_wider(names_from = call, values_from = n, values_fill = 0)
readr::write_csv(expression_summary, out_path("ExpressionSummary.csv"))
message("\n=== Expression status ===")
print(as.data.frame(expression_summary))

## ---------------------------------------------------------------------------
## 4. What the stability check costs
## ---------------------------------------------------------------------------

## Reported per group, because the check is repeated once per library and its
## severity depends on how much of the pooled depth each library carries. A
## group whose calls move more than its peers' is less resolved for reasons
## that are not biological, and any comparison involving it should say so.
gate_cost <- purrr::map_dfr(names(cell_stats), function(cl) {
  s <- cell_stats[[cl]]
  tibble::tibble(
    cell = cl, n_reps = s$n_reps[1], n_panel_genes = nrow(s),
    demoted_from_expressed     = sum(s$call_ungated == "expressed"     & s$call == "indeterminate"),
    demoted_from_not_expressed = sum(s$call_ungated == "not_expressed" & s$call == "indeterminate"),
    promoted                   = sum(s$call_ungated == "indeterminate" & s$call != "indeterminate"),
    total_moved                = sum(s$call_ungated != s$call),
    pct_moved                  = round(100 * sum(s$call_ungated != s$call) / nrow(s), 2))
}) %>% dplyr::arrange(dplyr::desc(pct_moved))
readr::write_csv(gate_cost, out_path("GateCost_ByCell.csv"))
message("\n=== Cost of the leave-one-out stability check, per group ===")
print(as.data.frame(gate_cost))
if (sum(gate_cost$promoted) > 0) {
  message("  WARNING: the check promoted ", sum(gate_cost$promoted), " call(s). ",
          "It is one-way by construction; investigate before proceeding.")
}

## ---------------------------------------------------------------------------
## 5. Assessability
## ---------------------------------------------------------------------------

contrasts_usable <- CONTRASTS[CONTRASTS$cell_a %in% groups$cell &
                              CONTRASTS$cell_b %in% groups$cell, , drop = FALSE]
if (nrow(contrasts_usable) < nrow(CONTRASTS)) {
  message("\nCONTRASTS DROPPED (a group is missing or too small): ",
          paste(setdiff(CONTRASTS$label, contrasts_usable$label), collapse = ", "))
}

## Classify one contrast. `a` is the test group, `b` the reference.
classify_contrast <- function(label, cell_a, cell_b, min_fold) {
  A <- cell_stats[[cell_a]]; B <- cell_stats[[cell_b]]

  ## The conservative fold difference: the smallest ratio compatible with both
  ## intervals. Greater than 1 exactly when the intervals do not overlap, which
  ## makes non-overlap and the fold threshold the same criterion. [4]
  fold_a_over_b <- A$lcl / B$ucl
  fold_b_over_a <- B$lcl / A$ucl
  bounded_fold  <- pmax(fold_a_over_b, fold_b_over_a)

  one_expressed <- xor(A$call == "expressed", B$call == "expressed")
  other_unresolved <- A$call == "indeterminate" | B$call == "indeterminate"

  tibble::tibble(
    contrast = label, gene_id = A$gene_id,
    cell_a = cell_a, cell_b = cell_b,
    call_a = A$call, call_b = B$call,
    cpm_a = A$pooled_cpm, cpm_b = B$pooled_cpm,
    lcl_a = A$lcl, ucl_a = A$ucl, lcl_b = B$lcl, ucl_b = B$ucl,
    intervals_separated = bounded_fold > 1,
    bounded_fold = round(bounded_fold, 3),
    point_fold = round(pmax(A$pooled_cpm, B$pooled_cpm) /
                       pmax(pmin(A$pooled_cpm, B$pooled_cpm), 1e-6), 2),
    higher_in = dplyr::case_when(A$pooled_cpm > B$pooled_cpm ~ cell_a,
                                 B$pooled_cpm > A$pooled_cpm ~ cell_b,
                                 TRUE ~ NA_character_),
    assessable = dplyr::case_when(
      A$call == "expressed" & B$call == "expressed"
        ~ "test quantitatively",
      (A$call == "expressed" & B$call == "not_expressed") |
      (A$call == "not_expressed" & B$call == "expressed")
        ~ "categorical (absent in one group)",
      A$call == "not_expressed" & B$call == "not_expressed"
        ~ "off in both",
      one_expressed & other_unresolved & bounded_fold >= min_fold
        ~ "categorical (interval separated)",
      one_expressed & other_unresolved
        ~ "cannot assess (one group expressed, intervals overlap)",
      TRUE
        ~ "cannot assess (neither group confidently expressed)"))
}

assessability <- purrr::pmap_dfr(
  list(contrasts_usable$label, contrasts_usable$cell_a, contrasts_usable$cell_b),
  classify_contrast, min_fold = CATEGORICAL_MIN_FOLD) %>%
  dplyr::inner_join(panel_tbl, by = "gene_id", relationship = "many-to-many") %>%
  dplyr::left_join(dplyr::select(contrasts_usable, contrast = label, axis, model),
                   by = "contrast")
readr::write_csv(assessability, out_path("Assessability.csv"))

assessability_summary <- assessability %>%
  dplyr::count(panel, contrast, assessable) %>%
  tidyr::pivot_wider(names_from = assessable, values_from = n, values_fill = 0)
readr::write_csv(assessability_summary, out_path("Assessability_Summary.csv"))
message("\n=== Assessability (categorical minimum fold ", CATEGORICAL_MIN_FOLD, ") ===")
print(as.data.frame(assessability_summary))

## ---------------------------------------------------------------------------
## 6. The categorical genes
## ---------------------------------------------------------------------------

categorical <- assessability %>%
  dplyr::filter(startsWith(assessable, "categorical")) %>%
  dplyr::mutate(direction = ifelse(higher_in == cell_a, "up in test", "up in reference")) %>%
  dplyr::select(contrast, axis, panel, symbol, gene_id, assessable, direction,
                higher_in, call_a, call_b, cpm_a, cpm_b,
                lcl_a, ucl_a, lcl_b, ucl_b, bounded_fold, point_fold) %>%
  dplyr::arrange(panel, contrast, dplyr::desc(bounded_fold))
readr::write_csv(categorical, out_path("CategoricalGenes.csv"))

## The table is one row per gene per PANEL, as everywhere else here. The totals
## are not: summing panel rows across panels counts each of the thirty shared
## genes once per panel it sits on. Distinct gene x contrast calls are counted
## separately, and both numbers are printed so a number quoted in the paper can
## be traced to the right one.
categorical_distinct <- categorical %>%
  dplyr::distinct(contrast, symbol, direction)

## A THIRD denominator, and the one the manuscript quotes. CONTRASTS carries the
## hinge comparison under two labels -- HINGE["A"] from the developmental model
## and HINGE["B"] from the opacity model -- but both compare the SAME two cells.
## A categorical call is model-free, so those two labels hold identical rows and
## counting them both counts one comparison twice. Everything below is over
## DISTINCT COMPARISONS, with the Model A label dropped.
hinge_a <- categorical_distinct %>% dplyr::filter(contrast == HINGE[["A"]])
hinge_b <- categorical_distinct %>% dplyr::filter(contrast == HINGE[["B"]])
if (!setequal(paste(hinge_a$symbol, hinge_a$direction),
              paste(hinge_b$symbol, hinge_b$direction))) {
  stop("The two hinge labels (", HINGE[["A"]], ", ", HINGE[["B"]], ") gave ",
       "different categorical calls. They compare the same cells, so this ",
       "cannot happen unless CONTRASTS or the classification logic is wrong.")
}
categorical_comparisons <- categorical_distinct %>%
  dplyr::filter(contrast != HINGE[["A"]])

message("\n=== Categorical differences: ", nrow(categorical_distinct),
        " gene x contrast call(s), from ", nrow(categorical),
        " panel-membership row(s) ===")
print(as.data.frame(categorical %>% dplyr::count(contrast, panel, assessable)))
message("  Rows above count panel membership and must not be summed across panels.")
message("\n  Direction, over distinct gene x contrast calls: ",
        sum(categorical_distinct$direction == "up in test"), " higher in the test group, ",
        sum(categorical_distinct$direction == "up in reference"), " higher in the reference group.")
message("  That total counts the hinge comparison twice, once per model label.")
message("\n  QUOTE THIS ONE. Over the ", dplyr::n_distinct(categorical_comparisons$contrast),
        " distinct comparisons (", HINGE[["A"]], " dropped as a duplicate of ",
        HINGE[["B"]], "): ", nrow(categorical_comparisons), " call(s), ",
        sum(categorical_comparisons$direction == "up in test"), " higher in the test group, ",
        sum(categorical_comparisons$direction == "up in reference"),
        " higher in the reference group.")

## ---------------------------------------------------------------------------
## 7. Sensitivity of the categorical class to its fold threshold
## ---------------------------------------------------------------------------

## Printed every run. The count at fold = 1 is bare interval non-overlap; the
## difference between that and the chosen threshold is what the threshold buys.
fold_sensitivity <- purrr::map_dfr(CATEGORICAL_FOLD_GRID, function(f) {
  a <- purrr::pmap_dfr(list(contrasts_usable$label, contrasts_usable$cell_a,
                            contrasts_usable$cell_b), classify_contrast, min_fold = f) %>%
    dplyr::inner_join(panel_tbl, by = "gene_id", relationship = "many-to-many")
  tibble::tibble(
    min_fold = f,
    interval_separated = sum(a$assessable == "categorical (interval separated)"),
    absent_in_one_group = sum(a$assessable == "categorical (absent in one group)"),
    categorical_total = sum(startsWith(a$assessable, "categorical")),
    ## Three denominators sit in this table and they are not interchangeable.
    ## The three counts above are PANEL-MEMBERSHIP ROWS: a gene on two panels
    ## is counted under each. `distinct_calls` is one per gene per contrast,
    ## which is the number quoted in the paper, and `unique_genes` is one per
    ## gene however many contrasts it appears in.
    distinct_calls = dplyr::n_distinct(
      paste(a$symbol, a$contrast)[startsWith(a$assessable, "categorical")]),
    ## As distinct_calls, but with the duplicated hinge label removed, so this
    ## is the count over distinct COMPARISONS. This is the figure to quote.
    distinct_comparison_calls = dplyr::n_distinct(
      paste(a$symbol, a$contrast)[startsWith(a$assessable, "categorical") &
                                  a$contrast != HINGE[["A"]]]),
    unique_genes = dplyr::n_distinct(a$symbol[startsWith(a$assessable, "categorical")]),
    test_quantitatively = sum(a$assessable == "test quantitatively"))
})
readr::write_csv(fold_sensitivity, out_path("Categorical_FoldThresholdSensitivity.csv"))
message("\n=== Categorical class against its fold threshold ===")
print(as.data.frame(fold_sensitivity))
message("  interval_separated, absent_in_one_group and categorical_total count")
message("  PANEL-MEMBERSHIP ROWS; distinct_calls counts gene x contrast; and")
message("  unique_genes counts genes. distinct_calls still carries the hinge")
message("  comparison under both of its model labels, so QUOTE")
message("  distinct_comparison_calls in the paper, not distinct_calls.")
message("  `test quantitatively` is identical at every threshold. The categorical")
message("  class draws only from genes that could never enter the differential")
message("  model, so no fold change or adjusted P value in this analysis depends")
message("  on the value chosen.")

if (length(unique(fold_sensitivity$test_quantitatively)) != 1) {
  stop("The quantitative test set changed with the categorical fold threshold. ",
       "That must never happen; the classification logic is wrong.")
}

message("\n", strrep("=", 74))
message("03_expression_status complete. Outputs in: ", OUT_ROOT)
message("  ExpressionCalls.csv    status, pooled counts and Poisson limits per group")
message("  GateCost_ByCell.csv    what the stability check costs each group")
message("  Assessability.csv      what statement each gene supports, per contrast")
message("  CategoricalGenes.csv   genes on in one group and off in the other")
message(strrep("=", 74))
