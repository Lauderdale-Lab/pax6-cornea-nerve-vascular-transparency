###############################################################################
## 04_differential.R -- differential expression, equivalence, omnibus and interactions
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
## 04_differential.R
##
## Fits the differential-expression models and reports, for every gene the
## previous script judged testable, one of four outcomes.
##
##   changed                   the effect differs from zero
##   no change detected        the effect is significantly SMALLER than the
##     (bounded)               equivalence bound, so an effect that size or
##                             larger is excluded
##   undetermined              neither test is significant. No effect was
##                             detected and none was excluded; this experiment
##                             does not distinguish the two
##   not testable              the model returned no P value for this gene
##
## Separating the second and third outcomes is the point of the design. A
## non-significant result is two different findings, and collapsing them is what
## lets a false negative pass as evidence of no effect.
##
## TWO MODELS, ONE PER BIOLOGICAL REGIME
##   A  development and genotype -- all transparent libraries, every age
##   B  opacity -- adult libraries only, three groups
## Both use design ~ batch + male_frac + cell_group, where cell_group is a
## single factor whose levels are the observed combinations of age, genotype
## and transparency. Each cell has its own mean, so no additivity between age
## and genotype is assumed.
##
## The regimes are fitted separately because the developing cornea and the
## adult cornea undergoing opacification have different variance structures. A
## combined fit gives near-identical effect estimates but tests the adult
## contrasts against variability contributed by the neonatal groups. Within a
## model every contrast shares a dispersion estimate and counts of differing
## genes are comparable; across models they are not, and are not compared.
##
## Adult transparent mutant versus wild type is estimable under both. It is
## fitted in both, reported from Model B -- an adult contrast whose dispersion
## comes from adult libraries alone -- and Model A's version is retained as a
## consistency check, flagged reportable = FALSE.
##
## ALSO FITTED
##   An omnibus likelihood-ratio test across the three adult groups, testing
##   any difference among wild type, transparent mutant and opaque mutant
##   without nominating a pair.
##
##   Interaction contrasts on the age axis: whether the genotype effect itself
##   changes between two ages, as a difference of differences. A claim that a
##   gene is "altered from onset" rather than later is an interaction claim,
##   and comparing two separately significant contrasts does not test it.
##   These are reported separately and are the least powered thing here.
##
## FOLD CHANGES ARE MAXIMUM-LIKELIHOOD ESTIMATES. No shrinkage is applied.
## Shrinking estimates toward zero would narrow the intervals that the
## equivalence test depends on and make claims of no change appear better
## supported than the data warrant. Rankings and figures use the same estimates
## as the tests, so no two numbers in this analysis come from different
## estimators.
##
## MULTIPLICITY. Adjusted P values are computed within panel and contrast, over
## the pre-specified testable set. A gene listed on more than one panel is
## corrected once, in the primary panel named in 00_config.R, and that value is
## shown wherever the gene appears. The panel each gene was corrected in is
## written out as `correction_panel`.
##
## STATISTICAL BASIS
##   [1] Negative-binomial generalised linear models with empirical-Bayes
##       dispersion shrinkage, and Wald tests on the fitted coefficients.
##       Love MI, Huber W, Anders S (2014) Genome Biology 15:550.
##   [2] Equivalence testing against the composite null |log2 fold change| >=
##       bound, via lfcThreshold with altHypothesis = "lessAbs". The P value is
##       the larger of two one-sided tests. Schuirmann DJ (1987) J Pharmacokinet
##       Biopharm 15:657-680; implemented in DESeq2 as above. Requires
##       betaPrior = FALSE, which is asserted at fit time.
##   [3] Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
##   [4] Independent filtering is deliberately not used for the reported
##       adjusted P values: the testable set is already an expression-based
##       prefilter chosen before any test was run, so filtering again on mean
##       count would select on the data twice. Bourgon R, Gentleman R, Huber W
##       (2010) PNAS 107:9546-9551.
##   [5] Genes with an influential observation by Cook's distance are reported
##       as not testable rather than tested. Group sizes here are below the
##       threshold at which DESeq2 replaces outlying counts, so the gene is
##       set aside instead. Cook RD (1977) Technometrics 19:15-18.
##   [6] Likelihood-ratio test of the full model against one lacking the group
##       term, for the adult omnibus comparison. Love et al. (2014), above.
##
## Inputs : 00_config.R, 01_load.R, 03_expression_status.R (or its Assessability.csv)
## Outputs: DE_Results.csv, DE_Summary.csv, AdultOmnibusLRT.csv,
##          InteractionResults.csv, HingeCheck.csv, TestabilityAudit.csv,
##          CooksFiltering.csv, SharedGenePrimaryPanel.csv
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(purrr)
  library(DESeq2); library(SummarizedExperiment)
})

if (!exists("count_mat")) stop("Source 00_config.R and 01_load.R first.")

if (!exists("assessability")) {
  f <- out_path("Assessability.csv")
  if (!file.exists(f)) stop("Assessability.csv not found. Run 03_expression_status.R first.")
  assessability <- readr::read_csv(f, show_col_types = FALSE)
}

## ---------------------------------------------------------------------------
## 1. Which panel each gene is corrected in
## ---------------------------------------------------------------------------

panel_membership <- panel_tbl %>% dplyr::count(symbol, name = "n_panels")
shared_symbols <- panel_membership$symbol[panel_membership$n_panels > 1]

## PANEL_PRIMARY may arrive as NULL if it was written c() rather than
## character(0). NULL[x] returns NULL, and every is.na() test on it is then
## length zero, which quietly defeats the checks below. Normalise first.
primary_map <- if (length(PANEL_PRIMARY)) {
  stats::setNames(as.character(PANEL_PRIMARY), names(PANEL_PRIMARY))
} else stats::setNames(character(0), character(0))

assigned <- primary_map[shared_symbols]
bad <- shared_symbols[!is.na(assigned) & !(assigned %in% names(PANEL_DIRS))]
if (length(bad)) {
  stop("PANEL_PRIMARY names a panel that does not exist, for: ",
       paste(bad, collapse = ", "), ". Valid panels: ",
       paste(names(PANEL_DIRS), collapse = ", "))
}
unassigned <- shared_symbols[is.na(assigned)]

if (length(unassigned) && identical(PANEL_PRIMARY_FALLBACK, "stop")) {
  stop(length(unassigned), " of ", length(shared_symbols), " gene(s) on more than ",
       "one panel have no entry in PANEL_PRIMARY, so the multiplicity frame is ",
       "not fully specified:\n  ", paste(unassigned, collapse = ", "),
       "\n\nRun 02_shared_gene_assignment.R, which writes ",
       "PanelPrimaryAssignment.R beside the scripts and applies it to the ",
       "session. To proceed on the precedence rule alone instead, set ",
       "PAX6_PANEL_PRIMARY_FALLBACK=rule.")
}

## Fallback: the precedence rule, identical to what 02b would have proposed.
## Recorded per gene rather than absorbed, because a gene assigned this way was
## never actually reviewed.
primary_of <- function(sym) {
  explicit <- unname(primary_map[sym])
  if (length(explicit) == 1 && !is.na(explicit)) return(explicit)
  in_panels <- unique(panel_tbl$panel[panel_tbl$symbol == sym])
  PANEL_PRECEDENCE[which(PANEL_PRECEDENCE %in% in_panels)[1]]
}

panel_primary <- panel_tbl %>%
  dplyr::distinct(symbol) %>%
  dplyr::mutate(
    n_panels = panel_membership$n_panels[match(symbol, panel_membership$symbol)],
    primary_panel = vapply(symbol, primary_of, character(1)),
    assignment = dplyr::case_when(
      n_panels == 1 ~ "only panel",
      !is.na(primary_map[symbol]) ~ "assigned in config",
      TRUE ~ "fallback: precedence rule"))

## The pasted block and the override table are two records of one decision. If
## they disagree, the block is stale -- someone changed an override and did not
## regenerate. Fail rather than silently use whichever was sourced last.
if (nrow(PANEL_PRIMARY_OVERRIDE)) {
  have <- unname(primary_map[PANEL_PRIMARY_OVERRIDE$symbol])
  clash <- which(is.na(have) | have != PANEL_PRIMARY_OVERRIDE$panel)
  if (length(clash) && !length(unassigned)) {
    stop("PANEL_PRIMARY disagrees with PANEL_PRIMARY_OVERRIDE for: ",
         paste(sprintf("%s (block says %s, override says %s)",
                       PANEL_PRIMARY_OVERRIDE$symbol[clash],
                       ifelse(is.na(have[clash]), "nothing", have[clash]),
                       PANEL_PRIMARY_OVERRIDE$panel[clash]), collapse = "; "),
         ". Re-run 02_shared_gene_assignment.R and paste the new block.")
  }
}
readr::write_csv(panel_primary %>% dplyr::filter(n_panels > 1),
                 out_path("SharedGenePrimaryPanel.csv"))

n_fallback <- sum(panel_primary$assignment == "fallback: precedence rule")
message("Shared genes: ", length(shared_symbols), " on more than one panel; ",
        sum(panel_primary$assignment == "assigned in config"), " assigned in config, ",
        n_fallback, " by fallback.")
if (n_fallback > 0) {
  message("  WARNING: ", n_fallback, " shared gene(s) fell back to the precedence ",
          "rule rather than being reviewed and recorded in PANEL_PRIMARY. ",
          "See SharedGenePrimaryPanel.csv.")
}

## ---------------------------------------------------------------------------
## 2. Model fitting
## ---------------------------------------------------------------------------

de_pool <- meta %>%
  dplyr::filter(tissue == TISSUE, laboratory %in% LABORATORY,
                !sample_id %in% DE_EXCLUDE)
if (!nrow(de_pool)) stop("No libraries remain for the differential models.")

build_cell_group <- function(df) {
  dplyr::case_when(
    df$age_group != "Adult" ~ paste0(df$age_group, "_", df$genotype),
    df$age_group == "Adult" & df$genotype == "WT" ~ "Adult_WT",
    df$age_group == "Adult" & df$transparency == "T" ~ "Adult_Sey_T",
    df$age_group == "Adult" & df$transparency == "O" ~ "Adult_Sey_O",
    TRUE ~ NA_character_)
}

## Diagnostics accumulate in an environment so they survive regardless of where
## this script is sourced.
diag_env <- new.env(parent = emptyenv())
diag_env$cooks <- list()

fit_regime <- function(regime) {
  samples <- if (identical(regime, "A")) {
    de_pool %>% dplyr::filter(transparency == "T", !is.na(age_group))
  } else {
    de_pool %>% dplyr::filter(age_group == "Adult", transparency %in% c("T", "O"))
  }
  samples <- samples %>%
    dplyr::mutate(cell_group = build_cell_group(.)) %>%
    dplyr::filter(!is.na(cell_group))

  col_data <- as.data.frame(samples %>% tibble::column_to_rownames("sample_id"))
  reference <- if (identical(regime, "A")) "P3_P4_WT" else "Adult_WT"
  col_data$cell_group <- stats::relevel(factor(col_data$cell_group), ref = reference)
  col_data$batch <- droplevels(factor(col_data$batch))

  terms <- c(if (nlevels(col_data$batch) > 1) "batch",
             if (USE_SEX_COVARIATE) "male_frac", "cell_group")
  formula_full <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- try(stats::model.matrix(formula_full, col_data), silent = TRUE)
  if (inherits(mm, "try-error") || qr(mm)$rank < ncol(mm)) {
    stop("Model ", regime, " design is not full rank: ", deparse(formula_full))
  }

  ## Genes are pre-filtered for the fit only. Expression status is decided in
  ## 03 on unfiltered counts, so this filter never removes a gene from a
  ## denominator -- only from a model that could not estimate it.
  cm <- count_mat[, rownames(col_data), drop = FALSE]
  cm <- cm[rowSums(cm >= 10) >= 2, , drop = FALSE]

  message("\n=== Model ", regime, ": ",
          if (identical(regime, "A")) "development and genotype (transparent, all ages)"
          else "opacity (adult only)", " ===")
  message("  ", ncol(cm), " libraries, ", nrow(cm), " genes after the count pre-filter")
  message("  design ", deparse(formula_full))
  print(table(col_data$cell_group))

  dds <- DESeq2::DESeq(
    DESeq2::DESeqDataSetFromMatrix(cm, col_data, design = formula_full), quiet = TRUE)
  beta_prior <- attr(dds, "betaPrior")
  if (!is.null(beta_prior) && isTRUE(beta_prior)) {
    stop("betaPrior = TRUE; the equivalence test [2] requires betaPrior = FALSE.")
  }

  cooks <- tryCatch({
    m <- as.matrix(SummarizedExperiment::assays(dds)[["cooks"]]); m[is.na(m)] <- 0; m
  }, error = function(e) NULL)

  list(dds = dds, col_data = col_data, cooks = cooks, regime = regime,
       levels = levels(col_data$cell_group), reference = reference)
}

## Records which library carried the largest Cook's distance for each gene the
## model declined to test. [5]
record_cooks <- function(fit, label, res) {
  na_rows <- which(is.na(res$pvalue))
  if (length(na_rows) && !is.null(fit$cooks)) {
    worst <- colnames(fit$cooks)[max.col(fit$cooks[na_rows, , drop = FALSE],
                                         ties.method = "first")]
    tb <- as.data.frame(table(worst), stringsAsFactors = FALSE)
    names(tb) <- c("sample_id", "n_genes_flagged")
  } else {
    tb <- data.frame(sample_id = NA_character_, n_genes_flagged = 0L,
                     stringsAsFactors = FALSE)
  }
  diag_env$cooks[[length(diag_env$cooks) + 1L]] <- tibble::tibble(
    model = fit$regime, contrast = label,
    n_genes_in_fit = nrow(res), n_pvalue_na = length(na_rows),
    sample_id = tb$sample_id, n_genes_flagged = as.integer(tb$n_genes_flagged))
}

## Wald test plus the matching equivalence test, for one contrast. [1] [2]
test_contrast <- function(fit, label, axis, contrast_arg) {
  res <- DESeq2::results(fit$dds, contrast = contrast_arg, alpha = FDR)
  eq  <- DESeq2::results(fit$dds, contrast = contrast_arg, alpha = FDR,
                         lfcThreshold = EQUIV_BOUND, altHypothesis = "lessAbs")
  stopifnot(identical(rownames(res), rownames(eq)))
  record_cooks(fit, label, res)
  message("  [", label, "]  ", sum(is.na(res$pvalue)),
          " gene(s) set aside by outlier filtering")
  tibble::tibble(
    model = fit$regime, contrast = label, axis = axis, gene_id = rownames(res),
    log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
    pvalue = res$pvalue, pvalue_equiv = eq$pvalue)
}

fit_A <- fit_regime("A")
fit_B <- fit_regime("B")
fits <- list(A = fit_A, B = fit_B)

pair_results <- purrr::pmap_dfr(
  list(CONTRASTS$model, CONTRASTS$label, CONTRASTS$axis,
       CONTRASTS$level_a, CONTRASTS$level_b),
  function(model, label, axis, level_a, level_b) {
    fit <- fits[[model]]
    if (!all(c(level_a, level_b) %in% fit$levels)) {
      message("  [", label, "] a group level is absent -- skipped"); return(NULL)
    }
    test_contrast(fit, label, axis, c("cell_group", level_a, level_b))
  })

## ---------------------------------------------------------------------------
## 3. Reported outcomes
## ---------------------------------------------------------------------------

testable <- assessability %>%
  dplyr::filter(assessable == "test quantitatively") %>%
  dplyr::select(panel, contrast, gene_id, symbol)

de_joined <- testable %>%
  dplyr::inner_join(pair_results, by = c("contrast", "gene_id"),
                    relationship = "many-to-many") %>%
  dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel),
                   by = "symbol")

## Correct once per gene, in its primary panel's family, then show that value
## wherever the gene appears. [3] [4]
corrections <- de_joined %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::group_by(correction_panel = panel, contrast) %>%
  dplyr::mutate(padj = stats::p.adjust(pvalue, method = "BH"),
                padj_equiv = stats::p.adjust(pvalue_equiv, method = "BH"),
                n_in_correction_family = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(contrast, gene_id, correction_panel, padj, padj_equiv,
                n_in_correction_family)

de <- de_joined %>%
  dplyr::left_join(corrections, by = c("contrast", "gene_id")) %>%
  dplyr::mutate(
    ci_low  = log2FoldChange - stats::qnorm(1 - FDR / 2) * lfcSE,
    ci_high = log2FoldChange + stats::qnorm(1 - FDR / 2) * lfcSE,
    largest_effect_not_excluded = round(pmax(abs(ci_low), abs(ci_high)), 3),
    outcome = dplyr::case_when(
      is.na(padj) & is.na(padj_equiv)            ~ "not testable",
      !is.na(padj) & padj < FDR                  ~ "changed",
      !is.na(padj_equiv) & padj_equiv < FDR      ~ "no change detected (bounded)",
      TRUE                                       ~ "undetermined"),
    ## Model A's version of the hinge contrast duplicates a Model B result on
    ## the same libraries. Kept for the check below; never quoted.
    reportable = contrast != HINGE[["A"]]) %>%
  dplyr::arrange(panel, contrast, padj)
readr::write_csv(de, out_path("DE_Results.csv"))

de_summary <- de %>%
  dplyr::count(panel, contrast, reportable, outcome) %>%
  tidyr::pivot_wider(names_from = outcome, values_from = n, values_fill = 0)
readr::write_csv(de_summary, out_path("DE_Summary.csv"))
message("\n=== Outcomes (FDR ", FDR, ", equivalence bound +/-", EQUIV_BOUND, " log2) ===")
print(as.data.frame(de_summary))
message("  `undetermined` is not evidence of no effect. The largest effect still")
message("  compatible with the data is given per gene in DE_Results.csv.")

## ---------------------------------------------------------------------------
## 4. Adult omnibus test [6]
## ---------------------------------------------------------------------------

## Any difference among wild type, transparent mutant and opaque mutant,
## without nominating a pair. Fitted on the adult regime, dropping only the
## group term so that batch and pool composition stay in the reduced model.
reduced_terms <- c(if (nlevels(fit_B$col_data$batch) > 1) "batch",
                   if (USE_SEX_COVARIATE) "male_frac")
formula_reduced <- if (length(reduced_terms))
  stats::as.formula(paste("~", paste(reduced_terms, collapse = " + "))) else ~ 1

dds_lrt <- DESeq2::DESeq(fit_B$dds, test = "LRT", reduced = formula_reduced, quiet = TRUE)
res_lrt <- DESeq2::results(dds_lrt, alpha = FDR)

omnibus <- assessability %>%
  dplyr::filter(contrast == "SeyO_vs_WT_Adult", assessable == "test quantitatively") %>%
  dplyr::select(panel, gene_id, symbol) %>%
  dplyr::inner_join(
    tibble::tibble(gene_id = rownames(res_lrt), lrt_pvalue = res_lrt$pvalue),
    by = "gene_id") %>%
  dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel), by = "symbol")

omnibus_padj <- omnibus %>%
  dplyr::filter(panel == primary_panel) %>%
  dplyr::group_by(panel, .drop = TRUE) %>%
  dplyr::mutate(lrt_padj = stats::p.adjust(lrt_pvalue, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_id, lrt_padj)

omnibus <- omnibus %>%
  dplyr::left_join(omnibus_padj, by = "gene_id") %>%
  dplyr::mutate(differs_among_groups = !is.na(lrt_padj) & lrt_padj < FDR) %>%
  dplyr::arrange(panel, lrt_padj)
readr::write_csv(omnibus, out_path("AdultOmnibusLRT.csv"))
message("\n=== Adult omnibus LRT: any difference among the three groups ===")
print(as.data.frame(omnibus %>% dplyr::count(panel, differs_among_groups)))
message("  Reduced model: ", deparse(formula_reduced),
        ". The omnibus set is the adult testable set; genes assessable only in ")
message("  some pairwise contrasts are not included.")

## ---------------------------------------------------------------------------
## 5. Interaction contrasts
## ---------------------------------------------------------------------------

## Whether the genotype effect itself differs between two ages, as a difference
## of differences. Estimable directly from the cell-group parameterisation:
## comparing two separately significant pairwise contrasts is NOT this test.
##
## Only the developmental regime has the factorial structure this needs. The
## adult regime has three groups and no second factor, so its "interaction" is
## already the SeyO-vs-SeyT contrast reported above.
INTERACTIONS <- tibble::tribble(
  ~label,                        ~a1,           ~b1,        ~a2,           ~b2,
  "GenotypeEffect_P15_vs_P3_P4", "P15_Sey",     "P15_WT",   "P3_P4_Sey",   "P3_P4_WT",
  "GenotypeEffect_Adult_vs_P15", "Adult_Sey_T", "Adult_WT", "P15_Sey",     "P15_WT",
  "GenotypeEffect_Adult_vs_P3_P4","Adult_Sey_T","Adult_WT", "P3_P4_Sey",   "P3_P4_WT")

## An indicator over resultsNames for one cell-group level. The reference level
## has no coefficient and is represented by a vector of zeros.
level_vector <- function(dds, level, reference) {
  rn <- DESeq2::resultsNames(dds)
  v <- numeric(length(rn))
  if (identical(level, reference)) return(v)
  nm <- paste0("cell_group_", level, "_vs_", reference)
  idx <- match(nm, rn)
  if (is.na(idx)) stop("No coefficient named ", nm, " in the fitted model.")
  v[idx] <- 1
  v
}

interaction_results <- purrr::pmap_dfr(
  list(INTERACTIONS$label, INTERACTIONS$a1, INTERACTIONS$b1,
       INTERACTIONS$a2, INTERACTIONS$b2),
  function(label, a1, b1, a2, b2) {
    if (!all(c(a1, b1, a2, b2) %in% fit_A$levels)) {
      message("  [", label, "] a group level is absent -- skipped"); return(NULL)
    }
    v <- level_vector(fit_A$dds, a1, fit_A$reference) -
         level_vector(fit_A$dds, b1, fit_A$reference) -
         level_vector(fit_A$dds, a2, fit_A$reference) +
         level_vector(fit_A$dds, b2, fit_A$reference)
    test_contrast(fit_A, label, "interaction", v)
  })

if (nrow(interaction_results)) {
  ## Restricted to genes testable in BOTH of the contrasts being compared: an
  ## interaction is only interpretable where each component effect is estimable.
  component_map <- list(
    GenotypeEffect_P15_vs_P3_P4    = c("Sey_vs_WT_P15", "Sey_vs_WT_P3_P4"),
    GenotypeEffect_Adult_vs_P15    = c("SeyT_vs_WT_Adult", "Sey_vs_WT_P15"),
    GenotypeEffect_Adult_vs_P3_P4  = c("SeyT_vs_WT_Adult", "Sey_vs_WT_P3_P4"))

  eligible <- purrr::imap_dfr(component_map, function(pair, label) {
    both <- assessability %>%
      dplyr::filter(contrast %in% pair, assessable == "test quantitatively") %>%
      dplyr::count(panel, gene_id, symbol) %>%
      dplyr::filter(n == 2) %>%
      dplyr::mutate(contrast = label) %>%
      dplyr::select(panel, contrast, gene_id, symbol)
    both
  })

  interactions <- eligible %>%
    dplyr::inner_join(interaction_results, by = c("contrast", "gene_id"),
                      relationship = "many-to-many") %>%
    dplyr::left_join(dplyr::select(panel_primary, symbol, primary_panel), by = "symbol")

  int_padj <- interactions %>%
    dplyr::filter(panel == primary_panel) %>%
    dplyr::group_by(panel, contrast) %>%
    dplyr::mutate(padj = stats::p.adjust(pvalue, method = "BH"),
                  padj_equiv = stats::p.adjust(pvalue_equiv, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::select(contrast, gene_id, padj, padj_equiv)

  interactions <- interactions %>%
    dplyr::left_join(int_padj, by = c("contrast", "gene_id")) %>%
    dplyr::mutate(
      ci_low  = log2FoldChange - stats::qnorm(1 - FDR / 2) * lfcSE,
      ci_high = log2FoldChange + stats::qnorm(1 - FDR / 2) * lfcSE,
      outcome = dplyr::case_when(
        is.na(padj) & is.na(padj_equiv)       ~ "not testable",
        !is.na(padj) & padj < FDR             ~ "genotype effect differs between ages",
        !is.na(padj_equiv) & padj_equiv < FDR ~ "no difference detected (bounded)",
        TRUE                                  ~ "undetermined")) %>%
    dplyr::arrange(panel, contrast, padj)
  readr::write_csv(interactions, out_path("InteractionResults.csv"))
  message("\n=== Interaction: does the genotype effect change with age? ===")
  print(as.data.frame(interactions %>% dplyr::count(panel, contrast, outcome)))
  message("  This is the test behind any claim that a gene is altered from onset")
  message("  rather than later. It is the least powered analysis here: a difference")
  message("  of differences carries the variance of both component contrasts, so")
  message("  `undetermined` should be expected to dominate and is not evidence")
  message("  that the genotype effect is constant across ages.")
}

## ---------------------------------------------------------------------------
## 6. Diagnostics
## ---------------------------------------------------------------------------

## Where genes are lost between being judged testable and receiving a P value.
testability <- testable %>%
  dplyr::count(panel, contrast, name = "n_testable") %>%
  dplyr::left_join(de %>% dplyr::count(panel, contrast, name = "n_in_output"),
                   by = c("panel", "contrast")) %>%
  dplyr::left_join(de %>% dplyr::filter(is.na(pvalue)) %>%
                     dplyr::count(panel, contrast, name = "n_no_pvalue"),
                   by = c("panel", "contrast")) %>%
  dplyr::mutate(across(c(n_in_output, n_no_pvalue), ~ dplyr::coalesce(., 0L)),
                n_lost_to_prefilter = n_testable - n_in_output,
                n_tested = n_in_output - n_no_pvalue)
readr::write_csv(testability, out_path("TestabilityAudit.csv"))
message("\n=== Testable -> fitted -> tested ===")
print(as.data.frame(testability))
if (sum(testability$n_lost_to_prefilter) > 0) {
  message("  ", sum(testability$n_lost_to_prefilter), " gene x contrast cell(s) were ",
          "judged testable and then dropped by the count pre-filter. The testable ",
          "counts quoted in the paper are then not the counts actually tested.")
}

cooks_table <- dplyr::bind_rows(diag_env$cooks)
readr::write_csv(cooks_table, out_path("CooksFiltering.csv"))
flagged <- cooks_table %>% dplyr::filter(!is.na(sample_id), n_genes_flagged > 0)
message("\n=== Outlier filtering by library ===")
if (nrow(flagged)) {
  print(as.data.frame(flagged %>% dplyr::group_by(sample_id) %>%
    dplyr::summarise(n_genes_flagged = sum(n_genes_flagged),
                     n_contrasts = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_genes_flagged))))
  message("  If one library dominates, that is a sample-quality question ",
          "separate from any exclusion decision already made.")
} else {
  message("  No gene was set aside by outlier filtering in any contrast.")
}

## The hinge: one comparison, two models, identical libraries. Effect estimates
## should agree closely; the count of significant genes may not, and that
## difference is the evidence that the two regimes differ in variance.
hinge <- de %>%
  dplyr::filter(contrast %in% HINGE) %>%
  dplyr::select(panel, contrast, gene_id, symbol, log2FoldChange, padj, outcome) %>%
  tidyr::pivot_wider(names_from = contrast,
                     values_from = c(log2FoldChange, padj, outcome))
lfc_a <- paste0("log2FoldChange_", HINGE[["A"]])
lfc_b <- paste0("log2FoldChange_", HINGE[["B"]])
out_a <- paste0("outcome_", HINGE[["A"]]); out_b <- paste0("outcome_", HINGE[["B"]])
if (all(c(lfc_a, lfc_b) %in% names(hinge))) {
  readr::write_csv(hinge, out_path("HingeCheck.csv"))
  message("\n=== Hinge: adult transparent vs wild type, Model A against Model B ===")
  print(as.data.frame(hinge %>% dplyr::group_by(panel) %>% dplyr::summarise(
    n = dplyr::n(),
    changed_model_A = sum(.data[[out_a]] == "changed", na.rm = TRUE),
    changed_model_B = sum(.data[[out_b]] == "changed", na.rm = TRUE),
    agree = sum(.data[[out_a]] == .data[[out_b]], na.rm = TRUE),
    r_log2FC = round(stats::cor(.data[[lfc_a]], .data[[lfc_b]],
                                use = "complete.obs"), 4), .groups = "drop")))
  message("  Report the Model B value for this comparison everywhere, including ",
          "in the developmental narrative.")
}

message("\n", strrep("=", 74))
message("04_differential complete. Outputs in: ", OUT_ROOT)
message(strrep("=", 74))
