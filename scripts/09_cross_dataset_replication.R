## ---------------------------------------------------------------------------
## 09_cross_dataset_replication.R
##
## James D. Lauderdale, PhD; Department of Cellular Biology,
## University of Georgia, Athens, GA 30602, USA
## Study: "Nerve and vascular abnormalities precede loss of corneal
##         transparency in Pax6-haploinsufficient mice"
## Study authors: Sneha K. Mohan, James D. Lauderdale
##
## Repository : https://github.com/Lauderdale-Lab/pax6-cornea-nerve-vascular-transparency
## Archived   : https://doi.org/10.5281/zenodo.22965379
## Licence    : MIT (see LICENSE)
## Contact    : James D. Lauderdale, jdlauder@uga.edu
##
## PURPOSE
##   Ask whether the three-state expression calls of this study reproduce in an
##   independent adult cornea dataset: GSE183742 (Duncan lab), which carries the
##   Pax6^tm1Pgr allele rather than Pax6^Sey-Neu.
##
##   TIER 1, the strong test -- WILD TYPE against WILD TYPE.
##     Adult wild-type cornea in both datasets. Requires no opacity information,
##     so the fact that GSE183742 does not report corneal transparency costs
##     nothing here. It tests the claim this paper actually makes about what the
##     normal adult cornea does and does not express, in a different strain,
##     laboratory and library preparation.
##
##   TIER 2, the weaker test -- THE GENOTYPE EFFECT.
##     Mutant versus wild type in each dataset, compared for direction. Weaker
##     for three reasons, all reported alongside the result: a different allele,
##     n = 3 per group, and unknown opacity status. Because the GSE183742
##     animals were probably fibrotic but this was never recorded, their
##     contrast conflates genotype with opacity and is compared against
##     SeyO_vs_WT_Adult, not SeyT_vs_WT_Adult.
##
## WHAT MAKES THE COMPARISON MEAN ANYTHING
##   [A] The two datasets are called by the SAME RULE. This script recomputes
##       status for a Lauderdale cell and STOPS unless it reproduces
##       ExpressionCalls.csv exactly. Without that check, any concordance could
##       be an artefact of two slightly different implementations.
##   [B] Concordance is measured against a NULL. Two datasets of the same tissue
##       agree on most genes by construction, so a raw agreement percentage is
##       uninterpretable. The null here is abundance-matched non-panel genes,
##       the same device used by 05_program_level.R.
##   [C] Genes that either dataset cannot resolve are reported as such, never
##       counted as agreement or disagreement. With n = 3 per group, expect the
##       indeterminate class to be large. That is the framework working.
##   [D] The exceptions are printed by name, and a DISSECTION-MARGIN CHECK
##       asks whether they are explained by the cut rather than the biology:
##       off-panel markers of limbus, angle, pigment and lens are compared
##       between the two wild types, and between GSE183742's own groups. The
##       GEO record does not state where the corneas were cut, and a rim of
##       limbus brings in vessels, meshwork and nerve that the central cornea
##       lacks -- exactly the genes on which the datasets were found to differ.
##
## STATISTICAL BASIS
##   [1] Garwood (1936); Ulm (1990) -- exact Poisson confidence limits.
##   [2] Pooling counts within a group makes the criterion invariant to the
##       number of libraries, which is what allows a 3-library dataset to be
##       compared with a 5-library one at all.
##   [3] Efron (1979) -- the leave-one-out stability check.
##   [4] Cohen (1960) -- kappa, agreement corrected for chance. Reported
##       alongside the permutation null, not instead of it.
##   [5] Schenker & Gentleman (2001) Am Stat 55:182-186 -- non-overlap of
##       intervals as a conservative criterion for a difference.
##
## USAGE
##   Sys.setenv(PAX6_DUNCAN_COUNTS = "<path to Duncan matrix>")
##   Sys.setenv(PAX6_DUNCAN_META   = "<path to Duncan metadata csv>")
##   source("09_cross_dataset_replication.R")
##
##   Requires 00_config.R, 01_load.R, 01c_contamination_indices.R and a
##   completed run of the pipeline (ExpressionCalls.csv must exist) for the
##   Lauderdale side.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

if (!exists("OUT_ROOT")) stop("Source 00_config.R first.")
if (!exists("count_mat")) stop("Source 01_load.R first.")
## The margin check uses the marker identifiers and the per-library rule that
## 01c defines, so that GSE183742 is judged on exactly the terms of Supp. Table 6.
if (!exists("ot_marker_ids") || !exists("ot_flag_library"))
  stop("Source 01c_contamination_indices.R first.")

X_OUT <- file.path(OUT_ROOT, "CrossDataset")
dir.create(X_OUT, showWarnings = FALSE, recursive = TRUE)

## Two different failures used to give the same message -- unset, and set to a
## path that does not exist -- which sends you looking for the wrong problem.
## They are now reported separately, and the resolved path is echoed so a
## leading-slash or working-directory mistake is visible immediately.
x_require_file <- function(value, var, what) {
  if (!nzchar(value))
    stop(var, " is not set in this R session. Run:\n",
         "  Sys.setenv(", var, " = normalizePath(\"<path to ", what,
         ">\", mustWork = TRUE))\n",
         "  Sys.setenv does NOT survive an R restart, so this must be re-run\n",
         "  in every new session before sourcing this script.")
  if (!file.exists(value))
    stop(var, " is set, but no file exists there:\n",
         "  value            : ", value, "\n",
         "  working directory: ", getwd(), "\n",
         "  A relative path is resolved against the working directory; a path\n",
         "  beginning with / is absolute and is NOT relative to the project.")
  normalizePath(value, mustWork = TRUE)
}

X_DUNCAN_COUNTS <- x_require_file(Sys.getenv("PAX6_DUNCAN_COUNTS", unset = ""),
                                  "PAX6_DUNCAN_COUNTS",
                                  "GSE183742 count matrix")
X_DUNCAN_META   <- x_require_file(Sys.getenv("PAX6_DUNCAN_META", unset = ""),
                                  "PAX6_DUNCAN_META",
                                  "GSE183742 sample metadata CSV")
message("GSE183742 inputs resolved:\n  counts  : ", X_DUNCAN_COUNTS,
        "\n  metadata: ", X_DUNCAN_META)

X_N_PERM <- as.integer(Sys.getenv("PAX6_X_N_PERM", unset = "2000"))
set.seed(SEED)

## ---------------------------------------------------------------------------
## The rule, restated here so both datasets are called by one implementation
## ---------------------------------------------------------------------------
## Deliberately NOT a second, subtly different copy: the acceptance test below
## refuses to continue unless this reproduces 03_expression_status.R exactly on
## Lauderdale data.

x_poisson_lcl <- function(k, L) qgamma(POISSON_ALPHA, k) / L * 1e6          # [1]
x_poisson_ucl <- function(k, L) qgamma(1 - POISSON_ALPHA, k + 1) / L * 1e6  # [1]

x_call_one <- function(k, L) {
  lcl <- x_poisson_lcl(k, L); ucl <- x_poisson_ucl(k, L)
  ifelse(lcl >= T_HI, "expressed",
         ifelse(ucl < T_LO, "not_expressed", "indeterminate"))
}

## Pooled call over a set of libraries, with the leave-one-out stability gate.
## Any call that moves when a single library is dropped becomes indeterminate.
##
## `lib` is the FULL-MATRIX library size per sample, and it is a required
## argument on purpose. Deriving the denominator from colSums(cm) would make it
## depend on which ROWS were passed in, so a call restricted to the panel would
## express every gene as a fraction of the panel total rather than of the
## library -- and the panel is exactly the thing that changes between groups.
## 03_expression_status.R uses sum(lib_sizes[ids]); this must match it.
x_status <- function(cm, cols, gate = STABILITY_GATE, lib) {                # [2]
  if (missing(lib)) stop("x_status(): pass lib = full-matrix library sizes.")
  stopifnot(all(cols %in% names(lib)))
  sub <- cm[, cols, drop = FALSE]
  k   <- rowSums(sub)
  Lj  <- lib[cols]
  L   <- sum(Lj)
  out <- x_call_one(k, L)
  if (gate && length(cols) > 1) {                                          # [3]
    for (j in seq_along(cols)) {
      kk <- k - sub[, j]
      LL <- L - Lj[j]
      moved <- x_call_one(kk, LL) != out
      out[moved] <- "indeterminate"
    }
  }
  data.frame(gene_id = rownames(cm), status = out,
             k = k, L = L, cpm = k / L * 1e6, stringsAsFactors = FALSE)
}

## Conservative bounded fold change: > 1 exactly when the two Poisson intervals
## do not overlap, so non-overlap and magnitude are ONE criterion.            [5]
x_bounded_fold <- function(ka, La, kb, Lb) {
  a_l <- x_poisson_lcl(ka, La); a_u <- x_poisson_ucl(ka, La)
  b_l <- x_poisson_lcl(kb, Lb); b_u <- x_poisson_ucl(kb, Lb)
  pmax(a_l / b_u, b_l / a_u)
}

## ---------------------------------------------------------------------------
## [A] ACCEPTANCE TEST -- does this implementation match 03_expression_status.R?
## ---------------------------------------------------------------------------
x_calls_file <- file.path(OUT_ROOT, "ExpressionCalls.csv")
if (!file.exists(x_calls_file))
  stop("ExpressionCalls.csv not found. Run the pipeline through 03 first.")

x_ref <- as.data.frame(data.table::fread(x_calls_file))
x_cellcol <- intersect(c("cell", "cell_group", "group"), names(x_ref))[1]
if (is.na(x_cellcol)) stop("Could not find the cell column in ExpressionCalls.csv")

## 03_expression_status.R writes this column as `call`, not `status`. Resolve it
## by name rather than assuming, because selecting a column that is not there
## fails inside merge() with a message that names neither the column nor the file.
x_statuscol <- intersect(c("call", "status"), names(x_ref))[1]
if (is.na(x_statuscol))
  stop("No status column found in ExpressionCalls.csv (looked for `call`, ",
       "`status`). Columns present: ", paste(names(x_ref), collapse = ", "))

x_test_cell <- "Adult:WT"
x_ref_cell  <- x_ref[x_ref[[x_cellcol]] == x_test_cell, c("gene_id", x_statuscol)]
names(x_ref_cell) <- c("gene_id", "status")
if (!nrow(x_ref_cell))
  stop("Cell '", x_test_cell, "' not present in ExpressionCalls.csv")

## ExpressionCalls.csv holds one row per gene per PANEL, and Adult:WT appears on
## both the age and the opacity axis, so one gene can carry several identical
## rows. Reduce to one row per gene before comparing; otherwise the merge is
## many-to-one and any disagreement is counted once per panel membership.
x_ref_cell <- unique(x_ref_cell)
if (anyDuplicated(x_ref_cell$gene_id))
  stop("ExpressionCalls.csv gives more than one DISTINCT status for the same ",
       "gene in cell '", x_test_cell, "'. That should be impossible.")

## The laboratory filter matters: 03_expression_status.R applies it, so
## omitting it here would build a different cell from the same metadata and the
## acceptance test would fail for the wrong reason.
## 03_expression_status.R also drops DETECTION_EXCLUDE from its calling pool.
## Omitting it here would build a different cell from the same metadata, and the
## acceptance test would then fail for a reason that has nothing to do with the rule.
x_excl <- if (exists("DETECTION_EXCLUDE")) DETECTION_EXCLUDE else character(0)
x_test_ids <- meta$sample_id[meta$tissue == "Cornea" &
                             meta$laboratory %in% LABORATORY &
                             meta$age_group == "Adult" &
                             meta$genotype == "WT" &
                             !meta$sample_id %in% x_excl]
x_test_ids <- intersect(x_test_ids, colnames(count_mat))
x_mine <- x_status(count_mat, x_test_ids, lib = lib_sizes)

x_cmp <- merge(x_mine[, c("gene_id", "status")], x_ref_cell,
               by = "gene_id", suffixes = c("_new", "_ref"))
x_cmp <- x_cmp[x_cmp$gene_id %in% panel_genes, ]
x_disagree <- sum(x_cmp$status_new != x_cmp$status_ref)

cat("\n=== Acceptance test: rule reproduces 03_expression_status.R ===\n")
cat(sprintf("  cell %s, %d libraries, %d panel genes compared, %d disagreements\n",
            x_test_cell, length(x_test_ids), nrow(x_cmp), x_disagree))
if (x_disagree > 0) {
  print(utils::head(x_cmp[x_cmp$status_new != x_cmp$status_ref, ], 10))
  stop("The status rule in this script does not reproduce 03_expression_status.R. ",
       "Any concordance measured with it would be uninterpretable. Fix before continuing.")
}
cat("  PASS -- both datasets will be called by an identical rule.\n")

## ---------------------------------------------------------------------------
## Load GSE183742
## ---------------------------------------------------------------------------
x_dm <- as.data.frame(data.table::fread(X_DUNCAN_COUNTS, skip = "Geneid"))
x_gene <- sub("\\.[0-9]+$", "", x_dm$Geneid)
x_dm <- as.matrix(x_dm[, setdiff(names(x_dm),
          c("Geneid", "Chr", "Start", "End", "Strand", "Length")), drop = FALSE])
rownames(x_dm) <- x_gene
colnames(x_dm) <- sub("_hisat2\\.sorted\\.bam$", "", basename(colnames(x_dm)))
x_dm <- x_dm[!duplicated(rownames(x_dm)), , drop = FALSE]

## Library sizes over the WHOLE matrix, taken before any gene subsetting below.
x_dlib <- colSums(x_dm)

x_meta <- as.data.frame(data.table::fread(X_DUNCAN_META))
names(x_meta) <- tolower(names(x_meta))
## The BAMs are named by SRR run, so that is the join key.
x_key <- intersect(c("srr_run", "run", "sample_id", "sampleid"), names(x_meta))[1]
if (is.na(x_key)) stop("No SRR/sample column found in the Duncan metadata.")
x_meta$join_id <- as.character(x_meta[[x_key]])

x_missing <- setdiff(colnames(x_dm), x_meta$join_id)
if (length(x_missing))
  stop("Duncan count columns with no metadata row: ", paste(x_missing, collapse = ", "),
       "\n  The BAMs are named by SRR accession; make that the join key.")

x_wt  <- x_meta$join_id[grepl("^wt$|wild", x_meta$genotype, ignore.case = TRUE)]
x_mut <- setdiff(x_meta$join_id, x_wt)
x_wt  <- intersect(x_wt,  colnames(x_dm))
x_mut <- intersect(x_mut, colnames(x_dm))

cat(sprintf("\nGSE183742: %d genes, %d libraries (%d wild type, %d mutant)\n",
            nrow(x_dm), ncol(x_dm), length(x_wt), length(x_mut)))
if (length(x_wt) < 2 || length(x_mut) < 2)
  stop("Need at least two libraries per group.")

## Common gene universe. Both datasets are counted against GENCODE vM25, so
## identifiers match once version suffixes are stripped; verified, not assumed.
x_common <- intersect(rownames(count_mat), rownames(x_dm))
cat(sprintf("Genes shared with this study: %d of %d and %d\n",
            length(x_common), nrow(count_mat), nrow(x_dm)))
if (length(x_common) < 20000)
  stop("Too few shared gene identifiers -- the two matrices are probably not ",
       "counted against the same annotation. Check before interpreting anything.")

x_panel <- intersect(panel_genes, x_common)
cat(sprintf("Panel genes available in both: %d of %d\n",
            length(x_panel), length(panel_genes)))

## ---------------------------------------------------------------------------
## TIER 1 -- wild type against wild type
## ---------------------------------------------------------------------------
x_laud_wt_ids <- meta$sample_id[meta$tissue == "Cornea" &
                                meta$laboratory %in% LABORATORY &
                                meta$age_group == "Adult" &
                                meta$genotype == "WT" &
                                !meta$sample_id %in% x_excl]
x_laud_wt_ids <- intersect(x_laud_wt_ids, colnames(count_mat))

x_a <- x_status(count_mat[x_common, , drop = FALSE], x_laud_wt_ids, lib = lib_sizes)
x_b <- x_status(x_dm[x_common, , drop = FALSE],      x_wt,          lib = x_dlib)

x_wt_tab <- merge(x_a[, c("gene_id", "status", "k", "L", "cpm")],
                  x_b[, c("gene_id", "status", "k", "L", "cpm")],
                  by = "gene_id", suffixes = c("_lauderdale", "_duncan"))
x_wt_tab$panel <- x_wt_tab$gene_id %in% x_panel

## Concordance is computed ONLY where both datasets resolved the gene. Genes
## indeterminate in either are counted separately and reported, never folded
## into agreement or disagreement.
x_concordance <- function(df) {
  res <- df[df$status_lauderdale != "indeterminate" &
            df$status_duncan     != "indeterminate", ]
  n_res <- nrow(res)
  n_agree <- sum(res$status_lauderdale == res$status_duncan)
  ## Cohen's kappa on the resolved subset                                    [4]
  tb <- table(factor(res$status_lauderdale, levels = c("expressed", "not_expressed")),
              factor(res$status_duncan,     levels = c("expressed", "not_expressed")))
  n <- sum(tb)
  po <- if (n > 0) sum(diag(tb)) / n else NA_real_
  pe <- if (n > 0) sum(rowSums(tb) * colSums(tb)) / n^2 else NA_real_
  kappa <- if (!is.na(pe) && pe < 1) (po - pe) / (1 - pe) else NA_real_
  list(n_total = nrow(df), n_resolved = n_res,
       n_indeterminate = nrow(df) - n_res,
       n_agree = n_agree,
       pct_agree = if (n_res > 0) 100 * n_agree / n_res else NA_real_,
       kappa = kappa, tab = tb)
}

x_panel_conc <- x_concordance(x_wt_tab[x_wt_tab$panel, ])
x_bg_conc    <- x_concordance(x_wt_tab[!x_wt_tab$panel, ])

cat("\n=== TIER 1: adult wild-type cornea, this study vs GSE183742 ===\n")
cat(sprintf("  Panel genes      : %d resolved in both of %d, %.1f%% agree, kappa %.3f\n",
            x_panel_conc$n_resolved, x_panel_conc$n_total,
            x_panel_conc$pct_agree, x_panel_conc$kappa))
cat(sprintf("  Non-panel genes  : %d resolved in both of %d, %.1f%% agree, kappa %.3f\n",
            x_bg_conc$n_resolved, x_bg_conc$n_total,
            x_bg_conc$pct_agree, x_bg_conc$kappa))
cat(sprintf("  Not resolved in one or both: %d panel, %d non-panel\n",
            x_panel_conc$n_indeterminate, x_bg_conc$n_indeterminate))
cat("\n  Status cross-tabulation, panel genes (rows this study, cols GSE183742):\n")
print(with(x_wt_tab[x_wt_tab$panel, ],
           table(this_study = status_lauderdale, GSE183742 = status_duncan)))

## ---------------------------------------------------------------------------
## [B] THE NULL -- what agreement would abundance-matched genes give?
## ---------------------------------------------------------------------------
## Two adult corneal datasets agree on most genes by construction. The panel is
## only interesting if it agrees MORE than comparable genes chosen at random.
## Matching is on abundance because agreement is easiest for genes that are very
## highly or very lowly expressed, and the panel is not a random draw on that axis.

x_bg_pool <- x_wt_tab[!x_wt_tab$panel, ]
x_bg_pool <- x_bg_pool[x_bg_pool$cpm_lauderdale > 0 | x_bg_pool$cpm_duncan > 0, ]
x_panel_rows <- x_wt_tab[x_wt_tab$panel, ]

## Twenty quantile bins is a request, not a guarantee. A large share of the
## background sits at exactly zero CPM in this study, so log10(cpm + 0.01) is
## -2 for all of them and several quantiles land on the same value. cut() then
## fails with "'breaks' are not unique". Collapsing tied breaks is the right
## answer rather than a workaround: genes that are all at the detection floor
## ARE one abundance stratum, and matching cannot distinguish within it.
x_abund_bg    <- log10(x_bg_pool$cpm_lauderdale + 0.01)
x_abund_panel <- log10(x_panel_rows$cpm_lauderdale + 0.01)
x_breaks <- unique(stats::quantile(x_abund_bg, probs = seq(0, 1, length.out = 21),
                                   na.rm = TRUE))
if (length(x_breaks) < 3)
  stop("Fewer than two usable abundance bins: the background is almost entirely ",
       "at one abundance. Matching would be meaningless.")

x_bins       <- cut(x_abund_bg,    breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
x_panel_bins <- cut(x_abund_panel, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)

cat(sprintf("\n  Abundance matching on %d bin(s); 20 were requested and tied\n",
            length(x_breaks) - 1L))
cat("  breaks at the detection floor were collapsed.\n")
x_unmatched <- sum(is.na(x_panel_bins))
if (x_unmatched)
  cat(sprintf("  %d panel gene(s) fall outside the background abundance range and\n  are counted in the observed agreement but cannot be matched in the null.\n",
              x_unmatched))

x_null <- numeric(X_N_PERM)
for (i in seq_len(X_N_PERM)) {
  idx <- integer(0)
  for (b in unique(x_panel_bins[!is.na(x_panel_bins)])) {
    need <- sum(x_panel_bins == b, na.rm = TRUE)
    pool <- which(x_bins == b)
    if (!length(pool)) next
    ## pool[sample.int(...)], never sample(pool, ...): when a bin holds exactly
    ## one gene, sample() treats the scalar as 1:n and draws the wrong indices.
    idx <- c(idx, pool[sample.int(length(pool), need, replace = length(pool) < need)])
  }
  x_null[i] <- x_concordance(x_bg_pool[idx, ])$pct_agree
}
x_null <- x_null[is.finite(x_null)]
x_p <- (1 + sum(x_null >= x_panel_conc$pct_agree)) / (1 + length(x_null))

cat("\n=== Panel concordance against an abundance-matched null ===\n")
cat(sprintf("  observed %.1f%%; null median %.1f%% (2.5-97.5%%: %.1f-%.1f); p = %.4g\n",
            x_panel_conc$pct_agree, median(x_null),
            quantile(x_null, 0.025), quantile(x_null, 0.975), x_p))
cat("  A panel that merely matches the null is still a REPLICATION of the calls;\n")
cat("  it simply is not more reproducible than comparable genes. Report both.\n")

## ---------------------------------------------------------------------------
## TIER 2 -- direction of the genotype effect
## ---------------------------------------------------------------------------
## GSE183742 does not report corneal transparency. Those animals were probably
## fibrotic, so their contrast conflates genotype with opacity and is compared
## against SeyO_vs_WT_Adult. This is stated in the output, not buried here.

x_laud_op_ids <- meta$sample_id[meta$tissue == "Cornea" &
                                meta$laboratory %in% LABORATORY &
                                meta$age_group == "Adult" &
                                meta$genotype == "Sey" &
                                !meta$sample_id %in% x_excl &
                                toupper(substr(as.character(meta$transparency), 1, 1)) == "O"]
x_laud_op_ids <- intersect(x_laud_op_ids, colnames(count_mat))

x_dir <- function(cm, ref_ids, test_ids, genes, lib) {
  r <- x_status(cm[genes, , drop = FALSE], ref_ids,  lib = lib)
  t <- x_status(cm[genes, , drop = FALSE], test_ids, lib = lib)
  bf <- x_bounded_fold(t$k, t$L, r$k, r$L)
  up <- t$cpm > r$cpm
  data.frame(gene_id = r$gene_id,
             direction = ifelse(bf < CATEGORICAL_MIN_FOLD, "unresolved",
                                ifelse(up, "up_in_mutant", "down_in_mutant")),
             bounded_fold = bf, stringsAsFactors = FALSE)
}

## Agreement statistics for a direction table. Three numbers are returned
## because no one of them can be read alone:
##   pct_agree  the raw proportion, which is inflated whenever one direction
##              dominates the margins -- as it does here
##   kappa      agreement corrected for the marginal composition             [4]
##   tab        the 2x2 itself, which is what lets a reader judge the other two
x_dir_agree_stats <- function(df) {
  lv  <- c("up_in_mutant", "down_in_mutant")
  res <- df[df$direction_lauderdale %in% lv & df$direction_duncan %in% lv, ]
  tb  <- table(factor(res$direction_lauderdale, levels = lv),
               factor(res$direction_duncan,     levels = lv))
  n     <- sum(tb)
  agree <- sum(diag(tb))
  po <- if (n > 0) agree / n else NA_real_
  pe <- if (n > 0) sum(rowSums(tb) * colSums(tb)) / n^2 else NA_real_
  list(n_total = nrow(df), n_resolved = n, n_agree = agree,
       pct_agree = if (n > 0) 100 * po else NA_real_,
       pct_chance = if (n > 0) 100 * pe else NA_real_,
       kappa = if (!is.na(pe) && pe < 1) (po - pe) / (1 - pe) else NA_real_,
       tab = tb)
}

if (length(x_laud_op_ids) >= 2) {
  ## Direction is computed for EVERY shared gene, not only the panel, on
  ## identical terms. Without the non-panel genes there is no null to judge the
  ## panel against, which is the same argument [B] makes for Tier 1.
  x_d_laud <- x_dir(count_mat, x_laud_wt_ids, x_laud_op_ids, x_common, lib = lib_sizes)
  x_d_dunc <- x_dir(x_dm,      x_wt,          x_mut,         x_common, lib = x_dlib)
  x_dir_all <- merge(x_d_laud, x_d_dunc, by = "gene_id",
                     suffixes = c("_lauderdale", "_duncan"))
  x_dir_all$panel <- x_dir_all$gene_id %in% x_panel
  x_dir_all <- merge(x_dir_all, x_wt_tab[, c("gene_id", "cpm_lauderdale")],
                     by = "gene_id", all.x = TRUE)
  x_dir_all$abund_bin <- cut(log10(x_dir_all$cpm_lauderdale + 0.01),
                             breaks = x_breaks, include.lowest = TRUE, labels = FALSE)

  x_dir_tab <- x_dir_all[x_dir_all$panel,
                         c("gene_id", "direction_lauderdale", "bounded_fold_lauderdale",
                           "direction_duncan", "bounded_fold_duncan")]
  x_obs <- x_dir_agree_stats(x_dir_all[x_dir_all$panel, ])

  ## Exact conditional test, both margins fixed. This asks whether the two
  ## datasets' direction calls are ASSOCIATED, which is the question, and it
  ## conditions away the marginal skew that makes the raw percentage look
  ## impressive. It still treats genes as independent -- see the note printed
  ## below -- so the permutation null is the number to prefer.
  x_dir_fisher <- if (sum(x_obs$tab) > 0)
    stats::fisher.test(x_obs$tab, alternative = "greater")$p.value else NA_real_

  ## Abundance-matched permutation null, the same device as Tier 1. Background
  ## genes are drawn to match the panel's abundance profile and scored the same
  ## way. This is the strongest of the three because the background carries a
  ## comparable correlation structure: non-panel genes in these same two
  ## contrasts are co-regulated too, so a coherent block moving together is
  ## present in the null as well as in the panel.
  x_dp <- x_dir_all[x_dir_all$panel, ]
  x_db <- x_dir_all[!x_dir_all$panel & !is.na(x_dir_all$abund_bin), ]
  x_dp_bins <- x_dp$abund_bin[!is.na(x_dp$abund_bin)]

  ## One draw: the panel's per-bin abundance profile, scaled by `mult`.
  x_draw <- function(mult) {
    idx <- integer(0)
    for (b in unique(x_dp_bins)) {
      need <- max(1L, round(sum(x_dp_bins == b) * mult))
      pool <- which(x_db$abund_bin == b)
      if (!length(pool)) next
      idx <- c(idx, pool[sample.int(length(pool), need, replace = length(pool) < need)])
    }
    idx
  }

  ## SIZE MATCHING, and why it is not optional. `resolved in both` requires a
  ## bounded fold >= the categorical threshold in BOTH datasets, which depends
  ## on effect size, not only on abundance. Panel genes clear it far more often
  ## than abundance-matched background, so a draw of the same SIZE resolves far
  ## fewer genes -- and a percentage computed on a dozen genes is both noisier
  ## and pushed toward 100% by the ceiling. Comparing it with a 47-gene
  ## percentage would not be like for like. The draw is therefore scaled until
  ## the null resolves a comparable NUMBER of comparisons.
  x_mults <- c(1, 2, 4, 8, 16, 32)
  x_med_n <- vapply(x_mults, function(m)
    stats::median(vapply(seq_len(200), function(i)
      x_dir_agree_stats(x_db[x_draw(m), ])$n_resolved, numeric(1))), numeric(1))
  x_mult <- x_mults[which.min(abs(x_med_n - x_obs$n_resolved))]

  x_dnull   <- rep(NA_real_,    X_N_PERM)
  x_dnull_n <- rep(NA_integer_, X_N_PERM)
  for (i in seq_len(X_N_PERM)) {
    st <- x_dir_agree_stats(x_db[x_draw(x_mult), ])
    x_dnull[i]   <- st$pct_agree
    x_dnull_n[i] <- st$n_resolved
  }
  x_dnull_ok <- x_dnull[is.finite(x_dnull)]
  x_dir_p <- if (length(x_dnull_ok) && is.finite(x_obs$pct_agree))
    (1 + sum(x_dnull_ok >= x_obs$pct_agree)) / (1 + length(x_dnull_ok)) else NA_real_

  cat("\n=== TIER 2: direction of the genotype effect, panel genes ===\n")
  cat(sprintf("  resolved in both: %d of %d; agree on direction: %d (%.1f%%)\n",
              x_obs$n_resolved, x_obs$n_total, x_obs$n_agree, x_obs$pct_agree))
  cat("\n  Direction cross-tabulation (rows this study, cols GSE183742):\n")
  print(x_obs$tab)
  cat(sprintf("\n  Agreement expected from the MARGINS ALONE: %.1f%%\n", x_obs$pct_chance))
  cat(sprintf("  Cohen's kappa (chance-corrected agreement)  : %.3f\n", x_obs$kappa))
  cat(sprintf("  Fisher exact, both margins fixed, one-sided : %.4g\n", x_dir_fisher))
  cat(sprintf("  Abundance-matched permutation null: observed %.1f%%, null median %.1f%%\n",
              x_obs$pct_agree, stats::median(x_dnull_ok)))
  cat(sprintf("    (2.5-97.5%%: %.1f-%.1f); p = %.4g\n",
              stats::quantile(x_dnull_ok, 0.025), stats::quantile(x_dnull_ok, 0.975),
              x_dir_p))
  cat(sprintf("    size-matched: draw scaled %gx so the null resolves a median of %d\n",
              x_mult, as.integer(stats::median(x_dnull_n, na.rm = TRUE))))
  cat(sprintf("    comparison(s) against the panel's %d. Unscaled draws resolve a median of %d.\n",
              x_obs$n_resolved, as.integer(x_med_n[x_mults == 1])))

  cat("\n  WHICH NUMBER TO QUOTE, and why not the obvious one.\n")
  cat("    A binomial test against 0.5 is NOT reported here, deliberately. It\n")
  cat("    asks whether the two datasets agree more often than a coin, but the\n")
  cat("    margins above are heavily skewed toward one direction, so two\n")
  cat("    UNRELATED sets of calls would already agree at the rate printed as\n")
  cat("    `expected from the margins alone`. Testing against 0.5 therefore\n")
  cat("    returns a p value that is too small by many orders of magnitude.\n")
  cat("    Report the counts, the 2x2, kappa, and the permutation p. Fisher is\n")
  cat("    the right fixed-margin test but still assumes genes are independent\n")
  cat("    of one another, which co-regulated panel members are not; only the\n")
  cat("    permutation null puts comparable dependence into the reference.\n")
  cat("    Note what the null itself shows: background genes that clear the\n")
  cat("    categorical threshold in BOTH datasets agree on direction almost\n")
  cat("    always. Direction agreement among strongly categorical genes is a\n")
  cat("    generic property, so a high panel figure is not panel-specific\n")
  cat("    evidence unless it exceeds this null.\n")

  cat("\n  CAVEATS, to be carried into any sentence quoting these numbers:\n")
  cat("    - different allele: Pax6^Sey-Neu/+ here, Pax6^tm1Pgr in GSE183742\n")
  cat("    - n = 3 per group in GSE183742\n")
  cat("    - corneal transparency was NOT REPORTED for GSE183742; the animals\n")
  cat("      were probably fibrotic, so the contrast is compared against\n")
  cat("      SeyO_vs_WT_Adult and cannot be assigned to genotype alone\n")
  cat("    - the genes resolved in both are selected for a LARGE categorical\n")
  cat("      difference in both datasets, which favours agreement; the\n")
  cat("      permutation null is drawn under the same selection\n")
  cat("    - the dissection boundary is NOT REPORTED for GSE183742; see the\n")
  cat("      dissection-margin check below for what the markers say\n")

  ## x_dir_tab itself is written in the Outputs section, once symbols are on it.
  write.csv(as.data.frame(x_obs$tab),
            file.path(X_OUT, "CrossDataset_Direction_CrossTab.csv"), row.names = FALSE)
  write.csv(data.frame(null_pct_agree = x_dnull, null_n_resolved = x_dnull_n),
            file.path(X_OUT, "CrossDataset_Direction_NullDistribution.csv"), row.names = FALSE)
} else {
  cat("\n[skip] TIER 2: fewer than two opaque adult mutant libraries found.\n")
  x_dir_tab <- NULL
}

## ---------------------------------------------------------------------------
## THE EXCEPTIONS, BY NAME
## ---------------------------------------------------------------------------
## Both tiers are summarised as counts, and the genes behind the counts are the
## part a reader asks about. They are printed here with symbols, and the
## per-gene outputs carry a symbol column, so that they never have to be
## recovered from Ensembl identifiers by hand. panel_tbl comes from 01_load.R.
x_sym <- unique(panel_tbl[, c("gene_id", "symbol")])
x_symbol_of <- function(ids) x_sym$symbol[match(ids, x_sym$gene_id)]

x_wt_tab$symbol <- x_symbol_of(x_wt_tab$gene_id)
x_t1_disagree <- x_wt_tab[x_wt_tab$panel &
                          x_wt_tab$status_lauderdale != "indeterminate" &
                          x_wt_tab$status_duncan     != "indeterminate" &
                          x_wt_tab$status_lauderdale != x_wt_tab$status_duncan, ]
cat("\n=== TIER 1 exceptions: panel genes resolved in both and called differently ===\n")
if (nrow(x_t1_disagree)) {
  print(data.frame(symbol = x_t1_disagree$symbol,
                   this_study = x_t1_disagree$status_lauderdale,
                   cpm_this_study = signif(x_t1_disagree$cpm_lauderdale, 3),
                   GSE183742 = x_t1_disagree$status_duncan,
                   cpm_GSE183742 = signif(x_t1_disagree$cpm_duncan, 3)),
        row.names = FALSE)
  x_t1_dir <- table(factor(paste(x_t1_disagree$status_lauderdale, "here,",
                                 x_t1_disagree$status_duncan, "there")))
  cat("  Direction of the exceptions:\n")
  print(x_t1_dir)
  cat("  Exceptions that all run one way -- absent here, present there -- are\n")
  cat("  what tissue ADDED to one dissection produces, and what a library-\n")
  cat("  chemistry difference does not. See the dissection-margin check below.\n")
} else cat("  none\n")

if (!is.null(x_dir_tab)) {
  x_dir_tab$symbol <- x_symbol_of(x_dir_tab$gene_id)
  x_dir_tab <- x_dir_tab[, c("symbol", setdiff(names(x_dir_tab), "symbol"))]
  x_lv <- c("up_in_mutant", "down_in_mutant")
  x_t2_disagree <- x_dir_tab[x_dir_tab$direction_lauderdale %in% x_lv &
                             x_dir_tab$direction_duncan %in% x_lv &
                             x_dir_tab$direction_lauderdale != x_dir_tab$direction_duncan, ]
  cat("\n=== TIER 2 exceptions: resolved in both, opposite direction ===\n")
  if (nrow(x_t2_disagree)) {
    x_t2_wt <- x_wt_tab[match(x_t2_disagree$gene_id, x_wt_tab$gene_id), ]
    print(data.frame(symbol = x_t2_disagree$symbol,
                     this_study = x_t2_disagree$direction_lauderdale,
                     bounded_fold_this_study = round(x_t2_disagree$bounded_fold_lauderdale, 2),
                     GSE183742 = x_t2_disagree$direction_duncan,
                     bounded_fold_GSE183742 = round(x_t2_disagree$bounded_fold_duncan, 2),
                     WT_cpm_this_study = signif(x_t2_wt$cpm_lauderdale, 3),
                     WT_cpm_GSE183742 = signif(x_t2_wt$cpm_duncan, 3)),
          row.names = FALSE)
    cat("  The two wild-type columns say whether a reversal is two opposite mutant\n")
    cat("  responses or one wild-type baseline that is present in one dataset and\n")
    cat("  absent in the other. Read them before calling it a biological difference.\n")
  } else cat("  none\n")
}

## ---------------------------------------------------------------------------
## DISSECTION-MARGIN CHECK
## ---------------------------------------------------------------------------
## The Tier 1 exceptions are endothelial and neural-crest genes present at a
## few CPM in the GSE183742 wild type and near zero here. That is the signature
## of a dissection that includes some limbus or iridocorneal angle, tissues
## that carry vessels, meshwork, pigment and nerve that the central cornea
## does not. The GEO record does not say where the corneas were cut. This
## block asks the count matrices directly, using marker genes chosen for the
## tissues a generous dissection would add and NOT for any result above: if
## limbal, meshwork, pigment or lens markers sit at a few CPM in one dataset's
## wild type and near zero in the other's, the dissections differ, and every
## cross-dataset sentence has to carry that.
##
## Four corneal markers are included as positive controls, two epithelial
## (Krt12, Aldh3a1) and two stromal (Kera, Angptl7). They should be high in
## every group. A control that is present in both wild types but MANY-FOLD
## lower in one says the two dissections took different proportions of
## epithelium and stroma, which matters as much as added tissue: a library
## with less stroma has less of everything the stroma carries. So the controls
## are judged on fold, not on status, with the tolerance X_CONTROL_MAX_FOLD.
## Angptl7 (CDT6) is a keratocyte gene abundant in corneal stroma; an earlier
## version listed it under trabecular meshwork, where it is also expressed,
## and it then flagged a stromal difference as an angle difference.
##
## THE MARKERS AND THE TOLERANCE COME FROM 00_config.R (OFFTARGET_MARKERS,
## OFFTARGET_CONTROL_MAX_FOLD), shared with 01c_contamination_indices.R, and
## the identifiers are the ones 01c resolved against the count matrix. This
## script kept its own list until 2026-09-24, and it had drifted from 01c's:
## it still counted Pmel as pigment (present at 12-55 CPM in every cornea) and
## Chi3l1 as angle (induced with inflammation in opaque cornea, so it read
## opacity as a dissection difference). Both are now excluded, with the reason
## recorded in the config table, and the lens and retina sets are 01c's.
X_CONTROL_MAX_FOLD <- OFFTARGET_CONTROL_MAX_FOLD
X_MARGIN_MARKERS <- OFFTARGET_MARKERS[OFFTARGET_MARKERS$margin &
                                      is.na(OFFTARGET_MARKERS$excluded),
                                      c("symbol", "compartment", "family")]
X_STROMAL_CONTROLS    <- X_MARGIN_MARKERS$symbol[X_MARGIN_MARKERS$family %in% "stromal_control"]
X_EPITHELIAL_CONTROLS <- X_MARGIN_MARKERS$symbol[X_MARGIN_MARKERS$family %in% "epithelial_control"]

if (nrow(X_MARGIN_MARKERS)) {
  X_MARGIN_MARKERS$gene_id <- unname(ot_marker_ids[X_MARGIN_MARKERS$symbol])
  x_unmapped <- X_MARGIN_MARKERS$symbol[is.na(X_MARGIN_MARKERS$gene_id) |
                                        !X_MARGIN_MARKERS$gene_id %in% x_common]
  if (length(x_unmapped)) {
    cat("\n  Margin markers not mapped or not in both matrices, dropped: ",
        paste(x_unmapped, collapse = ", "), "\n")
  }
  x_mm <- X_MARGIN_MARKERS[!X_MARGIN_MARKERS$symbol %in% x_unmapped, ]
  x_mm$on_panel <- x_mm$gene_id %in% panel_genes

  ## The Lauderdale transparent adult mutant, for the same reason the opaque
  ## set is used in Tier 2: a margin difference between GSE183742's own groups
  ## would confound their contrast, and the same question is asked of ours.
  x_laud_tr_ids <- meta$sample_id[meta$tissue == "Cornea" &
                                  meta$laboratory %in% LABORATORY &
                                  meta$age_group == "Adult" &
                                  meta$genotype == "Sey" &
                                  !meta$sample_id %in% x_excl &
                                  toupper(substr(as.character(meta$transparency), 1, 1)) == "T"]
  x_laud_tr_ids <- intersect(x_laud_tr_ids, colnames(count_mat))

  x_margin_group <- function(cm, ids, lib, label) {
    if (length(ids) < 2) return(NULL)
    s <- x_status(cm[x_mm$gene_id, , drop = FALSE], ids, lib = lib)
    data.frame(gene_id = s$gene_id, group = label, status = s$status,
               cpm = s$cpm, k = s$k, L = s$L, stringsAsFactors = FALSE)
  }
  x_margin_long <- rbind(
    x_margin_group(count_mat, x_laud_wt_ids, lib_sizes, "this_study_WT"),
    x_margin_group(count_mat, x_laud_tr_ids, lib_sizes, "this_study_Sey_transparent"),
    x_margin_group(count_mat, x_laud_op_ids, lib_sizes, "this_study_Sey_opaque"),
    x_margin_group(x_dm, x_wt,  x_dlib, "GSE183742_WT"),
    x_margin_group(x_dm, x_mut, x_dlib, "GSE183742_mutant"))

  x_margin <- x_mm
  for (g in unique(x_margin_long$group)) {
    sub <- x_margin_long[x_margin_long$group == g, ]
    x_margin[[paste0("cpm_", g)]]    <- signif(sub$cpm[match(x_margin$gene_id, sub$gene_id)], 3)
    x_margin[[paste0("status_", g)]] <- sub$status[match(x_margin$gene_id, sub$gene_id)]
  }
  ## Conservative fold between the two wild types, from the Poisson limits, as
  ## everywhere else in this script.                                        [5]
  x_a <- x_margin_long[x_margin_long$group == "this_study_WT", ]
  x_b <- x_margin_long[x_margin_long$group == "GSE183742_WT", ]
  x_a <- x_a[match(x_margin$gene_id, x_a$gene_id), ]
  x_b <- x_b[match(x_margin$gene_id, x_b$gene_id), ]
  x_margin$bounded_fold_WT <- round(x_bounded_fold(x_b$k, x_b$L, x_a$k, x_a$L), 2)
  ## A bounded fold below 1 means the two Poisson intervals overlap, and then
  ## neither wild type is "higher" in any sense this script would defend.
  x_margin$higher_in_WT <- ifelse(x_margin$bounded_fold_WT < 1, "intervals overlap",
                                  ifelse(x_b$cpm > x_a$cpm, "GSE183742", "this study"))

  cat("\n=== DISSECTION-MARGIN CHECK: off-panel markers of tissue a generous cut adds ===\n")
  print(x_margin[, c("symbol", "compartment", "on_panel",
                     "cpm_this_study_WT", "cpm_GSE183742_WT",
                     "bounded_fold_WT", "higher_in_WT")], row.names = FALSE)

  x_ctrl <- x_margin[x_margin$compartment == "central cornea (control)", ]
  x_test <- x_margin[x_margin$compartment != "central cornea (control)", ]
  x_flag <- x_test[x_test$status_GSE183742_WT == "expressed" &
                   x_test$status_this_study_WT != "expressed" &
                   x_test$bounded_fold_WT >= CATEGORICAL_MIN_FOLD &
                   x_test$higher_in_WT == "GSE183742", ]
  x_flag_rev <- x_test[x_test$status_this_study_WT == "expressed" &
                       x_test$status_GSE183742_WT != "expressed" &
                       x_test$bounded_fold_WT >= CATEGORICAL_MIN_FOLD &
                       x_test$higher_in_WT == "this study", ]
  x_ctrl_ok <- x_ctrl$status_this_study_WT == "expressed" &
               x_ctrl$status_GSE183742_WT == "expressed" &
               x_ctrl$bounded_fold_WT < X_CONTROL_MAX_FOLD
  cat(sprintf("\n  Controls expressed in both wild types AND within %g-fold: %d of %d\n",
              X_CONTROL_MAX_FOLD, sum(x_ctrl_ok), nrow(x_ctrl)))
  if (any(!x_ctrl_ok)) {
    x_bad <- x_ctrl[!x_ctrl_ok, ]
    cat("  CONTROLS OUTSIDE TOLERANCE -- the two wild types are not the same tissue mix:\n")
    print(data.frame(symbol = x_bad$symbol,
                     cpm_this_study_WT = x_bad$cpm_this_study_WT,
                     cpm_GSE183742_WT = x_bad$cpm_GSE183742_WT,
                     bounded_fold_WT = x_bad$bounded_fold_WT,
                     higher_in = x_bad$higher_in_WT), row.names = FALSE)
    x_bad_stromal <- intersect(x_bad$symbol, X_STROMAL_CONTROLS)
    x_bad_epi     <- intersect(x_bad$symbol, X_EPITHELIAL_CONTROLS)
    if (length(x_bad_stromal) && !length(x_bad_epi))
      cat("  Both stromal controls are out of tolerance and both epithelial controls\n",
          "  are within it: the libraries differ in their STROMAL share, not in\n",
          "  whether they are cornea. Stromal genes are then depleted in the dataset\n",
          "  with the lower Kera/Angptl7 before any biology is invoked.\n", sep = "")
  }
  cat(sprintf("  Margin markers expressed in GSE183742 wild type, not here, >= %.0f-fold: %d of %d",
              CATEGORICAL_MIN_FOLD, nrow(x_flag), nrow(x_test)))
  if (nrow(x_flag)) cat(" (", paste(x_flag$symbol, collapse = ", "), ")", sep = "")
  cat("\n")
  cat(sprintf("  Margin markers expressed here, not in GSE183742 wild type, >= %.0f-fold: %d",
              CATEGORICAL_MIN_FOLD, nrow(x_flag_rev)))
  if (nrow(x_flag_rev)) cat(" (", paste(x_flag_rev$symbol, collapse = ", "), ")", sep = "")
  cat("\n")

  ## Does the margin differ between the two groups of ONE dataset? If it does,
  ## that contrast carries a dissection difference as well as genotype (and
  ## opacity), and any comparison built on it inherits that. The question is
  ## asked of GSE183742's groups, because Tier 2 uses their contrast, and of
  ## this study's own adult groups, because the transparent-mutant-vs-WT and
  ## opaque-vs-transparent comparisons are what the paper rests on. A pigment,
  ## angle or lens marker that rises with genotype or with opacity is either
  ## tissue the dissection took (adherent iris, angle) or a change in the cornea
  ## itself; the count matrix cannot tell those apart, and a sentence about a
  ## gene of that tissue has to say so.
  x_within_check <- function(group_ref, group_test, label) {
    x_r <- x_margin_long[x_margin_long$group == group_ref, ]
    x_t <- x_margin_long[x_margin_long$group == group_test, ]
    if (!nrow(x_r) || !nrow(x_t)) return(invisible(NULL))
    x_r <- x_r[match(x_margin$gene_id, x_r$gene_id), ]
    x_t <- x_t[match(x_margin$gene_id, x_t$gene_id), ]
    bf  <- x_bounded_fold(x_t$k, x_t$L, x_r$k, x_r$L)
    tab <- data.frame(symbol = x_margin$symbol, compartment = x_margin$compartment,
                      cpm_ref = signif(x_r$cpm, 3), cpm_test = signif(x_t$cpm, 3),
                      bounded_fold = round(bf, 2),
                      higher_in = ifelse(bf < 1, "intervals overlap",
                                         ifelse(x_t$cpm > x_r$cpm, group_test, group_ref)),
                      stringsAsFactors = FALSE)
    names(tab)[3:4] <- paste0("cpm_", c(group_ref, group_test))
    flag <- tab[bf >= CATEGORICAL_MIN_FOLD &
                (x_r$status == "expressed" | x_t$status == "expressed") &
                x_margin$compartment != "central cornea (control)", ]
    cat(sprintf("\n  %s: margin markers differing >= %.0f-fold: %d of %d",
                label, CATEGORICAL_MIN_FOLD, nrow(flag),
                sum(x_margin$compartment != "central cornea (control)")))
    if (nrow(flag)) {
      cat("\n"); print(flag, row.names = FALSE)
    } else cat(" (none)\n")
    x_margin[[paste0("bounded_fold_", group_test, "_vs_", group_ref)]] <<- round(bf, 2)
    x_margin[[paste0("higher_in_", group_test, "_vs_", group_ref)]]   <<- tab$higher_in
    invisible(flag)
  }
  x_within_flag <- x_within_check("GSE183742_WT", "GSE183742_mutant",
                                  "GSE183742 mutant vs wild type (the Tier 2 contrast)")
  x_within_check("this_study_WT", "this_study_Sey_transparent",
                 "This study, transparent adult mutant vs wild type")
  x_within_check("this_study_Sey_transparent", "this_study_Sey_opaque",
                 "This study, opaque vs transparent adult mutant")
  if (!is.null(x_within_flag) && nrow(x_within_flag) &&
      exists("x_t2_disagree") && nrow(x_t2_disagree)) {
    cat("\n  The Tier 2 reversal(s) above -- ", paste(x_t2_disagree$symbol, collapse = ", "),
        " -- should be read against this list: a gene of a tissue whose marker\n",
        "  differs between GSE183742's own groups is a dissection difference until\n",
        "  shown otherwise, not a difference between the alleles.\n", sep = "")
  }

  cat("\n  HOW TO READ THIS. Markers of limbus, angle, pigment or lens at a few CPM\n")
  cat("  in one wild type and near zero in the other mean the two dissections\n")
  cat("  took different amounts of peripheral tissue, and every cross-dataset\n")
  cat("  disagreement in a gene of those tissues is then explained before any\n")
  cat("  biology is invoked. The same markers rising between two groups of ONE\n")
  cat("  dataset add a dissection term to that contrast that the count matrix\n")
  cat("  cannot separate from the biology. Neither can be excluded from the GEO\n")
  cat("  record, which does not state the cut.\n")

  write.csv(x_margin, file.path(X_OUT, "CrossDataset_DissectionMargin.csv"), row.names = FALSE)

  ## PER LIBRARY, ON THE TERMS OF SUPP. TABLE 6. The block above compares
  ## pooled groups; this one judges each GSE183742 library by the rule 01c
  ## applies to ours, against the reference 01c uses for an adult cornea (this
  ## study's adult wild type and transparent mutant). The same flags, the same
  ## thresholds and the same stroma-poor criterion, so a sentence comparing the
  ## two datasets' contamination rests on one rule rather than two.
  x_ot_ref  <- ot_marker_cpm(count_mat, lib_sizes, c(x_laud_wt_ids, x_laud_tr_ids))
  x_ot_dcpm <- ot_marker_cpm(x_dm, x_dlib, c(x_wt, x_mut))
  x_ot_duncan <- do.call(rbind, lapply(colnames(x_ot_dcpm), function(s) {
    cbind(data.frame(library_id = s,
                     group = if (s %in% x_wt) "GSE183742_WT" else "GSE183742_mutant",
                     stringsAsFactors = FALSE),
          as.data.frame(ot_indices(x_ot_dcpm[, s, drop = FALSE])[, -1]),
          as.data.frame(ot_flag_library(x_ot_dcpm[, s], x_ot_ref)))
  }))
  cat(sprintf("\n=== GSE183742 per library, by the Supp. Table 6 rule (reference: this study's adult WT + transparent Sey, n = %d) ===\n",
              ncol(x_ot_ref)))
  print(x_ot_duncan[, c("library_id", "group", "Kera", "Rho", "pigment",
                        "contamination_flags")], row.names = FALSE)
  write.csv(x_ot_duncan, file.path(X_OUT, "CrossDataset_GSE183742_OffTarget.csv"),
            row.names = FALSE)
} else {
  cat("\n[skip] DISSECTION-MARGIN CHECK: no margin markers defined in OFFTARGET_MARKERS.\n")
}

## ---------------------------------------------------------------------------
## Outputs
## ---------------------------------------------------------------------------
write.csv(x_wt_tab[x_wt_tab$panel, c("symbol", setdiff(names(x_wt_tab), "symbol"))],
          file.path(X_OUT, "CrossDataset_WT_StatusConcordance.csv"), row.names = FALSE)
write.csv(as.data.frame(x_panel_conc$tab),
          file.path(X_OUT, "CrossDataset_WT_CrossTab.csv"), row.names = FALSE)
write.csv(data.frame(null_pct_agree = x_null),
          file.path(X_OUT, "CrossDataset_NullDistribution.csv"), row.names = FALSE)
if (!is.null(x_dir_tab)) {
  write.csv(x_dir_tab, file.path(X_OUT, "CrossDataset_Direction.csv"), row.names = FALSE)
}

cat("\n==========================================================\n")
cat("09_cross_dataset_replication complete. Outputs in:\n  ", X_OUT, "\n", sep = "")
cat("  CrossDataset_WT_StatusConcordance.csv   per-gene, both datasets, with symbols\n")
cat("  CrossDataset_WT_CrossTab.csv            the 2x2 on resolved genes\n")
cat("  CrossDataset_NullDistribution.csv       abundance-matched null\n")
if (!is.null(x_dir_tab)) {
  cat("  CrossDataset_Direction.csv              tier 2, per gene, with symbols\n")
  cat("  CrossDataset_Direction_CrossTab.csv     tier 2, the 2x2 on resolved genes\n")
  cat("  CrossDataset_Direction_NullDistribution.csv  tier 2, abundance-matched null\n")
}
if (file.exists(file.path(X_OUT, "CrossDataset_DissectionMargin.csv"))) {
  cat("  CrossDataset_DissectionMargin.csv       off-panel margin markers, every group\n")
  cat("  CrossDataset_GSE183742_OffTarget.csv    GSE183742 per library, Supp. Table 6 rule\n")
}
cat("\nWhat this analysis CANNOT say: that a gene behaves identically in the\n")
cat("two alleles. It says whether the same expressed/not-expressed calls are\n")
cat("recovered in independent adult cornea, which is a statement about the\n")
cat("robustness of the calls, not about equivalence of the mutants.\n")
cat("==========================================================\n")
