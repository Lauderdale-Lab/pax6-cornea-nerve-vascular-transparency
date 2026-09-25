###############################################################################
## 05_program_level.R -- abundance-matched program-level permutation tests
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
## Archived   : <ZENODO_DOI>
## Licence    : MIT (see LICENSE)
## Contact    : James D. Lauderdale, jdlauder@uga.edu
##
## Run the pipeline with run_all.R. Scripts are numbered in execution order and
## share one R session by design; see run_all.R for why.
###############################################################################

###############################################################################
## 05_program_level.R
##
## Asks whether a whole curated program is more, or less, of it confidently
## expressed in one group than another -- a question about the panel rather than
## about any gene in it.
##
## WHY THIS TEST EXISTS. Every per-gene result in this analysis depends on
## thresholds, on a stability check, and on how many libraries happened to fall
## in each group. Counting how many panel genes are confidently expressed in one
## group and not the other inherits all of that. A group whose calls are less
## resolved -- because its libraries are less complex, or because the stability
## check bit harder there -- will appear to express fewer genes for reasons that
## have nothing to do with biology.
##
## The control is to ask the same question of genes that are NOT in the panel
## but sit at the same abundances, drawn from the same libraries and scored by
## exactly the same rule. If the difference is an artefact of how these
## particular libraries were sequenced, the matched control genes show it too,
## and the null centres away from zero. What is reported is the EXCESS over that
## null, not the raw difference. [1] [2]
##
## This is the one instrument in the analysis that no per-gene fragility
## touches, and it is where a negative result can be given a magnitude: the
## permutation distribution says directly how large an excess would have been
## detected, so "no program-level difference" is reported with that limit
## attached rather than as bare absence.
##
## Read that limit against the panel's own size. A small panel produces a narrow
## null and therefore detects a SMALLER excess in absolute gene counts, not a
## larger one; as a proportion of the panel the sensitivity of all three is
## similar. A negative from the 25-gene panel is not the weak statement it might
## appear to be.
##
## WHAT IS COMPARED. For each panel and contrast, two counts:
##   delta_not         panel genes confidently SILENT in group A minus those in B
##   delta_expressed   panel genes confidently EXPRESSED in A minus those in B
## Both use the same gated calls as everything else, so the cost of the
## stability check is present in the observed value AND in the null, and
## therefore cancels.
##
## THE CONTROL UNIVERSE is autosomal, non-panel, and detected somewhere. Sex
## chromosomes are excluded because the libraries are mixed-sex pools of varying
## composition; the mitochondrial genome because its abundance does not behave
## like a nuclear gene's. Genes with no counts anywhere are excluded here only:
## they carry no information and would inflate the silent side of every drawn
## set. They remain in the panel denominators, where a gene never detected is a
## legitimate result.
##
## LIMITATION TO STATE IN THE PAPER. This is a competitive test: the panel is
## compared against other genes rather than against a null of no change within
## itself. Panel genes are co-regulated, and positive correlation between genes
## in a set inflates competitive tests relative to their nominal level. [3] The
## permutation preserves the correlation among the CONTROL genes drawn, but it
## cannot reproduce the panel's own internal correlation, so the P values here
## should be read as ordering the panels and contrasts rather than as exact.
##
## STATISTICAL BASIS
##   [1] Competitive gene-set testing, and the distinction from a self-contained
##       null. Goeman JJ, Buhlmann P (2007) Bioinformatics 23:980-987.
##   [2] Matching the control set on abundance, because detection probability
##       depends on it and an unmatched control would measure that instead of
##       the effect. The same reasoning as length and count bias corrections in
##       gene-set analysis: Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010)
##       Genome Biology 11:R14.
##   [3] Inter-gene correlation inflates competitive set tests. Wu D, Smyth GK
##       (2012) Nucleic Acids Research 40:e133.
##   [4] Permutation P values are computed as (1 + r) / (n + 1) and can never be
##       zero; a value at that floor means only "smaller than the resolution of
##       this many permutations". Phipson B, Smyth GK (2010) Statistical
##       Applications in Genetics and Molecular Biology 9:Article 39.
##   [5] Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300, applied once
##       across every test in this script.
##
## Inputs : 00_config.R, 01_load.R, 03_expression_status.R -- run in that order
##          in one session. This script needs transcriptome-wide expression
##          status, which is not written to disk.
## Outputs: ProgramLevel.csv
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(purrr)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

if (!exists("cell_status_genomewide")) {
  stop("Transcriptome-wide expression status not found. Source 00_config.R, ",
       "01_load.R and 03_expression_status.R in that order in one session; ",
       "this analysis scores all ", nrow(count_mat), " genes, not just the ",
       "panels, and that is not written to disk.")
}

## ---------------------------------------------------------------------------
## 1. Control universe
## ---------------------------------------------------------------------------

chromosome <- try(suppressWarnings(suppressMessages(
  AnnotationDbi::mapIds(org.Mm.eg.db, keys = rownames(count_mat), column = "CHR",
                        keytype = "ENSEMBL", multiVals = "first"))), silent = TRUE)

if (inherits(chromosome, "try-error")) {
  message("NOTE: chromosome annotation unavailable; the control universe will ",
          "include sex-linked and mitochondrial genes. Pool sex composition ",
          "varies between libraries, so this weakens the control.")
  autosomal <- rep(TRUE, nrow(count_mat))
} else {
  autosomal <- !(unname(chromosome[rownames(count_mat)]) %in% c("X", "Y", "MT"))
  autosomal[is.na(autosomal)] <- TRUE   # unannotated genes are kept
}

detected_anywhere <- rowSums(count_mat) > 0
universe <- rownames(count_mat)[autosomal & detected_anywhere &
                                !rownames(count_mat) %in% panel_genes]

message("Control universe: ", length(universe), " genes (",
        sum(!autosomal), " sex-linked or mitochondrial, ",
        sum(!detected_anywhere), " never detected, and ", length(panel_genes),
        " panel genes excluded).")

## Gene id to row position, built once. Every lookup in the permutation loop is
## then integer subsetting rather than a hash lookup on a 55,000-element name
## vector, which is what made an earlier version appear to hang at 20,000
## permutations.
gene_position <- stats::setNames(seq_len(nrow(count_mat)), rownames(count_mat))
universe_idx <- unname(gene_position[universe])

## ---------------------------------------------------------------------------
## 2. The test
## ---------------------------------------------------------------------------

## The hinge contrast appears twice in CONTRASTS, once per differential model.
## This test depends only on the expression calls of the two groups, which are
## identical, so it is run once under the Model B label.
prog_contrasts <- CONTRASTS[CONTRASTS$label != HINGE[["A"]] &
                            CONTRASTS$cell_a %in% names(cell_status_genomewide) &
                            CONTRASTS$cell_b %in% names(cell_status_genomewide), ,
                            drop = FALSE]
panels <- names(PANEL_DIRS)
n_blocks <- nrow(prog_contrasts) * length(panels)
block <- 0L

set.seed(SEED)

program <- purrr::pmap_dfr(
  list(prog_contrasts$label, prog_contrasts$cell_a,
       prog_contrasts$cell_b, prog_contrasts$axis),
  function(label, cell_a, cell_b, axis) {

    ## Abundance is measured on the two groups pooled, so the matching variable
    ## is neutral between them and cannot itself encode the difference.
    ids <- unique(c(groups$ids[[match(cell_a, groups$cell)]],
                    groups$ids[[match(cell_b, groups$cell)]]))
    abundance <- log10(rowSums(count_mat[, ids, drop = FALSE]) /
                       sum(lib_sizes[ids]) * 1e6 + 1e-3)

    A <- cell_status_genomewide[[cell_a]][rownames(count_mat)]
    B <- cell_status_genomewide[[cell_b]][rownames(count_mat)]
    a_silent <- unname(A == "not_expressed"); a_expressed <- unname(A == "expressed")
    b_silent <- unname(B == "not_expressed"); b_expressed <- unname(B == "expressed")

    purrr::map_dfr(panels, function(pl) {
      block <<- block + 1L
      panel_idx <- unname(gene_position[panel_tbl$gene_id[panel_tbl$panel == pl]])

      breaks <- unique(stats::quantile(abundance[panel_idx],
                                       probs = seq(0, 1, length.out = N_BINS + 1)))
      if (length(breaks) < 3) {
        message("  [", pl, " / ", label, "] too few distinct abundance bins -- skipped")
        return(NULL)
      }
      panel_bin <- cut(abundance[panel_idx], breaks = breaks, include.lowest = TRUE)
      need <- as.integer(table(panel_bin))
      pools <- split(universe_idx,
                     cut(abundance[universe_idx], breaks = breaks, include.lowest = TRUE))

      ## Adequacy of the control pool, computed once. These are properties of
      ## the binning, not of the permutations.
      pool_n <- vapply(pools, function(z) if (is.null(z)) 0L else length(z), integer(1))
      bins_empty      <- sum(need > 0L & pool_n == 0L)
      bins_undersized <- sum(need > 0L & pool_n > 0L & pool_n < need)

      observed_not <- sum(a_silent[panel_idx]) - sum(b_silent[panel_idx])
      observed_exp <- sum(a_expressed[panel_idx]) - sum(b_expressed[panel_idx])

      message(sprintf("  [%d/%d] %s / %s ... %d permutations",
                      block, n_blocks, pl, label, N_PERM))

      null <- vapply(seq_len(N_PERM), function(i) {
        drawn <- unlist(lapply(seq_along(need), function(j) {
          nb <- need[j]
          if (nb == 0L) return(integer(0))
          pool <- pools[[j]]
          if (is.null(pool) || !length(pool)) return(integer(0))
          ## Index the pool explicitly. sample(pool, ...) would treat a pool of
          ## length one as the number 1:pool and draw arbitrary row positions.
          if (length(pool) < nb) pool[sample.int(length(pool), nb, replace = TRUE)]
          else                   pool[sample.int(length(pool), nb)]
        }), use.names = FALSE)
        c(sum(a_silent[drawn])    - sum(b_silent[drawn]),
          sum(a_expressed[drawn]) - sum(b_expressed[drawn]))
      }, numeric(2))

      ## Two-sided empirical P about the null's own centre. [4]
      empirical_p <- function(observed, draws) {
        centre <- mean(draws)
        (1 + sum(abs(draws - centre) >= abs(observed - centre))) / (N_PERM + 1)
      }
      ## The smallest excess this test would have flagged, from the null's own
      ## 95% interval. This is what gives a negative result a magnitude.
      detection_limit <- function(draws) {
        q <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
        round(max(abs(q - mean(draws))), 2)
      }

      tibble::tibble(
        panel = pl, contrast = label, axis = axis,
        observed_delta_silent    = observed_not,
        null_mean_silent         = round(mean(null[1, ]), 2),
        null_sd_silent           = round(stats::sd(null[1, ]), 2),
        excess_silent            = round(observed_not - mean(null[1, ]), 2),
        detectable_silent        = detection_limit(null[1, ]),
        p_silent                 = empirical_p(observed_not, null[1, ]),
        observed_delta_expressed = observed_exp,
        null_mean_expressed      = round(mean(null[2, ]), 2),
        null_sd_expressed        = round(stats::sd(null[2, ]), 2),
        excess_expressed         = round(observed_exp - mean(null[2, ]), 2),
        detectable_expressed     = detection_limit(null[2, ]),
        p_expressed              = empirical_p(observed_exp, null[2, ]),
        n_panel_genes = length(panel_idx), n_bins = length(need),
        n_bins_empty = bins_empty, n_bins_undersized = bins_undersized)
    })
  })

## ---------------------------------------------------------------------------
## 3. Multiplicity and outcomes
## ---------------------------------------------------------------------------

## Every test in this script is one family: two counts for each panel and
## contrast, all addressing the same question. Reported separately they would
## invite exactly the over-reading this design exists to prevent. The two counts
## are strongly dependent, which makes a joint correction conservative rather
## than anti-conservative. [5]
p_all <- c(program$p_silent, program$p_expressed)
padj_all <- stats::p.adjust(p_all, method = "BH")
n_tests <- length(p_all)
p_floor <- 1 / (N_PERM + 1)

program <- program %>%
  dplyr::mutate(
    padj_silent    = padj_all[seq_len(nrow(program))],
    padj_expressed = padj_all[nrow(program) + seq_len(nrow(program))],
    ## Three outcomes, not two. A count whose excess reaches its own detection
    ## limit but whose adjusted P does not clear FDR is not a negative: the
    ## sentence "no difference detected (an excess of X would have been)" would
    ## be self-contradictory, because an excess of X WAS observed. Such a
    ## comparison sits at the unadjusted 5% boundary and is lost to the joint
    ## correction; report it as unresolved, with both numbers.
    outcome_silent = dplyr::case_when(
      padj_silent < FDR ~ "differs",
      abs(excess_silent) >= detectable_silent ~
        paste0("unresolved (excess ", excess_silent, " reaches the detectable ",
               detectable_silent, " but padj ", signif(padj_silent, 2), " after correction)"),
      TRUE ~ paste0("no difference detected (an excess of ", detectable_silent,
                    " genes would have been)")),
    outcome_expressed = dplyr::case_when(
      padj_expressed < FDR ~ "differs",
      abs(excess_expressed) >= detectable_expressed ~
        paste0("unresolved (excess ", excess_expressed, " reaches the detectable ",
               detectable_expressed, " but padj ", signif(padj_expressed, 2), " after correction)"),
      TRUE ~ paste0("no difference detected (an excess of ", detectable_expressed,
                    " genes would have been)")),
    at_permutation_floor = p_silent <= p_floor | p_expressed <= p_floor) %>%
  dplyr::relocate(padj_silent, .after = p_silent) %>%
  dplyr::relocate(padj_expressed, .after = p_expressed)

readr::write_csv(program, out_path("ProgramLevel.csv"))

## This file was written as ProgrammeLevel.csv until the spelling was
## standardised. Remove the old one so the output directory cannot hold two
## copies of the same result under different names.
obsolete_csv <- file.path(OUT_ROOT, "ProgrammeLevel.csv")
if (file.exists(obsolete_csv)) {
  file.remove(obsolete_csv)
  message("Removed obsolete ProgrammeLevel.csv (now written as ProgramLevel.csv).")
}

message("\n=== Program-level, abundance-matched (", N_PERM, " permutations) ===")
print(as.data.frame(program %>% dplyr::select(
  panel, contrast, excess_silent, padj_silent, excess_expressed, padj_expressed)))

n_survive <- sum(program$padj_silent < FDR | program$padj_expressed < FDR)
message("\n  ", n_survive, " of ", nrow(program), " panel x contrast comparisons ",
        "reach significance on at least one count,")
message("  after correcting all ", n_tests, " tests in this script together.")

message("\n=== Non-significant comparisons, with the excess each would have detected ===")
negatives <- program %>%
  dplyr::filter(padj_silent >= FDR & padj_expressed >= FDR) %>%
  dplyr::mutate(status = dplyr::case_when(
    abs(excess_silent) >= detectable_silent &
      abs(excess_expressed) >= detectable_expressed ~ "UNRESOLVED (both counts)",
    abs(excess_silent) >= detectable_silent    ~ "UNRESOLVED (silent count)",
    abs(excess_expressed) >= detectable_expressed ~ "UNRESOLVED (expressed count)",
    TRUE ~ "bounded negative")) %>%
  dplyr::select(panel, contrast, excess_silent, detectable_silent,
                excess_expressed, detectable_expressed, status)
print(as.data.frame(negatives))
n_unres <- sum(negatives$status != "bounded negative")
if (n_unres) {
  message("\n  ", n_unres, " of these ", nrow(negatives), " are NOT bounded negatives: the",
          " observed excess reaches the")
  message("  comparison's own detection limit, and the result is non-significant only")
  message("  after the joint correction. Quote those as unresolved, with the excess and")
  message("  the limit both stated, never as `no difference`.")
}
message("  For a bounded negative, read `detectable` as the sensitivity of that comparison: the smallest")
message("  excess it would have flagged. Compare it with the panel size, not")
message("  across panels. A small panel gives a NARROW null, so it detects a")
message("  smaller excess in absolute genes -- roughly 1 gene in 25 here against")
message("  4 in 92 and 5 in 177, which is about 4 to 5 per cent of the panel in")
message("  every case. These negatives are comparably sensitive, and each should")
message("  be quoted with its own limit rather than hedged as underpowered.")

if (any(program$at_permutation_floor)) {
  message("\n  ", sum(program$at_permutation_floor), " comparison(s) sit at the ",
          "permutation floor of ", signif(p_floor, 3), ".")
  message("  Those P values are not distinguishable from one another and must not")
  message("  be ranked. Raise PAX6_N_PERM if the ordering among them matters.")
}
if (sum(program$n_bins_undersized) > 0) {
  message("\n  ", sum(program$n_bins_undersized), " abundance bin(s) held fewer ",
          "control genes than the panel needed, so genes were drawn with")
  message("  replacement there. That correlates the draws and narrows the null,")
  message("  making those P values anti-conservative. Lower PAX6_N_BINS if this")
  message("  is more than a handful.")
}
if (sum(program$n_bins_empty) > 0) {
  message("\n  ", sum(program$n_bins_empty), " abundance bin(s) had no matched ",
          "control gene at all, so the drawn sets are smaller than the panel")
  message("  in those comparisons.")
}

message("\n", strrep("=", 74))
message("05_program_level complete. Outputs in: ", OUT_ROOT)
message("  Report the EXCESS over the null, never the raw difference in counts.")
message("  The raw difference contains whatever these libraries do to detection;")
message("  the excess is what remains after the matched controls absorb it.")
message(strrep("=", 74))
