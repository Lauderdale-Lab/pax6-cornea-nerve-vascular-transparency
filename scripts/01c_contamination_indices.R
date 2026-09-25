###############################################################################
## 01c_contamination_indices.R -- off-target tissue in each library (Supp. Table 6)
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
## 01c_contamination_indices.R
##
## Asks of every library how much tissue OTHER than cornea it carries, and
## whether its share of corneal stroma and epithelium matches its peers.
##
## WHY THIS EXISTS. A dissected cornea can bring adjacent tissue with it: iris
## stuck to an opaque cornea, a rim of angle or limbus, a fragment of lens or
## retina. Each carries genes the central cornea does not express, and a
## difference between groups in how much was carried looks, in a count matrix,
## exactly like a difference in expression. The reverse problem is composition:
## a library with less stroma has less of every stromal gene, for reasons that
## have nothing to do with genotype. The independent adult dataset used in 09
## (GSE183742) was found to have both problems. This script asks the same
## questions of this study's own 48 libraries, by a stated rule, so the answer
## for our data is on the record and reproducible rather than asserted.
##
## At dissection, iris was seen adhering to the opaque regions of opaque mutant
## corneas, so pigment signal in those libraries is EXPECTED. The table reports
## it; it does not exclude on it.
##
## WHAT IS COMPUTED, per library, in CPM against the full-matrix library size:
##   an index per marker family (the summed CPM of its genes), for pigment,
##   retina, lens, angle and conjunctiva; single-gene columns for the markers
##   quoted in the supplement; and stromal and epithelial control indices.
##
## WHAT IS FLAGGED, for corneal libraries only:
##   added tissue    a family index far above the same-age reference, carried by
##                   more than one of its genes (rule in 00_config.R)
##   composition     stromal or epithelial controls more than
##                   OFFTARGET_CONTROL_MAX_FOLD below the same-age reference
## Trigeminal libraries get indices but no flags: none of these families is a
## meaningful contaminant of ganglion.
##
## THE REFERENCE is the other corneal libraries of the same age and laboratory,
## excluding opaque mutants and excluding the library being judged. Including
## the library itself would pull the reference toward it and understate its
## own excess; with five libraries per group that is not a small effect.
##
## NO LIBRARY IS EXCLUDED on this basis. Like 01b, this script reports and
## changes nothing downstream. What a flag means for a given gene is argued in
## the text, with this table cited.
##
## The marker table and every threshold live in 00_config.R and are shared with
## the margin check in 09_cross_dataset_replication.R, which also calls the
## functions defined here so GSE183742's libraries are judged by the same rule.
##
## STATISTICAL BASIS
##   Descriptive only. Medians rather than means for the reference, because a
##   single contaminated library in a reference group would otherwise move it.
##
## Inputs : 00_config.R, 01_load.R
## Outputs: Contamination_Indices.csv          every library, every column
##          SuppTable6_Contamination.csv       the published subset
## Leaves in the session for 09: ot_marker_ids, ot_marker_cpm(),
##          ot_flag_library()
##
## Every top-level object here is prefixed ot_ / OT_. The scripts share one
## global environment and a bare name would overwrite an upstream object.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

if (!exists("count_mat")) stop("Source 00_config.R and 01_load.R first.")
if (!exists("OFFTARGET_MARKERS")) {
  stop("OFFTARGET_MARKERS is not defined. This 00_config.R predates ",
       "01c_contamination_indices.R; use the matching version.")
}

## ---------------------------------------------------------------------------
## 1. Marker symbols to count-matrix identifiers
## ---------------------------------------------------------------------------

## Symbols are resolved against the count matrix itself: every row id is mapped
## to its symbol, and a marker takes the FIRST row carrying that symbol. That
## is how the table was first built, so the identifiers are unchanged; what is
## new is that a missing or ambiguous symbol is named instead of stopping on an
## anonymous stopifnot().
ot_used <- OFFTARGET_MARKERS[is.na(OFFTARGET_MARKERS$excluded), , drop = FALSE]

ot_row_symbol <- suppressMessages(AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = rownames(count_mat), column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"))
ot_ids_by_symbol <- split(names(ot_row_symbol), unname(ot_row_symbol))

## A symbol the row-wise mapping does not produce can still be reachable the
## other way round: org.Mm.eg.db's ENSEMBL -> SYMBOL "first" and SYMBOL ->
## ENSEMBL do not always agree for genes with several identifiers (Muc4 is
## the known case). Try that direction before giving up, keeping only ids that
## are actually rows of the count matrix.
ot_missing <- setdiff(ot_used$symbol, names(ot_ids_by_symbol))
if (length(ot_missing)) {
  ot_back <- suppressMessages(AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = ot_missing, column = "ENSEMBL",
    keytype = "SYMBOL", multiVals = "list"))
  for (s in ot_missing) {
    hit <- intersect(ot_back[[s]], rownames(count_mat))
    if (length(hit)) {
      ot_ids_by_symbol[[s]] <- hit
      message("NOTE: ", s, " resolved by symbol -> Ensembl lookup (", hit[1], ").")
    }
  }
  ot_missing <- setdiff(ot_used$symbol, names(ot_ids_by_symbol))
}
## Still missing. A gene that feeds an index or a flag is essential, and its
## absence stops the run. A gene shown only in 09's margin table is dropped
## with a note, which is what 09 did with an unmapped marker before.
if (length(ot_missing)) {
  ot_essential <- intersect(ot_missing, ot_used$symbol[!is.na(ot_used$family)])
  if (length(ot_essential)) {
    stop("Off-target marker symbol(s) used in an index not found among the ",
         "count-matrix genes: ", paste(ot_essential, collapse = ", "),
         ". Correct the symbol in OFFTARGET_MARKERS (00_config.R) or mark the ",
         "gene excluded with a reason.")
  }
  message("NOTE: margin-table marker(s) not found in the count matrix, dropped: ",
          paste(ot_missing, collapse = ", "))
  ot_used <- ot_used[!ot_used$symbol %in% ot_missing, , drop = FALSE]
}
ot_ambiguous <- ot_used$symbol[lengths(ot_ids_by_symbol[ot_used$symbol]) > 1]
if (length(ot_ambiguous)) {
  message("NOTE: symbol(s) carried by more than one gene id; the first in matrix ",
          "order is used: ", paste(ot_ambiguous, collapse = ", "))
}
ot_marker_ids <- vapply(ot_used$symbol, function(s) ot_ids_by_symbol[[s]][1],
                        character(1))

## ---------------------------------------------------------------------------
## 2. The rule, as functions -- 09 calls these for GSE183742
## ---------------------------------------------------------------------------

ot_family_genes <- function(fam) ot_used$symbol[ot_used$family %in% fam]

## Marker CPM, rows = marker symbols, columns = libraries. `lib` must be the
## FULL-matrix library size: a denominator taken from the marker rows alone
## would express each marker as a share of the markers.
ot_marker_cpm <- function(cm, lib, ids) {
  miss <- setdiff(ot_marker_ids, rownames(cm))
  if (length(miss)) {
    stop("Marker gene id(s) absent from this count matrix: ",
         paste(names(ot_marker_ids)[ot_marker_ids %in% miss], collapse = ", "))
  }
  stopifnot(all(ids %in% colnames(cm)), all(ids %in% names(lib)))
  m <- cm[ot_marker_ids, ids, drop = FALSE]
  m <- sweep(m, 2, lib[ids], "/") * 1e6
  rownames(m) <- names(ot_marker_ids)
  m
}

## Indices for a set of libraries, from a marker CPM matrix.
OT_INDEX_FAMILIES <- c("pigment", "retina", "lens", "angle", "conj")
ot_indices <- function(cpm) {
  out <- tibble::tibble(library_id = colnames(cpm))
  for (k in OT_INDEX_FAMILIES) out[[k]] <- colSums(cpm[ot_family_genes(k), , drop = FALSE])
  for (g in OFFTARGET_SINGLE_GENES) out[[g]] <- as.numeric(cpm[g, ])
  out$stromal_index    <- colSums(cpm[ot_family_genes("stromal_control"), , drop = FALSE])
  out$epithelial_index <- colSums(cpm[ot_family_genes("epithelial_control"), , drop = FALSE])
  ## Rounded for the table only; every flag is computed from unrounded CPM.
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], signif, digits = 3)
  out
}

## Judge ONE library against a reference set. `cpm_lib` is a named vector over
## marker symbols; `cpm_ref` a marker x library matrix that must not contain the
## library being judged. Returns one row: the fold for each family and the
## composition ratios, and the flag text.
OT_FLAG_LABEL <- c(pigment = "pigment cells = iris", retina = "retina",
                   lens = "lens", angle = "angle")
ot_flag_library <- function(cpm_lib, cpm_ref) {
  if (ncol(cpm_ref) < 2) {
    return(tibble::tibble(reference_n = ncol(cpm_ref),
                          contamination_flags = "not evaluated (reference n < 2)"))
  }
  row <- tibble::tibble(reference_n = ncol(cpm_ref))
  flags <- character(0)

  for (k in OFFTARGET_FLAG_FAMILIES) {
    gs <- ot_family_genes(k)
    fold <- sum(cpm_lib[gs]) /
      max(stats::median(colSums(cpm_ref[gs, , drop = FALSE])), OFFTARGET_INDEX_FLOOR)
    gene_ref <- apply(cpm_ref[gs, , drop = FALSE], 1, stats::median)
    hits <- gs[cpm_lib[gs] > OFFTARGET_GENE_FOLD * pmax(gene_ref, OFFTARGET_GENE_FLOOR)]
    row[[paste0("fold_", k)]] <- round(fold, 2)
    if (fold > OFFTARGET_FOLD && length(hits) >= min(2L, length(gs))) {
      flags <- c(flags, sprintf("%s %.0fx (%s)", OT_FLAG_LABEL[[k]], fold,
                                paste(hits, collapse = ", ")))
    } else if (fold > OFFTARGET_TRACE_FOLD) {
      flags <- c(flags, sprintf("%s trace %.0fx", k, fold))
    }
  }

  ## Composition: how many fold BELOW the reference each control gene sits.
  ## Judged per gene, as 09 judges the controls, and flagged only when EVERY
  ## control of the layer is beyond tolerance -- the smallest shortfall is the
  ## one reported. Summing the controls first would let one abundant gene hide
  ## the other, and a single gene falling alone is that gene, not the layer.
  for (ctl in c("stromal", "epithelial")) {
    gs <- ot_family_genes(paste0(ctl, "_control"))
    ref_med <- apply(cpm_ref[gs, , drop = FALSE], 1, stats::median)
    shortfall <- min(ref_med / pmax(cpm_lib[gs], 1e-3))
    row[[paste0(ctl, "_shortfall")]] <- round(shortfall, 2)
    if (shortfall > OFFTARGET_CONTROL_MAX_FOLD) {
      flags <- c(flags, sprintf("%s-poor (%s all >= %.0fx below reference; %s)",
                                if (ctl == "stromal") "stroma" else "epithelium",
                                paste(gs, collapse = ", "), shortfall,
                                paste(sprintf("%s %.0f CPM", gs, cpm_lib[gs]),
                                      collapse = ", ")))
    }
  }
  row$contamination_flags <- paste(flags, collapse = "; ")
  row
}

## ---------------------------------------------------------------------------
## 3. This study's libraries
## ---------------------------------------------------------------------------

ot_meta <- meta %>%
  dplyr::filter(laboratory %in% LABORATORY, sample_id %in% colnames(count_mat))
ot_cpm <- ot_marker_cpm(count_mat, lib_sizes, ot_meta$sample_id)

ot_is_cornea <- ot_meta$tissue %in% TISSUE & !is.na(ot_meta$age_group)
ot_ref_ok    <- ot_is_cornea & !(ot_meta$group3 %in% "Sey_O")

ot_flags <- dplyr::bind_rows(lapply(seq_len(nrow(ot_meta)), function(j) {
  id <- ot_meta$sample_id[j]
  if (!ot_is_cornea[j]) {
    return(tibble::tibble(library_id = id, evaluated = FALSE, reference_n = NA_integer_,
                          contamination_flags = NA_character_))
  }
  ref <- ot_meta$sample_id[ot_ref_ok & ot_meta$age_group == ot_meta$age_group[j] &
                           ot_meta$sample_id != id]
  dplyr::bind_cols(tibble::tibble(library_id = id, evaluated = TRUE),
                   ot_flag_library(ot_cpm[, id], ot_cpm[, ref, drop = FALSE]))
}))

ot_table <- ot_meta %>%
  dplyr::transmute(library_id = sample_id, tissue, age_group, genotype, group3, batch) %>%
  dplyr::left_join(ot_indices(ot_cpm), by = "library_id") %>%
  dplyr::left_join(ot_flags, by = "library_id") %>%
  dplyr::arrange(tissue, factor(age_group, levels = AGE_LEVELS), group3, library_id)

readr::write_csv(ot_table, out_path("Contamination_Indices.csv"))

## The published subset. Column names say which gene or family each is.
ot_supp <- ot_table %>%
  dplyr::transmute(library_id, pigment_index = pigment, Rho, Cryaa,
                   angle_Myoc = angle, Kera, Krt12, contamination_flags)
readr::write_csv(ot_supp, out_path("SuppTable6_Contamination.csv"))

## ---------------------------------------------------------------------------
## 4. Report
## ---------------------------------------------------------------------------

message("=== Off-target tissue indices: ", sum(ot_is_cornea), " corneal libraries evaluated, ",
        sum(!ot_is_cornea), " non-corneal reported without flags ===")
message("Reference: same age and laboratory, opaque mutants and the library itself excluded.")

ot_flagged <- ot_table %>% dplyr::filter(evaluated, nzchar(contamination_flags))
if (nrow(ot_flagged)) {
  message("\nFlagged libraries:")
  print(as.data.frame(ot_flagged %>%
    dplyr::select(library_id, age_group, group3, contamination_flags)), row.names = FALSE)
} else {
  message("\nNo corneal library is flagged.")
}

message("\nFlags by group (libraries with any non-trace flag / libraries):")
print(as.data.frame(ot_table %>% dplyr::filter(evaluated) %>%
  dplyr::group_by(age_group, group3) %>%
  dplyr::summarise(n = dplyr::n(),
                   flagged = sum(grepl("[0-9]x \\(|-poor", contamination_flags)),
                   .groups = "drop")), row.names = FALSE)

## The stroma-poor criterion changed when this step moved into the pipeline.
## The old absolute rule is printed beside the new one so the change is
## auditable rather than silent.
ot_old_rule <- ot_table$library_id[ot_table$evaluated & ot_table$Kera < 100]
ot_new_rule <- ot_table$library_id[ot_table$evaluated &
                                   grepl("stroma-poor", ot_table$contamination_flags)]
message("\nStroma-poor, current rule (stromal controls > ", OFFTARGET_CONTROL_MAX_FOLD,
        "x below same-age reference): ",
        if (length(ot_new_rule)) paste(ot_new_rule, collapse = ", ") else "none")
message("Stroma-poor, superseded rule (Kera < 100 CPM, any age)              : ",
        if (length(ot_old_rule)) paste(ot_old_rule, collapse = ", ") else "none")
if (!setequal(ot_old_rule, ot_new_rule)) {
  message("  The two rules disagree. The current one is age-matched; say so if an")
  message("  earlier Supp. Table 6 used the old one.")
}

message("\n", strrep("=", 74))
message("Methods sentence, generated from this run:")
message("  Off-target tissue was assessed per library from summed CPM of marker")
message("  genes for pigmented tissue, retina, lens and iridocorneal angle, relative")
message("  to the median of same-age corneal libraries excluding opaque mutants and")
message("  the library itself; a family was flagged above ", OFFTARGET_FOLD,
        "-fold when carried by")
message("  at least two of its genes. Stromal and epithelial share were judged from")
message("  Kera/Angptl7 and Krt12/Aldh3a1, flagged beyond ", OFFTARGET_CONTROL_MAX_FOLD,
        "-fold below reference.")
ot_n_strong <- sum(grepl("[0-9]x \\(|-poor", ot_flagged$contamination_flags))
message("  ", ot_n_strong, " of ", sum(ot_is_cornea),
        " corneal libraries carried a flag (", nrow(ot_flagged) - ot_n_strong,
        " more at trace level only). No library was excluded.")
message(strrep("=", 74))
