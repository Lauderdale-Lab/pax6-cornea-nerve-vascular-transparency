###############################################################################
## 02_shared_gene_assignment.R -- primary panel for genes on more than one panel
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
## 02_shared_gene_assignment.R
##
## Works out which curated-panel genes appear on more than one panel, assigns
## each one a single primary panel by the stated precedence rule, and writes the
## assignment out in two forms: a table for the supplementary material, and a
## block of R ready to paste into PANEL_PRIMARY in 00_config.R.
##
## WHY THIS EXISTS. Thirty gene symbols sit on more than one of the three
## curated panels -- twenty-nine of the ninety-two axon-guidance genes are also
## vascular/lymphatic genes. If each panel's adjusted P values were computed
## independently, such a gene would receive two of them for a single test and
## could be reported as significant in one panel's table and not the other.
## Assigning each gene one primary panel and correcting it once removes that,
## at the cost of a choice that has to be made and defended.
##
## THE CHOICE IS CONSEQUENTIAL. The correction families differ in size by
## roughly tenfold -- of the order of fifteen tested genes for myelination
## against a hundred and thirty for vascular/lymphatic -- so the panel a shared
## gene is corrected in changes the Benjamini-Hochberg threshold it faces.
## Assigning a gene to a smaller panel is the more permissive choice. The rule
## must therefore be fixed on curation grounds and recorded before anyone looks
## at which way it moves a call.
##
## THE RULE, from PANEL_PRECEDENCE in 00_config.R: a gene on the axon-guidance
## panel is corrected there; otherwise a gene on the myelination panel is
## corrected there; otherwise vascular/lymphatic. Guidance and myelination are
## narrow, mechanism-specific curations; vascular/lymphatic is a broad annotated
## survey, and a semaphorin that also carries a vascular annotation is a
## guidance gene first.
##
## Two genes are claimed by both of the first two clauses and are resolved by
## the order rather than by the clauses themselves. This script names them
## explicitly in its output so the decision cannot pass unnoticed.
##
## PANEL_PRIMARY_OVERRIDE in 00_config.R takes precedence over the rule for
## individual genes decided on biological grounds. Each override carries its
## reason, and the reason travels into both outputs. An override that moves a
## gene into a LARGER panel subjects it to a stricter threshold and costs power;
## the script reports the direction so that cannot be glossed over.
##
## This script makes no assignment permanent by itself. It proposes; the
## assignment becomes real only when the block it writes is pasted into
## 00_config.R, where a reader can see all thirty at once.
##
## Inputs : 00_config.R, 01_load.R
## Outputs: SharedGeneAssignment.csv, PanelPrimaryAssignment.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(purrr); library(tidyr)
})

if (!exists("panel_tbl")) stop("Source 00_config.R and 01_load.R first.")
if (!all(PANEL_PRECEDENCE %in% names(PANEL_DIRS))) {
  stop("PANEL_PRECEDENCE names a panel that does not exist. Valid panels: ",
       paste(names(PANEL_DIRS), collapse = ", "))
}
if (!setequal(PANEL_PRECEDENCE, names(PANEL_DIRS))) {
  stop("PANEL_PRECEDENCE must list every panel exactly once, so that no gene ",
       "can be left unassigned. Missing: ",
       paste(setdiff(names(PANEL_DIRS), PANEL_PRECEDENCE), collapse = ", "))
}

## ---------------------------------------------------------------------------
## 1. Membership
## ---------------------------------------------------------------------------

membership <- panel_tbl %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    n_panels  = dplyr::n_distinct(panel),
    panels    = paste(sort(unique(panel)), collapse = "; "),
    gene_ids  = paste(sort(unique(gene_id)), collapse = "; "),
    .groups   = "drop")

shared <- membership %>% dplyr::filter(n_panels > 1)

message("Panel membership: ", nrow(membership), " unique symbols across ",
        length(PANEL_DIRS), " panels; ", nrow(shared), " on more than one.")
if (!nrow(shared)) {
  message("No shared genes. PANEL_PRIMARY can stay empty."); return(invisible(NULL))
}

## A symbol mapping to more than one Ensembl id would make "the gene" ambiguous
## and the assignment meaningless. Fail rather than pick one.
ambiguous <- panel_tbl %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(n_ids = dplyr::n_distinct(gene_id), .groups = "drop") %>%
  dplyr::filter(n_ids > 1)
if (nrow(ambiguous)) {
  stop("These symbols map to more than one gene id, so a per-symbol assignment ",
       "is not well defined: ", paste(ambiguous$symbol, collapse = ", "))
}

## ---------------------------------------------------------------------------
## 2. Apply the precedence rule
## ---------------------------------------------------------------------------

## Overrides are validated against membership: naming a panel a gene is not on
## would silently produce an empty correction family for it.
if (nrow(PANEL_PRIMARY_OVERRIDE)) {
  ov <- PANEL_PRIMARY_OVERRIDE
  unknown_panel <- ov$symbol[!ov$panel %in% names(PANEL_DIRS)]
  if (length(unknown_panel)) {
    stop("PANEL_PRIMARY_OVERRIDE names a panel that does not exist, for: ",
         paste(unknown_panel, collapse = ", "))
  }
  not_shared <- setdiff(ov$symbol, shared$symbol)
  if (length(not_shared)) {
    stop("PANEL_PRIMARY_OVERRIDE names gene(s) that are not on more than one ",
         "panel, so there is nothing to decide: ",
         paste(not_shared, collapse = ", "))
  }
  wrong <- purrr::map2_lgl(ov$symbol, ov$panel, function(s, p) {
    !p %in% trimws(strsplit(shared$panels[shared$symbol == s], ";")[[1]])
  })
  if (any(wrong)) {
    stop("PANEL_PRIMARY_OVERRIDE assigns gene(s) to a panel they are not on: ",
         paste(sprintf("%s -> %s", ov$symbol[wrong], ov$panel[wrong]),
               collapse = "; "))
  }
  if (anyDuplicated(ov$symbol)) {
    stop("PANEL_PRIMARY_OVERRIDE lists a gene more than once: ",
         paste(unique(ov$symbol[duplicated(ov$symbol)]), collapse = ", "))
  }
}

override_for <- function(symbol) {
  if (!nrow(PANEL_PRIMARY_OVERRIDE)) return(NA_character_)
  i <- match(symbol, PANEL_PRIMARY_OVERRIDE$symbol)
  if (is.na(i)) NA_character_ else PANEL_PRIMARY_OVERRIDE$panel[i]
}
reason_for <- function(symbol) {
  if (!nrow(PANEL_PRIMARY_OVERRIDE)) return(NA_character_)
  i <- match(symbol, PANEL_PRIMARY_OVERRIDE$symbol)
  if (is.na(i)) NA_character_ else PANEL_PRIMARY_OVERRIDE$reason[i]
}

assign_one <- function(symbol, panels_string) {
  ov <- override_for(symbol)
  if (!is.na(ov)) return(ov)
  in_panels <- trimws(strsplit(panels_string, ";")[[1]])
  PANEL_PRECEDENCE[which(PANEL_PRECEDENCE %in% in_panels)[1]]
}

## A gene is "decided by precedence" when more than one clause of the rule
## claims it -- that is, when it sits on two panels that both rank above the
## rest. Those are the assignments a reader will want justified.
decided_by_order <- function(panels_string) {
  in_panels <- trimws(strsplit(panels_string, ";")[[1]])
  ranks <- which(PANEL_PRECEDENCE %in% in_panels)
  length(ranks) > 1 && ranks[2] <= 2
}

assignment <- shared %>%
  dplyr::mutate(
    primary_panel = purrr::map2_chr(symbol, panels, assign_one),
    override_reason = vapply(symbol, reason_for, character(1)),
    corrected_elsewhere = purrr::map2_chr(panels, primary_panel, function(p, prim) {
      others <- setdiff(trimws(strsplit(p, ";")[[1]]), prim)
      if (length(others)) paste(others, collapse = "; ") else NA_character_
    }),
    resolution = dplyr::case_when(
      !is.na(override_reason)                     ~ "overridden on biological grounds",
      vapply(panels, decided_by_order, logical(1)) ~ "resolved by precedence order",
      TRUE                                         ~ "unambiguous under the rule")) %>%
  dplyr::arrange(primary_panel, dplyr::desc(n_panels), symbol)

readr::write_csv(
  assignment %>% dplyr::select(symbol, gene_ids, n_panels, panels,
                               primary_panel, corrected_elsewhere, resolution,
                               override_reason),
  out_path("SharedGeneAssignment.csv"))

message("\n=== Assignment by primary panel ===")
print(as.data.frame(assignment %>% dplyr::count(primary_panel, resolution)))

overridden <- assignment %>% dplyr::filter(resolution == "overridden on biological grounds")
if (nrow(overridden)) {
  message("\n*** OVERRIDDEN ON BIOLOGICAL GROUNDS, against the precedence rule.")
  ## Loop variable named defensively: this script is sourced into the global
  ## environment, so a bare `i` here would overwrite a caller's `i`.
  for (ov_i in seq_len(nrow(overridden))) {
    would_be <- PANEL_PRECEDENCE[which(PANEL_PRECEDENCE %in%
      trimws(strsplit(overridden$panels[ov_i], ";")[[1]]))[1]]
    direction <- if (PANEL_SIZES[[overridden$primary_panel[ov_i]]] >
                     PANEL_SIZES[[would_be]]) {
      "into a LARGER panel: stricter threshold, costs power"
    } else if (PANEL_SIZES[[overridden$primary_panel[ov_i]]] <
               PANEL_SIZES[[would_be]]) {
      "into a SMALLER panel: more permissive threshold -- justify carefully"
    } else "panels of equal size"
    message("    ", overridden$symbol[ov_i], ": rule said ", would_be,
            ", assigned to ", overridden$primary_panel[ov_i])
    message("      ", direction)
    message("      reason: ", overridden$override_reason[ov_i])
  }
}

by_order <- assignment %>% dplyr::filter(resolution == "resolved by precedence order")
if (nrow(by_order)) {
  message("\n*** DECIDED BY THE ORDER OF PANEL_PRECEDENCE, NOT BY THE RULE ITSELF.")
  message("    These sit on two panels that both rank above the rest, so the")
  message("    assignment follows from the order alone and needs a stated reason")
  message("    in Methods:")
  print(as.data.frame(by_order %>%
    dplyr::select(symbol, panels, primary_panel)))
}

## ---------------------------------------------------------------------------
## 3. Impact on the correction families
## ---------------------------------------------------------------------------

## What each panel gains and loses. A panel that loses genes gets a smaller
## Benjamini-Hochberg family and therefore a more permissive threshold for the
## genes it keeps, so this table belongs in the record even though the effect is
## a consequence of removing duplicate tests rather than a choice in itself.
family_change <- purrr::map_dfr(names(PANEL_DIRS), function(p) {
  on_panel <- panel_tbl$symbol[panel_tbl$panel == p]
  keeps <- sum(assignment$primary_panel[assignment$symbol %in% on_panel] == p) +
           sum(!on_panel %in% assignment$symbol)
  tibble::tibble(panel = p,
                 genes_on_panel = length(on_panel),
                 shared_on_panel = sum(on_panel %in% assignment$symbol),
                 shared_kept = sum(assignment$primary_panel[assignment$symbol %in% on_panel] == p),
                 corrected_in_this_panel = keeps)
})
message("\n=== Correction families after assignment ===")
print(as.data.frame(family_change))
message("  `corrected_in_this_panel` is the panel-level universe. The family for")
message("  any one contrast is its intersection with that contrast's testable set,")
message("  which 04_differential.R reports as n_in_correction_family.")

total_before <- nrow(panel_tbl)
total_after  <- sum(family_change$corrected_in_this_panel)
message("\n  Gene-panel rows before assignment: ", total_before,
        "; distinct genes corrected once: ", total_after,
        " (", total_before - total_after, " duplicate tests removed).")

## ---------------------------------------------------------------------------
## 4. The block to paste into 00_config.R
## ---------------------------------------------------------------------------

width <- max(nchar(assignment$symbol)) + 2
lines <- c(
  "## PanelPrimaryAssignment.R -- GENERATED by 02_shared_gene_assignment.R.",
  "## Do not edit by hand. Change PANEL_PRECEDENCE or PANEL_PRIMARY_OVERRIDE in",
  "## 00_config.R and re-run 02; this file is rewritten from them.",
  paste0("## Written ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "."),
  "",
  "PANEL_PRIMARY <- c(",
  sprintf('  %-*s = "%s"%s%s',
          width, paste0('"', assignment$symbol, '"'),
          assignment$primary_panel,
          c(rep(",", nrow(assignment) - 1), ""),
          ifelse(is.na(assignment$override_reason), "",
                 "   # biological override, see 00_config.R")),
  ")")
## Written beside the scripts, where 00_config.R sources it, and a copy kept
## with the outputs so a result set carries the assignment that produced it.
writeLines(lines, file.path(SCRIPT_DIR, "PanelPrimaryAssignment.R"))
writeLines(lines, out_path("PanelPrimaryAssignment.R"))

## Also set it now, so that a single pass of run_all.R uses this assignment
## rather than whatever 00_config.R loaded before this script regenerated it.
PANEL_PRIMARY <- stats::setNames(assignment$primary_panel, assignment$symbol)

message("\n", strrep("=", 74))
message(nrow(assignment), " assignment(s) written to PanelPrimaryAssignment.R and")
message("applied to this session. 00_config.R sources that file on every")
message("subsequent run, so there is no step to remember and nothing to paste.")
message("")
message("  SharedGeneAssignment.csv    the derivation, for the supplement")
message("  PanelPrimaryAssignment.R    the assignment itself, generated")
message(strrep("=", 74))
