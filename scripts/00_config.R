###############################################################################
## 00_config.R -- every analysis parameter, in one place
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
## 00_config.R
##
## Every parameter of the PAX6 cornea RNA-seq analysis, in one file.
##
## Nothing here computes a result. This is the single place a reader can look
## to see each choice the analysis makes, what it is set to, and why. No later
## script defines a threshold, a path or an exclusion of its own.
##
## Any value may be overridden by an environment variable of the same name
## prefixed PAX6_, which is how sensitivity analyses are run without editing
## code. Every run prints its resolved configuration.
###############################################################################

## --------------------------------------------------------------------------
## Paths
## --------------------------------------------------------------------------

PROJECT_NAME <- Sys.getenv("PAX6_PROJECT_NAME", unset = "Combined_2025_2026")

## Where the analysis scripts live. Used to locate PanelPrimaryAssignment.R.
## Resolution order: PAX6_SCRIPT_DIR if set; otherwise scripts/ under the
## working directory when it exists (a clone of the repository, run from its
## root); otherwise the working directory itself.
SCRIPT_DIR <- local({
  explicit <- Sys.getenv("PAX6_SCRIPT_DIR", unset = "")
  if (nzchar(explicit)) explicit
  else if (file.exists(file.path("scripts", "00_config.R"))) "scripts"
  else getwd()
})
SCRIPT_DIR <- normalizePath(SCRIPT_DIR, mustWork = FALSE)

resolve_project_root <- function() {
  explicit <- Sys.getenv("PAX6_PROJECT_ROOT", unset = NA_character_)
  if (!is.na(explicit) && nzchar(explicit)) {
    return(normalizePath(path.expand(explicit), mustWork = FALSE))
  }
  ## No machine-specific fallbacks: the project root is the working directory
  ## unless PAX6_PROJECT_ROOT says otherwise, so a run always analyses the
  ## copy of the data it was started from.
  normalizePath(getwd(), mustWork = FALSE)
}

PROJECT_ROOT <- resolve_project_root()

## Every output of a run is written under results/. Set PAX6_OUT_SUBDIR to
## write a run to results/<name> instead, so that two runs (a sensitivity
## analysis, say) never overwrite each other and can be compared directly.
OUT_SUBDIR <- Sys.getenv("PAX6_OUT_SUBDIR", unset = "")
OUT_ROOT <- if (nzchar(OUT_SUBDIR)) {
  file.path(PROJECT_ROOT, "results", OUT_SUBDIR)
} else file.path(PROJECT_ROOT, "results")

GENE_PANEL_DIR <- Sys.getenv("PAX6_GENE_PANEL_DIR",
                             unset = file.path(PROJECT_ROOT, "Gene_Panels_and_Reference_Lists"))

FILE_METADATA <- file.path(PROJECT_ROOT, "Metadata_and_Sample_Guide", PROJECT_NAME,
                           paste0(PROJECT_NAME, "_sample_metadata.csv"))
FILE_COUNTS   <- file.path(PROJECT_ROOT, "RNAseq_Quantification_Matrices", PROJECT_NAME,
                           "gene_counts_featureCounts_all48.txt")

## Curated panels. 00b_panel_universes.R builds each panel from its curated
## CSV and writes it as a master list in PANEL_MASTER_DIR; every later script
## reads the master lists, never the CSVs. The master lists are committed with
## the repository, so a rebuilt panel is checked against the published one.
## The expected sizes are asserted at load time: a master list that has
## silently changed size invalidates every denominator in the paper.
PANEL_MASTER_DIR <- file.path(GENE_PANEL_DIR, "master_lists")
PANEL_MASTER_FILES <- c("Axon guidance"      = "AxonGuidance_Master_List.csv",
                        "Myelination"        = "Myelination_Master_List.csv",
                        "Vascular/Lymphatic" = "VascularLymphatic_Master_List.csv")
PANEL_SIZES <- c("Axon guidance" = 92L, "Myelination" = 25L, "Vascular/Lymphatic" = 177L)

## Curated source files, read by 00b_panel_universes.R from GENE_PANEL_DIR. The
## axon-guidance panel is built from TWO files and cannot be reconstructed from
## the curated list alone; see that script's header.
PANEL_SOURCE_FILES <- list(
  "Axon guidance" = c(
    curated  = "Curated_CANONICAL_guidance_genes_with_ENSEMBL_release112_updated.csv",
    evidence = "cornea_axon_guidance_genes.csv"),
  "Myelination" = c(
    curated  = "Curated_Myelination_Genes.csv"),
  "Vascular/Lymphatic" = c(
    curated  = "Curated_Vascular_Lymphatic_Genes_Mouse_annotated.csv"))

## Refuse to write a panel whose size has changed. Every denominator in the
## paper depends on these three numbers.
PANEL_STRICT <- as.logical(Sys.getenv("PAX6_PANEL_STRICT", "TRUE"))

## Optional: a GTF to take gene symbols from. Unset uses org.Mm.eg.db. The
## axon-guidance union is matched by symbol, so the panel-size check above is
## what confirms the chosen source agrees with the curated lists.
GTF_FILE <- Sys.getenv("PAX6_GTF_FILE", unset = "")

## --------------------------------------------------------------------------
## Sample handling
## --------------------------------------------------------------------------

as_list <- function(x) { v <- trimws(strsplit(x, ",")[[1]]); v[nzchar(v)] }

## No library is excluded from either the expression calls or the differential
## models. A detection exclusion of A_CW3 was applied in earlier versions of
## this analysis; it was introduced to protect a detection rule requiring a
## non-zero count in every replicate, which this analysis no longer uses, and
## was removed once a matched-depth survey of all 48 libraries showed A_CW3 to
## be third-worst on that rule's own metric behind two libraries that were
## always retained. See the analytical development record.
DETECTION_EXCLUDE <- as_list(Sys.getenv("PAX6_DETECTION_EXCLUDE", ""))
DE_EXCLUDE        <- as_list(Sys.getenv("PAX6_DE_EXCLUDE", ""))

## Libraries entering the analysis at all.
TISSUE     <- "Cornea"
## Which laboratories' libraries to analyse. Comma-separated for several.
## The filter exists to exclude other groups' samples that share a metadata
## file, not to hard-code one lab: GSE183742 (Duncan) is analysed with
## PAX6_LABORATORY=Duncan, and further datasets (Lagali) will follow.
LABORATORY <- trimws(strsplit(
  Sys.getenv("PAX6_LABORATORY", unset = "Lauderdale"), ",")[[1]])

## A group smaller than this is not called. Three is the minimum at which the
## leave-one-out stability check below has any meaning.
MIN_GROUP_N <- as.integer(Sys.getenv("PAX6_MIN_GROUP_N", "3"))

## --------------------------------------------------------------------------
## Expression status
## --------------------------------------------------------------------------

## Counts-per-million thresholds on the exact Poisson limits, not on the point
## estimate. T_HI is not chosen by optimising any downstream quantity -- doing
## so is circular, because coverage rises whenever the band narrows. It is set
## from sequencing depth: at 40-50 M assigned reads, 1 CPM is about 40 reads.
T_HI <- as.numeric(Sys.getenv("PAX6_T_HI", "1"))    # lower limit must clear this to call ON
T_LO <- as.numeric(Sys.getenv("PAX6_T_LO", "0.2"))  # upper limit must fall below to call OFF

## Confidence level for the Poisson limits (two one-sided 5% tails).
POISSON_ALPHA <- as.numeric(Sys.getenv("PAX6_POISSON_ALPHA", "0.05"))

## Leave-one-out stability. When TRUE, any status that changes when a single
## library is removed is reported as indeterminate. The ungated status is
## written alongside so the cost of the check is always measurable.
STABILITY_GATE <- as.logical(Sys.getenv("PAX6_STABILITY_GATE", "TRUE"))

## --------------------------------------------------------------------------
## Categorical differences (genes on in one group and off in the other)
## --------------------------------------------------------------------------

## A gene confidently expressed in one group and confidently silent in the
## other is a categorical difference that needs no model. A gene confidently
## expressed in one group whose counterpart is indeterminate is NOT
## automatically comparable -- but when the two Poisson intervals do not
## overlap, the data do support a difference, and the smallest fold difference
## compatible with both intervals is a defensible statement.
##
## CATEGORICAL_MIN_FOLD is applied to that conservative bound (the lower limit
## of the higher group divided by the upper limit of the lower group), not to
## the ratio of point estimates. A value of 1 is bare non-overlap; the default
## of 2 requires the data to exclude anything smaller than a doubling.
##
## Requiring one side to be confidently EXPRESSED is essential. Two silent or
## two near-silent groups can also have non-overlapping intervals, and calling
## those categorical would be meaningless.
CATEGORICAL_MIN_FOLD <- as.numeric(Sys.getenv("PAX6_CATEGORICAL_MIN_FOLD", "2"))

## Thresholds reported in the sensitivity table, so the effect of the choice
## above is visible in every run rather than assumed.
CATEGORICAL_FOLD_GRID <- c(1, 1.5, 2, 3, 5)

## --------------------------------------------------------------------------
## Differential expression
## --------------------------------------------------------------------------

FDR <- as.numeric(Sys.getenv("PAX6_FDR", "0.05"))

## Equivalence bound in log2 units for the second, formal test of no change.
## 0.585 is 1.5-fold. A two-fold floor is arguably the wrong expectation for a
## haploinsufficiency model whose primary lesion is itself a halving of dosage.
EQUIV_BOUND <- as.numeric(Sys.getenv("PAX6_EQUIV_BOUND", "0.585"))

## Adjusted P values are computed within panel and contrast, with each gene
## corrected exactly once -- see PANEL_PRIMARY below for how genes on more than
## one panel are handled.

## Libraries are mixed-sex pools, so sex is a continuous composition variable
## and can never be a per-sample factor.
USE_SEX_COVARIATE <- as.logical(Sys.getenv("PAX6_USE_SEX_COVARIATE", "TRUE"))

## Trigeminal ganglion (10_trigeminal.R). Two libraries per group, one batch:
## the models carry no batch term, and by default no male-fraction covariate,
## which would spend most of the remaining residual degrees of freedom. Each
## library's male fraction is written to Trigeminal_SampleTable.csv instead.
TRIGEMINAL_TISSUE <- "Trigeminal"
TRIGEMINAL_SEX_COVARIATE <- as.logical(Sys.getenv("PAX6_TRIGEMINAL_SEX_COVARIATE", "FALSE"))

## Thirty gene symbols appear on more than one curated panel. Correcting each
## panel separately would give such a gene two adjusted P values for one test,
## and it can then be significant in one panel's table and not the other. Each
## shared gene is therefore assigned ONE primary panel and corrected once, in
## that panel's family; the resulting adjusted P value is shown in every panel
## the gene belongs to.
##
## The assignment is consequential and must be made on curation grounds BEFORE
## looking at which way it moves any call. Family sizes differ by an order of
## magnitude -- roughly 15-17 tested genes for myelination against 111-132 for
## vascular/lymphatic -- so the panel a shared gene is corrected in changes the
## Benjamini-Hochberg threshold it faces. Assigning a gene to the smaller panel
## is the more permissive choice, and a reviewer is entitled to ask why each
## gene sits where it does. Record the reason with the assignment.
##
## The assignment rule adopted for this study, in precedence order: a gene on
## the axon-guidance panel is corrected there; otherwise a gene on the
## myelination panel is corrected there; otherwise vascular/lymphatic. The
## rationale is that the guidance and myelination panels are narrow, mechanism-
## specific curations, while the vascular/lymphatic panel is a broad annotated
## survey -- semaphorins, ephrins and netrins are guidance molecules that are
## also annotated in vascular biology, not the other way round.
##
## The order matters for exactly two genes, which sit on the guidance AND
## myelination panels and so are claimed by both clauses: Gas1 and Sema4d.
## Both were curated as guidance genes first and added to the myelination panel
## later, so guidance-first is applied. Say so in Methods rather than leaving it
## to the reader to infer from the order of a vector.
PANEL_PRECEDENCE <- c("Axon guidance", "Myelination", "Vascular/Lymphatic")

## Genes whose primary panel is decided on biological grounds rather than by
## the precedence rule. An override always wins, and always carries a reason,
## because the whole justification for this scheme is that the assignment was
## made from biology rather than from what it does to a P value.
##
## Note the direction of travel for each one. Moving a gene from a small panel
## to a large one subjects it to a stricter Benjamini-Hochberg threshold, so an
## override in that direction costs power and cannot be read as favourable.
PANEL_PRIMARY_OVERRIDE <- data.frame(
  symbol = c("Notch1"),
  panel  = c("Vascular/Lymphatic"),
  reason = c(paste("Core angiogenic receptor: Dll4-Notch1 signalling specifies",
                   "tip and stalk cell identity during sprouting. Its role in",
                   "myelination, as a brake on Schwann cell differentiation, is",
                   "secondary. Moves the gene from the smallest correction",
                   "family to the largest.")),
  stringsAsFactors = FALSE)

## The resulting assignment, gene by gene. It is NOT written here by hand.
## 02_shared_gene_assignment.R derives it and writes PanelPrimaryAssignment.R
## beside this file; that file is sourced below if it exists, and it also lists
## all thirty explicitly so a reader sees every assignment without re-deriving
## it.
##
## Keeping it in a separate file is deliberate. It was previously pasted into
## this one, which meant that editing the configuration for any other reason
## silently discarded the assignment and left the multiplicity frame
## unspecified. A value that is derived should not live in a file that is
## hand-edited.
##
## Declared as a zero-length CHARACTER vector rather than c(), which is NULL and
## silently defeats every length check that depends on it.
PANEL_PRIMARY <- character(0)
PANEL_PRIMARY_FILE <- file.path(SCRIPT_DIR, "PanelPrimaryAssignment.R")
if (file.exists(PANEL_PRIMARY_FILE)) source(PANEL_PRIMARY_FILE, local = FALSE)

## What to do about a shared gene with no entry above.
##   "stop"  refuse to run. The default: an unassigned gene means the
##           multiplicity frame is not fully specified.
##   "rule"  fall back to PANEL_PRECEDENCE below. Deterministic and identical to
##          what 02_shared_gene_assignment.R would have proposed, but it hides
##          the fact that the assignment was never reviewed. Use knowingly.
PANEL_PRIMARY_FALLBACK <- Sys.getenv("PAX6_PANEL_PRIMARY_FALLBACK", "stop")
stopifnot(PANEL_PRIMARY_FALLBACK %in% c("stop", "rule"))

## --------------------------------------------------------------------------
## Library quality control
## --------------------------------------------------------------------------

## Detected-gene count rises with sequencing depth, so raw detection cannot
## separate a shallow library from a low-complexity one. Every library is
## downsampled to a common depth; what remains after depth is removed is
## complexity. Rarefaction is used HERE ONLY and never for a differential test.
RAREFY_DEPTH <- as.numeric(Sys.getenv("PAX6_RAREFY_DEPTH", NA))  # NA = the smallest library
RAREFY_REPS  <- as.integer(Sys.getenv("PAX6_RAREFY_REPS", "3"))

## Each library is judged within its own tissue and age stratum, because
## detection genuinely differs by both. Robust z uses the median and MAD.
DETECT_MAD_CUTOFF  <- as.numeric(Sys.getenv("PAX6_DETECT_MAD_CUTOFF", "3"))
DETECT_SOFT_RATIO  <- as.numeric(Sys.getenv("PAX6_DETECT_SOFT_RATIO", "0.90"))

## --------------------------------------------------------------------------
## Off-target tissue in a library (01c_contamination_indices.R, and the
## dissection-margin check in 09_cross_dataset_replication.R)
## --------------------------------------------------------------------------

## ONE marker table, read by both scripts. They previously kept separate lists,
## and the lists had drifted: the margin check still counted Pmel as a pigment
## marker and Chi3l1 as an angle marker after both had been shown to mislead in
## this tissue. Two scripts must not judge the same tissue with different genes.
##
##   family    the per-library index the gene contributes to in 01c. NA means
##             the gene is shown in the 09 margin table only.
##   margin    TRUE: listed in 09's dissection-margin table.
##   excluded  why a gene that looks like a marker is NOT used. Such genes stay
##             in the table so the exclusion is on the record, not silent.
##
## The conjunctival index is computed and reported but never FLAGGED:
## conjunctivalisation of the Pax6 mutant cornea is a known phenotype of the
## tissue itself, not a dissection artefact, so a high value is biology.
OFFTARGET_MARKERS <- data.frame(
  symbol = c(
    "Krt12", "Aldh3a1", "Kera", "Angptl7",
    "Krt15", "Krt19", "Muc4", "Krt13", "Krt4", "Muc5ac",
    "Myoc", "Mgp", "Chi3l1",
    "Tyr", "Tyrp1", "Dct", "Mlana", "Pmel",
    "Vwf", "Cldn5",
    "Cryaa", "Cryba1", "Crybb2", "Crygd", "Mip", "Bfsp2", "Cryab",
    "Rho", "Gnat1", "Sag", "Pde6b", "Rcvrn", "Crx", "Nrl", "Rpe65",
    "Optc"),
  compartment = c(
    rep("central cornea (control)", 4),
    rep("limbal / conjunctival epithelium", 6),
    rep("trabecular meshwork / angle", 3),
    rep("pigmented tissue: iris, ciliary body, limbal melanocytes", 5),
    rep("vascular endothelium: limbal arcade", 2),
    rep("lens", 7),
    rep("retina / RPE", 8),
    "ciliary body / vitreous"),
  family = c(
    "epithelial_control", "epithelial_control", "stromal_control", "stromal_control",
    NA, NA, NA, "conj", "conj", "conj",
    "angle", NA, NA,
    "pigment", "pigment", "pigment", "pigment", NA,
    NA, NA,
    "lens", "lens", "lens", "lens", "lens", "lens", NA,
    rep("retina", 8),
    NA),
  margin = c(
    rep(TRUE, 12), FALSE,
    rep(TRUE, 4), FALSE,
    TRUE, TRUE,
    rep(TRUE, 6), FALSE,
    rep(TRUE, 8),
    TRUE),
  excluded = c(
    rep(NA, 12),
    "induced with inflammation in opaque Sey cornea; would read opacity as an angle difference",
    rep(NA, 4),
    "12-55 CPM in every corneal library; not specific to adherent pigmented tissue",
    rep(NA, 8),
    "expressed by the cornea itself",
    rep(NA, 9)),
  stringsAsFactors = FALSE)

## Genes also reported individually, one column each, in the per-library table.
## (Myoc is not repeated here: it is the whole of the angle index.)
OFFTARGET_SINGLE_GENES <- c("Optc", "Rho", "Cryaa", "Kera", "Krt12")

## The indices that can raise a flag, and the rule. A library is flagged for a
## family when its summed index exceeds OFFTARGET_FOLD times the median of its
## reference libraries AND at least two of the family's genes (one, for a
## one-gene family) each exceed OFFTARGET_GENE_FOLD times their own reference
## median. The second clause stops a single noisy gene carrying a family.
## Between OFFTARGET_TRACE_FOLD and OFFTARGET_FOLD it is reported as trace.
## The floors keep a reference median of zero from making any count a flag.
##
## The reference is the other corneal libraries of the same age, same
## laboratory, EXCLUDING opaque mutants (where adherent iris is expected and
## would raise the bar) and EXCLUDING the library being judged.
OFFTARGET_FLAG_FAMILIES <- c("pigment", "retina", "lens", "angle")
OFFTARGET_FOLD       <- as.numeric(Sys.getenv("PAX6_OFFTARGET_FOLD", "10"))
OFFTARGET_TRACE_FOLD <- as.numeric(Sys.getenv("PAX6_OFFTARGET_TRACE_FOLD", "5"))
OFFTARGET_GENE_FOLD  <- as.numeric(Sys.getenv("PAX6_OFFTARGET_GENE_FOLD", "5"))
OFFTARGET_INDEX_FLOOR <- 0.5   # CPM, floor on a reference index median
OFFTARGET_GENE_FLOOR  <- 0.2   # CPM, floor on a reference gene median

## Tissue COMPOSITION, as opposed to added tissue. A library whose stromal (or
## epithelial) control markers ALL sit more than this many fold below its
## reference took a different share of that layer. This is the ONE criterion for
## "stroma-poor" in the pipeline: 01c applies it per library against same-age
## references, and 09 applies it between the two wild types. It replaces an
## absolute rule (Kera below 100 CPM) that ignored age.
## PAX6_X_CONTROL_MAX_FOLD is the name 09 used before this was shared.
OFFTARGET_CONTROL_MAX_FOLD <- as.numeric(Sys.getenv(
  "PAX6_OFFTARGET_CONTROL_MAX_FOLD",
  Sys.getenv("PAX6_X_CONTROL_MAX_FOLD", "4")))

## Sample PCA in 07_concordance.R: a descriptive QC figure only. The number of
## most-variable genes it uses (500 is the DESeq2 plotPCA convention).
PCA_NTOP <- as.integer(Sys.getenv("PAX6_PCA_NTOP", "500"))

## --------------------------------------------------------------------------
## Program-level test
## --------------------------------------------------------------------------

N_PERM <- as.integer(Sys.getenv("PAX6_N_PERM", "20000"))
N_BINS <- as.integer(Sys.getenv("PAX6_N_BINS", "20"))
SEED   <- as.integer(Sys.getenv("PAX6_SEED", "123"))

## --------------------------------------------------------------------------
## Design
## --------------------------------------------------------------------------

AGE_GROUP_MAP <- c(P3 = "P3_P4", P4 = "P3_P4", P15 = "P15", Adult = "Adult")
AGE_LEVELS    <- c("P3_P4", "P15", "Adult")

## Two models, one per biological regime. The developing cornea and the adult
## cornea undergoing opacification have different variance structures; a pooled
## fit makes the adult contrasts pay for neonatal variability. Within a model
## every contrast shares a dispersion estimate and gene counts are comparable;
## across models they are not, and are not compared.
##
## Adult transparent mutant versus wild type is estimable under both models. It
## is fitted in both and reported from Model B, which estimates its dispersion
## from adult libraries alone. Model A's version is a consistency check.
CONTRASTS <- data.frame(
  model    = c("A",               "A",             "A",
               "B",                "B",                "B"),
  axis     = c("age",             "age",           "age",
               "opacity",          "opacity",          "opacity"),
  label    = c("Sey_vs_WT_P3_P4", "Sey_vs_WT_P15", "Sey_vs_WT_Adult",
               "SeyT_vs_WT_Adult", "SeyO_vs_WT_Adult", "SeyO_vs_SeyT_Adult"),
  cell_a   = c("P3_P4:Sey",       "P15:Sey",       "Adult:Sey_T",
               "Adult:Sey_T",      "Adult:Sey_O",      "Adult:Sey_O"),
  cell_b   = c("P3_P4:WT",        "P15:WT",        "Adult:WT",
               "Adult:WT",         "Adult:WT",         "Adult:Sey_T"),
  level_a  = c("P3_P4_Sey",       "P15_Sey",       "Adult_Sey_T",
               "Adult_Sey_T",      "Adult_Sey_O",      "Adult_Sey_O"),
  level_b  = c("P3_P4_WT",        "P15_WT",        "Adult_WT",
               "Adult_WT",         "Adult_WT",         "Adult_Sey_T"),
  stringsAsFactors = FALSE)

## The one comparison fitted twice. Report the Model B label.
HINGE <- c(A = "Sey_vs_WT_Adult", B = "SeyT_vs_WT_Adult")

## --------------------------------------------------------------------------

config_report <- function() {
  message(strrep("=", 74))
  message("PAX6 cornea RNA-seq -- resolved configuration")
  message(strrep("=", 74))
  message("  project root         : ", PROJECT_ROOT)
  message("  outputs              : ", OUT_ROOT)
  message("  detection exclusions : ",
          if (length(DETECTION_EXCLUDE)) paste(DETECTION_EXCLUDE, collapse = ", ") else "(none)")
  message("  DE exclusions        : ",
          if (length(DE_EXCLUDE)) paste(DE_EXCLUDE, collapse = ", ") else "(none)")
  message("  expression band      : ON at lower limit >= ", T_HI,
          " CPM, OFF at upper limit < ", T_LO, " CPM")
  message("  stability gate       : ", STABILITY_GATE)
  message("  categorical min fold : ", CATEGORICAL_MIN_FOLD, " (conservative bound)")
  message("  FDR                  : ", FDR,
          "   equivalence bound +/-", EQUIV_BOUND, " log2")
  message("  multiplicity frame   : within panel x contrast, each gene corrected once")
  message("  shared-gene rule     : ", paste(PANEL_PRECEDENCE, collapse = " > "),
          if (nrow(PANEL_PRIMARY_OVERRIDE))
            paste0("  (", nrow(PANEL_PRIMARY_OVERRIDE), " biological override(s): ",
                   paste(PANEL_PRIMARY_OVERRIDE$symbol, collapse = ", "), ")") else "")
  message("  shared-gene fallback : ", PANEL_PRIMARY_FALLBACK,
          " (", length(PANEL_PRIMARY), " assignment(s) loaded from ",
          if (file.exists(PANEL_PRIMARY_FILE)) basename(PANEL_PRIMARY_FILE)
          else "nowhere yet -- 02 will write it", ")")
  message("  sex covariate        : ", USE_SEX_COVARIATE)
  message("  QC matched depth     : ",
          if (is.na(RAREFY_DEPTH)) "smallest library" else format(RAREFY_DEPTH, big.mark = ","),
          "   ", RAREFY_REPS, " replicate(s)")
  message("  off-target flag      : index > ", OFFTARGET_FOLD, "x same-age reference (trace > ",
          OFFTARGET_TRACE_FOLD, "x); composition: controls > ",
          OFFTARGET_CONTROL_MAX_FOLD, "x below reference")
  message("  permutations / bins  : ", N_PERM, " / ", N_BINS, "   seed ", SEED)
  message(strrep("=", 74))
}

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
out_path <- function(...) file.path(OUT_ROOT, ...)
