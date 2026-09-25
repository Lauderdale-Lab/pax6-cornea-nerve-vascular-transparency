###############################################################################
## 01_load.R -- metadata, counts, pool sex composition, panels
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
## 01_load.R
##
## The shared input layer. Reads sample metadata, the count matrix and the
## curated panel universes, and derives the pooled-library sex composition.
## Sourced by every analysis script so that all of them see identical inputs.
##
## Defines: meta, count_mat, lib_sizes, panel_tbl, panel_genes
##
## STATISTICAL BASIS
##   [1] Counts from featureCounts. Liao Y, Smyth GK, Shi W (2014)
##       Bioinformatics 30:923-930.
##   [2] Libraries are pools of animals of both sexes, confirmed by libraries
##       carrying both Y-chromosome reads and Xist above background. Sex is
##       therefore a continuous composition variable, not a per-sample factor,
##       and is carried as the pool's male read fraction.
###############################################################################

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(stringr)
  library(janitor); library(data.table); library(purrr)
})

if (!exists("PROJECT_ROOT")) stop("Source 00_config.R before 01_load.R.")

## --------------------------------------------------------------------------
## Sample metadata
## --------------------------------------------------------------------------

if (!file.exists(FILE_METADATA)) stop("Metadata not found: ", FILE_METADATA)
meta_raw <- readr::read_csv(FILE_METADATA, show_col_types = FALSE) %>% janitor::clean_names()

id_col <- intersect(c("sample_id", "sampleid", "display_label"), names(meta_raw))[1]
if (is.na(id_col)) {
  stop("No sample identifier column in the metadata. Columns present: ",
       paste(names(meta_raw), collapse = ", "))
}
## Both sample pools filter on laboratory. If the column is absent the
## comparison yields NA for every row and silently discards the whole dataset,
## so its absence is fatal rather than defaulted.
if (!"laboratory" %in% names(meta_raw)) {
  stop("Metadata has no `laboratory` column, which both sample pools filter on. ",
       "Columns present: ", paste(names(meta_raw), collapse = ", "))
}

blank_to_empty <- function(x) dplyr::coalesce(as.character(x), "")

meta <- meta_raw %>%
  dplyr::mutate(
    sample_id = as.character(.data[[id_col]]),
    genotype = dplyr::case_when(
      stringr::str_to_upper(blank_to_empty(genotype)) %in% c("WT", "WILD-TYPE", "WILDTYPE") ~ "WT",
      stringr::str_detect(stringr::str_to_upper(blank_to_empty(genotype)), "SEY") ~ "Sey",
      TRUE ~ as.character(genotype)),
    tissue = dplyr::case_when(
      stringr::str_starts(stringr::str_to_upper(blank_to_empty(tissue)), "CORNEA") ~ "Cornea",
      stringr::str_starts(stringr::str_to_upper(blank_to_empty(tissue)), "TRIGEMINAL") ~ "Trigeminal",
      TRUE ~ as.character(tissue)),
    transparency = dplyr::case_when(
      stringr::str_to_upper(blank_to_empty(transparency)) %in% c("T", "TRANSPARENT", "CLEAR") ~ "T",
      stringr::str_to_upper(blank_to_empty(transparency)) %in% c("O", "OPAQUE", "SCAR") ~ "O",
      TRUE ~ NA_character_),
    age_group = unname(AGE_GROUP_MAP[age]),
    batch = factor(as.character(batch)),
    ## The three-level adult grouping used by the opacity model.
    group3 = dplyr::case_when(
      genotype == "WT" ~ "WT",
      genotype == "Sey" & transparency == "T" ~ "Sey_T",
      genotype == "Sey" & transparency == "O" ~ "Sey_O",
      TRUE ~ NA_character_)) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

unmapped_age <- meta$sample_id[is.na(meta$age_group)]
if (length(unmapped_age)) {
  message("NOTE: ", length(unmapped_age), " sample(s) have an age outside ",
          paste(names(AGE_GROUP_MAP), collapse = "/"), " and will not enter any ",
          "age group: ", paste(unmapped_age, collapse = ", "))
}
undefined_group <- meta$sample_id[meta$tissue == TISSUE & meta$genotype == "Sey" &
                                    is.na(meta$transparency)]
if (length(undefined_group)) {
  message("NOTE: ", length(undefined_group), " mutant corneal sample(s) have no ",
          "transparency call and cannot be assigned to a group: ",
          paste(undefined_group, collapse = ", "))
}

## --------------------------------------------------------------------------
## Count matrix
## --------------------------------------------------------------------------

if (!file.exists(FILE_COUNTS)) stop("Count matrix not found: ", FILE_COUNTS)

## featureCounts writes the aligner's file paths as column names.
strip_alignment_suffix <- function(x) {
  x <- basename(x)
  for (pat in c("_hisat2(\\.sorted)?\\.bam$", "_hisat2$",
                "\\.Aligned\\.sortedByCoord\\.out$", "\\.Aligned\\.out$",
                "\\.sorted\\.bam$", "\\.(bam|sam|cram)$")) {
    x <- base::sub(pat, "", x, ignore.case = TRUE)
  }
  trimws(x)
}
match_sample_id <- function(x, valid_ids) {
  if (x %in% valid_ids) return(x)
  hits <- valid_ids[startsWith(x, paste0(valid_ids, "_"))]
  if (length(hits) == 1) hits else x
}

annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
header <- data.table::fread(FILE_COUNTS, nrows = 0, header = TRUE, sep = "\t")
sample_cols <- setdiff(names(header), annotation_cols)

counts_dt <- data.table::fread(FILE_COUNTS, sep = "\t", comment.char = "#",
                               select = c("Geneid", sample_cols))
data.table::setnames(counts_dt, "Geneid", "gene_id")
data.table::setnames(counts_dt, sample_cols, strip_alignment_suffix(sample_cols))
counts_dt[, gene_id := base::sub("\\.[0-9]+$", "", gene_id)]
counts_dt <- unique(counts_dt, by = "gene_id")

count_mat <- as.matrix(counts_dt[, -1])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- counts_dt[["gene_id"]]
colnames(count_mat) <- vapply(colnames(count_mat), match_sample_id, character(1),
                              valid_ids = unique(meta$sample_id))
rm(counts_dt); invisible(gc())

unmatched <- setdiff(colnames(count_mat), meta$sample_id)
if (length(unmatched)) {
  stop("Count matrix columns with no metadata row: ", paste(unmatched, collapse = ", "))
}

## Genes with no counts anywhere are KEPT. A gene never detected is a valid
## "not expressed" result and must stay in the panel denominator. Such genes
## are excluded only from the control universe of the program-level test,
## where they would carry no information.
meta <- meta %>% dplyr::filter(sample_id %in% colnames(count_mat))
lib_sizes <- colSums(count_mat)

## --------------------------------------------------------------------------
## Pooled-library sex composition [2]
## --------------------------------------------------------------------------

## Ddx3y, Eif2s3y, Uty, Kdm5d.
Y_GENE_IDS <- c("ENSMUSG00000069045", "ENSMUSG00000069049",
                "ENSMUSG00000068457", "ENSMUSG00000056673")
y_present <- intersect(Y_GENE_IDS, rownames(count_mat))
if (!length(y_present)) stop("No Y-chromosome marker genes found in the count matrix.")
y_cpm <- colSums(t(t(count_mat[y_present, , drop = FALSE]) / lib_sizes) * 1e6)

## Scaled within tissue and age group, because Y-gene abundance differs by
## both. The result is a relative composition index on [0, 1], not an absolute
## male proportion, and is used only as a covariate.
meta <- meta %>%
  dplyr::mutate(y_cpm = as.numeric(y_cpm[sample_id])) %>%
  dplyr::group_by(tissue, age_group) %>%
  dplyr::mutate(male_frac = if (max(y_cpm, na.rm = TRUE) > 0)
                              y_cpm / max(y_cpm, na.rm = TRUE) else 0) %>%
  dplyr::ungroup()

## --------------------------------------------------------------------------
## Curated panels
## --------------------------------------------------------------------------

missing_panels <- vapply(PANEL_DIRS, function(d)
  !file.exists(file.path(PAPER_ROOT, d, "GenesOfInterest_Master_List.csv")), logical(1))
if (any(missing_panels)) {
  stop("Panel master list(s) not found:\n  ",
       paste(file.path(PAPER_ROOT, PANEL_DIRS[missing_panels],
                       "GenesOfInterest_Master_List.csv"), collapse = "\n  "),
       "\nRun 02_panel_universes.R first. On a cloud-synced volume the file may ",
       "be present but evicted; download it before re-running.")
}

panel_tbl <- purrr::imap(PANEL_DIRS, function(subdir, label) {
  d <- readr::read_csv(file.path(PAPER_ROOT, subdir, "GenesOfInterest_Master_List.csv"),
                       show_col_types = FALSE) %>% janitor::clean_names()
  tibble::tibble(panel = label,
                 gene_id = base::sub("\\.[0-9]+$", "", d$gene_id),
                 symbol = d$symbol) %>%
    dplyr::filter(gene_id %in% rownames(count_mat)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
}) %>% dplyr::bind_rows()

for (nm in names(PANEL_SIZES)) {
  got <- sum(panel_tbl$panel == nm)
  if (got != PANEL_SIZES[[nm]]) {
    stop(nm, ": master list holds ", got, " genes, expected ", PANEL_SIZES[[nm]],
         ". Every denominator in the analysis depends on this. Re-run ",
         "02_panel_universes.R.")
  }
}

panel_genes <- unique(panel_tbl$gene_id)

message("Loaded ", ncol(count_mat), " libraries, ", nrow(count_mat), " genes.")
message("Panels: ", paste(sprintf("%s %d", names(PANEL_SIZES),
                                  as.integer(PANEL_SIZES[names(PANEL_SIZES)])),
                          collapse = " | "),
        "  (", length(panel_genes), " unique genes across panels)")
