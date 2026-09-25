###############################################################################
## 00b_panel_universes.R -- builds the three curated gene panels
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
## 00b_panel_universes.R
##
## Writes a master list for each of the three curated panels.
## These three files define the gene universes for the entire analysis, and
## every denominator in the paper is one of their sizes. Nothing downstream can
## run until they exist, which is why this step comes before the data are loaded.
##
## THE THREE PANELS ARE BUILT DIFFERENTLY, and the differences are not
## incidental. They are reproduced here exactly as originally constructed.
##
##   AXON GUIDANCE   a FULL OUTER JOIN of the curated list, which carries mouse
##                   Ensembl IDs, with a cornea-evidence list whose entries
##                   carry HUMAN Ensembl IDs and are matched to mouse IDs by
##                   UPPERCASED SYMBOL against a lookup built over the whole
##                   counts matrix. The panel is therefore a union of two
##                   sources and CANNOT be rebuilt from the curated file alone.
##                   known_cornea_expressed means present in the evidence list.
##
##   MYELINATION     the curated list, deduplicated on gene ID, first row kept.
##                   known_cornea_expressed requires the evidence field to match
##                   a citation pattern, NOT merely to be non-empty: most rows
##                   record where the identifier came from (GO, Reactome, MGI)
##                   rather than any corneal observation, and treating those as
##                   evidence would overstate the benchmark.
##
##   VASCULAR        the curated list, deduplicated on gene ID keeping the
##                   FIRST-LISTED category as primary, with the full
##                   multi-category membership preserved in all_categories_raw.
##                   known_cornea_expressed requires a non-empty evidence field.
##
## SAFETY. Nothing is written unless every panel reproduces its expected size.
## If a master list already exists it is compared gene by gene first and any
## difference is reported. A panel that has silently changed size would move
## every count in the paper without anything else failing.
##
## SYMBOL LOOKUP. Set GTF_FILE in 00_config.R to take symbols from a GTF;
## otherwise org.Mm.eg.db is used. The axon-guidance union depends on this
## lookup, so the panel-size check is what confirms the chosen source agrees
## with the curated lists.
##
## Inputs : 00_config.R, and the curated CSVs named in PANEL_SOURCE_FILES
## Outputs: <PANEL_MASTER_DIR>/<panel>_Master_List.csv  (x3)
##          <OUT_ROOT>/PanelUniverses_BuildReport.txt
###############################################################################

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(stringr)
  library(janitor); library(data.table); library(rlang)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

if (!exists("PANEL_SOURCE_FILES")) stop("Source 00_config.R first.")

log_lines <- character(0)
say <- function(...) { t <- paste0(...); message(t); log_lines <<- c(log_lines, t) }

say("=== Panel universes ===")
say("Gene panel directory : ", GENE_PANEL_DIR)
say("Written to           : ", PANEL_MASTER_DIR)

strip_version <- function(x) base::sub("\\.[0-9]+$", "", as.character(x))

panel_source <- function(panel, which) {
  nm <- PANEL_SOURCE_FILES[[panel]][[which]]
  p <- file.path(GENE_PANEL_DIR, nm)
  if (!file.exists(p)) {
    stop(nm, " not found in ", GENE_PANEL_DIR, "\nCSVs present:\n  ",
         paste(list.files(GENE_PANEL_DIR, pattern = "\\.csv$"), collapse = "\n  "))
  }
  p
}

## Curated files were prepared at different times and their column names differ.
## Resolution is by candidate list rather than by position, so a reordered or
## renamed column fails loudly instead of silently supplying the wrong field.
find_col <- function(df, variants, what, required = TRUE) {
  hit <- intersect(variants, names(df))
  if (!length(hit)) {
    if (!required) return(NA_character_)
    stop("No ", what, " column. Columns present: ", paste(names(df), collapse = ", "))
  }
  hit[1]
}
V_ID  <- c("ensembl_id", "gene_id", "ensembl_gene_id", "geneid", "ensembl")
V_SYM <- c("mouse_symbol", "symbol", "gene_symbol", "gene", "mgi_symbol")
V_CAT <- c("category", "gene_family", "family", "class")
V_SUB <- c("subcategory_or_pathway", "subcategory", "pathway", "subclass")
V_EVI <- c("evidence_source", "evidence", "source", "publication_source")

## ---------------------------------------------------------------------------
## Symbol lookup over the counts matrix
## ---------------------------------------------------------------------------

if (!file.exists(FILE_COUNTS)) stop("Count matrix not found: ", FILE_COUNTS)
gene_ids_all <- unique(strip_version(
  data.table::fread(FILE_COUNTS, sep = "\t", comment.char = "#",
                    select = "Geneid")[["Geneid"]]))
say("Counts matrix gene universe: ", length(gene_ids_all), " genes")

if (nzchar(GTF_FILE) && file.exists(GTF_FILE)) {
  say("Symbol lookup source: GTF (", basename(GTF_FILE), ")")
  gtf <- data.table::fread(GTF_FILE, sep = "\t", header = FALSE, showProgress = FALSE)
  attr9 <- gtf[[9]][gtf[[3]] == "gene"]
  symbol_lookup <- tibble::tibble(
    gene_id = strip_version(base::sub('.*gene_id "([^"]+)".*', "\\1", attr9)),
    symbol  = base::sub('.*gene_name "([^"]+)".*', "\\1", attr9)) %>%
    dplyr::filter(!is.na(symbol), nzchar(symbol))
} else {
  say("Symbol lookup source: org.Mm.eg.db")
  sym <- suppressMessages(AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = gene_ids_all, column = "SYMBOL",
    keytype = "ENSEMBL", multiVals = "first"))
  symbol_lookup <- tibble::tibble(gene_id = names(sym), symbol = unname(sym)) %>%
    dplyr::filter(!is.na(symbol), nzchar(symbol))
}
symbol_lookup <- symbol_lookup %>%
  dplyr::mutate(symbol_upper = toupper(symbol)) %>%
  dplyr::distinct(symbol_upper, .keep_all = TRUE)
say("Symbol lookup entries: ", nrow(symbol_lookup))

## ---------------------------------------------------------------------------
## Axon guidance -- curated list UNION cornea-evidence list
## ---------------------------------------------------------------------------

build_axon <- function() {
  f_curated  <- panel_source("Axon guidance", "curated")
  f_evidence <- panel_source("Axon guidance", "evidence")
  say("  curated : ", basename(f_curated))
  say("  evidence: ", basename(f_evidence))

  curated_raw  <- readr::read_csv(f_curated, show_col_types = FALSE) %>% janitor::clean_names()
  evidence_raw <- readr::read_csv(f_evidence, show_col_types = FALSE) %>% janitor::clean_names()

  curated <- curated_raw %>%
    dplyr::transmute(
      gene_id         = strip_version(.data[[find_col(curated_raw, V_ID, "gene id")]]),
      symbol_curated  = .data[[find_col(curated_raw, V_SYM, "symbol")]],
      category        = .data[[find_col(curated_raw, V_CAT, "category")]],
      subcategory     = .data[[find_col(curated_raw, V_SUB, "subcategory")]],
      evidence_source = .data[[find_col(curated_raw, V_EVI, "evidence", FALSE)]]) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)

  ## The evidence list carries human identifiers, so the join to mouse genes is
  ## by symbol. This is the step that makes the panel a union.
  evidence <- evidence_raw %>%
    dplyr::transmute(
      symbol_evidence  = .data[[find_col(evidence_raw, c("gene_symbol", V_SYM), "symbol")]],
      human_ensembl_id = .data[[find_col(evidence_raw, c("ensembl_gene_id", V_ID), "human id")]],
      cell_type        = .data[[find_col(evidence_raw, c("corneal_cell_type", "cell_type"),
                                         "cell type", FALSE)]],
      evidence         = .data[[find_col(evidence_raw, c("evidence"), "evidence", FALSE)]],
      source           = .data[[find_col(evidence_raw, c("publication_source", "source"),
                                         "source", FALSE)]]) %>%
    dplyr::mutate(symbol_upper = toupper(symbol_evidence)) %>%
    dplyr::left_join(symbol_lookup, by = "symbol_upper") %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)

  unmatched <- evidence$symbol_evidence[is.na(evidence$gene_id)]
  if (length(unmatched)) {
    say("  NOTE: ", length(unmatched), " evidence-list gene(s) unmatched by symbol: ",
        paste(unmatched, collapse = ", "))
  }

  curated %>%
    dplyr::full_join(
      evidence %>% dplyr::select(gene_id, symbol_evidence, human_ensembl_id,
                                 cell_type, evidence, source),
      by = "gene_id") %>%
    dplyr::mutate(symbol = dplyr::coalesce(symbol_curated, symbol_evidence),
                  known_cornea_expressed = !is.na(symbol_evidence)) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::arrange(category, symbol)
}

## ---------------------------------------------------------------------------
## Myelination -- benchmark requires a citation, not a non-empty field
## ---------------------------------------------------------------------------

build_myelination <- function() {
  f <- panel_source("Myelination", "curated")
  say("  curated : ", basename(f))
  raw <- readr::read_csv(f, show_col_types = FALSE) %>% janitor::clean_names()
  ev <- find_col(raw, V_EVI, "evidence", FALSE)

  out <- raw %>%
    dplyr::rename(gene_id     = !!rlang::sym(find_col(raw, V_ID,  "gene id")),
                  symbol      = !!rlang::sym(find_col(raw, V_SYM, "symbol")),
                  category    = !!rlang::sym(find_col(raw, V_CAT, "category")),
                  subcategory = !!rlang::sym(find_col(raw, V_SUB, "subcategory")))
  if (!is.na(ev) && ev %in% names(out)) {
    out <- dplyr::rename(out, evidence_source = !!rlang::sym(ev))
  } else out$evidence_source <- NA_character_

  n_before <- nrow(out)
  out <- out %>%
    dplyr::mutate(gene_id = strip_version(gene_id)) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::mutate(known_cornea_expressed = !is.na(evidence_source) &
                    stringr::str_detect(evidence_source,
                                        stringr::regex("et al", ignore_case = TRUE))) %>%
    dplyr::arrange(category, symbol)
  if (n_before != nrow(out)) {
    say("  NOTE: ", n_before - nrow(out), " duplicate gene id row(s) collapsed")
  }
  out
}

## ---------------------------------------------------------------------------
## Vascular / lymphatic -- first-listed category primary, membership preserved
## ---------------------------------------------------------------------------

build_vascular <- function() {
  f <- panel_source("Vascular/Lymphatic", "curated")
  say("  curated : ", basename(f))
  raw <- readr::read_csv(f, show_col_types = FALSE) %>% janitor::clean_names()
  ev <- find_col(raw, V_EVI, "evidence", FALSE)

  out <- raw %>%
    dplyr::rename(gene_id     = !!rlang::sym(find_col(raw, V_ID,  "gene id")),
                  symbol      = !!rlang::sym(find_col(raw, V_SYM, "symbol")),
                  category    = !!rlang::sym(find_col(raw, V_CAT, "category")),
                  subcategory = !!rlang::sym(find_col(raw, V_SUB, "subcategory")))
  if (!is.na(ev) && ev %in% names(out)) {
    out <- dplyr::rename(out, evidence_source = !!rlang::sym(ev))
  } else out$evidence_source <- NA_character_

  out <- out %>%
    dplyr::mutate(gene_id = strip_version(gene_id)) %>%
    dplyr::filter(!is.na(gene_id), nzchar(gene_id))
  say("  ", nrow(out), " rows, ", dplyr::n_distinct(out$gene_id), " unique genes (",
      nrow(out) - dplyr::n_distinct(out$gene_id), " multi-category duplicates collapsed)")

  all_categories <- out %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(all_categories_raw = paste(
      unique(paste0(category, " (", subcategory, ")")), collapse = "; "),
      .groups = "drop")

  out %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::left_join(all_categories, by = "gene_id") %>%
    dplyr::mutate(known_cornea_expressed = !is.na(evidence_source) &
                    nzchar(trimws(dplyr::coalesce(evidence_source, "")))) %>%
    dplyr::arrange(category, symbol)
}

## ---------------------------------------------------------------------------
## Build, verify, write
## ---------------------------------------------------------------------------

BUILDERS <- list("Axon guidance" = build_axon, "Myelination" = build_myelination,
                 "Vascular/Lymphatic" = build_vascular)

built <- list(); problems <- character(0)

for (label in names(BUILDERS)) {
  say(""); say("--- ", label)
  genes <- BUILDERS[[label]]()
  expected <- PANEL_SIZES[[label]]
  say("  built ", nrow(genes), " genes (expected ", expected, ")")
  if (nrow(genes) != expected) {
    problems <- c(problems, sprintf("%s: built %d, expected %d",
                                    label, nrow(genes), expected))
  }

  out_dir <- PANEL_MASTER_DIR
  out_csv <- file.path(out_dir, PANEL_MASTER_FILES[[label]])
  if (file.exists(out_csv)) {
    existing <- readr::read_csv(out_csv, show_col_types = FALSE) %>% janitor::clean_names()
    existing_ids <- strip_version(existing$gene_id)
    added <- setdiff(genes$gene_id, existing_ids)
    lost  <- setdiff(existing_ids, genes$gene_id)
    if (!length(added) && !length(lost)) {
      say("  identical to the existing master list")
    } else {
      say("  DIFFERS from the existing master list: ", length(added), " added, ",
          length(lost), " lost")
      if (length(added)) say("    added: ", paste(utils::head(added, 12), collapse = ", "))
      if (length(lost))  say("    lost : ", paste(utils::head(lost, 12), collapse = ", "))
      problems <- c(problems, sprintf("%s: gene set differs from the existing master list",
                                      label))
    }
  } else say("  no existing master list to compare against")

  built[[label]] <- list(genes = genes, dir = out_dir, csv = out_csv)
}

say("")
if (length(problems)) {
  say("PROBLEMS:")
  for (q in problems) say("  - ", q)
  if (isTRUE(PANEL_STRICT)) {
    writeLines(log_lines, out_path("PanelUniverses_BuildReport.txt"))
    stop("Panel universes do not match expectations -- NOTHING WRITTEN. ",
         "Every denominator in the analysis depends on these sizes. Investigate, ",
         "or set PAX6_PANEL_STRICT=FALSE to write anyway and record why.")
  }
  say("PANEL_STRICT disabled -- writing despite the problems above.")
}

for (label in names(built)) {
  b <- built[[label]]
  dir.create(b$dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(b$genes, b$csv)
  say("wrote ", b$csv, " (", nrow(b$genes), " genes)")
}

writeLines(log_lines, out_path("PanelUniverses_BuildReport.txt"))
message("\nBuild report: ", out_path("PanelUniverses_BuildReport.txt"))
