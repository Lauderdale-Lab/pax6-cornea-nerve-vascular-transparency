###############################################################################
## 08_figures.R -- the manuscript's RNA-seq figures
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
###############################################################################
##
## THE FIGURE IS BUILT AROUND ONE CONSTRAINT: wild-type corneas do not become
## opaque under these conditions, so genotype and opacity are not two levels of
## one factor. They are two different comparisons that share only the adult
## transparent mutant, and the figure keeps them in separate blocks throughout.
##
##   LEFT BLOCK -- genotype in the transparent cornea
##     Sey versus wild type at P3/P4, at P15, and in the adult while the mutant
##     cornea is still transparent. Every cornea in this block is transparent,
##     in both genotypes, so any difference is attributable to genotype and not
##     to disease state.
##
##   RIGHT BLOCK -- the transition to opacity
##     Opaque versus transparent adult mutant. Both groups are Pax6 Sey-Neu/+,
##     abbreviated "Sey" in the figure labels, so any
##     difference is attributable to the loss of transparency and not to
##     genotype. There is no wild-type counterpart to this comparison and the
##     figure never implies one.
##
## THE OPAQUE-MUTANT-VERSUS-WILD-TYPE COMPARISON IS DELIBERATELY ABSENT. It is
## the sum of both axes and cannot separate them, which is precisely what this
## figure exists to do. It is available as a supplementary panel (set
## INCLUDE_CONFOUNDED = TRUE) and should be described as a combined comparison
## wherever it appears.
##
## THIS SCRIPT DRAWS THREE FIGURES, not one, and they are deliberately not
## combined. Panel B alone needs about 0.15 inch per gene row, so the
## program-level figure is already taller than a printed page before a heatmap
## is added; a single figure holding everything cannot be printed legibly at
## this gene count. The three answer three different questions:
##
##   FIG_DEVELOPMENT  what the normal cornea does with age, and where the
##                    mutant fails to follow. Genes whose confident expression
##                    status changes across a developmental interval, within
##                    genotype, grouped by which genotypes they cross in. This
##                    is the counterpart to the imaging in Figures 1 and 2, and
##                    it needs no fold change and no adjusted P value.
##   FIG_PROGRAM   A  does the program as a whole differ, once library
##                    detection is controlled for? Excess over an
##                    abundance-matched null, with the smallest excess that
##                    comparison could have detected.
##                 B  which individual genes cross between expressed and
##                    silent? A status difference, needing no model and no
##                    fold change.
##   FIG_NERVE        among genes expressed in both groups, which change and
##                    by how much -- axon guidance and myelination.
##   FIG_VASCULAR     the same, for the vascular/lymphatic panel. Identical
##                    construction; separate only because the two together
##                    exceed a printed page.
##
## The heatmaps carry NO panel letter. A letter is only correct inside a
## figure that has the others; a lone "C" reads as a panel lost in production.
##
## FIGURE NUMBERS ARE PARAMETERS, not string literals. They have changed twice
## already, and a journal can change them again. Set them in one place below.
## THE CONFOUNDED RUN. Setting PAX6_FIG7_CONFOUNDED=TRUE adds the opaque-mutant-
## versus-wild-type comparison, which differs from wild type in genotype and in
## transparency at once. That run is a SEPARATE, SUPPLEMENTARY FIGURE and is
## treated as one throughout: it writes only the program-level figure, under
## FigureS3_Confounded by default, with its own source data and its own caption
## numbers, and it touches none of the main figures. It is safe to run in any
## order relative to the main pass.
##
## OPEN SYMBOLS. Figures 7 and 8B draw status calls, which come from pooled
## Poisson intervals and carry no replicate variance. 07_concordance.R asks
## whether each call is also significant, in the same direction, under a
## replicate-aware negative-binomial model fitted to every panel gene. A call
## that is not is still a status call and is still drawn, but as an OPEN
## symbol, so the figure carries its own caveat. The concordance files are
## therefore required inputs here, and 07 must run before this script.
##
## Inputs : 00_config.R (for OUT_ROOT), and from that directory
##          ProgramLevel.csv, CategoricalGenes.csv, DE_Results.csv,
##          Developmental_StatusCrossings.csv,
##          Concordance/Concordance_AgeAxis.csv,
##          Concordance/Concordance_GenotypeAxis.csv
## Outputs: <FIG_DEVELOPMENT>.pdf/.png, <FIG_PROGRAM>.pdf/.png,
##          <FIG_NERVE>.pdf/.png, <FIG_VASCULAR>.pdf/.png, one
##          <FIG_*>_SourceData.csv per figure, and Figure_CaptionNumbers.txt
##          covering all four
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

if (!exists("OUT_ROOT")) stop("Source 00_config.R first.")
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required to assemble the three panels. ",
       "Install it with install.packages(\"patchwork\").")
}
library(patchwork)

## ---------------------------------------------------------------------------
## Figure parameters
## ---------------------------------------------------------------------------

## Printed at the top of every run and written into the caption file, so that
## a figure set can be traced to the version of this script that drew it. A
## layout change that leaves every number unchanged is otherwise invisible in
## the log, and two copies of this script with different legends have already
## been mistaken for one another.
FIGURES_VERSION <- "2026-09-19c (stacked legends; open symbols; S3 invocation)"
message("08_figures.R version ", FIGURES_VERSION)

## Include the opaque-mutant-versus-wild-type comparison. FALSE for the main
## figure: it confounds the two axes. TRUE produces the supplementary version.
INCLUDE_CONFOUNDED <- as.logical(Sys.getenv("PAX6_FIG7_CONFOUNDED", "FALSE"))

## Panel C can show every gene that changes in any block, or only those that
## change in at least this many comparisons. 1 shows all; 2 shows the
## reproducible core and is easier to read at journal size.
PANEL_C_MIN_CONTRASTS <- as.integer(Sys.getenv("PAX6_FIG7_MIN_CONTRASTS", "1"))

## Panel C at 103 genes is unreadable at journal width. The main figure carries
## the two nerve panels, which is what this study is about; the vascular panel
## is written as a companion figure with identical construction. Set
## PAX6_FIG7_PANELS="all" to put every panel in one figure instead.
## Resolved below, once PANEL_COLOUR exists.
PANEL_C_PANELS <- Sys.getenv("PAX6_FIG7_PANELS", "nerve")

## Output geometry. DEFINED HERE, with the other parameters, because the
## subtitle wrapping below depends on the width the figure is actually saved
## at. Reading them further down -- after build_panel_c has already been
## called -- is the same forward-reference defect as the PANEL_COLOUR one:
## it appears to work whenever a value survives in the global environment
## from a previous run, and fails in a clean session.
##
## The heatmaps have only four columns, so at full text width each tile is
## about ten times wider than it is tall and the panel reads as bars rather
## than cells. They are drawn narrower and with slightly taller rows; Figure 7
## keeps the full width because its panels use it.
HEATMAP_WIDTH   <- as.numeric(Sys.getenv("PAX6_FIG_HEATMAP_WIDTH", "5.0"))
DEV_WIDTH       <- as.numeric(Sys.getenv("PAX6_FIG_DEV_WIDTH", "6.5"))
FIGURE_WIDTH    <- as.numeric(Sys.getenv("PAX6_FIG_WIDTH", "7.2"))

## Figure numbers. Every output file name is built from these strings, so
## renumbering the manuscript is a few edits rather than a search for "Figure8"
## across the script.
##
## THE CONFOUNDED RUN IS A DIFFERENT FIGURE AND TAKES DIFFERENT NAMES BY
## DEFAULT. Adding the combined comparison changes panel A's y limits, panel
## B's gene list and ordering, and the heatmap colour limits, so nothing it
## produces is a variant of the main figures and none of it may be written over
## them. In that mode this script writes the program-level figure ONLY, under a
## supplementary name, and leaves Figures 7, 9 and 10 untouched on disk.
FIG_DEVELOPMENT <- Sys.getenv("PAX6_FIG_DEVELOPMENT", "Figure7")
FIG_PROGRAM  <- Sys.getenv("PAX6_FIG_PROGRAM",
                           if (INCLUDE_CONFOUNDED) "FigureS3_Confounded" else "Figure8")
FIG_NERVE    <- Sys.getenv("PAX6_FIG_NERVE",    "Figure9")
FIG_VASCULAR <- Sys.getenv("PAX6_FIG_VASCULAR", "Figure10")

## The caption numbers follow the figures they describe. This was the one
## output whose name did not, so a confounded run silently replaced the main
## run's numbers with numbers that belong to a different figure.
CAPTION_FILE <- Sys.getenv(
  "PAX6_FIG_CAPTION_FILE",
  if (INCLUDE_CONFOUNDED) paste0(FIG_PROGRAM, "_CaptionNumbers.txt")
  else "Figure_CaptionNumbers.txt")

## Subtitle point size, needed by the wrap calculation as well as the theme.
SUBTITLE_PT <- 7

## ggplot2 DOES NOT WRAP subtitle text. Each line is drawn as one run and is
## clipped at the device edge with no warning, so a note that fits at one
## width silently loses its ending at another. The first version of the
## grey/blank note ended "...the gene was not co" at 5.0 in for exactly this
## reason. The text is therefore wrapped here against the width the figure is
## saved at, and the panel height is derived from the number of lines the wrap
## produced rather than being set independently.
##
## 0.5 em per character is the usual average for mixed-case Helvetica; the 0.9
## factor is deliberate slack, because a clipped line leaves no trace in the
## output and cannot be caught by any check the script runs on itself.
wrap_subtitle <- function(txt, width_in, pt = SUBTITLE_PT, margin_in = 0.15) {
  n_chars <- floor(0.9 * (width_in - margin_in) * 72.27 / (pt * 0.5))
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  paste(vapply(lines,
               function(s) paste(strwrap(s, width = n_chars), collapse = "\n"),
               character(1), USE.NAMES = FALSE),
        collapse = "\n")
}

n_text_lines <- function(s) length(strsplit(s, "\n", fixed = TRUE)[[1]])

## Written as separate sentences on separate source lines so that the wrap has
## natural break points and so that shortening one note cannot run it into the
## next.
PANEL_C_SUBTITLE <- wrap_subtitle(paste(
  "Genes altered while transparent sit above genes that change only with the transition",
  "Grey: tested in that comparison but not significant",
  "Blank: not testable there -- the gene was not confidently expressed in both groups",
  sep = "\n"), HEATMAP_WIDTH)

TRANSPARENT <- c("Sey_vs_WT_P3_P4", "Sey_vs_WT_P15", "SeyT_vs_WT_Adult")
OPACITY     <- "SeyO_vs_SeyT_Adult"
CONFOUNDED  <- "SeyO_vs_WT_Adult"

CONTRAST_LABEL <- c(
  Sey_vs_WT_P3_P4    = "P3/P4",
  Sey_vs_WT_P15      = "P15",
  SeyT_vs_WT_Adult   = "Adult",
  SeyO_vs_SeyT_Adult = "Opaque vs\ntransparent",
  SeyO_vs_WT_Adult   = "Opaque Sey\nvs wild type")

BLOCK_LABEL <- c(
  transparent = "Transparent cornea\nSey vs wild type",
  opacity     = "Transition to opacity\nwithin adult Sey",
  confounded  = "Combined\nnot separable")

## Row labels for panel B. Two short lines survive rotation; the nested
## block-and-contrast strips did not.
STRIP_B <- c(
  "P3/P4"                    = "Sey vs WT\nP3/P4",
  "P15"                      = "Sey vs WT\nP15",
  "Adult"                    = "Sey vs WT\nAdult",
  "Opaque vs\ntransparent"   = "Within Sey\nOpaque vs transparent",
  "Opaque Sey\nvs wild type" = "Combined\nOpaque Sey vs WT")

PANEL_COLOUR <- c("Axon guidance"      = "#0072B2",
                  "Myelination"        = "#009E73",
                  "Vascular/Lymphatic" = "#D55E00")

## Which panels panel C draws. Placed HERE, not with the other parameters,
## because it reads PANEL_COLOUR: a forward reference works only when a stale
## copy happens to be left in the global environment from an earlier run, and
## fails the first time the script is sourced into a clean session.
## Braces are required: at top level R closes the expression at the newline and
## a bare `else` on the next line is a parse error.
panel_c_keep <- if (identical(PANEL_C_PANELS, "all")) {
  names(PANEL_COLOUR)
} else {
  c("Axon guidance", "Myelination")
}
panel_c_rest <- setdiff(names(PANEL_COLOUR), panel_c_keep)

fig_contrasts <- c(TRANSPARENT, OPACITY, if (INCLUDE_CONFOUNDED) CONFOUNDED)

block_of <- function(contrast) {
  dplyr::case_when(contrast %in% TRANSPARENT ~ "transparent",
                   contrast == OPACITY       ~ "opacity",
                   TRUE                      ~ "confounded")
}

as_factors <- function(df) {
  df %>% dplyr::mutate(
    block    = factor(BLOCK_LABEL[block_of(contrast)], levels = unname(BLOCK_LABEL)),
    contrast = factor(CONTRAST_LABEL[contrast],
                      levels = unname(CONTRAST_LABEL[fig_contrasts])),
    panel    = factor(panel, levels = names(PANEL_COLOUR)))
}

read_out <- function(f) {
  p <- file.path(OUT_ROOT, f)
  if (!file.exists(p)) stop(f, " not found in ", OUT_ROOT, ". Run run_all.R first.")
  readr::read_csv(p, show_col_types = FALSE)
}

## The concordance verdicts. Read here, once, so that both figures that draw
## status calls take the same definition of "confirmed". A missing file is an
## error, not a fallback to all-filled symbols: a figure that silently drew
## every call as confirmed would be worse than no figure.
read_concordance <- function(f) {
  p <- file.path(OUT_ROOT, "Concordance", f)
  if (!file.exists(p)) {
    stop("Concordance/", f, " not found. 07_concordance.R must run before this ",
         "script; Figures 7 and 8B mark calls it does not confirm.")
  }
  readr::read_csv(p, show_col_types = FALSE)
}
conc_age  <- read_concordance("Concordance_AgeAxis.csv") %>%
  dplyr::transmute(genotype, step, gene_id, confirmed = verdict == "confirmed")
conc_geno <- read_concordance("Concordance_GenotypeAxis.csv") %>%
  dplyr::filter(route == "categorical") %>%
  dplyr::transmute(contrast, gene_id, confirmed = verdict == "significant, same direction")

SUPPORT_LABEL <- c(`TRUE`  = "Confirmed by replicate-aware model",
                   `FALSE` = "Status level only")
support_factor <- function(confirmed) {
  factor(SUPPORT_LABEL[as.character(confirmed)], levels = unname(SUPPORT_LABEL))
}

base_theme <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey94", colour = NA),
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(0.4, "cm"),
    plot.title = ggplot2::element_text(face = "bold", size = 10))

## ---------------------------------------------------------------------------
## A. Program level
## ---------------------------------------------------------------------------

## The grey band is the null's own 95% interval: the region in which an excess
## would not have been distinguishable from the abundance-matched controls. A
## point inside the band is a bounded negative, not an absence of information,
## and the band's width is the statement of what was excluded.
prog <- read_out("ProgramLevel.csv") %>%
  dplyr::filter(contrast %in% fig_contrasts) %>%
  dplyr::select(panel, contrast,
                excess_silent, detectable_silent, padj_silent,
                excess_expressed, detectable_expressed, padj_expressed) %>%
  tidyr::pivot_longer(
    -c(panel, contrast),
    names_to = c(".value", "metric"),
    names_pattern = "(excess|detectable|padj)_(silent|expressed)") %>%
  dplyr::mutate(metric = factor(
    dplyr::recode(metric, expressed = "Confidently\nexpressed",
                  silent = "Confidently\nsilent"),
    levels = c("Confidently\nexpressed", "Confidently\nsilent")),
    significant = !is.na(padj) & padj < FDR) %>%
  as_factors()

panel_a <- ggplot2::ggplot(prog, ggplot2::aes(contrast, excess, colour = panel)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = -detectable, ymax = detectable),
    position = ggplot2::position_dodge(width = 0.6),
    colour = "grey75", linewidth = 2.2, alpha = 0.55) +
  ggplot2::geom_point(ggplot2::aes(shape = significant),
                      position = ggplot2::position_dodge(width = 0.6), size = 2.2) +
  ggplot2::facet_grid(metric ~ block, scales = "free_x", space = "free_x") +
  ggplot2::scale_colour_manual(values = PANEL_COLOUR, name = NULL) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1), name = NULL,
    labels = c(`TRUE` = "adjusted P < 0.05", `FALSE` = "not significant")) +
  ## Every point and every null band must be inside the panel; a clipped point
  ## silently understates the largest effect in the figure.
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = 0.12),
    limits = range(c(prog$excess, prog$detectable, -prog$detectable)) * 1.1) +
  ggplot2::labs(x = NULL, y = "Excess over abundance-matched null (genes)",
                title = "A   Program level",
                ## The bar is this comparison's own unadjusted sensitivity;
                ## the filled points are significant after correcting all 30
                ## tests in 05 together. A point can therefore sit just outside
                ## its bar and still be open, which is not an error and needs
                ## saying on the panel as well as in the caption.
                subtitle = paste("Grey bar: smallest excess this comparison could",
                                 "have detected, before multiplicity correction")) +
  base_theme +
  ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7, colour = "grey30"))

## ---------------------------------------------------------------------------
## B. Status changes
## ---------------------------------------------------------------------------

## Genes confidently expressed in one group and either confidently silent or
## separated by non-overlapping Poisson intervals in the other. The x axis is
## the smallest fold difference compatible with both intervals, so every point
## is a lower bound rather than a point estimate.
## Genes are ordered once, globally, by the largest separation they achieve in
## any comparison; free y scales then show each block only its own genes.
categ_symbol_order <- function(df) {
  df %>% dplyr::group_by(symbol) %>%
    dplyr::summarise(m = max(bounded_fold), .groups = "drop") %>%
    dplyr::arrange(m) %>% dplyr::pull(symbol)
}

categ <- read_out("CategoricalGenes.csv") %>%
  dplyr::filter(contrast %in% fig_contrasts) %>%
  dplyr::distinct(contrast, symbol, .keep_all = TRUE) %>%
  ## Verdict joined on the raw contrast label, before as_factors() relabels it.
  dplyr::left_join(conc_geno, by = c("contrast", "gene_id")) %>%
  as_factors() %>%
  dplyr::mutate(
    direction = factor(ifelse(direction == "up in test",
                              "Higher in Sey / opaque", "Higher in reference"),
                       levels = c("Higher in Sey / opaque", "Higher in reference")),
    symbol = factor(symbol, levels = categ_symbol_order(.)),
    ## One strip per row, short enough not to be truncated when rotated. The
    ## nested block/contrast strips were clipped at every width tried.
    strip = factor(STRIP_B[as.character(contrast)], levels = unname(STRIP_B)),
    support  = support_factor(confirmed),
    fill_key = ifelse(confirmed, as.character(panel), "open"))
if (anyNA(categ$confirmed)) {
  stop(sum(is.na(categ$confirmed)), " categorical call(s) have no concordance ",
       "verdict. Concordance_GenotypeAxis.csv is from a different run than ",
       "CategoricalGenes.csv; re-run 03 and 07.")
}
n_b_confirmed <- sum(categ$confirmed)

panel_b <- ggplot2::ggplot(categ,
    ggplot2::aes(bounded_fold, symbol, colour = panel, fill = fill_key,
                 shape = direction, size = support)) +
  ggplot2::geom_vline(xintercept = CATEGORICAL_MIN_FOLD,
                      linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_point(stroke = 0.4) +
  ## Horizontal strip text is what makes this work. With space = "free_y" a
  ## strip box is only as tall as its panel, so ROTATED text in the three-gene
  ## facets was clipped at both ends. Horizontal text needs strip WIDTH, which
  ## is shared across all facets, and two lines of height, which even the
  ## shortest facet has.
  ggplot2::facet_grid(rows = ggplot2::vars(strip),
                      scales = "free_y", space = "free_y") +
  ggplot2::scale_x_log10() +
  ggplot2::scale_colour_manual(values = PANEL_COLOUR, name = NULL) +
  ggplot2::scale_fill_manual(values = c(PANEL_COLOUR, open = "white"), guide = "none") +
  ## Fillable circle and triangle, so an unconfirmed call can be drawn open.
  ggplot2::scale_shape_manual(values = c(21, 24), name = NULL) +
  ggplot2::scale_size_manual(values = c(1.8, 1.8), name = NULL, drop = FALSE) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(order = 1),
    shape = ggplot2::guide_legend(
      order = 2, override.aes = list(colour = "grey25", fill = "grey25", size = 1.8)),
    size = ggplot2::guide_legend(
      order = 3, override.aes = list(shape = 21, colour = "grey25",
                                     fill = c("grey25", "white"), size = 1.8))) +
  ggplot2::labs(x = "Minimum fold difference (95% interval bound)", y = NULL,
                title = "B   Genes crossing between expressed and silent") +
  base_theme +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                 strip.text.y = ggplot2::element_text(angle = 0, size = 7,
                                                      hjust = 0),
                 ## As for Figure 7: stacked legend groups, height paid below.
                 legend.box = "vertical",
                 legend.box.just = "left",
                 legend.spacing.y = ggplot2::unit(0, "pt"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))

## ---------------------------------------------------------------------------
## C. Quantitative changes among genes expressed in both groups
## ---------------------------------------------------------------------------

de <- read_out("DE_Results.csv") %>%
  dplyr::filter(reportable, contrast %in% fig_contrasts)

changed <- de %>% dplyr::filter(outcome == "changed")

## Ordering carries the message: genes already altered while the cornea is
## transparent are placed above genes that change only with the transition.
gene_order <- changed %>%
  dplyr::group_by(panel, symbol) %>%
  dplyr::summarise(
    n_contrasts = dplyr::n_distinct(contrast),
    in_transparent = any(contrast %in% TRANSPARENT),
    first_change = dplyr::first(contrast[order(match(contrast, fig_contrasts))]),
    .groups = "drop") %>%
  dplyr::filter(n_contrasts >= PANEL_C_MIN_CONTRASTS) %>%
  dplyr::arrange(panel, dplyr::desc(in_transparent),
                 match(first_change, fig_contrasts), symbol)

build_panel_c <- function(keep, title) {
  ord <- gene_order %>% dplyr::filter(panel %in% keep)
  if (!nrow(ord)) return(NULL)
  h <- de %>%
    dplyr::semi_join(ord, by = c("panel", "symbol")) %>%
    dplyr::filter(panel %in% keep) %>%
    dplyr::mutate(
      shown = ifelse(outcome == "changed", log2FoldChange, NA_real_),
      symbol = factor(symbol, levels = rev(unique(ord$symbol)))) %>%
    as_factors()
  l <- max(abs(h$shown), na.rm = TRUE)
  ggplot2::ggplot(h, ggplot2::aes(contrast, symbol, fill = shown)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::facet_grid(panel ~ block, scales = "free", space = "free") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-l, l), na.value = "grey92",
      name = expression(log[2]~fold~change)) +
    ggplot2::labs(
      x = NULL, y = NULL, title = title,
      ## The grey/blank distinction is not self-evident and is easy to misread
      ## as "no data". It is stated on the panel as well as in the caption,
      ## pre-wrapped to the output width (see wrap_subtitle).
      subtitle = PANEL_C_SUBTITLE) +
    base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6),
                   panel.spacing.x = ggplot2::unit(2, "pt"),
                   panel.grid = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = SUBTITLE_PT,
                                                         colour = "grey30"))
}

heat <- de %>%
  dplyr::semi_join(gene_order, by = c("panel", "symbol")) %>%
  dplyr::mutate(
    shown = ifelse(outcome == "changed", log2FoldChange, NA_real_),
    symbol = factor(symbol, levels = rev(unique(gene_order$symbol)))) %>%
  as_factors()

if (!nrow(gene_order)) {
  stop("No gene meets PANEL_C_MIN_CONTRASTS = ", PANEL_C_MIN_CONTRASTS,
       "; lower it or check DE_Results.csv.")
}
## No panel letter: this heatmap is a figure in its own right.
panel_c <- build_panel_c(panel_c_keep, "Genes expressed in both groups that change")

## ---------------------------------------------------------------------------
## D. Developmental status crossings -- its own figure
## ---------------------------------------------------------------------------

## The other three figures compare genotypes or disease states. This one asks
## what happens with age WITHIN a genotype, which is the axis the imaging in
## Figures 1-2 is about and the only axis on which the wild-type cornea can be
## seen doing anything at all.
##
## A crossing is a change in confident expression status between two ages. It
## carries no model, no fold change and no adjusted P value, which is why it can
## be shown for genes the differential model cannot test.

DEV_STEP_LABEL <- c(P15_vs_P3_P4   = "P3/P4\nto P15",
                    Adult_vs_P15   = "P15\nto adult",
                    Adult_vs_P3_P4 = "P3/P4\nto adult")

DEV_GENOTYPE_LABEL <- c(WT = "Wild type", Sey = "Sey")

DEV_PATTERN_LABEL <- c(
  WT   = "Crosses in\nwild type only",
  both = "Crosses in\nboth genotypes",
  Sey  = "Crosses in\nSey only")

DEV_DIRECTION_LABEL <- c("off with age" = "Turns off with age",
                         "on with age"  = "Turns on with age")

## The fourth line is not decoration. An empty cell is the most persuasive
## thing in this figure and the easiest to over-read: a genotype with no point
## is a genotype in which no crossing was DETECTED, which is not a tested
## absence of change. The tested version of that claim is the genotype
## contrast, and the caption should cite it.
##
## The third line says "compares the first and last ages directly" because
## that is what the column is: its own P3/P4-versus-adult comparison under the
## same rule as the other two, NOT the union of them. A gene can therefore
## cross within an interval and show nothing across the span -- a transient
## peak that returns to where it began (Mpz), or a dip at P15 that leaves the
## first and last ages too close to separate (Mmp3 in wild type). An earlier
## wording, "the two intervals combined", made those rows look like errors.
## The line is kept under 100 characters so that it stays on one line at
## DEV_WIDTH and the figure height, which follows the line count, is unchanged.
DEV_SUBTITLE <- wrap_subtitle(paste(
  "Each point is a gene whose confident expression status changes across that interval",
  "Rows are grouped by which genotypes the gene crosses in, over any interval",
  "The third column compares the first and last ages directly: a restatement, not independent evidence",
  "A blank cell means no crossing was detected there, not a tested absence of change",
  sep = "\n"), DEV_WIDTH)

cross <- read_out("Developmental_StatusCrossings.csv")

## One row per gene per panel is the pipeline's convention, so a gene on two
## panels arrives twice. A figure has one row and one colour per gene, so the
## primary-panel assignment decides which. Written by 06; if it is absent the
## file predates that change and the figure would silently draw shared genes
## twice, in two colours.
if (!"primary_panel" %in% names(cross)) {
  stop("Developmental_StatusCrossings.csv has no primary_panel column. ",
       "Re-run 06_developmental.R, which now writes it, then re-run this script.")
}
cross <- cross %>% dplyr::filter(panel == primary_panel)

## Attach the verdict while genotype and step are still the raw codes the
## concordance file uses. Every crossing must have one; a crossing the
## concordance step did not see means the two files come from different runs.
cross <- cross %>% dplyr::left_join(conc_age, by = c("genotype", "step", "gene_id"))
if (anyNA(cross$confirmed)) {
  stop(sum(is.na(cross$confirmed)), " crossing(s) have no concordance verdict. ",
       "Concordance_AgeAxis.csv is from a different run than ",
       "Developmental_StatusCrossings.csv; re-run 06 and 07.")
}

## Which genotypes a gene crosses in, over ANY interval. Deliberately not taken
## from the whole-series column alone: some crossings appear only there, and a
## few genes cross only within an interval and not across the span.
dev_pattern <- cross %>%
  dplyr::distinct(symbol, genotype) %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    pattern = if (dplyr::n() > 1) "both" else if (genotype[1] == "WT") "WT" else "Sey",
    .groups = "drop")

## Ordered by pattern, then program, then by the largest separation the gene
## achieves anywhere -- so the strongest sit at the top of each run.
dev_order <- cross %>%
  dplyr::left_join(dev_pattern, by = "symbol") %>%
  dplyr::group_by(pattern, primary_panel, symbol) %>%
  dplyr::summarise(m = max(bounded_fold), .groups = "drop") %>%
  dplyr::arrange(factor(pattern, levels = names(DEV_PATTERN_LABEL)),
                 factor(primary_panel, levels = names(PANEL_COLOUR)),
                 dplyr::desc(m))

dev <- cross %>%
  dplyr::left_join(dev_pattern, by = "symbol") %>%
  dplyr::mutate(
    symbol    = factor(symbol, levels = rev(dev_order$symbol)),
    step      = factor(DEV_STEP_LABEL[step], levels = unname(DEV_STEP_LABEL)),
    genotype  = factor(DEV_GENOTYPE_LABEL[genotype], levels = unname(DEV_GENOTYPE_LABEL)),
    pattern   = factor(DEV_PATTERN_LABEL[pattern], levels = unname(DEV_PATTERN_LABEL)),
    panel     = factor(primary_panel, levels = names(PANEL_COLOUR)),
    direction = factor(DEV_DIRECTION_LABEL[direction],
                       levels = unname(DEV_DIRECTION_LABEL)),
    support   = support_factor(confirmed),
    ## The interior of a symbol is the panel colour when the call is confirmed
    ## and white when it is not; the outline is always the panel colour. One
    ## fill scale holds both, and its legend is suppressed because the colour
    ## scale already names the panels.
    fill_key  = ifelse(confirmed, primary_panel, "open"))

n_dev <- dplyr::n_distinct(dev$symbol)
n_dev_confirmed <- sum(dev$confirmed)

panel_dev <- ggplot2::ggplot(
    dev, ggplot2::aes(step, symbol, colour = panel, fill = fill_key,
                      shape = direction, size = support)) +
  ## The whole-series column compares the same libraries as the two beside it,
  ## first age against last, so it is not new evidence. The rule sets it apart
  ## on the panel; the subtitle says what it is.
  ggplot2::geom_vline(xintercept = 2.5, linetype = "dashed",
                      colour = "grey70", linewidth = 0.3) +
  ggplot2::geom_point(stroke = 0.4) +
  ggplot2::facet_grid(pattern ~ genotype, scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = PANEL_COLOUR, name = NULL) +
  ggplot2::scale_fill_manual(values = c(PANEL_COLOUR, open = "white"), guide = "none") +
  ## 25 and 24 are the down and up triangles that take an interior from `fill`
  ## and an outline from `colour`. Shape carries direction alone.
  ggplot2::scale_shape_manual(values = c(25, 24), name = NULL) +
  ## `size` carries no information -- both levels are the same size. It is
  ## mapped only so that a legend exists for the filled/open distinction, and
  ## its keys are overridden to show a filled and an open triangle.
  ggplot2::scale_size_manual(values = c(1.7, 1.7), name = NULL, drop = FALSE) +
  ggplot2::guides(
    shape = ggplot2::guide_legend(
      order = 2, override.aes = list(colour = "grey25", fill = "grey25", size = 1.7)),
    size = ggplot2::guide_legend(
      order = 3, override.aes = list(shape = 24, colour = "grey25",
                                     fill = c("grey25", "white"), size = 1.7)),
    colour = ggplot2::guide_legend(order = 1)) +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Genes whose expression status changes with age",
                subtitle = DEV_SUBTITLE) +
  base_theme +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6),
                 panel.grid.major.y = ggplot2::element_line(colour = "grey93",
                                                            linewidth = 0.2),
                 panel.grid.major.x = ggplot2::element_blank(),
                 strip.text.y = ggplot2::element_text(angle = 0, size = 7),
                 plot.subtitle = ggplot2::element_text(size = SUBTITLE_PT,
                                                       colour = "grey30"),
                 ## Three legend groups on one row overrun the figure width and
                 ## are clipped at both ends without warning. Stacked instead;
                 ## LEGEND_ROWS_EXTRA below pays for the height.
                 legend.box = "vertical",
                 legend.box.just = "left",
                 legend.spacing.y = ggplot2::unit(0, "pt"),
                 legend.margin = ggplot2::margin(0, 0, 0, 0))

## ---------------------------------------------------------------------------
## Assemble and save
## ---------------------------------------------------------------------------

n_b <- nrow(categ)
n_c <- sum(gene_order$panel %in% panel_c_keep)
n_v <- sum(gene_order$panel %in% panel_c_rest)

## Fixed vertical cost of a heatmap: title, axis, legend and one line per line
## of wrapped subtitle. Tying it to the wrap means the two cannot drift apart
## when a note is reworded or the width is changed.
HEATMAP_FIXED <- 1.5 + 0.20 * n_text_lines(PANEL_C_SUBTITLE)

## The developmental figure needs a little more per row than the heatmaps: it
## is points on a grid rather than adjacent tiles, so the rows have to be
## separable by eye rather than merely distinguishable in colour.
## Two extra legend rows (three groups stacked instead of one row) in the two
## figures that draw status calls. About 0.2 in per row at 9 pt.
LEGEND_ROWS_EXTRA <- 0.4

DEV_FIXED <- 1.5 + 0.20 * n_text_lines(DEV_SUBTITLE) + LEGEND_ROWS_EXTRA
h_d <- DEV_FIXED + max(3, n_dev * 0.14)

h_a <- 2.8
h_b <- n_b * 0.15 + LEGEND_ROWS_EXTRA
h_c <- HEATMAP_FIXED + max(3, n_c * 0.13)
h_v <- HEATMAP_FIXED + max(3, n_v * 0.13)

## cairo_pdf is preferred because it embeds fonts properly, but capabilities()
## reporting cairo as available is not a promise that it will draw: it can be
## compiled in and still fail on a particular plot or font. PAX6_FIG_PDF_DEVICE
## forces the base device instead, which is the first thing to try if a PDF
## goes missing while its PNG is written.
PDF_DEVICE_NAME <- Sys.getenv("PAX6_FIG_PDF_DEVICE",
                              if (capabilities("cairo")) "cairo_pdf" else "pdf")
pdf_device <- switch(PDF_DEVICE_NAME,
                     cairo_pdf = grDevices::cairo_pdf,
                     pdf       = grDevices::pdf,
                     stop("PAX6_FIG_PDF_DEVICE must be cairo_pdf or pdf, not ",
                          PDF_DEVICE_NAME))
message("PDF device: ", PDF_DEVICE_NAME)

## Every file this script produces is written and then CHECKED. ggsave can
## return without error and leave nothing behind -- a graphics device that
## fails at draw time, or a filesystem that accepts the open and loses the
## write, both look like success from R. A figure missing from the output
## directory is the kind of thing that is noticed at submission, so it is
## caught here instead.
wrote_file <- function(path, what) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(what, " was not written: ", path,
         "\n  ggsave returned without error but the file is missing or empty.",
         "\n  If this is the PDF, try PAX6_FIG_PDF_DEVICE=pdf to use the base",
         "\n  device instead of cairo_pdf, and check free space and sync status",
         "\n  on the output directory.")
  }
  invisible(path)
}

## One place where a figure is written, so the pdf/png pair can never drift
## apart in size and no figure can be saved without its source data.
save_figure <- function(plot, name, width, height, source_data) {
  pdf_path <- file.path(OUT_ROOT, paste0(name, ".pdf"))
  png_path <- file.path(OUT_ROOT, paste0(name, ".png"))
  csv_path <- file.path(OUT_ROOT, paste0(name, "_SourceData.csv"))

  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  units = "in", device = pdf_device)
  wrote_file(pdf_path, paste(name, "PDF"))

  ggplot2::ggsave(png_path, plot, width = width, height = height,
                  units = "in", dpi = 400)
  wrote_file(png_path, paste(name, "PNG"))

  ## Source data is written PER FIGURE because that is how journals ask for
  ## it, and because one combined file named after one of three figures is
  ## wrong whichever name it takes.
  readr::write_csv(source_data, csv_path)
  wrote_file(csv_path, paste(name, "source data"))
  invisible(NULL)
}

## The developmental figure does not use fig_contrasts, so a confounded run
## would rewrite it byte for byte -- but rewriting a main figure during a
## supplementary run moves its timestamp and makes the output directory lie
## about which run produced what. It is skipped instead.
if (!INCLUDE_CONFOUNDED) {
  save_figure(
    panel_dev, FIG_DEVELOPMENT, DEV_WIDTH, h_d,
    dev %>% dplyr::transmute(genotype, interval = step, panel, symbol, direction,
                             assessable, cpm_earlier, cpm_later, bounded_fold,
                             confirmed_by_replicate_model = confirmed))
}

figure_program <- panel_a / panel_b +
  patchwork::plot_layout(heights = c(h_a, h_b))

save_figure(
  figure_program, FIG_PROGRAM, FIGURE_WIDTH, h_a + h_b,
  dplyr::bind_rows(
    prog %>% dplyr::transmute(figure_panel = "A", panel, contrast, metric,
                              value = excess, limit = detectable, padj),
    categ %>% dplyr::transmute(figure_panel = "B", panel, contrast, symbol,
                               value = bounded_fold, direction,
                               confirmed_by_replicate_model = confirmed)))

## THE HEATMAPS ARE NOT WRITTEN IN CONFOUNDED MODE. They would carry the same
## row order and the same titles as Figures 9 and 10 but a different colour
## scale, because the fill limit is taken from the largest effect present. Two
## near-identical heatmaps on incomparable scales is worse than one, so the
## combined comparison is shown at program level only.
if (!INCLUDE_CONFOUNDED) {
  save_figure(
    panel_c, FIG_NERVE, HEATMAP_WIDTH, h_c,
    heat %>% dplyr::filter(!is.na(shown), panel %in% panel_c_keep) %>%
      dplyr::transmute(panel, contrast, symbol, value = shown, padj))

  ## The vascular panel is a figure of the same standing as the nerve one, not
  ## a companion to it: it carries six of the seven significant program-level
  ## comparisons. Identical construction, separate file only because the two
  ## together exceed a printed page.
  if (length(panel_c_rest)) {
    panel_v <- build_panel_c(
      panel_c_rest,
      paste("Genes expressed in both groups that change:",
            paste(panel_c_rest, collapse = ", ")))
    save_figure(
      panel_v, FIG_VASCULAR, HEATMAP_WIDTH, h_v,
      heat %>% dplyr::filter(!is.na(shown), panel %in% panel_c_rest) %>%
        dplyr::transmute(panel, contrast, symbol, value = shown, padj))
  }
}

## ---------------------------------------------------------------------------
## The numbers the captions need
## ---------------------------------------------------------------------------

## One file, sectioned by figure, because the captions are written in one
## sitting and the cross-figure totals (how the 103 entries split between the
## two heatmaps) are only checkable when they sit together.
## Both are consumed only by the main-run sections of the caption file below.
axis_split <- changed %>%
  dplyr::group_by(panel, symbol) %>%
  dplyr::summarise(in_transparent = any(contrast %in% TRANSPARENT),
                   in_opacity = any(contrast == OPACITY), .groups = "drop")

dev_counts <- dev %>%
  dplyr::distinct(symbol, pattern, panel) %>%
  dplyr::count(pattern, panel)

caption_numbers <- c(
  if (INCLUDE_CONFOUNDED)
    "Numbers for the SUPPLEMENTARY figure: the combined comparison"
  else
    "Numbers for the RNA-seq figure captions",
  paste("08_figures.R version", FIGURES_VERSION),
  strrep("=", 60), "",
  if (INCLUDE_CONFOUNDED) c(
    "This run includes the combined opaque-mutant-versus-wild-type comparison,",
    "which differs from wild type in genotype AND in transparency at once. Only",
    "the program-level figure is written; the heatmaps and the developmental",
    "figure are left as the main run produced them.",
    "",
    "The transparent / both / transition split is NOT reported here. A gene that",
    "changes only in the combined comparison belongs to neither block, so those",
    "three counts would not add up to the total and would misdescribe the figure.",
    ""),
  if (!INCLUDE_CONFOUNDED) c(
    sprintf("%s  -- developmental status crossings", FIG_DEVELOPMENT),
    sprintf("  %d distinct crossings over %d genes, in %d genotype x interval blocks.",
            nrow(dev), n_dev, dplyr::n_distinct(paste(dev$genotype, dev$step))),
    sprintf("  %d turn off with age, %d turn on.",
            sum(dev$direction == DEV_DIRECTION_LABEL[["off with age"]]),
            sum(dev$direction == DEV_DIRECTION_LABEL[["on with age"]])),
    sprintf("  %d of %d confirmed by the replicate-aware model (filled); %d drawn open.",
            n_dev_confirmed, nrow(dev), nrow(dev) - n_dev_confirmed),
    sprintf("  Open crossings: %s.",
            paste(sort(unique(as.character(dev$symbol[!dev$confirmed]))), collapse = ", ")),
    paste0("  ", apply(dev_counts, 1, function(r)
      sprintf("%s / %s: %s", gsub("\n", " ", r[["pattern"]]), r[["panel"]], r[["n"]]))),
    "  Counts are distinct crossings, resolved to the primary panel, so no gene",
    "  is counted twice. Panel-membership totals in 06's console output are larger.",
    ""),
  sprintf("%s  -- program level and status changes", FIG_PROGRAM),
  sprintf("  Panel A: %d panel x comparison tests shown; %d reach significance.",
          nrow(prog), sum(prog$significant)),
  sprintf("  Panel B: %d gene x comparison status changes, %d unique genes.",
          nrow(categ), dplyr::n_distinct(categ$symbol)),
  sprintf("    transparent block: %d   opacity block: %d",
          sum(categ$block == BLOCK_LABEL[["transparent"]]),
          sum(categ$block == BLOCK_LABEL[["opacity"]])),
  sprintf("    %d of %d confirmed by the replicate-aware model (filled); %d drawn open: %s.",
          n_b_confirmed, nrow(categ), nrow(categ) - n_b_confirmed,
          paste(sort(unique(as.character(categ$symbol[!categ$confirmed]))), collapse = ", ")),
  ## The combined block is the whole point of this run, so its size and its
  ## exclusive membership are stated rather than left to be derived from the
  ## total minus the other two.
  if (INCLUDE_CONFOUNDED) c(
    sprintf("    combined block: %d",
            sum(categ$block == BLOCK_LABEL[["confounded"]])),
    sprintf("    of which %d gene(s) appear in no other comparison: %s",
            length(setdiff(categ$symbol[categ$block == BLOCK_LABEL[["confounded"]]],
                           categ$symbol[categ$block != BLOCK_LABEL[["confounded"]]])),
            paste(sort(setdiff(
              categ$symbol[categ$block == BLOCK_LABEL[["confounded"]]],
              categ$symbol[categ$block != BLOCK_LABEL[["confounded"]]])),
              collapse = ", ")),
    sprintf("    %d of %d combined-block calls are higher in the opaque mutant.",
            sum(categ$block == BLOCK_LABEL[["confounded"]] &
                categ$direction == "Higher in Sey / opaque"),
            sum(categ$block == BLOCK_LABEL[["confounded"]]))),
  "",
  if (!INCLUDE_CONFOUNDED) c(
    sprintf("%s  -- %s", FIG_NERVE, paste(panel_c_keep, collapse = " + ")),
    sprintf("  %d gene x panel entries shown.", n_c),
    "",
    sprintf("%s  -- %s", FIG_VASCULAR, paste(panel_c_rest, collapse = ", ")),
    sprintf("  %d gene x panel entries shown.", n_v),
    "",
    "Across both heatmaps:",
    sprintf("  changed while transparent only  : %d",
            sum(axis_split$in_transparent & !axis_split$in_opacity)),
    sprintf("  changed in both blocks          : %d",
            sum(axis_split$in_transparent & axis_split$in_opacity)),
    sprintf("  changed with the transition only: %d",
            sum(!axis_split$in_transparent & axis_split$in_opacity)),
    ""),
  sprintf("FDR %.2f; equivalence bound +/-%.3f log2; categorical minimum fold %.1f.",
          FDR, EQUIV_BOUND, CATEGORICAL_MIN_FOLD),
  "Group sizes: P3/P4 wild type 4, Sey 5; P15 5 and 5; adult wild type 5,",
  "transparent Sey 5, opaque Sey 5.",
  "",
  "Wild-type corneas do not become opaque under these conditions, so the right",
  "block has no wild-type counterpart and none is implied.")
writeLines(caption_numbers, file.path(OUT_ROOT, CAPTION_FILE))

message("\n", paste(caption_numbers, collapse = "\n"))

## ---------------------------------------------------------------------------
## Remove outputs written under names this script no longer uses
## ---------------------------------------------------------------------------

## Earlier versions wrote the figures as Figure7/Figure8 and the vascular
## panel as a supplement. Left in place they would sit beside the current
## files, and there would be no way to tell from the directory which run
## produced which. Only names this script itself once used are listed; nothing
## is removed by pattern.
## NOTE Figure7.pdf/.png are NOT listed: Figure 7 is now the developmental
## figure, so those names are current outputs. `current` below would protect
## them anyway, but leaving them out of the list makes the intent explicit.
obsolete <- c("Figure7_supplement.pdf", "Figure7_supplement.png",
              "Figure8_supplement.pdf", "Figure8_supplement.png",
              "Figure7_SourceData.csv", "Figure7_CaptionNumbers.txt")
current <- c(paste0(rep(c(FIG_DEVELOPMENT, FIG_PROGRAM, FIG_NERVE, FIG_VASCULAR),
                        each = 3),
                    c(".pdf", ".png", "_SourceData.csv")),
             CAPTION_FILE, "Figure_CaptionNumbers.txt")
## A supplementary run does not tidy the main run's output directory. It knows
## only about its own figure, and deleting anything else from here would be a
## side effect of asking for one extra panel.
if (!INCLUDE_CONFOUNDED) {
  for (f in setdiff(obsolete, current)) {
    p <- file.path(OUT_ROOT, f)
    if (file.exists(p)) {
      file.remove(p)
      message("Removed ", f, ", written under an earlier figure numbering.")
    }
  }
}

message("\nWritten to ", OUT_ROOT, ":")
if (INCLUDE_CONFOUNDED) {
  message("  ", FIG_PROGRAM, ".pdf / .png    SUPPLEMENTARY: includes the combined")
  message("                                 opaque-mutant-versus-wild-type comparison")
  message("  ", FIG_PROGRAM, "_SourceData.csv")
  message("  ", CAPTION_FILE)
  message("\nNOT written in this mode, and left as the main run produced them:")
  message("  ", FIG_DEVELOPMENT, ", ", FIG_NERVE, ", ", FIG_VASCULAR,
          " and Figure_CaptionNumbers.txt")
  message("Re-run with PAX6_FIG7_CONFOUNDED unset to regenerate the main figures.")
} else {
  message("  ", FIG_DEVELOPMENT, ".pdf / .png    developmental status crossings")
  message("  ", FIG_PROGRAM, ".pdf / .png    panels A and B (program level, status changes)")
  message("  ", FIG_NERVE, ".pdf / .png    quantitative changes: ",
          paste(panel_c_keep, collapse = " + "))
  if (length(panel_c_rest)) {
    message("  ", FIG_VASCULAR, ".pdf / .png    quantitative changes: ",
            paste(panel_c_rest, collapse = ", "))
  }
  message("  one <figure>_SourceData.csv per figure, and ", CAPTION_FILE)
}
