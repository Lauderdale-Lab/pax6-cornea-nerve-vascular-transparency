###############################################################################
## 01b_qc_library_depth.R -- library complexity at matched sequencing depth
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
## 01b_qc_library_depth.R
##
## Measures how sensitively each library detects genes, with sequencing depth
## removed from the comparison, and reports what that means for the groups the
## analysis pools over.
##
## WHY DEPTH HAS TO BE REMOVED. The number of genes detected rises with
## sequencing depth, so a raw detection count cannot distinguish a library that
## was sequenced shallowly from one whose complexity is genuinely lower. Every
## library is downsampled to a common depth by multinomial subsampling; what
## remains after depth is held constant is complexity. The depth series then
## separates the two cases directly: a merely shallow library is still climbing
## at matched depth, while a complexity-limited one has already flattened.
##
## Rarefaction is used HERE ONLY. It is not appropriate for differential testing
## -- discarding reads to equalise depth throws away information that the count
## model uses properly -- and no differential result in this analysis touches a
## rarefied count.
##
## WHY THIS MATTERS FOR THIS PIPELINE. Expression status is decided from counts
## POOLED across the libraries in a group. Pooling is deliberate: it makes the
## criterion depend on total counts and total depth rather than on how many
## libraries those were split across. But it also means a library contributes to
## a group's calls in proportion to its share of the pooled reads, so a library
## that is both large and low in complexity can shift a group's calls on its own.
## Section 3 reports each library's share of its group's pooled depth alongside
## its complexity, which is the pair of numbers that would reveal such a case.
##
## Whether any group's calls actually depend on one library is a separate
## question, answered directly by the leave-one-out check in
## 03_expression_status.R and reported in GateCost_ByCell.csv.
##
## NO LIBRARY IS EXCLUDED FROM THIS ANALYSIS. This script exists to state that
## with evidence rather than by assertion, and to make the criterion available
## for anyone who wants to apply a different one. It changes nothing.
##
## STATISTICAL BASIS
##   [1] Rarefaction by multinomial subsampling, used to compare detection at
##       matched depth and for no other purpose. On why it is inadmissible for
##       differential abundance: McMurdie PJ, Holmes S (2014) PLoS
##       Computational Biology 10:e1003531.
##   [2] Effective number of expressed genes, exp(Shannon entropy): complexity
##       in one number, and unlike a detection count it depends on no threshold.
##   [3] Median and median absolute deviation, within tissue and age stratum.
##       Detection genuinely differs by both, so a library is compared with its
##       own stratum and never with the whole set.
##
## Inputs : 00_config.R, 01_load.R
## Outputs: QC_LibraryDepth_PerSample.csv, QC_LibraryDepth_Series.csv,
##          QC_LibraryDepth_PoolShare.csv, QC_LibraryDepth_Rarefied.png,
##          QC_LibraryDepth_Series.png
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr)
  library(purrr); library(ggplot2)
})

if (!exists("count_mat")) stop("Source 00_config.R and 01_load.R first.")

qc_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]
qc_lib <- colSums(qc_mat)
message("Libraries: ", ncol(qc_mat), "   genes with any counts: ", nrow(qc_mat))
message("Library sizes (M): min ", round(min(qc_lib) / 1e6, 1),
        " (", names(which.min(qc_lib)), "), max ", round(max(qc_lib) / 1e6, 1),
        " (", names(which.max(qc_lib)), ")")

## ---------------------------------------------------------------------------
## 1. Detection at matched depth
## ---------------------------------------------------------------------------

set.seed(SEED)
subsample <- function(counts, depth) {
  as.integer(stats::rmultinom(1, size = depth, prob = counts / sum(counts)))
}

depth_target <- if (is.na(RAREFY_DEPTH)) min(qc_lib) else min(RAREFY_DEPTH, min(qc_lib))
message("\nDownsampling every library to ", format(round(depth_target), big.mark = ","),
        " assigned reads (", RAREFY_REPS, " replicate(s), seed ", SEED, ")")

per_sample <- purrr::map_dfr(colnames(qc_mat), function(s) {
  counts <- qc_mat[, s]
  reps <- vapply(seq_len(RAREFY_REPS), function(i) {
    d <- subsample(counts, depth_target)
    p <- d[d > 0] / depth_target
    c(detected = sum(d > 0), detected_ge5 = sum(d >= 5), detected_ge10 = sum(d >= 10),
      effective_genes = exp(-sum(p * log(p))))          # [2]
  }, numeric(4))
  m <- rowMeans(reps)
  tibble::tibble(
    sample_id = s, lib_size = qc_lib[[s]],
    detected_raw = sum(counts > 0),
    detected_matched = m[["detected"]],
    detected_matched_ge5 = m[["detected_ge5"]],
    detected_matched_ge10 = m[["detected_ge10"]],
    effective_genes = m[["effective_genes"]],
    pct_reads_top100 = 100 * sum(sort(counts, decreasing = TRUE)[1:100]) / sum(counts))
})

## [3] Judged within tissue and age stratum.
per_sample <- per_sample %>%
  dplyr::inner_join(dplyr::select(meta, sample_id, tissue, genotype, age_group,
                                  transparency, batch), by = "sample_id") %>%
  dplyr::group_by(tissue, age_group) %>%
  dplyr::mutate(
    stratum_n = dplyr::n(),
    stratum_median = stats::median(detected_matched),
    stratum_mad = stats::mad(detected_matched),
    ratio_to_median = detected_matched / stratum_median,
    robust_z = ifelse(stratum_mad > 0,
                      (detected_matched - stratum_median) / stratum_mad, NA_real_),
    flag_hard = !is.na(robust_z) & robust_z < -DETECT_MAD_CUTOFF,
    flag_soft = ratio_to_median < DETECT_SOFT_RATIO) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ratio_to_median)
readr::write_csv(per_sample, out_path("QC_LibraryDepth_PerSample.csv"))

message("\n=== Detection at matched depth, least sensitive first ===")
print(as.data.frame(per_sample %>% dplyr::transmute(
  sample_id, tissue, age_group, genotype,
  lib_M = round(lib_size / 1e6, 1),
  detected_raw, detected = round(detected_matched),
  effective_genes = round(effective_genes),
  top100_pct = round(pct_reads_top100, 1),
  ratio = round(ratio_to_median, 3), robust_z = round(robust_z, 2),
  flag = dplyr::case_when(flag_hard ~ "HARD", flag_soft ~ "soft", TRUE ~ ""))))

message("\nStratum sizes (MAD is unstable below about n = 5; read robust_z accordingly):")
print(as.data.frame(per_sample %>% dplyr::distinct(
  tissue, age_group, stratum_n, stratum_median = round(stratum_median))))

hard <- per_sample$sample_id[per_sample$flag_hard]
soft <- per_sample$sample_id[per_sample$flag_soft & !per_sample$flag_hard]
message("\nHARD flags (robust z < -", DETECT_MAD_CUTOFF, "): ",
        if (length(hard)) paste(hard, collapse = ", ") else "none")
message("soft flags (< ", DETECT_SOFT_RATIO * 100, "% of stratum median): ",
        if (length(soft)) paste(soft, collapse = ", ") else "none")
message("\n  Two metrics are reported because they can disagree. A library whose")
message("  deficit is confined to rare transcripts loses detected genes while its")
message("  effective gene number stays mid-range; that pattern is degraded or")
message("  low-input RNA rather than a failed library. Read both columns.")

## ---------------------------------------------------------------------------
## 2. Depth series: shallow, or complexity-limited?
## ---------------------------------------------------------------------------

depths <- sort(unique(c(5e6, 10e6, 20e6, 30e6, depth_target)))
depths <- depths[depths <= min(qc_lib)]
message("\nDepth series at: ", paste(round(depths / 1e6, 1), collapse = ", "), " M reads")

series <- purrr::map_dfr(colnames(qc_mat), function(s) {
  counts <- qc_mat[, s]
  tibble::tibble(sample_id = s, depth = depths,
                 detected = vapply(depths, function(d) sum(subsample(counts, d) > 0),
                                   numeric(1)))
}) %>%
  dplyr::inner_join(dplyr::select(meta, sample_id, tissue, genotype, age_group),
                    by = "sample_id")
readr::write_csv(series, out_path("QC_LibraryDepth_Series.csv"))

## ---------------------------------------------------------------------------
## 3. Each library's weight in the group it is pooled with
## ---------------------------------------------------------------------------

## Expression status is decided on pooled counts, so a library influences its
## group's calls in proportion to its share of the pooled reads. A library that
## carries an unusually large share AND detects unusually poorly is the case
## worth knowing about: it would pull the group's pooled CPM down for exactly
## the low-abundance genes whose status is least certain.
pool_share <- meta %>%
  dplyr::filter(tissue == TISSUE, laboratory %in% LABORATORY, !is.na(age_group)) %>%
  dplyr::mutate(pool = dplyr::case_when(
    age_group != "Adult" ~ paste0(age_group, ":", genotype),
    genotype == "WT" ~ "Adult:WT",
    transparency == "T" ~ "Adult:Sey_T",
    transparency == "O" ~ "Adult:Sey_O",
    TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(pool)) %>%
  dplyr::select(sample_id, pool) %>%
  dplyr::inner_join(dplyr::select(per_sample, sample_id, lib_size, detected_matched,
                                  effective_genes, ratio_to_median), by = "sample_id") %>%
  dplyr::group_by(pool) %>%
  dplyr::mutate(
    n_in_pool = dplyr::n(),
    pct_of_pooled_reads = round(100 * lib_size / sum(lib_size), 1),
    even_share_pct = round(100 / dplyr::n(), 1),
    weight_ratio = round(pct_of_pooled_reads / even_share_pct, 2)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(weight_ratio))
readr::write_csv(pool_share, out_path("QC_LibraryDepth_PoolShare.csv"))

message("\n=== Weight in the pooled group, against detection sensitivity ===")
print(as.data.frame(pool_share %>% dplyr::transmute(
  pool, sample_id, n_in_pool, pct_of_pooled_reads, even_share_pct, weight_ratio,
  detection_ratio = round(ratio_to_median, 3))))

concerning <- pool_share %>% dplyr::filter(weight_ratio > 1.2, ratio_to_median < 0.95)
if (nrow(concerning)) {
  message("\n  CARRIES MORE THAN ITS SHARE OF THE POOL AND DETECTS BELOW ITS STRATUM:")
  print(as.data.frame(concerning %>% dplyr::select(pool, sample_id, weight_ratio,
                                                   ratio_to_median)))
  message("  Check these against GateCost_ByCell.csv from 03_expression_status.R.")
  message("  If the leave-one-out check moves that group's calls more than its")
  message("  peers', the two observations agree and belong together in Methods.")
} else {
  message("\n  No library both carries more than its even share of a pooled group")
  message("  and detects below its stratum. Pooling is not dominated by any one")
  message("  library.")
}

## ---------------------------------------------------------------------------
## 4. Figures
## ---------------------------------------------------------------------------

genotype_colours <- c(WT = "#0072B2", Sey = "#D55E00")

p1 <- ggplot2::ggplot(per_sample, ggplot2::aes(x = age_group, y = detected_matched,
                                               colour = genotype, label = sample_id)) +
  ggplot2::geom_point(size = 2.5,
                      position = ggplot2::position_jitter(width = 0.12, height = 0)) +
  ggplot2::geom_text(size = 2.4, vjust = -0.8, show.legend = FALSE) +
  ggplot2::facet_wrap(~ tissue, scales = "free") +
  ggplot2::scale_colour_manual(values = genotype_colours) +
  ggplot2::labs(x = NULL,
                y = paste0("Genes detected at ", round(depth_target / 1e6, 1), "M reads"),
                title = "Detection sensitivity at matched depth") +
  ggplot2::theme_bw()
ggplot2::ggsave(out_path("QC_LibraryDepth_Rarefied.png"), p1,
                width = 10, height = 5.5, dpi = 200)

p2 <- ggplot2::ggplot(series, ggplot2::aes(depth / 1e6, detected, group = sample_id,
                                           colour = genotype)) +
  ggplot2::geom_line(alpha = 0.6) + ggplot2::geom_point(size = 1) +
  ggplot2::facet_grid(tissue ~ age_group, scales = "free_y") +
  ggplot2::scale_colour_manual(values = genotype_colours) +
  ggplot2::labs(x = "Assigned reads (M)", y = "Genes detected",
                title = "Rarefaction curves: a low curve that has flattened is complexity-limited") +
  ggplot2::theme_bw()
ggplot2::ggsave(out_path("QC_LibraryDepth_Series.png"), p2,
                width = 11, height = 6, dpi = 200)

## ---------------------------------------------------------------------------

message("\n", strrep("=", 74))
message("Methods sentence, generated from this run:")
message("  All ", ncol(qc_mat), " libraries were downsampled to ",
        format(round(depth_target), big.mark = ","), " assigned reads (",
        RAREFY_REPS, " replicate(s)) and detection was compared within tissue and")
message("  age stratum. ",
        if (length(hard)) paste0(length(hard), " library(ies) fell more than ",
                                 DETECT_MAD_CUTOFF, " median absolute deviations below ",
                                 "the stratum median.")
        else paste0("No library fell more than ", DETECT_MAD_CUTOFF,
                    " median absolute deviations below the stratum median."))
message("  No library was excluded from any analysis.")
message("")
message("  The criterion was applied to every library and is reported whichever")
message("  way it falls. It was not chosen to capture a library already decided")
message("  upon -- which is the only version of this check a reviewer can trust.")
message(strrep("=", 74))
