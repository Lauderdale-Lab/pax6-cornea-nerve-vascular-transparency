###############################################################################
## Pax6_levels_JCI_models.R
##
## Reports the log2 fold change, SE and P value of ONE gene (default Pax6) from
## the JCI Insight pipeline's own DESeq2 fits, for every genotype and opacity
## contrast the pipeline defines:
##   Model A (all transparent cornea libraries): Sey vs WT at P3/P4, P15, Adult
##   Model B (adult only): Sey-T vs WT, Sey-O vs WT, Sey-O vs Sey-T
## Design in both: ~ batch + male_frac + cell_group (see 04_differential.R).
##
## It runs scripts 00, 01, 02, 03 and 04 unchanged, in one R session as
## run_all.R does, and then reads the fitted models. Nothing in the repository
## is modified; outputs go to results/<PAX6_OUT_SUBDIR> (default Pax6_levels).
##
## Run from the repository root, with
## RNAseq_Quantification_Matrices/GEO_GSE348661/gene_counts_featureCounts_all48.txt
## in place (MD5 108f41293abf5c827980d2d997feb4aa):
##
##   Rscript scripts/Pax6_levels_JCI_models.R                # Pax6
##   PAX6_TARGET_GENE=ENSMUSG00000026728 Rscript scripts/Pax6_levels_JCI_models.R
##
## Takes about 3-5 minutes (the two DESeq2 fits dominate).
## Writes, under results/<subdir>/Gene_levels/:
##   <gene>_DE_by_contrast.csv   log2FC, SE, P, genome-wide BH padj, per contrast
##   <gene>_normalised_counts.csv  size-factor-normalised counts per library
##   <gene>_CPM_all_libraries.csv  CPM for every library in the matrix
##
## The padj column is the DESeq2 adjusted P over ALL genes tested in that
## model. It is NOT the panel-restricted padj used in DE_Results.csv, which is
## corrected within a gene panel; Pax6 is in no panel.
###############################################################################

TARGET_ID <- Sys.getenv("PAX6_TARGET_GENE", unset = "ENSMUSG00000027168")  # Pax6
TARGET_ID <- sub("\\.[0-9]+$", "", TARGET_ID)

SD <- if (nzchar(Sys.getenv("PAX6_SCRIPT_DIR"))) Sys.getenv("PAX6_SCRIPT_DIR") else
  file.path(getwd(), "scripts")
Sys.setenv(PAX6_SCRIPT_DIR = SD)
if (!nzchar(Sys.getenv("PAX6_OUT_SUBDIR"))) Sys.setenv(PAX6_OUT_SUBDIR = "Pax6_levels")

for (f in c("00_config.R", "01_load.R", "02_shared_gene_assignment.R",
            "03_expression_status.R", "04_differential.R")) {
  message("\n#### ", f)
  source(file.path(SD, f), local = FALSE, echo = FALSE)
}

if (!(TARGET_ID %in% rownames(count_mat))) stop(TARGET_ID, " is not in the count matrix.")

rows <- list()
for (i in seq_len(nrow(CONTRASTS))) {
  fit <- fits[[CONTRASTS$model[i]]]
  la <- CONTRASTS$level_a[i]; lb <- CONTRASTS$level_b[i]
  if (!all(c(la, lb) %in% fit$levels)) next
  r <- DESeq2::results(fit$dds, contrast = c("cell_group", la, lb), alpha = FDR)
  if (!(TARGET_ID %in% rownames(r))) {
    rows[[length(rows) + 1L]] <- data.frame(contrast = CONTRASTS$label[i],
      model = CONTRASTS$model[i], note = "gene removed by the count pre-filter in this model")
    next
  }
  rows[[length(rows) + 1L]] <- data.frame(
    contrast = CONTRASTS$label[i], model = CONTRASTS$model[i],
    baseMean = r[TARGET_ID, "baseMean"], log2FC = r[TARGET_ID, "log2FoldChange"],
    fold = 2^r[TARGET_ID, "log2FoldChange"], lfcSE = r[TARGET_ID, "lfcSE"],
    pvalue = r[TARGET_ID, "pvalue"], padj_all_genes = r[TARGET_ID, "padj"], note = "")
}
res <- dplyr::bind_rows(rows)

outdir <- file.path(OUT_ROOT, "Gene_levels"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(outdir, paste0(TARGET_ID, "_DE_by_contrast.csv")), row.names = FALSE)

nl <- dplyr::bind_rows(lapply(names(fits), function(m) {
  v <- DESeq2::counts(fits[[m]]$dds, normalized = TRUE)
  if (!(TARGET_ID %in% rownames(v))) return(NULL)
  data.frame(model = m, sample_id = colnames(v), norm_count = as.numeric(v[TARGET_ID, ]),
             cell_group = as.character(fits[[m]]$col_data[colnames(v), "cell_group"]))
}))
write.csv(nl, file.path(outdir, paste0(TARGET_ID, "_normalised_counts.csv")), row.names = FALSE)

cpm <- count_mat[TARGET_ID, ] / colSums(count_mat) * 1e6
write.csv(data.frame(sample_id = names(cpm), cpm = as.numeric(cpm)),
          file.path(outdir, paste0(TARGET_ID, "_CPM_all_libraries.csv")), row.names = FALSE)

message("\n", TARGET_ID, " in the JCI Insight models:")
print(res, digits = 4, row.names = FALSE)
message("\nWritten to ", outdir)
