# 10_contamination_indices.R -- off-target tissue indices per library (Supp. Table 6)
# Reproduces the contamination columns of RNAseq_Library_Metadata_Supp_Table_4_6.
# Same marker families as the 09 margin check; run after 01_load.R (needs `counts`, a
# gene x library integer matrix with version-stripped Ensembl IDs as rownames).
suppressPackageStartupMessages({ library(org.Mm.eg.db); library(AnnotationDbi) })

CONTAM_SETS <- list(
  pigment = c("Tyr", "Tyrp1", "Dct", "Mlana"),   # melanocyte-lineage markers; = adherent iris in cornea. Pmel excluded (12-55 CPM in every cornea)
  retina = c("Rho", "Gnat1", "Sag", "Pde6b", "Rcvrn", "Crx", "Nrl", "Rpe65"),
  lens   = c("Cryaa", "Cryba1", "Crybb2", "Crygd", "Mip", "Bfsp2"),   # Cryab excluded: corneal expression
  angle  = c("Myoc"),                                                 # Chi3l1 excluded: inflammation in Sey_O
  conj   = c("Krt13", "Krt4", "Muc5ac")
)
CONTAM_SINGLE   <- c("Optc", "Rho", "Cryaa", "Kera", "Krt12")
CONTAM_FOLD     <- 10   # index vs median of same-age cornea libraries excluding Sey_O
CONTAM_TRACE    <- 5
CONTAM_GENEFOLD <- 5    # >= 2 constituent genes must each exceed this (angle: single gene)
STROMA_POOR_CPM <- 100  # Kera below this = stroma-poor (criterion applied to Duncan)

contamination_indices <- function(counts, meta) {
  sym <- AnnotationDbi::mapIds(org.Mm.eg.db, rownames(counts), "SYMBOL", "ENSEMBL")
  cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
  get <- function(genes) { i <- match(genes, sym); stopifnot(!anyNA(i)); cpm[i, , drop = FALSE] }
  idx <- as.data.frame(sapply(CONTAM_SETS, function(gs) colSums(get(gs))))
  for (g in CONTAM_SINGLE) idx[[g]] <- as.numeric(get(g))
  idx$library_id <- colnames(counts)
  m <- meta[match(idx$library_id, meta$sample_id), ]
  idx$contamination_flags <- ""
  for (i in which(m$tissue == "cornea")) {
    ref <- which(m$tissue == "cornea" & m$age_group == m$age_group[i] & m$group3 != "Sey_O")
    f <- character()
    for (k in c("pigment", "retina", "lens", "angle")) {
      fold <- idx[[k]][i] / max(median(idx[[k]][ref]), 0.5)
      gs <- CONTAM_SETS[[k]]; gc <- get(gs)
      hits <- gs[gc[, i] > CONTAM_GENEFOLD * pmax(apply(gc[, ref, drop = FALSE], 1, median), 0.2)]
      if (fold > CONTAM_FOLD && (k == "angle" || length(hits) >= 2))
        f <- c(f, sprintf("%s %.0fx (%s)", if (k == "pigment") "pigment cells = iris" else k, fold, paste(hits, collapse = ", ")))
      else if (fold > CONTAM_TRACE) f <- c(f, sprintf("%s trace %.0fx", k, fold))
    }
    if (idx$Kera[i] < STROMA_POOR_CPM) f <- c(f, sprintf("stroma-poor (Kera %.0f CPM)", idx$Kera[i]))
    idx$contamination_flags[i] <- paste(f, collapse = "; ")
  }
  idx
}
# Published Supp. Table 6 subset: pigment, Rho, Cryaa, angle (Myoc), Kera, Krt12, contamination_flags
# contam <- contamination_indices(counts, meta)
# write.csv(contam, file.path(OUT_DIR, "Contamination_Indices.csv"), row.names = FALSE)
