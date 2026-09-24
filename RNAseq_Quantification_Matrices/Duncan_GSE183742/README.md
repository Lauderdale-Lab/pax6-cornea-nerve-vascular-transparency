# Count matrix for GSE183742 (not stored in this repository)

The optional cross-dataset check (`scripts/09_cross_dataset_replication.R`)
expects one file in this folder:

    gene_counts_featureCounts.txt

It is NOT a file distributed by GEO. It was produced by this study from the raw
reads of GEO series GSE183742 (Krishnan A, Faranda AP, Novo SG, Wang Y,
Duncan MK; SRA SRP336260; six adult cornea libraries), using the same
alignment and counting pipeline as this study's own data (`upstream/`), with
strandedness set to reverse (`-s 2`) because those libraries are
ribo-depleted and stranded.

**Where to get it:** download it from the Zenodo record `<ZENODO_DATA_DOI>` and
save it here under the name above, or regenerate it from SRA SRP336260 with the
scripts in `upstream/`.

**Check it is the right file:** its MD5 checksum should be

    <MD5 of gene_counts_featureCounts.txt>

Then point the pipeline at it before running:

    Sys.setenv(PAX6_DUNCAN_COUNTS = "RNAseq_Quantification_Matrices/Duncan_GSE183742/gene_counts_featureCounts.txt",
               PAX6_DUNCAN_META   = "Metadata_and_Sample_Guide/Duncan_GSE183742/GSE183742_sample_metadata.csv")

Without these two settings, `run_all.R` skips the cross-dataset check and
everything else runs as normal.
