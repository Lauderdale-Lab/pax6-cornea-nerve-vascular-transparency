# Count matrix for GSE183742

The optional cross-dataset check (`scripts/09_cross_dataset_replication.R`)
reads one file in this folder, which is included in this repository:

    gene_counts_featureCounts.txt

It is NOT a file distributed by GEO. It was produced by this study from the raw
reads of GEO series GSE183742 (Krishnan A, Faranda AP, Novo SG, Wang Y,
Duncan MK; SRA SRP336260; six adult cornea libraries), using the same
alignment and counting pipeline as this study's own data ([pax6-cornea-trigeminal-rnaseq-upstream](https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream)), with
strandedness set to reverse (`-s 2`) because those libraries are
ribo-depleted and stranded.

**To regenerate it** from SRA SRP336260, use the scripts in [pax6-cornea-trigeminal-rnaseq-upstream](https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream).

**Check it is the right file:** its MD5 checksum should be

    bdb3e653f08241a9951322a5ea73b795

Then point the pipeline at it before running:

    Sys.setenv(PAX6_DUNCAN_COUNTS = "RNAseq_Quantification_Matrices/Duncan_GSE183742/gene_counts_featureCounts.txt",
               PAX6_DUNCAN_META   = "Metadata_and_Sample_Guide/Duncan_GSE183742/GSE183742_sample_metadata.csv")

Without these two settings, `run_all.R` skips the cross-dataset check and
everything else runs as normal.
