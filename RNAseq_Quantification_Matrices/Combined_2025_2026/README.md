# Count matrix for this study (not stored in this repository)

The pipeline expects one file in this folder:

    gene_counts_featureCounts_all48.txt

It is the featureCounts gene-level count matrix for all 48 libraries of this
study (34 cornea, 14 trigeminal ganglion), counted against GENCODE vM25 on
GRCm38. The scripts that produced it from the raw reads are in [pax6-cornea-trigeminal-rnaseq-upstream](https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream).

**Where to get it:** download it from the processed-data files of GEO series
`<ACCESSION>` and save it here under the name above, unchanged.

**Check it is the right file:** its MD5 checksum should be

    67735b701252702ff4d1be50366e30e6

(`tools::md5sum("gene_counts_featureCounts_all48.txt")` in R, or
`md5 gene_counts_featureCounts_all48.txt` on macOS). `run_all.R` records the
checksum of the file it used in `RunProvenance.txt` on every run.

Count files are kept out of git (see `.gitignore`); data belong in a data
repository, code in this one.
