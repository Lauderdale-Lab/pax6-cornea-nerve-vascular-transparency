# Count matrix for this study (not stored in this repository)

The pipeline expects one file in this folder:

    gene_counts_featureCounts_all48.txt

It is the gene-level fragment-count matrix for all 48 libraries of this study
(34 cornea, 14 trigeminal ganglion), counted with featureCounts against
GENCODE vM25 on GRCm38. It was produced from the raw reads by run
`reprocess_2026-09` of
[pax6-cornea-trigeminal-rnaseq-upstream](https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream)
(commit 9182fd5; release v1.0.0).

**Where to get it:** download
`GSE348661_gene_fragment_counts_raw_GRCm38_GENCODEvM25.tsv.gz` from the
supplementary files of GEO series
[GSE348661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348661),
then decompress it into this folder under the name above:

    gunzip -c GSE348661_gene_fragment_counts_raw_GRCm38_GENCODEvM25.tsv.gz > gene_counts_featureCounts_all48.txt

**Check it is the right file:**

| File | MD5 |
|---|---|
| `GSE348661_gene_fragment_counts_raw_GRCm38_GENCODEvM25.tsv.gz` (as downloaded) | `260c993816591f5a6e0d55e525c24e9a` |
| `gene_counts_featureCounts_all48.txt` (after decompressing) | `108f41293abf5c827980d2d997feb4aa` |

(`tools::md5sum("gene_counts_featureCounts_all48.txt")` in R, or
`md5 gene_counts_featureCounts_all48.txt` on macOS.) `run_all.R` records the
checksum of the file it used in `RunProvenance.txt` on every run.

The file is tab-delimited: a `gene_id` column of versioned Ensembl gene IDs
(55,401 genes), then one column per library named by sample ID, holding raw
integer fragment counts. The loader also accepts the raw featureCounts output
of the upstream pipeline (the same counts, with coordinate columns and BAM
paths as column names).

Count files are kept out of git (see `.gitignore`); data belong in a data
repository, code in this one.
