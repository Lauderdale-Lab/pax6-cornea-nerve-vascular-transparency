# References

Methods cited in the analysis pipeline for **Nerve and vascular abnormalities precede loss of corneal transparency in *Pax6*-haploinsufficient mice** (Sneha K. Mohan, James D. Lauderdale; Department of Cellular
Biology, University of Georgia, Athens, GA 30602, USA).

Each script carries a numbered
`STATISTICAL BASIS` block in its header, and inline `[n]` markers at the line
that implements the method. The numbering is local to each script; this file is
the shared list.

> Volumes and page numbers should be verified against the primary sources before
> submission. Software versions are not listed here because they vary by run —
> `RunProvenance.txt` records the exact versions used for any given set of
> results.

---

## Counts and quantification

**Kim D, Paggi JM, Park C, Bennett C, Salzberg SL** (2019) Graph-based genome
alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*
37:907–915.

**Liao Y, Smyth GK, Shi W** (2014) featureCounts: an efficient general purpose
program for assigning sequence reads to genomic features. *Bioinformatics*
30:923–930.

**Anders S, Huber W** (2010) Differential expression analysis for sequence count
data. *Genome Biology* 11:R106.
— Median-of-ratios size factors, used for the differential models. Expression
status is computed on counts-per-million instead, because it asks about absolute
detection rather than relative change.

## Expression status

**Garwood F** (1936) Fiducial limits for the Poisson distribution. *Biometrika*
28:437–442.
— Exact confidence limits for a Poisson rate. Used in `03_expression_status.R`
and, for developmental status crossings, in `06_developmental.R`.

**Ulm K** (1990) A simple method to calculate the confidence interval of a
standardized mortality ratio. *American Journal of Epidemiology* 131:373–375.
— The gamma-quantile form of the Garwood limits as implemented here:
lower = `qgamma(α, k)`, upper = `qgamma(1−α, k+1)`.

**Efron B** (1979) Bootstrap methods: another look at the jackknife. *Annals of
Statistics* 7:1–26.
— Basis for the leave-one-out stability check. A status that changes when a
single library is removed is reported as indeterminate.

**Schenker N, Gentleman JF** (2001) On judging the significance of differences by
examining the overlap between confidence intervals. *The American Statistician*
55:182–186.
— Non-overlap of two intervals is a *conservative* criterion for a difference,
more conservative than a direct test at the same level. This is what makes it
suitable for the interval-separated categorical class, which asserts a
difference without fitting a model.

## Differential expression

**Love MI, Huber W, Anders S** (2014) Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
— Negative-binomial generalised linear models, empirical-Bayes dispersion
shrinkage, Wald tests, and the likelihood-ratio test used for the adult omnibus
comparison.

**Schuirmann DJ** (1987) A comparison of the two one-sided tests procedure and
the power approach for assessing the equivalence of average bioavailability.
*Journal of Pharmacokinetics and Biopharmaceutics* 15:657–680.
— The equivalence test behind every "no change detected (bounded)" outcome,
implemented in DESeq2 as `lfcThreshold` with `altHypothesis = "lessAbs"`.

**Cook RD** (1977) Detection of influential observation in linear regression.
*Technometrics* 19:15–18.
— Outlier flagging. Group sizes here are below the threshold at which DESeq2
replaces outlying counts, so an affected gene is reported as not testable rather
than tested.

**Benjamini Y, Hochberg Y** (1995) Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *Journal of the Royal
Statistical Society B* 57:289–300.

**Bourgon R, Gentleman R, Huber W** (2010) Independent filtering increases
detection power for high-throughput experiments. *PNAS* 107:9546–9551.
— Cited for why independent filtering is **not** used for the reported adjusted
P values. The tested set is already an expression-based prefilter chosen before
any test was run, so filtering again on mean count would select on the data
twice.

## Program-level tests

**Goeman JJ, Bühlmann P** (2007) Analyzing gene expression data in terms of gene
sets: methodological issues. *Bioinformatics* 23:980–987.
— Competitive versus self-contained gene-set nulls. The program-level test
here is competitive: the panel is compared against matched non-panel genes.

**Young MD, Wakefield MJ, Smyth GK, Oshlack A** (2010) Gene ontology analysis for
RNA-seq: accounting for selection bias. *Genome Biology* 11:R14.
— The reasoning behind matching the control set on abundance. Detection
probability depends on abundance, so an unmatched control would measure that
rather than the effect.

**Wu D, Smyth GK** (2012) Camera: a competitive gene set test accounting for
inter-gene correlation. *Nucleic Acids Research* 40:e133.
— Positive correlation between genes in a set inflates competitive tests
relative to their nominal level. Panel genes are co-regulated, so the
program-level P values order the panels and contrasts rather than being exact.
This belongs in the paper's limitations.

**Phipson B, Smyth GK** (2010) Permutation P-values should never be zero:
calculating exact P-values when permutations are randomly drawn. *Statistical
Applications in Genetics and Molecular Biology* 9:Article 39.
— The `(1 + r) / (n + 1)` form, and the reason a value at that floor means only
"smaller than the resolution of this many permutations" and must not be ranked
against other floor values.

## Concordance with a replicate-aware analysis

The status calls are checked in `07_concordance.R` against the negative-binomial
models above, applied to every panel gene without the expression-status gate
(Love et al. 2014; Benjamini & Hochberg 1995, over the ungated family), and
against a second engine:

**Lund SP, Nelson D, McCarthy DJ, Smyth GK** (2012) Detecting differential
expression in RNA-sequence data using quasi-likelihood with shrunken dispersion
estimates. *Statistical Applications in Genetics and Molecular Biology*
11:Article 8.
— The quasi-likelihood F-test on negative-binomial GLMs, edgeR's `glmQLFit` /
`glmQLFTest`. Used only as an independent implementation against which the
DESeq2 calls are compared; no reported number comes from it.

**Chen Y, Lun ATL, Smyth GK** (2016) From reads to genes to pathways:
differential expression analysis of RNA-Seq experiments using Rsubread and the
edgeR quasi-likelihood pipeline. *F1000Research* 5:1438.
— The workflow followed for the second engine.

**Robinson MD, Oshlack A** (2010) A scaling normalization method for
differential expression analysis of RNA-seq data. *Genome Biology* 11:R25.
— TMM library-size normalisation, edgeR's default, used in the second engine
only. The DESeq2 models use median-of-ratios (Anders & Huber 2010, above).

Garwood (1936) and Schenker & Gentleman (2001), above, are restated locally in
`07_concordance.R` for the threshold sweep, with an acceptance test that the
local rule reproduces `03_expression_status.R` exactly.

## Cross-dataset replication

`09_cross_dataset_replication.R` calls the independent dataset GEO GSE183742 by
the same rule (Garwood 1936; Ulm 1990; Efron 1979; Schenker & Gentleman 2001,
above) and compares it against an abundance-matched permutation null (Phipson &
Smyth 2010, above).

**Cohen J** (1960) A coefficient of agreement for nominal scales. *Educational
and Psychological Measurement* 20:37–46.
— Kappa, agreement between the two datasets' calls corrected for the agreement
expected from their marginal composition. Reported beside the permutation null,
not in place of it.

**Krishnan A, Faranda AP, Novo SG, Wang Y, Duncan MK** (2024) Bioinformatic
analysis of aniridia related keratopathy. Gene Expression Omnibus, NCBI,
accession GSE183742 (BioProject PRJNA761786; SRA SRP336260). Submitted 8 Sep
2021, public 14 Sep 2024.
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183742>
— The independent adult cornea dataset (*Pax6*<sup>tm1Pgr</sup>/+; three wild-type
and three heterozygous corneas, Illumina NovaSeq 6000), re-quantified here with
the same pipeline as this study. The GEO record lists no associated
publication, so the dataset is cited by its accession.

## Design and quality control

**Nygaard V, Rødland EA, Hovig E** (2016) Methods that remove batch effects while
retaining group differences may lead to exaggerated confidence in downstream
analyses. *Biostatistics* 17:29–39.
— Why batch enters every design as a covariate rather than being removed from
the counts beforehand.

**Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, Geman D,
Baggerly K, Irizarry RA** (2010) Tackling the widespread and critical impact of
batch effects in high-throughput data. *Nature Reviews Genetics* 11:733–739.

**McMurdie PJ, Holmes S** (2014) Waste not, want not: why rarefying microbiome
data is inadmissible. *PLoS Computational Biology* 10:e1003531.
— Cited because it argues *against* rarefying for differential abundance.
Rarefaction is used here only to compare library complexity at matched depth in
the quality-control step, never for any differential test.

**Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK** (2015) limma
powers differential expression analyses for RNA-sequencing and microarray
studies. *Nucleic Acids Research* 43:e47.
— `removeBatchEffect`, used only to draw the batch-removed view of the sample
PCA in `07_concordance.R`. Nothing is computed from the adjusted values, which
is the use Nygaard et al. (2016, above) leave open: display, not inference.

The sample PCA uses DESeq2's variance-stabilising transformation (Anders &
Huber 2010; Love et al. 2014, above) with `blind = TRUE`, on the 500 most
variable genes. The off-target tissue screen in `01c_contamination_indices.R`
is descriptive (medians of marker CPM against same-age reference libraries)
and needs no methods citation.

---

## Reference genome and annotation

Alignment: mouse genome assembly **GRCm38**, HISAT2 SNP- and transcript-aware
index (`genome_snp_tran`).

Quantification: **GENCODE vM25**, the last GENCODE release built on GRCm38.
Counted at the `exon` feature level summarised to `gene_id`.

Panel curation: the curated gene lists carry identifiers from **Ensembl release
112**, which is GRCm39-era. ENSMUSG identifiers are stable across the two
assemblies, and `00b_panel_universes.R` verifies every panel gene against the
counted gene universe on each run, asserting the expected panel sizes.

Gene symbol and chromosome mapping: Bioconductor `org.Mm.eg.db`.

Record the exact package versions from `RunProvenance.txt` in the Methods.
