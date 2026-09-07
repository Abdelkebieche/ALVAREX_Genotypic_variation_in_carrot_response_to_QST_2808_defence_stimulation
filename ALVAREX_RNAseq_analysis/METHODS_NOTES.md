# Methods notes

This file records the key methodological details needed to keep the public code and manuscript description aligned.

## edgeR

Primary design:

```text
~ 0 + Genotype + Replicate
```

The Protected-vs-NonProtected contrast is the equally weighted mean of Protected genotype coefficients minus the equally weighted mean of NonProtected genotype coefficients.

Significance criteria:

```text
FDR < 0.05
|log2FC| > 1
```

## Leave-one-genotype-out sensitivity analysis

The implemented procedure retains the edgeR fit from the complete dataset.

At each iteration, one genotype is omitted from the corresponding class-mean contrast and the Protected-vs-NonProtected contrast is recomputed.

Samples are not removed and the edgeR model is not refitted at each iteration.

The manuscript should therefore describe this as a **leave-one-genotype-out contrast sensitivity analysis**.

## WGCNA

The supplied code uses:

```text
DESeq2 VST
blind = FALSE
design = ~ Protection + Replicate
replicate correction with limma::removeBatchEffect
MAD >= 60th percentile
signed network
bicor
soft-threshold power = 8 in the final analysis
deepSplit = 3
minModuleSize = 50
mergeCutHeight = 0.30
kME >= 0.80 for candidate hubs
```

MAD >= 60th percentile retains approximately the top 40% most variable genes.

## topGO

The analysis uses Fisher's exact test with the topology-aware `weight01` algorithm separately for BP, MF and CC.

The reference universe contains genes retained/tested in the differential-expression analysis that also have GO annotation.

Interpretation is based on:

```text
weight01 P < 0.05
```

with at least four genes from the tested query set assigned to the term.

GO enrichment is treated as exploratory rather than as an independent statistical validation.
