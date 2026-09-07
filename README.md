# ALVAREX RNA-seq analysis workflow

This repository contains the scripts used for RNA-seq read processing, gene-level quantification, differential-expression analysis, leave-one-genotype-out contrast sensitivity analysis, weighted gene co-expression network analysis (WGCNA), Gene Ontology enrichment, and candidate-gene prioritisation.

The repository is organised so that the full workflow can be followed from raw paired-end RNA-seq reads to the downstream statistical analyses used in the study.

## Repository structure

```text
ALVAREX_RNAseq_analysis/
├── README.md
├── METHODS_NOTES.md
├── .gitignore
│
├── scripts/
│   ├── 01_preprocessing/
│   │   ├── 01_build_STAR_index.sh
│   │   └── 02_rnaseq_pipeline.sh
│   │
│   └── 02_downstream/
│       ├── 03_edgeR_DEG_LOO_topGO.R
│       ├── 04_WGCNA.R
│       ├── 05_integrate_LOO_WGCNA.R
│       └── 06_integrate_functional_annotations.R
│
└── archive/
    └── legacy_edgeR_v4_alternative.R
```

The script stored under `archive/` is an older/alternative implementation and is not part of the primary workflow unless explicitly required to reproduce an earlier analysis.

---

## Security and portability

Absolute paths from the original HPC infrastructure are **not included** in this public repository.

Cluster usernames, personal email addresses, internal storage locations, and project-specific filesystem paths have also been removed.

Before running the scripts, insert the paths and SLURM settings corresponding to your own computing environment in the sections labelled:

```bash
# USER CONFIGURATION
```

or:

```r
## USER CONFIGURATION
```

This separation allows the analysis parameters to remain public and reproducible without exposing infrastructure-specific information.

---

# 1. Input data

The RNA-seq workflow assumes paired-end sequencing data.

The preprocessing script expects one directory per sample, for example:

```text
raw_data/
├── sample01/
│   ├── sample01_1.fq.gz
│   └── sample01_2.fq.gz
├── sample02/
│   ├── sample02_1.fq.gz
│   └── sample02_2.fq.gz
└── sample03/
    ├── sample03_1.fq.gz
    └── sample03_2.fq.gz
```

The workflow also accepts common alternatives such as:

```text
sample_R1.fq.gz / sample_R2.fq.gz
sample_1.fastq.gz / sample_2.fastq.gz
sample_R1.fastq.gz / sample_R2.fastq.gz
```

---

# 2. STAR genome index

Script:

```text
scripts/01_preprocessing/01_build_STAR_index.sh
```

Before running it, insert the paths to:

```bash
GENOME_FASTA="/path/to/reference/genome.fasta"
GTF_FILE="/path/to/reference/annotation.gtf"
STAR_INDEX_DIR="/path/to/reference/STAR_index"
```

Also provide the appropriate SLURM account and partition for your cluster:

```bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH -p YOUR_PARTITION
```

The script generates the STAR genome index required for alignment.

Submit with:

```bash
sbatch scripts/01_preprocessing/01_build_STAR_index.sh
```

The STAR index normally needs to be generated only once for a given genome/annotation configuration.

---

# 3. Read processing and gene-level quantification

Script:

```text
scripts/01_preprocessing/02_rnaseq_pipeline.sh
```

Insert your paths in the `USER CONFIGURATION` section:

```bash
DATA_DIR="/path/to/raw_data"
GENOME_DIR="/path/to/reference/STAR_index"
GTF_FILE="/path/to/reference/annotation.gtf"
RESULTS_ROOT="/path/to/results"
```

The script performs:

```text
Raw paired-end FASTQ
        |
        v
     FastQC
        |
        v
      fastp
        |
        v
FastQC after trimming
        |
        v
       STAR
        |
        v
coordinate-sorted BAM
        |
        +----> samtools index
        |
        +----> samtools flagstat
        |
        v
  featureCounts
        |
        v
gene-level count matrix
```

MultiQC summaries are generated for quality-control results.

Software used in this part of the workflow includes:

```text
FastQC      0.11.9
fastp       0.23.1
STAR        2.7.9a
samtools    1.15.1
MultiQC     1.13
Subread / featureCounts
```

The featureCounts command is configured for paired-end, reverse-stranded libraries:

```bash
-p -s 2
```

If a different library preparation protocol is used, the strandedness parameter must be adapted accordingly.

---

# 4. Differential-expression analysis

Script:

```text
scripts/02_downstream/03_edgeR_DEG_LOO_topGO.R
```

Insert:

```r
counts_file   <- "/path/to/input/Count.csv"
meta_file     <- "/path/to/input/MetaData.csv"
go_annot_file <- "/path/to/annotation/gene_to_GO.csv"
main_dir      <- "/path/to/results/edgeR_LOO_topGO"
```

The primary differential-expression analysis uses edgeR and includes:

```text
filterByExpr
TMMwsp normalisation
design: ~ 0 + Genotype + Replicate
robust dispersion estimation
robust quasi-likelihood fitting
Protected vs NonProtected mean contrast
FDR < 0.05
|log2FC| > 1
```

Positive log2 fold changes correspond to higher expression in the Protected genotype class.

---

# 5. Leave-one-genotype-out contrast sensitivity analysis

The sensitivity analysis is implemented in:

```text
03_edgeR_DEG_LOO_topGO.R
```

The edgeR model fitted to the complete dataset is retained.

At each iteration, one genotype coefficient is omitted from the calculation of its response-class mean contrast, and the Protected-vs-NonProtected contrast is recalculated using the remaining genotype coefficients.

Therefore, this analysis is a **leave-one-genotype-out contrast sensitivity analysis**. It does not remove the genotype samples and refit the dispersion/model at every iteration.

The procedure is used to identify class-associated differential-expression signals that are less dependent on any single genotype.

---

# 6. WGCNA

Script:

```text
scripts/02_downstream/04_WGCNA.R
```

Insert:

```r
counts_file   <- "/path/to/input/Count.csv"
meta_file     <- "/path/to/input/MetaData.csv"
go_annot_file <- "/path/to/annotation/gene_to_GO.csv"
base_out      <- "/path/to/results/WGCNA"
```

The network analysis includes:

```text
DESeq2 VST
blind = FALSE
design = ~ Protection + Replicate
replicate-effect correction with limma::removeBatchEffect
MAD filtering
signed WGCNA
bicor correlation
scale-free topology target R² >= 0.80
soft-threshold power = 8 in the final analysis
deepSplit = 3
minimum module size = 50
merge cut height = 0.30
```

Genes with MAD values at or above the 60th percentile are retained, corresponding approximately to the 40% most variable genes.

Module-response associations are analysed using mixed models of the form:

```text
module eigengene ~ response class + (1 | genotype)
```

Module membership (kME) is calculated using biweight midcorrelation.

Genes with:

```text
kME >= 0.80
```

are considered candidate hub genes.

---

# 7. Integration of differential expression, sensitivity and WGCNA

Script:

```text
scripts/02_downstream/05_integrate_LOO_WGCNA.R
```

Insert the directories containing the outputs from the edgeR and WGCNA analyses:

```r
loo_dir   <- "/path/to/results/edgeR_LOO_topGO/02_DEG"
wgcna_dir <- "/path/to/results/WGCNA/WGCNA_Protection_ALL_YYYYMMDD_HHMMSS"
```

This script integrates differential-expression robustness and network structure to prioritise candidate genes and exports tables suitable for downstream network visualisation.

---

# 8. Functional annotation integration

Script:

```text
scripts/02_downstream/06_integrate_functional_annotations.R
```

Insert the appropriate paths to:

```r
annotation_file   <- "/path/to/annotation/DH13M14_functional_annotation.txt"
hub_nodes_file    <- "/path/to/results/WGCNA/integration/Hub_network_nodes.csv"
tier1_both_file   <- "/path/to/results/edgeR_LOO_topGO/02_DEG/DEG_Tier1_Both_Sides.csv"
tier1_both_up_file <- "/path/to/results/edgeR_LOO_topGO/02_DEG/Tier_Split_By_Direction/Tier1_both_Up.csv"
out_dir           <- "/path/to/results/functional_annotation_integration"
```

---

# 9. Gene Ontology enrichment

GO enrichment is implemented in the primary edgeR/topGO workflow.

Analyses are performed separately for:

```text
Biological Process (BP)
Molecular Function (MF)
Cellular Component (CC)
```

using:

```text
topGO
Fisher's exact test
weight01 algorithm
```

The reference universe contains genes retained/tested in the differential-expression analysis for which GO annotations are available.

GO terms are retained for interpretation when:

```text
weight01 P < 0.05
```

and at least four genes from the tested query set are assigned to the term.

Because selection is based on the unadjusted topology-aware `weight01` P-value, GO enrichment is interpreted as an exploratory functional summary.

---

# 10. Main downstream software versions

The analyses were performed in R under Bioconductor using the following versions reported for the study:

```text
R          4.3.2
edgeR      4.2.2
WGCNA      1.73
DESeq2     1.44.0
lme4       1.1-37
lmerTest   3.1-3
topGO      2.56.0
```

The scripts also save `sessionInfo()` where applicable, allowing the exact package environment of an analysis run to be archived.

---

# 11. Full workflow

```text
Reference FASTA + GTF
        |
        v
STAR genome index
        |
        v
Raw paired-end FASTQ
        |
        v
FastQC -> fastp -> FastQC
        |
        v
STAR alignment
        |
        v
sorted BAM
        |
        v
featureCounts
        |
        v
gene-count matrix
        |
        +--------------------------+
        |                          |
        v                          v
      edgeR                      WGCNA
        |                          |
        v                          |
LOO contrast sensitivity          |
        |                          |
        +------------+-------------+
                     |
                     v
          candidate prioritisation
                     |
                     v
             functional annotation
                     |
                     v
                GO enrichment
```

---

## Reproducibility

The repository contains analysis code but does not distribute raw sequencing data or private HPC filesystem information.

To reproduce the workflow:

1. obtain the raw RNA-seq data;
2. obtain the DH13M14 reference genome and corresponding gene annotation;
3. configure the paths in each script;
4. generate the STAR index;
5. run the RNA-seq preprocessing pipeline;
6. use the resulting gene-count matrix for the downstream R analyses.

Analysis parameters are retained in the scripts so that the computational workflow can be inspected and reproduced independently of the original HPC environment.
