#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH --job-name=RNAseq_STAR_FC
#SBATCH --cpus-per-task=75
#SBATCH --mem=300G
#SBATCH -p YOUR_PARTITION
#SBATCH -t 4-0:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=RNAseq_STAR_FC_%j.out
#SBATCH --error=RNAseq_STAR_FC_%j.err

set -euo pipefail

############################################
# 1) Load required modules
############################################

module load fastqc/0.11.9
module load fastp/0.23.1
module load star/2.7.9a
module load samtools/1.15.1
module load subread/2.0.6
module load multiqc/1.13

############################################
# 2) USER CONFIGURATION
#    Replace ONLY the paths below.
############################################

# Directory containing one folder per sample.
# Expected example:
# DATA_DIR/
# ├── sample01/
# │   ├── sample01_1.fq.gz
# │   └── sample01_2.fq.gz
# └── sample02/
#     ├── sample02_1.fq.gz
#     └── sample02_2.fq.gz
DATA_DIR="/path/to/raw_data"

# STAR index directory generated with 01_build_STAR_index.sh
GENOME_DIR="/path/to/reference/STAR_index"

# GTF annotation file used for the analysis
GTF_FILE="/path/to/reference/annotation.gtf"

# Root directory where all results will be written
RESULTS_ROOT="/path/to/results"

############################################
# 3) Output directories
############################################

QC_DIR="$RESULTS_ROOT/01_fastqc_multiqc"
TRIM_DIR="$RESULTS_ROOT/02_trimming"
MAP_DIR="$RESULTS_ROOT/03_alignment"
COUNT_DIR="$RESULTS_ROOT/04_counts"
LOG_DIR="$RESULTS_ROOT/LOGS"

FASTQC_RAW="$QC_DIR/FastQC_Raw"
FASTQC_CLEAN="$QC_DIR/FastQC_Cleaned"
MULTIQC_RAW="$QC_DIR/MultiQC_Raw"
MULTIQC_CLEAN="$QC_DIR/MultiQC_Cleaned"

TRIM_OUT="$TRIM_DIR/fastp"
STAR_OUT="$MAP_DIR/STAR"
FC_OUT="$COUNT_DIR/featureCounts"

mkdir -p \
  "$FASTQC_RAW" "$FASTQC_CLEAN" \
  "$MULTIQC_RAW" "$MULTIQC_CLEAN" \
  "$TRIM_OUT" "$STAR_OUT" "$FC_OUT" "$LOG_DIR"

############################################
# 4) Input checks
############################################

echo "Checking input files and directories..."

[[ -d "$DATA_DIR" ]] || { echo "ERROR: DATA_DIR not found: $DATA_DIR"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: GENOME_DIR not found: $GENOME_DIR"; exit 1; }
[[ -f "$GTF_FILE" ]] || { echo "ERROR: GTF file not found: $GTF_FILE"; exit 1; }

if [[ ! -f "$GENOME_DIR/SA" ]]; then
  echo "ERROR: STAR index not detected in: $GENOME_DIR"
  echo "Run 01_build_STAR_index.sh first or correct GENOME_DIR."
  exit 1
fi

############################################
# 5) Threads
############################################

CPUS="${SLURM_CPUS_PER_TASK:-20}"

FASTQC_T=$(( CPUS < 20 ? CPUS : 20 ))
FASTP_T=$(( CPUS < 20 ? CPUS : 20 ))
STAR_T="$CPUS"
SAM_T=$(( CPUS < 8 ? CPUS : 8 ))
FC_T=$(( CPUS < 20 ? CPUS : 20 ))

echo "Threads: FastQC=$FASTQC_T | fastp=$FASTP_T | STAR=$STAR_T | samtools=$SAM_T | featureCounts=$FC_T"

############################################
# 6) Per-sample processing
############################################

shopt -s nullglob

SAMPLE_DIRS=( "$DATA_DIR"/*/ )

if [[ ${#SAMPLE_DIRS[@]} -eq 0 ]]; then
  echo "ERROR: No sample directory found in: $DATA_DIR"
  exit 1
fi

for sample_dir in "${SAMPLE_DIRS[@]}"; do

  sample="$(basename "$sample_dir")"

  BAM_DONE="$STAR_OUT/$sample/${sample}_Aligned.sortedByCoord.out.bam"

  if [[ -f "$BAM_DONE" ]]; then
    echo "Already completed: $sample -> skipping alignment workflow"
    continue
  fi

  echo
  echo "========================================"
  echo "Sample: $sample"
  echo "========================================"

  ##########################################
  # 6.1 Detect paired-end FASTQ files
  ##########################################

  R1_raw="$sample_dir/${sample}_1.fq.gz"
  R2_raw="$sample_dir/${sample}_2.fq.gz"

  [[ -f "$R1_raw" ]] || R1_raw="$sample_dir/${sample}_1.fastq.gz"
  [[ -f "$R2_raw" ]] || R2_raw="$sample_dir/${sample}_2.fastq.gz"

  [[ -f "$R1_raw" ]] || R1_raw="$sample_dir/${sample}_R1.fastq.gz"
  [[ -f "$R2_raw" ]] || R2_raw="$sample_dir/${sample}_R2.fastq.gz"

  [[ -f "$R1_raw" ]] || R1_raw="$sample_dir/${sample}_R1.fq.gz"
  [[ -f "$R2_raw" ]] || R2_raw="$sample_dir/${sample}_R2.fq.gz"

  if [[ ! -f "$R1_raw" || ! -f "$R2_raw" ]]; then
    echo "WARNING: paired FASTQ files missing for $sample -> skipping"
    continue
  fi

  mkdir -p \
    "$FASTQC_RAW/$sample" \
    "$FASTQC_CLEAN/$sample" \
    "$TRIM_OUT/$sample" \
    "$STAR_OUT/$sample" \
    "$LOG_DIR/$sample"

  ##########################################
  # 6.2 FastQC on raw reads
  ##########################################

  echo "FastQC: raw reads"

  fastqc \
    -t "$FASTQC_T" \
    -q \
    --outdir "$FASTQC_RAW/$sample" \
    "$R1_raw" "$R2_raw" \
    > "$LOG_DIR/$sample/fastqc_raw.log" 2>&1

  ##########################################
  # 6.3 fastp trimming
  #     Default fastp filtering parameters
  ##########################################

  echo "fastp: trimming/filtering"

  R1_clean="$TRIM_OUT/$sample/${sample}_clean_1.fq.gz"
  R2_clean="$TRIM_OUT/$sample/${sample}_clean_2.fq.gz"

  fastp \
    -i "$R1_raw" \
    -I "$R2_raw" \
    -o "$R1_clean" \
    -O "$R2_clean" \
    -w "$FASTP_T" \
    --html "$TRIM_OUT/$sample/${sample}_fastp.html" \
    --json "$TRIM_OUT/$sample/${sample}_fastp.json" \
    > "$LOG_DIR/$sample/fastp.log" 2>&1

  ##########################################
  # 6.4 FastQC on cleaned reads
  ##########################################

  echo "FastQC: cleaned reads"

  fastqc \
    -t "$FASTQC_T" \
    -q \
    --outdir "$FASTQC_CLEAN/$sample" \
    "$R1_clean" "$R2_clean" \
    > "$LOG_DIR/$sample/fastqc_clean.log" 2>&1

  ##########################################
  # 6.5 STAR alignment
  ##########################################

  echo "STAR: alignment"

  STAR_PREFIX="$STAR_OUT/$sample/${sample}_"

  STAR \
    --genomeDir "$GENOME_DIR" \
    --runThreadN "$STAR_T" \
    --readFilesCommand zcat \
    --readFilesIn "$R1_clean" "$R2_clean" \
    --outFileNamePrefix "$STAR_PREFIX" \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile "$GTF_FILE" \
    --alignIntronMax 60000 \
    --alignMatesGapMax 10000 \
    --outFilterMultimapNmax 10 \
    --limitBAMsortRAM 20000000000 \
    --outSAMstrandField intronMotif \
    > "$LOG_DIR/$sample/STAR.log" 2>&1

  BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"

  [[ -f "$BAM" ]] || {
    echo "ERROR: STAR BAM not generated for $sample"
    exit 1
  }

  ##########################################
  # 6.6 BAM index and mapping summary
  ##########################################

  echo "samtools: BAM index and flagstat"

  samtools index -@ "$SAM_T" "$BAM"

  samtools flagstat \
    -@ "$SAM_T" \
    "$BAM" \
    > "$LOG_DIR/$sample/flagstat.txt"

done

############################################
# 7) MultiQC
############################################

echo
echo "Generating MultiQC reports..."

multiqc "$FASTQC_RAW" \
  --dirs \
  -o "$MULTIQC_RAW" \
  > "$LOG_DIR/multiqc_raw.log" 2>&1

multiqc "$FASTQC_CLEAN" \
  --dirs \
  -o "$MULTIQC_CLEAN" \
  > "$LOG_DIR/multiqc_clean.log" 2>&1

############################################
# 8) featureCounts global count matrix
############################################

echo
echo "featureCounts: generating global gene-count matrix..."

BAMS=( "$STAR_OUT"/*/*Aligned.sortedByCoord.out.bam )

if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "ERROR: No STAR BAM files found in: $STAR_OUT"
  exit 1
fi

featureCounts \
  -T "$FC_T" \
  -p \
  -t exon \
  -g gene_id \
  -s 2 \
  -a "$GTF_FILE" \
  -o "$FC_OUT/gene_counts_featureCounts.txt" \
  "${BAMS[@]}" \
  > "$LOG_DIR/featureCounts.log" 2>&1

echo
echo "Pipeline completed successfully."
echo "Counts: $FC_OUT/gene_counts_featureCounts.txt"
echo "Logs:   $LOG_DIR"
