#!/bin/bash
#SBATCH -A YOUR_SLURM_ACCOUNT
#SBATCH --job-name=STAR_index
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH -p YOUR_PARTITION
#SBATCH -t 0-2:00
#SBATCH --output=STAR_index_%j.out
#SBATCH --error=STAR_index_%j.err

set -euo pipefail

############################################
# 1) Load required module
############################################
module load star/2.7.9a

############################################
# 2) USER CONFIGURATION
#    Replace ONLY the paths/values below.
############################################

# Reference genome FASTA
GENOME_FASTA="/path/to/reference/genome.fasta"

# Genome annotation GTF
GTF_FILE="/path/to/reference/annotation.gtf"

# Directory where the STAR index will be created
STAR_INDEX_DIR="/path/to/reference/STAR_index"

# Read length used for RNA-seq sequencing.
# Example: 100 for 100-bp reads, 150 for 150-bp reads.
READ_LENGTH=100

# Parameter used in the original analysis.
# Adjust only if required for your genome size.
GENOME_SA_INDEX_NBASES=13

############################################
# 3) Checks and parameters
############################################

[[ -f "$GENOME_FASTA" ]] || { echo "ERROR: FASTA not found: $GENOME_FASTA"; exit 1; }
[[ -f "$GTF_FILE" ]] || { echo "ERROR: GTF not found: $GTF_FILE"; exit 1; }

mkdir -p "$STAR_INDEX_DIR"

CPUS="${SLURM_CPUS_PER_TASK:-24}"
SJDB_OVERHANG=$(( READ_LENGTH - 1 ))

echo "Building STAR genome index"
echo "Genome FASTA : $GENOME_FASTA"
echo "Annotation   : $GTF_FILE"
echo "Index output : $STAR_INDEX_DIR"
echo "Threads      : $CPUS"
echo "sjdbOverhang : $SJDB_OVERHANG"

############################################
# 4) STAR genome index generation
############################################

STAR \
  --runThreadN "$CPUS" \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX_DIR" \
  --genomeFastaFiles "$GENOME_FASTA" \
  --sjdbGTFfile "$GTF_FILE" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --genomeSAindexNbases "$GENOME_SA_INDEX_NBASES"

echo "STAR genome index successfully generated."
