#!/bin/bash
#SBATCH --job-name=kallisto_quant
#SBATCH --output=kallisto_%A_%a.out
#SBATCH --error=kallisto_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=02:30:00
#SBATCH --array=0-11  # For 12 samples
#SBATCH --qos=short

# Define paths (adjust as needed)
KALLISTO="/storage/gge/Joan/bin/kallisto"
INDEX="/storage/gge/Joan/fastq_memoria/Transcriptome/Saccharomyces_cerevisiae.idx"
TRIMMED_DIR="/storage/gge/Joan/new_fastq/trimmed_reads"
OUTPUT_DIR="/storage/gge/Joan/new_fastq/kallisto_quant"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Read sample names from trimmed.fofn (each line contains a unique sample name)
readarray -t samples < /storage/gge/Joan/new_fastq/trimmed.fofn
sample="${samples[$SLURM_ARRAY_TASK_ID]}"

echo "Processing sample: $sample"

# Construct full paths for paired-end FASTQ files using the base sample name
FASTQ_R1="${TRIMMED_DIR}/${sample}_1.fastq.gz"
FASTQ_R2="${TRIMMED_DIR}/${sample}_2.fastq.gz"

# Check if both input files exist
if [ ! -f "$FASTQ_R1" ] || [ ! -f "$FASTQ_R2" ]; then
    echo "Error: One or both input FASTQ files do not exist:"
    echo "  $FASTQ_R1"
    echo "  $FASTQ_R2"
    exit 1
fi

# Run Kallisto quantification
$KALLISTO quant \
    -i "$INDEX" \
    -t 8 \
    -o "$OUTPUT_DIR/$sample" \
    "$FASTQ_R1" "$FASTQ_R2"

echo "Finished processing sample: $sample"

