#!/bin/bash
#SBATCH --job-name=fastqs_multiqc
#SBATCH --output=fastqs_multiqc_%A_%a.out
#SBATCH --error=fastqs_multiqc_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --array=0-47  # Exactly 48 FastQ files

# Load FastQC module (if available)
module load fastqc 2>/dev/null || echo "FastQC module not found!"

# Activate Conda environment (if needed)
source activate mapper 2>/dev/null || echo "Conda environment 'mapper' not found!"

pwd; date  # Print working directory & timestamp

# Read FASTQ file paths into an array from fastq.fofn
readarray -t myarray < fastq.fofn
fastq=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "Processing: $fastq"

# Define output directories
FASTQC_DIR="/storage/gge/Joan/fastq_memoria/fastqc_reports"
MULTIQC_DIR="/storage/gge/Joan/fastq_memoria/multiqc_reports"

mkdir -p $FASTQC_DIR
mkdir -p $MULTIQC_DIR

# Run FastQC
fastqc -t 4 -o $FASTQC_DIR/ $fastq

# Run MultiQC only when the last job finishes
if [[ $SLURM_ARRAY_TASK_ID -eq 0 ]]; then
    echo "Running MultiQC..."
    
    # Ensure MultiQC is installed
    export PYTHONPATH=/storage/gge/Joan/bin/:$PYTHONPATH
    export PATH=/storage/gge/Joan/bin/bin/:$PATH
    
    multiqc $FASTQC_DIR -o $MULTIQC_DIR
fi

echo "✅ FastQC completed for $fastq"

