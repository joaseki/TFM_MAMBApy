#!/bin/bash
#SBATCH --job-name=fp
#SBATCH --output=fp_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --array=0-23

# Load the fastp module
module load fastp/0.24  # Use the latest version available

pwd; date  # Print working directory & timestamp

# Read file paths from pe.fofn
readarray -t myarray < pe.fofn
fastq=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "Processing: $fastq"

# Define input and output files
r1="_1.fastq.gz"
r2="_2.fastq.gz"
reports_dir="/storage/gge/Joan/fastq_memoria/trimmed_reads"
mkdir -p $reports_dir

# Run fastp
fastp --detect_adapter_for_pe \
      -i ${fastq}${r1} -I ${fastq}${r2} \
      -o ${reports_dir}/${fastq}${r1} -O ${reports_dir}/${fastq}${r2} \
      -w 4 -j ${reports_dir}/${fastq}.json

echo "Finished processing: $fastq"

