#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=sra_download.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --qos=short

module load sra-tools  # Load module if required

cat SRR_Acc_List.txt | xargs -n 1 -P 4 fasterq-dump --split-files --gzip
