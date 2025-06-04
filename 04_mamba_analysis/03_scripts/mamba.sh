#!/bin/bash
#SBATCH --job-name=mambapy_experiment
#SBATCH --output=mambapy_experiment.out
#SBATCH --error=mambapy_experiment.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=24:00:00

# Load any modules you need (if your cluster uses modules)
# module load python/3.8

# Create and activate virtual environment
python3 -m venv myenv
source myenv/bin/activate

# Install your dependencies, if not already present
# pip install -r requirements.txt

# Run MAMBApy with your inputs
MAMBApy iMM904_YPD.pkl \
         DEG_matrix_all_contrasts.csv \
         -m mets_full_experiment.csv \
         -i matrix_model.txt \
         -f results/MAMBA \
         -s

# Deactivate the environment when done
deactivate
