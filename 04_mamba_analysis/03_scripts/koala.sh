#!/bin/bash
# ===============================================
# Script: annotate_reactions_with_KOs.sh
# Description: 
#   Extract genes from iMM904 metabolic model,
#   retrieve their protein sequences, run KofamKOALA,
#   and generate KO annotations.
# Author: Joan Serrano
# ===============================================

# Exit on any error
set -e

# === Step 1: Extract gene-reaction pairs from the iMM904 model ===
echo "Step 1: Extracting genes from iMM904 model..."

python << EOF
from cobra.io import read_sbml_model
import pandas as pd

# Load model
model = read_sbml_model("iMM904.xml")

# Collect reaction-gene pairs
data = []
for rxn in model.reactions:
    genes = [g.id for g in rxn.genes]
    for gene in genes:
        data.append((rxn.id, gene))

# Save to CSV
df = pd.DataFrame(data, columns=["ReactionID", "Gene"])
df.to_csv("reaction_gene_map.csv", index=False)
EOF

# === Step 2: Download and decompress yeast protein sequences ===
echo "Step 2: Downloading yeast protein FASTA..."
wget -N https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz
gunzip -f orf_trans_all.fasta.gz

# === Step 3: Filter protein sequences for genes in the model ===
echo "Step 3: Filtering sequences for model genes..."
cut -d',' -f2 reaction_gene_map.csv | tail -n +2 > genes.txt
seqkit grep -f genes.txt orf_trans_all.fasta > filtered_proteins.fasta

# === Step 4: Download and set up KofamKOALA ===
echo "Step 4: Setting up KofamKOALA..."
mkdir -p kofamkoala
cd kofamkoala
wget -N https://github.com/takaram/kofamkoala/archive/refs/heads/master.zip
unzip -o master.zip
mv kofamkoala-master/* . && rmdir kofamkoala-master
wget -N ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
wget -N ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
gunzip -f ko_list.gz
tar -xzf profiles.tar.gz
cd ..

# === Step 5: Run KO annotation ===
echo "Step 5: Running KO annotation with KofamKOALA..."
./kofamkoala/exec_annotation \
  -o kofam_output.tsv \
  -p kofamkoala/profiles/ \
  -k kofamkoala/ko_list \
  --cpu 4 -f mapper \
  filtered_proteins.fasta

echo "KofamKOALA annotation complete. Output saved to kofam_output.tsv"
