# Evaluation of MAMBApy for Metabolic Modelling Using Multi-Omic Data: Application in the mip6Δ Mutant of Saccharomyces cerevisiae 

This repository contains all scripts and outputs used in the master's thesis project focused on integrating transcriptomic and metabolic data in *Saccharomyces cerevisiae* under thermal stress, with particular attention to the role of **Mip6**. The analysis combines RNA-seq processing, differential expression analysis, genome-scale metabolic modeling with **MAMBApy**, and benchmarking against other methods like **iMAT** and **GIMME**.

## Repository structure

```
01_preprocessing/         # Raw data QC, trimming and quantification with kallisto
02_DESeq2/                # Differential expression analysis using DESeq2
03_benchmarking/          # Benchmark of MAMBA vs iMAT and GIMME on external datasets
04_mamba_analysis/        # Flux analysis and functional interpretation with MAMBApy
```

### 01_preprocessing/
- **scripts/**: SLURM batch scripts for data preprocessing.
- **00_raw_qc/** and **01_trimmed_qc/**: FastQC + MultiQC outputs.
- **02_kallisto_quant/**: Quantification outputs for all samples.

### 02_DESeq2/
- **scripts/**: RMarkdown file for differential analysis.
- **01_counts/**: Normalised counts and DEG matrix.
- **02_comparisons/**: Volcano plots, MA plots, and GO/KEGG enrichment per contrast.

###  03_benchmarking/
- **01_inputs/**: Expression matrices from published datasets.
- **02_scripts/**: Python/R scripts for benchmarking GIMME, iMAT, MAMBA.
- **03_results/**: Predicted fluxes and plots (error/correlation).
- **00_environment/**: Conda environment file for reproducing the analysis.

###  04_mamba_analysis/
- **01_inputs/**: Model files, normalised expression, metabolite data.
- **02_outputs/**: MAMBA results, enrichment analyses, and final figures.
- **03_scripts/**: RMarkdown and bash scripts for metabolic modeling and analysis.

## Requirements

- Python ≥ 3.8 with `cobra`, `pandas`, `numpy`, `scikit-learn`
- R ≥ 4.2 with `DESeq2`, `clusterProfiler`, `GSVA`, `ComplexHeatmap`
- MAMBApy: [https://github.com/alexpan00/MAMBApy](https://github.com/alexpan00/MAMBApy)
- Troppo (for benchmarking): [https://github.com/BioSystemsUM/troppo](https://github.com/BioSystemsUM/troppo)

##  License

This repository is shared for academic purposes and reproducibility of the TFM analyses. Contact the author for any reuse beyond this scope.

---

##  Author

Joan Serrano Quílez  
MSc in Bioinformatics and Biostatistics  

