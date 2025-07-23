# SomaScan-pQTL
VIKING SomaScan Proteogenomics Analysis Pipeline

This repository contains code used to process, analyze, and interpret data from the VIKING cohort using SomaLogic's SomaScan v4.1 proteomics platform. The analysis includes pre-GWAS quality control, cis/trans pQTL mapping, colocalization, causal inference through Mendelian Randomization, and pQTL replication efforts.

# Project Overview

The pipeline is designed for:
- SomaScan v4.1 data QC and preprocessing
- GWAS and post-GWAS analyses of protein traits
- Flagging aptamers with potential technical artifacts
- Annotating associations with UniProt, gene symbols, transcription start site, and functional context
- Colocalization with external GWAS datasets
- Forward and reverse two-sample Mendelian Randomization
- Replication against other proteogenomics and eQTL studies

# Main requirements

Python ≥ 3.8
R ≥ 4.0
PLINK 1.90

# Citation
If using any of the code herein, please cite:

Kuliesius, J. et al. Efficient candidate drug target discovery through proteogenomics in a Scottish cohort. 2025. doi: https://doi.org/10.1101/2024.01.29.24301808

