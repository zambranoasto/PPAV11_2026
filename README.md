# Reproducibility-driven discovery and systematic benchmarking reveal a robust cerebrospinal fluid proteomic signature in Alzheimer’s disease

## Introduction
The repository contains the code for bioinformatic analysis performed in the paper "Reproducibility-driven discovery and systematic benchmarking reveal a robust cerebrospinal fluid proteomic signature in Alzheimer’s disease". 

This study aimed to identify a reproducible protein biomarker signature using publicly available CSF proteomic datasets and to perform a head-to-head comparison of our signature against published proteomic panels to rigorously assess the stability of their diagnostic and prognostic performance across independent cohorts.

## Content
1. Boxplots for protein abundance 
2. Principal component analysis 
3. ROC curve analysis for diagnosis: cross-validation model for SomaScan-derived data 
4. ROC curve analysis for diagnosis: independent training and testing for LC-MS/MS-derived data 
5. Hazard models for prognosis analysis
6. LIMMA analysis for correlation of protein abundance and clinical variables 
7. Functional enrichment and pathway analysis 

## Data
The proteomics datasets used in this study are available at:  
1. ADNI: https://adni.loni.usc.edu 
2. Bader et al. (2020): Supplemental Information https://doi.org/10.15252/msb.20199356 
3. Bangs et al. (2025): https://www.synapse.org/Synapse:syn65461849
4. Johnson et al. (2020): https://www.synapse.org/Synapse:syn20821165
5. PPMI: https://www.ppmi-info.org/ 
6. Tao et al. (2024): Supplemental Information https://doi.org/10.1016/j.xinn.2023.100544 

## Instructions
### Prerequisites 
The code was tested using R (v4.3.0) on Windows, but it should be compatible with later versions of R and other operating systems. Non-specific hardware is required to run these scripts; a conventional laptop is sufficient.

The following packages need to be installed before running the analyses. The code was tested using these specific versions, but later versions should be compatible. Installation typically takes a few seconds per package on a standard machine. tidyverse (v2.0.0), ggpubr (v0.6.0), rstatix (v0.7.2), clusterProfiler (v4.14.6), org.Hs.eg.db (v3.20.0), ReactomePA (v1.50.0), circlize (v0.4.16), RColorBrewer (v1.1.3), limma (v3.62.2), pheatmap (v1.0.12), viridis (v0.6.5), dplyr (v1.1.4), ggplot2 (v3.3.3), ggfortify (v0.4.17), pROC (v1.18.5), caret (v7.0.1), tibble (v3.2.1), readr (v2.1.5), survival (v3.7.0), tidyr (v1.3.1), GOplot (v1.0.2).

### Setup and execution
Before running the scripts, the working directory containing the input data must be specified at the beginning of each R script. Input dataset links and processing steps are described in DATA_README.md.

The scripts can be run independently, but are intended to be executed in the following order:
1. boxplots.R
2. pca.R
3. roc_cross_validation.R
4. roc_independent_train_test.R
5. prognosis_analysis.R
6. limma_analysis.R
7. enrichment_analysis.R

Note: The expected run time for each analysis ranges from a few seconds to less than 4 minutes on a conventional laptop.
