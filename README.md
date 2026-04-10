# A reproducibility-based cerebrospinal fluid proteomics signature outperforms published biomarkers for Alzheimer’s disease diagnosis and prognosis

## Introduction
The repository contains the code for bioinformatic analysis used in the paper "A reproducibility-based cerebrospinal fluid proteomics signature outperforms published biomarkers for Alzheimer’s disease diagnosis and prognosis". 

This study aimed to identify a reproducible protein biomarker signature using publicly available CSF proteomic datasets and to perform a head-to-head comparison of our signature against published proteomic panels to rigorously assess the stability of their diagnostic and prognostic performance across independent cohorts.

## Content
1. Boxplots for protein abundance - boxplots.R
2. Principal component analysis - pca.R
3. ROC curve analysis for diagnosis: cross-validation model for SomaScan-derived data - roc_cross_validation.R
4. ROC curve analysis for diagnosis: independent training and testing for LC-MS/MS-derived data - roc_independent_train_test.R
5. Hazard models for prognosis - prognosis_analysis.R
6. LIMMA analysis for correlation of protein abundance and clinical variables - limma_analysis.R
7. Functional enrichment and pathway analysis - enrichment_analysis.R

## Data
The proteomics data used in this study are available at:  
1. ADNI: https://adni.loni.usc.edu 
2. Bader et al. (2020): Supplemental Information https://doi.org/10.15252/msb.20199356 
3. Bangs et al. (2025): https://www.synapse.org/Synapse:syn65461849
4. Johnson et al. (2020): https://www.synapse.org/Synapse:syn20821165/wiki/596086
5. PPMI: https://www.ppmi-info.org/ 
6. Tao et al. (2024): Supplemental Information https://doi.org/10.1016/j.xinn.2023.100544 

## Instructions
The code was tested using R (v4.3.0) on Windows, but it should be compatible with later versions of R and other operating systems. 

Before running the scripts, the working directory containing the input data must be specified at the beginning of each R script. Input datasets links and preprocessing steps are described in DATA_README.md.

The scripts can be run independently, but are intended to be executed in the following order:
1. boxplots.R
2. pca.R
3. roc_cross_validation.R
4. roc_independent_train_test.R
5.prognosis_analysis.R
6. limma_analysis.R
7. enrichment_analysis.R

