# Introduction
This file describes how the proteomic and clinical datasets were acquired, processed, and organized for the bioinformatic analyses performed in the paper "Reproducibility-driven discovery and systematic benchmarking reveal a robust cerebrospinal fluid proteomic signature in Alzheimer’s disease".

# Boxplot analysis 
The following datasets are used as input for the script boxplots.R.

### Analysis on ADNI database
- **Original files:**
  1. Proteomics: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
  2. Clinical and demographic data: "ADNIMERGE_21Apr2025.csv" obtained from ADNI (https://adni.loni.usc.edu).
  3. Analyte names: "ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Processing:** Protein abundance log₂ transformation and inclusion of protein names, subject diagnosis, and sex information. 
- **Stratification:** Different stratifications were used for boxplot analysis. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/pTau ratio; n = 703), and clinical (established by ADNI based on clinical assessment) diagnosis (2) A-T- CU vs A+T+ MCI (n = 350), (3) A-T- CU vs A+T+ dementia (n = 265).
- **Input CSVs include:** RID, diagnosis (0 = control, 1 = group of interest), and protein log₂ abundance values per subject.

### Bangs et al. (2025) cohort
- **Original file:** "EmorySubCSF_483CleanDat.xml". Obtained from https://www.synapse.org/ (Synapse ID: syn65461849).
- **Processing:** Inclusion of subject sex and ethnicity information.
- **Stratification:** Different stratifications were used for boxplot analysis. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/tTau ratio; n = 431); biological (based on Aβ/tTau ratio) and clinical (based on MoCA) diagnosis (2) A-T- CU vs A+T+ MCI (n = 242), and (3) A-T- CU vs A+T+ dementia (n = 273).
- **Input CSVs include:** ID, diagnosis (0 = control, 1 = group of interest), and protein log₂ abundance values per subject.

### Johnson et al. (2020) cohort
- **Original file:** "2b.unregressed_Batch-corrected_cleanDat_Cohort2.csv". Obtained from https://www.synapse.org/ (Synapse ID: syn20821165).
- - **Processing:** Inclusion of subject sex information.
- **Stratification:** Biological stratification based on A−T− and A+T+.
- **Input CSVs include:** ID, diagnosis (0 = control, 1 = group of interest), and protein log₂ abundance values per subject.

# ROC curve cross-validation analysis
The following datasets are used as input for the script roc_cross_validation.R and pca.R.

### Analysis on ADNI database
- **Original files:**
  1. Proteomics: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
  2. Clinical and demographic data: "ADNIMERGE_21Apr2025.csv" obtained from ADNI (https://adni.loni.usc.edu).
  3. Analyte names: "ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Processing:** Protein abundance log₂ transformation, inclusion of protein names, and subject diagnosis.
- **Stratification:** Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/pTau ratio; n = 703), (2) A-T- vs A+T- (based on AB42, n = 388), (3) AV45- vs AV45+ (n = 488); biological (based on Aβ/pTau ratio) and clinical (established by ADNI based on clinical assessment) diagnosis (4) A-T- CU vs A+T+ MCI (n = 350), (5) A-T- CU vs A+T+ dementia (n = 265); and differential diagnosis (6) A-T- CU vs A-T- MCI (n = 276).
- **Input CSVs include:** RID, diagnosis (0 = control, 1 = group of interest), and protein log₂ abundance values per subject.

### Analysis on PPMI database
- **Original file:** "ppmi_project_151.xml" obtained from PPMI (https://www.ppmi-info.org/).
- **Processing:** Conversion of protein names to official gene symbols.
- **Stratification:** Subjects with PD included those classified as sporadic PD or LRRK2-associated PD. Control subjects corresponded to individuals classified as Healthy Controls.
- **Input CSV includes:** ID, group (0 = control, 1 = PD), and protein log₂ abundance values per subject.

# ROC curve independent training testing analysis
The following datasets are used as input for the script roc_independent_train_test.R.

## Training datasets

### Bader et al. (2020) cohort
- **Original file:** "NAME". Obtained from the supplementary information of https://doi.org/10.15252/msb.20199356
- **Processing:** Selection of the Sweden cohort. Data preprocessing performed in MetaboAnalyst included batch correction, exclusion of proteins with > 80% missing values, quantile normalization, and log₂ transformation. Z-score transformation was subsequently performed in R.

### Tao et al. (2024) cohort
- **Original file:** "Table S1 identified Protein groups in CSF.xml". Obtained from the supplementary information of https://doi.org/10.1016/j.xinn.2023.100544
- **Processing:** Data preprocessing in MetaboAnalyst included batch correction, exclusion of proteins with > 80% missing values, quantile normalization, and log₂ transformation. Z-score transformation was subsequently performed in R.

## Testing datasets

### Bangs et al. (2025) cohort
- **Original file:** "EmorySubCSF_483CleanDat.xml". Obtained from https://www.synapse.org/ (Synapse ID: syn65461849).
- **Stratification:** Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/tTau ratio; n = 431); biological (based on Aβ/tTau ratio) and clinical (based on MoCA) diagnosis (2) A-T- CU vs A+T+ MCI (n = 242), (3) A-T- CU vs A+T+ dementia (n = 273); and differential diagnosis (5) A-T- CU vs A-T- MCI (n = 213).
- **Processing:** Z-score transformation performed in R.

### Johnson et al. (2020) cohort
- **Original file:** "2b.unregressed_Batch-corrected_cleanDat_Cohort2". Obtained from https://www.synapse.org/ (Synapse ID: syn20821165).
- **Stratification:** Selection of A−T− and A+T+ cognitively impaired subjects.
- **Processing:** Z-score transformation performed in R.

- **Input CSVs:** The final input datasets included the two training cohorts and one validation cohort. Each file contained ID, group (0 = control, 1 = group of interest), and protein log₂ z-score abundance values per subject.

# Prognosis analysis
The following datasets are used as input for the script prognosis_analysis.R.

### Acquisition of signature scores 
- **Original files:**
  1. Proteomics: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
  2. Clinical and demographic data: "ADNIMERGE_21Apr2025.csv" obtained from ADNI (https://adni.loni.usc.edu).
  3. Analyte names: "ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Processing:** Protein abundance log₂ transformation, inclusion of protein names, and subject diagnosis.
- **Stratification:** Before analysis, the data were split into two separate CSV files (1) A-T- CU vs A+T+ MCI (n = 350), and (2) A+T+ MCI vs A+T+ dementia (n = 413). 
- **Input CSV includes:** RID, group (0 = control, 1 = PD), and protein log₂ abundance values per subject.

### Prognosis analysis
- **Input files:** Output files obtained from "Acquisition of signature scores"
- **Processing:** Inclusion of clinical data (age, sex, diagnosis, time of visit, and change in diagnosis) from "ADNIMERGE_21Apr2025". 
- **Stratification:** Only participants with longitudinal follow-up were included in the analysis. (1) CU individuals to analyze the transition from CU to MCI (n = 143), and (2) A+T+ individuals (defined using the Aβ/pTau ratio) with baseline MCI to evaluate the transition from MCI to dementia (n = 174). 
- **Input CSVs include:** RID, time_followup, status, age, sex, and signature scores.

# LIMMA analysis
The following datasets are used as input for the script limma_analysis.R.

## Expression proteomics matrix
- **Original file:** "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Stratification:** Only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408).
- **Processing:** protein abundance log2 transformation and inclusion of protein names.
- **Input CSV includes:** RID and protein log₂ abundance values for PPAV11 and core biomarkers per subject.

## Clinical data
- **Original file:** "ADNIMERGE_21Apr2025.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Stratification:** Only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408), together with the clinical variables of interest.
- **Input CSV includes:** RID, Age, Sex, PTEDUCAT, FDG, PIB, AV45, APOE4, ABETA, TAU, PTAU, pTau/AB, CDRSB, ADAS11, ADAS13, ADASQ4, MMSE, RAVLT_immediate, RAVLT_learning, RAVLT_forgetting, RAVLT_perc_forgetting, LDELTOTAL, DIGITSCOR, TRABSCOR, FAQ, MOCA, EcogPtMem, EcogPtLang, EcogPtVisspat, EcogPtPlan, EcogPtOrgan, EcogPtDivatt, EcogPtTotal, EcogSPMem, EcogSPLang, EcogSPVisspat, EcogSPPlan, EcogSPOrgan, EcogSPDivatt, EcogSPTotal, Ventricles, Hippocampus, WholeBrain, Entorhinal, Fusiform, MidTemp, ICV.


