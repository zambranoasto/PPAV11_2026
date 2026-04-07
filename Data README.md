# ROC curve cross-validation

### Analysis on ADNI database
- **Original file:** "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Processing:** Protein abundance log₂ transformation, inclusion of protein names, and subject diagnosis.
- **Stratification:** Different stratifications were used to evaluate ROC performance: Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/pTau ratio; n = 703), (2) A-T- vs A+T- (based on AB42, n = 388), (3) AV45- vs AV45+ (n = 488); biological (based on Aβ/pTau ratio) and clinical (stablished by ADNI based on clinical assessment) diagnosis (4) A-T- CU vs A+T+ MCI (n = 350), (5) A-T- CU vs A+T+ dementia (n = 265); and differential diagnosis (6) A-T- CU vs A-T- MCI (n = 276).
- **Input CSVs include:** RID, diagnosis (0 = control, 1 = group of interest), and protein log₂ abundance values per subject.

### Analysis on PPMI database
- **Original file:** "ppmi_project_151.xml" obtained from PPMI (https://www.ppmi-info.org/).
- **Processing:** Conversion of protein names to official gene symbols.
- **Stratification:** Subjects with PD included those classified as sporadic PD or LRRK2-associated PD. Control subjects corresponded to individuals classified as Healthy Controls.
- **Input CSV includes:** RID, group (0 = control, 1 = PD), and protein log₂ abundance values per subject.

# ROC curve independent training testing

## Training datasets

### Bader et al. (2020) cohort
- **Original file:** "NAME". Obtained from the supplementary information of https://doi.org/10.15252/msb.20199356
- **Processing:** Selection of the Sweden cohort. Data preprocessing performed in MetaboAnalyst included batch correction, exclusion of proteins with >80% missing values, quantile normalization, and log₂ transformation. Z-score transformation was subsequently performed in R.

### Tao et al. (2024) cohort
- **Original file:** "Table S1 identified Protein groups in CSF.xml". Obtained from the supplementary information of https://doi.org/10.1016/j.xinn.2023.100544
- **Processing:** Data preprocessing in MetaboAnalyst included batch correction, exclusion of proteins with >80% missing values, quantile normalization, and log₂ transformation. Z-score transformation was subsequently performed in R.

## Testing datasets

### Bangs et al. (2025) cohort
- **Original file:** "EmorySubCSF_483CleanDat.xml". Obtained from https://www.synapse.org/ (Synapse ID: syn65461849).
- **Processing:** Z-score transformation performed in R.
- **Stratification:** Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/tTau ratio; n = 431); biological (based on Aβ/tTau ratio) and clinical (based on MoCA) diagnosis (2) A-T- CU vs A+T+ MCI (n = 242), (3) A-T- CU vs A+T+ dementia (n = 273); and differential diagnosis (5) A-T- CU vs A-T- MCI (n = 213).
- **Processing:** Z-score transformation performed in R.

### Johnson et al. (2020) cohort
- **Original file:** "2b.unregressed_Batch-corrected_cleanDat_Cohort2". Obtained from https://www.synapse.org/ (Synapse ID: syn20821165).
- **Stratification:** Selection of A−T− and A+T+ cognitively impaired subjects.
- **Processing:** Z-score transformation performed in R.

- **Input CSVs:** The final input datasets included the two training cohorts and one validation cohort. Each file contained subject ID, group (0 = control, 1 = group of interest), and protein log₂ z-score abundance values per subject.


# Prognosis analysis

- **Original file:** "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Processing:** Protein abundance log₂ transformation, inclusion of protein names, and subject diagnosis.
- **Stratification:** After downloading the dataset and prior to analysis, the data were split into two separate CSV files to evaluate disease-stage transitions. (1) One dataset included cognitively unimpaired (CU) individuals with longitudinal follow-up to analyze the transition from CU to MCI (n = 143). (2) The second dataset included A+T+ individuals (defined using the Aβ/pTau ratio) with baseline MCI and longitudinal follow-up to evaluate the transition from MCI to dementia (n = 174). 
- **Input CSVs include:** RID, time_followup, status, age, sex, and signature values.  Signature values were derived using `ROC_curve_cross_validation.R` based on CU vs A+T+ MCI (n = 350) and A+T+ MCI vs A+T+ dementia (n = 357) datasets.

# LIMMA analysis

## Expression proteomics matrix
- **Original file:** "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Stratification:** Only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408).
- **Processing:** protein abundance log2 transformation and inclusion of protein names.
- **Input CSV includes:** RID and protein log₂ abundance values for PPAV11 and core biomarkers per subject.

## Clinical data
- **Original file:** "ADNIMERGE_21Apr2025" obtained from ADNI (https://adni.loni.usc.edu).
- **Stratification:** Only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408), together with the clinical variables of interest.
- **Input CSV includes:** RID, Age, Sex, PTEDUCAT, FDG, PIB, AV45, APOE4, ABETA, TAU, PTAU, pTau/AB, CDRSB, ADAS11, ADAS13, ADASQ4, MMSE, RAVLT_immediate, RAVLT_learning, RAVLT_forgetting, RAVLT_perc_forgetting, LDELTOTAL, DIGITSCOR, TRABSCOR, FAQ, MOCA, EcogPtMem, EcogPtLang, EcogPtVisspat, EcogPtPlan, EcogPtOrgan, EcogPtDivatt, EcogPtTotal, EcogSPMem, EcogSPLang, EcogSPVisspat, EcogSPPlan, EcogSPOrgan, EcogSPDivatt, EcogSPTotal, Ventricles, Hippocampus, WholeBrain, Entorhinal, Fusiform, MidTemp, ICV.


