# ROC curve cross-validation

### Analysis on ADNI database
- **Original file:** "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- **Processing:** Protein abundance log₂ transformation, inclusion of protein names, and subject diagnosis.
- **Stratification:** Different stratifications were used to evaluate ROC performance:

  Based on biological diagnosis:
  1. A−T− vs A+T+ (based on Aβ/pTau ratio; n = 703)
  2. A−T− vs A+T− (based on Aβ42; n = 388)
  3. AV45− vs AV45+ (n = 488)

  Based on biological (Aβ/pTau ratio) and clinical diagnosis (established by ADNI clinical assessment):
  4. A−T− CU vs A+T+ MCI (n = 350)
  5. A−T− CU vs A+T+ dementia (n = 265)

  Differential diagnosis:
  6. A−T− CU vs A−T− MCI (n = 276)

- **Input CSVs include:** RID, group (0 = control, 1 = group of interest), and protein log₂ abundance values per subject.





# ROC curve cross-validation
Analysis on ADNI database:
- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- Processing: protein abundance log2 transformation, inclusion of protein names, and diagnosis per subject.
- Stratification: Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/pTau ratio; n = 703), (2) A-T- vs A+T- (based on AB42, n = 388), (3) AV45- vs AV45+ (n = 488); biological (based on Aβ/pTau ratio) and clinical (stablished by ADNI based on clinical assessment) diagnosis (4) A-T- CU vs A+T+ MCI (n = 350), (5) A-T- CU vs A+T+ dementia (n = 265); and differential diagnosis (6) A-T- CU vs A-T- MCI (n = 276).
- Input CSVs include: RID, group (0-control or 1-group of interest), and protein log2 abundance values per subject.

Analysis on PPMI database:
- Original file: "ppmi_project_151.xlm" obtained from PPMI (https://www.ppmi-info.org/).
- Processing: Conversion for protein names into official gene names. 
- Stratification: Subjects with PD were those classified with sporadic PD and LRRK2; control subjects were the ones classified as Healthy Controls.
- Input CSV includes: RID, group (0-control or 1-PD), and protein log2 abundance values per subject.

# ROC curve independent training testing
Training 

Bader et al. (2020) cohort 
- Original file: "". Obtained from supplemental information in https://doi.org/10.15252/msb.20199356 
- Processing: Selection of the "Sweden" cohort. Processing performed in Metaboanalyst included batch correction, exclusion of proteins with > 80% missing values, quantile normalization, and log₂ transformation. Z-score transformation in R. 

Tao et al. (2024) cohort 
- Original file: "Table S1 identified Protein groups in CSF.xlm". Obtained from supplemental information in https://doi.org/10.1016/j.xinn.2023.100544 
- Processing: Processing performed in Metaboanalyst included batch correction, exclusion of proteins with > 80% missing values, quantile normalization, and log₂ transformation. Z-score transformation in R. 

Testing

Bangs et al. (2025) cohort 
- Original file: "EmorySubCSF_483CleanDat.xlm". Obtained from https://www.synapse.org/ (Synapse ID: syn65461849).
- Stratification: Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/tTau ratio; n = 431); biological (based on Aβ/tTau ratio) and clinical (based on MoCA) diagnosis (2) A-T- CU vs A+T+ MCI (n = 242), (3) A-T- CU vs A+T+ dementia (n = 273); and differential diagnosis (5) A-T- CU vs A-T- MCI (n = 213).
- Processing: Z-score transformation in R. 

Johnson et al. (2020) cohort
- Original file: "2b.unregressed_Batch-corrected_cleanDat_Cohort2". Obtained from https://www.synapse.org/ (Synapse ID: syn20821165)
- Stratification: Selection of A-T- and A+T+ CI subjects.
- - Processing: Z-score transformation in R. 

Input CSVs included the two training datasets and one validation dataset. Each document included ID, group (0-control or 1-group of interest), and protein log₂ z-score abundance values per subject. 

# Prognosis analysis
- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- Processing: protein abundance log2 transformation, inclusion of protein names, and diagnosis per subject. 
- Stratification: After downloading the dataset and prior to analysis, the data were split into two separate CSV files to evaluate disease-stage transitions. One dataset included cognitively unimpaired (CU) individuals with longitudinal follow-up to analyze the transition from CU to MCI (n = 143). The second dataset included A+T+ individuals (defined using the Aβ/pTau ratio) with baseline MCI and longitudinal follow-up to evaluate the transition from MCI to dementia (n = 174). 
- Input CSVs include: RID, time_followup, status, age, sex, and signature values. Values were acquired based on ROC_curve_cross_validation.R using CU vs A+T+ MCI (n = 350), and A+T+ MCI vs A+T+ dementia (n = 265) datasets.

# LIMMA analysis
Expression proteomics matrix: 

- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from ADNI (https://adni.loni.usc.edu).
- Stratification: For analysis, only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408).
- Input CSV includes: RID and protein log2 abundance values of PPAV11 and core biomarkers per subject.

Clinical data: 

- Original file: "ADNIMERGE_21Apr2025" obtained from ADNI (https://adni.loni.usc.edu).
- Processing: protein abundance log2 transformation, inclusion of protein names, and diagnosis per subject. 
- Stratification: For analysis, only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408), and clinical variables of interest.
- Input CSV includes: RID,	Age,	Sex,	PTEDUCAT,	FDG,	PIB,	AV45,	APOE4,	ABETA,	TAU,	PTAU,	pTau/AB,	CDRSB,	ADAS11,	ADAS13,	ADASQ4,	MMSE,	RAVLT_immediate,	RAVLT_learning,	RAVLT_forgetting,	RAVLT_perc_forgetting,	LDELTOTAL,	DIGITSCOR,	TRABSCOR,	FAQ,	MOCA,	EcogPtMem,	EcogPtLang,	EcogPtVisspat,	EcogPtPlan,	EcogPtOrgan,	EcogPtDivatt,	EcogPtTotal,	EcogSPMem,	EcogSPLang,	EcogSPVisspat,	EcogSPPlan,	EcogSPOrgan,	EcogSPDivatt,	EcogSPTotal, Ventricles,	Hippocampus,	WholeBrain,	Entorhinal,	Fusiform,	MidTemp,	ICV.
