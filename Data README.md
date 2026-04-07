# Prognosis analysis
- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from the ADNI database.
- Processing: protein abundance log2 transformation, inclusion of protein names, and diagnosis per subject. 
- Stratification: After downloading the dataset and prior to analysis, the data were split into two separate CSV files to evaluate disease-stage transitions. One dataset included cognitively unimpaired (CU) individuals with longitudinal follow-up to analyze the transition from CU to MCI (n = 143). The second dataset included A+T+ individuals (defined using the Aβ/pTau ratio) with baseline MCI and longitudinal follow-up to evaluate the transition from MCI to dementia (n = 174). 
- Input CSVs include: RID, time_followup, status, age, sex, and signature values. Values were acquired based on ROC_curve_cross_validation.R using CU vs A+T+ MCI (n = 350), and A+T+ MCI vs A+T+ dementia (n = 265) datasets.

# LIMMA analysis
Expression proteomics matrix: 

- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from the ADNI database.
- Stratification: For analysis, only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408).
- Input CSV includes: RID and protein log2 abundance values of PPAV11 and core biomarkers per subject.

Clinical data: 

- Original file: "ADNIMERGE_21Apr2025" obtained from the ADNI database.
- Processing: protein abundance log2 transformation, inclusion of protein names, and diagnosis per subject. 
- Stratification: For analysis, only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408), and clinical variables of interest.
- Input CSV includes: RID,	Age,	Sex,	PTEDUCAT,	FDG,	PIB,	AV45,	APOE4,	ABETA,	TAU,	PTAU,	pTau/AB,	CDRSB,	ADAS11,	ADAS13,	ADASQ4,	MMSE,	RAVLT_immediate,	RAVLT_learning,	RAVLT_forgetting,	RAVLT_perc_forgetting,	LDELTOTAL,	DIGITSCOR,	TRABSCOR,	FAQ,	MOCA,	EcogPtMem,	EcogPtLang,	EcogPtVisspat,	EcogPtPlan,	EcogPtOrgan,	EcogPtDivatt,	EcogPtTotal,	EcogSPMem,	EcogSPLang,	EcogSPVisspat,	EcogSPPlan,	EcogSPOrgan,	EcogSPDivatt,	EcogSPTotal, Ventricles,	Hippocampus,	WholeBrain,	Entorhinal,	Fusiform,	MidTemp,	ICV.

# ROC curve cross-validation
Analysis on ADNI database:
- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from the ADNI database.
- Processing: protein abundance log2 transformation, inclusion of protein names, and diagnosis per subject.
- Stratification: Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/pTau ratio; n = 703), (2) A-T- vs A+T- (based on AB42, n = 388), (3) AV45- vs AV45+ (n = 488); biological (based on Aβ/pTau ratio) and clinical (stablished by ADNI based on clinical assessment) diagnosis (4) A-T- CU vs A+T+ MCI (n = 350), (5) A-T- CU vs A+T+ dementia (n = 265); and differential diagnosis (6) A-T- CU vs A-T- MCI (n = 276).
- Input CSVs include: RID, group (0-control or 1-group of interest), and protein log2 abundance values per subject.

Analysis on PPMI database:
- Original file: "ppmi_project_151.xlm" obtained from the PPMI database.
- Processing: Conversion for protein names into official gene names. 
- Stratification: Subjects with PD were those classified with sporadic PD and LRRK2; control subjects were the ones classified as Healthy Controls.
- Input CSV includes: RID, group (0-control or 1-PD), and protein log2 abundance values per subject.

# ROC curve independent training testing
Training 
Bader et al. (2020) cohort 
- Original file:
- Processing: Selection of the "Sweden" cohort. Processing performed in Metaboanalyst included quantile normalization and log2 transformation.

Tao et al. (2024) cohort 
- Original file: "Table S1 identified Protein groups in CSF.xlm". Obtained from supplemental information in 
- Processing: Processing performed in Metaboanalyst included quantile normalization and log2 transformation.

Testing
Bangs et al. (2025) cohort 
- Original file: "EmorySubCSF_483CleanDat.xlm". Obtained from https://www.synapse.org/ (Synapse ID: syn65461849).
- Stratification: Different stratifications were used for ROC performance. Based on biological diagnosis (1) A-T- vs A+T+ (based on Aβ/tTau ratio; n = 431); biological (based on Aβ/tTau ratio) and clinical (based on MoCA) diagnosis (2) A-T- CU vs A+T+ MCI (n = 242), (3) A-T- CU vs A+T+ dementia (n = 273); and differential diagnosis (5) A-T- CU vs A-T- MCI (n = 213).

Johnson et al. (2020) cohort

Input CSVs included the two training datasets and one validation dataset. Each document included ID, group (0-control or 1-group of interest), and protein log2 abundance values per subject. 
