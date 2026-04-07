# Prognosis analysis
Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from the ADNI database.

After downloading the dataset and prior to analysis, the data were split into two separate CSV files to evaluate disease-stage transitions. One dataset included cognitively unimpaired (CU) individuals with longitudinal follow-up to analyze the transition from CU to MCI (n = 143). The second dataset included A+T+ individuals (defined using the Aβ/pTau ratio) with baseline MCI and longitudinal follow-up to evaluate the transition from MCI to dementia (n = 174). 

CSVs include: RID, time_followup, status, age, sex, and signature values. Values were acquired based on ROC_curve_cross_validation.R using CU vs A+T+ MCI (n = 350), and A+T+ MCI vs A+T+ dementia (n = 399) datasets 

# LIMMA analysis
Expression proteomics matrix

- Original file: "CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620-2.csv" obtained from the ADNI database.
- For analysis, only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408)
- CSV includes RID and protein values of PPAV11 and core biomarkers per subject 

Clinical data

- Original file: "ADNIMERGE_21Apr2025" obtained from the ADNI database. 
- For analysis, only A+T+ individuals (defined using the Aβ/pTau ratio) were selected (n = 408), and clinical variables of interest
- CSV includes: RID,	Age,	Sex,	PTEDUCAT,	FDG,	PIB,	AV45,	APOE4,	ABETA,	TAU,	PTAU,	pTau/AB,	CDRSB,	ADAS11,	ADAS13,	ADASQ4,	MMSE,	RAVLT_immediate,	RAVLT_learning,	RAVLT_forgetting,	RAVLT_perc_forgetting,	LDELTOTAL,	DIGITSCOR,	TRABSCOR,	FAQ,	MOCA,	EcogPtMem,	EcogPtLang,	EcogPtVisspat,	EcogPtPlan,	EcogPtOrgan,	EcogPtDivatt,	EcogPtTotal,	EcogSPMem,	EcogSPLang,	EcogSPVisspat,	EcogSPPlan,	EcogSPOrgan,	EcogSPDivatt,	EcogSPTotal, Ventricles,	Hippocampus,	WholeBrain,	Entorhinal,	Fusiform,	MidTemp,	ICV

# ROC curve cross-validation

# ROC curve independent training testing
