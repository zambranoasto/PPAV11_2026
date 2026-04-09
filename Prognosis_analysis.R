# Prognosis analysis of CSF proteomic biomarker signatures (Hazard models)
# Input data described in Data README.md

# Acquire signature scores
# Import required libraries
library(pROC)
library(caret)
library(dplyr)
library(tibble)
library(purrr)
library(readr)

# Load and prepare dataset
data <- read_csv("Path/to/file.csv",
                 col_types = cols())
colnames(data) <- trimws(colnames(data))
data$Group <- factor(data$Group)
levels(data$Group) <- c("0", "1")  # Positive level = "Case"

# Define biomarker signatures
signatures <- list(
  Bader = c("MAPT", "PKM", "YWHAZ", "ALDOC", "IMPA1"),
  Wang = c("SMOC1", "MAPT", "GFAP", "SUCLG2", "PRDX3", "NTN1"),
  PPAV11 = c("CHI3L1", "CYCS", "DDAH1", "GDA", "LRRC4B", "NPTX2", "PKM", "SMOC1", "SPON1", "YWHAG", "YWHAZ"),
  Sathe = c("NPTX2", "PKM", "YWHAG"),
  Tao = c("PLTP", "C16orf89", "CHRD", "FRZB", "MGP", "ART3", "COL18A1", "PCDHGC5", "CXCL16", "B3GALNT1", "SHBG", "SCRG1", "MGAT2", "GALNT7", "MAN2A2", "GM2A", "IGHM", "MAN1C1", "PCSK1N"),
  Shen = c("NPTX2", "SMOC1", "GFAP", "SMOC2", "PEA15", "TNFRSF1B"),
  Ali = c("TMOD2", "SMOC1", "YWHAG", "TMED2", "LRRN1", "HHIP", "CD248", "PPP3R1", "PPP1R1A", "MIF"),
  VZ = c("ALDOA", "PKM", "LDHB"),
  Liu = c("AT1B1", "SRGN", "PRDX3"),
  Campo = c("MMP10", "ABL1", "SDC4", "ITGB2", "CLEC5A", "TREM1", "SPON2", "THBD"),
  Guo1 = c("ACHE", "YWHAG", "PCSK1", "MMP10", "IRF1"),
  Guo2 = c("YWHAG", "SMOC1", "PIGR", "TMOD2"),
  Hou = c("YWHAZ", "LDHA", "PKM", "CHI3L1", "BASP1", "SMOC1", "ENO2", "PRDX1", "VSTM2A", "F2", "SOD1", "FN1"),
  Yun  = c("YWHAZ", "EPHA5", "CABIN1", "SST", "PPIA", "CNTN5", "GFAP", "POMC", "RSPO1", "TPT1", "MINDY2", "PRDX6")
)

# Define function to calculate scores 
calculate_scores_per_subject <- function(signature_data, k = 5, rep = 10) {
  scores <- list()
  set.seed(123)
  for (i in 1:rep) {
    folds <- createFolds(signature_data$Group, k = k, list = TRUE)
    for (j in 1:k) {
      test_idx <- folds[[j]]
      train <- signature_data[-test_idx, ]
      test <- signature_data[test_idx, ]
      
      if (length(unique(test$Group)) < 2) next
      model <- tryCatch(glm(Group ~ ., data = train, family = "binomial"), error = function(e) NULL)
      if (is.null(model)) next
      probs_control <- predict(model, newdata = test, type = "response")
      probs <- 1 - probs_control  # Probability of "Case"
      scores[[length(scores) + 1]] <- tibble(
        RID = test$RID,
        Group = test$Group,
        Probability = probs
      )
    }
  }
  bind_rows(scores)
}

signature_results <- list()

for (signature_name in names(signatures)) {
  genes <- intersect(signatures[[signature_name]], colnames(data))
  
  if (length(genes) == 0) {
    message("Signature ", signature_name, " has no genes in the dataset.")
    next
  }
  signature_data <- data %>%
    select(RID, Group, all_of(genes)) %>%
    na.omit()
  test_predictions <- calculate_scores_per_subject(signature_data)
  average_score <- test_predictions %>%
    group_by(RID, Group) %>%
    summarise(!!paste0("score_", signature_name) := mean(Probability), .groups = "drop")
  signature_results[[signature_name]] <- average_score
}

# Combine scores into a single dataframe
final_scores <- reduce(signature_results, full_join, by = c("RID", "Group"))

# Save combined scores to CSV
write_csv(final_scores, "output.csv")

# Prognosis analysis: all signatures comparison
# Import required libraries
library(survival)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)

# Define biomarker signatures
signatures <- c("score_PPAV11", "score_Bader", "score_Wang", "score_Sathe",
                "score_Tao", "score_Shen", "score_Ali", "score_VZ",
                "score_Liu", "score_Campo", "score_Guo1", "score_Guo2",
                "score_Hou", "score_Yun", "pTauAB", "tTauAB", "ABETA", "TAU", "PTAU")
risk_logic <- tibble(
  signature = signatures,
  clinical_risk = c("Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low",
                    "Low", "Low", "Low", "Low", "High", "High", "Low", "High", "High")
)

# Load and prepare dataset
dat <- read_csv("Path/to/file.csv") # Use the output csv from last step 
dat <- dat %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))


# Define Cox binary function
  fit_cox_clinical <- function(var) {
  cutoff <- median(dat[[var]], na.rm = TRUE)
  risk <- risk_logic %>% filter(signature == var) %>% pull(clinical_risk)
  dat_bin <- dat %>%
    mutate(signature_bin = ifelse(.data[[var]] <= cutoff, "Low", "High")) %>%
    mutate(signature_bin = factor(signature_bin, levels = if (risk == "Low") c("High", "Low") else c("Low", "High"))) %>%
    drop_na(signature_bin, time_followup, status, age, sex)
  if (nrow(dat_bin) == 0) return(NULL)
  fit <- coxph(Surv(time_followup, status) ~ signature_bin + age + sex, data = dat_bin)
  s <- summary(fit)
  hr <- s$coefficients[1, "exp(coef)"]
  lower <- s$conf.int[1, "lower .95"]
  upper <- s$conf.int[1, "upper .95"]
  pval <- s$coefficients[1, "Pr(>|z|)"]
  interpretation <- if (risk == "Low") {
    paste("Risk ↑ if", var, "is low")
  } else {
    paste("Risk ↑ if", var, "is high")
  }
  
  # Print results to console
  cat("\n==============================\n")
  cat("Signature:", var, "\n")
  cat("HR:", round(hr, 3), " [", round(lower, 3), "-", round(upper, 3), "]\n")
  cat("p-value:", signif(pval, 3), "\n")
  cat("Clinical interpretation:", interpretation, "\n")
  cat("==============================\n")
  data.frame(
    signature = var,
    HR_risk_group_vs_reference = hr,
    lower_CI = lower,
    upper_CI = upper,
    p_value = pval,
    interpretation = interpretation,
    clinical_risk = risk
  )
}

# Apply function to all signatures 
clinical_results <- lapply(signatures, fit_cox_clinical)
clinical_results <- do.call(rbind, clinical_results[!sapply(clinical_results, is.null)])

# Order results by HR
clinical_results <- clinical_results %>%
  arrange(desc(HR_risk_group_vs_reference)) %>%
  mutate(signature = factor(signature, levels = signature))

# Generate forest plot 
ggplot(clinical_results, aes(x = signature, y = HR_risk_group_vs_reference, color = clinical_risk)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Low" = "darkred", "High" = "darkblue")) +
  labs(
    title = "Hazard Ratio per signature according to clinical logic",
    subtitle = "HR > 1 indicates higher risk in clinically relevant group",
    y = "Hazard Ratio",
    x = "Signature",
    color = "Clinical risk"
  ) +
  geom_text(aes(label = paste0("p=", signif(p_value, 3))), hjust = -0.2, size = 4)


# Prognosis analysis: single signature (Kaplan-Meier)
# Import required libraries
library(survival)
library(survminer)
library(dplyr)
library(readr)

# Load and prepare dataset
dat <- read_csv("Path/to/file.csv")
dat <- dat %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    score_PPAV11_bin = ifelse(score_PPAV11 <= median(score_PPAV11, na.rm = TRUE), "Low", "High"),
    score_PPAV11_bin = factor(score_PPAV11_bin, levels = c("Low", "High"))
  ) %>%
  drop_na(score_PPAV11_bin, time_followup, status, age, sex)

# Fit Cox proportional hazards model
cox_ppav11 <- coxph(Surv(time_followup, status) ~ score_PPAV11_bin + age + sex, data = dat)
summary(cox_ppav11)

# Fit Kaplan-Meier survival curve                                  
fit_km <- survfit(Surv(time_followup, status) ~ score_PPAV11_bin, data = dat)

# Plot Kaplan-Meier curve 
ggsurvplot(
  fit_km,
  data = dat,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#D55E00", "#0072B2"),  # Low in red, High in blue
  legend.title = "PPAV11 group",
  legend.labs = c("Low (low values)", "High (high values)"),
  title = "Survival according to PPAV11 signature (Low vs High)",
  xlab = "Months of follow-up",
  ylab = "Survival probability",
  ggtheme = theme_minimal(base_size = 14),
  risk.table.height = 0.25
)

