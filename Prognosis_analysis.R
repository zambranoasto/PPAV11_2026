# Prognosis analysis of CSF proteomic biomarker signatures (Hazard models)
# Data for input is described in Data README.md

# All signatures comparison
# Required libraries
library(survival)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)

# Biomarker signatures
signatures <- c("score_PPAV11", "score_Bader", "score_Wang", "score_Sathe",
                "score_Tao", "score_Shen", "score_Ali", "score_VZ",
                "score_Liu", "score_Campo", "score_Guo1", "score_Guo2",
                "score_Hou", "score_Yun", "pTauAB", "tTauAB", "ABETA", "TAU", "PTAU")
risk_logic <- tibble(
  signature = signatures,
  clinical_risk = c("Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low",
                    "Low", "Low", "Low", "Low", "High", "High", "Low", "High", "High")  # ABETA corrected to "Low"
)

# Load data
dat <- read_csv("Path/to/file.csv") 
dat <- dat %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")))

# Cox binary function 
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
  
  # Print to console 
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

# Order by HR 
clinical_results <- clinical_results %>%
  arrange(desc(HR_risk_group_vs_reference)) %>%
  mutate(signature = factor(signature, levels = signature))

# Forest plot 
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


# One signature Kaplan-Meier                              
# Required libraries
library(survival)
library(survminer)
library(dplyr)
library(readr)

# Load data
dat <- read_csv("C:/Users/Username/OneDrive/Name.csv")

# Transformations 
dat <- dat %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    score_PPAV11_bin = ifelse(score_PPAV11 <= median(score_PPAV11, na.rm = TRUE), "Low", "High"),
    score_PPAV11_bin = factor(score_PPAV11_bin, levels = c("Low", "High"))
  ) %>%
  drop_na(score_PPAV11_bin, time_followup, status, age, sex)

# Cox model 
cox_ppav11 <- coxph(Surv(time_followup, status) ~ score_PPAV11_bin + age + sex, data = dat)
summary(cox_ppav11)

# Kaplan-Meier curve
fit_km <- survfit(Surv(time_followup, status) ~ score_PPAV11_bin, data = dat)

# Plot
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

