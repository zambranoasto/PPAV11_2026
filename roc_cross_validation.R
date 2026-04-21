# ROC curve analysis of CSF proteomic biomarker signatures
# Cross-validation logistic regression model, optimal threshold determined by Youden index
# Input data described in DATA_README.md

# Import required libraries
library(pROC)
library(ggplot2)
library(caret)
library(dplyr)
library(tibble)
library(purrr)
library(RColorBrewer)

# Load and prepare dataset
data <- read.csv("Path/to/file.csv",
                 header = TRUE, stringsAsFactors = FALSE)
colnames(data) <- trimws(colnames(data))
data$Group <- factor(data$Group)
levels(data$Group) <- c("Control", "Case")

# Define biomarker signatures
signature_list <- list(
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

# Assign colors to signatures
signature_colors <- c(
  PPAV11 = "olivedrab3", PPA5 = "darkorchid4", Bader = "blue4", Wang = "darkorchid",
  Sathe = "aquamarine3", Tao = "maroon4", Shen = "lightpink3", Ali = "cyan3",
  VZ = "darkseagreen2", Liu = "plum3", Campo = "violet", Guo1 = "deeppink2",
  Guo2 = "gold1", Hou = "turquoise4", Yun = "salmon3"
)

# Filter genes present in dataset
filter_genes <- function(signature) intersect(signature, colnames(data))
signatures <- lapply(signature_list, filter_genes)
print(signatures)

# Define function to compute ROC curves
calculate_roc_cv <- function(genes, signature_name, k = 5, rep = 10) {
  signature_data <- data[, c(genes, "Group")]
  signature_data <- na.omit(signature_data)
  if (nrow(signature_data) < 7 || length(unique(signature_data$Group)) < 2) return(NULL)
  all_rocs <- list()
  aucs <- c()
  all_predictions <- list()
  all_betas <- list()
  sensitivities <- c()
  specificities <- c()
  set.seed(123)

  for (i in 1:rep) {
    folds <- createFolds(signature_data$Group, k = k, list = TRUE)
    for (j in 1:k) {
      test_idx <- folds[[j]]
      train <- signature_data[-test_idx, ]
      test <- signature_data[test_idx, ]
      if (length(unique(test$Group)) < 2) next
      model <- tryCatch(
        glm(Group ~ ., data = train, family = "binomial"),
        error = function(e) NULL
      )

      if (is.null(model)) next
      probs <- predict(model, newdata = test, type = "response")
      roc_cv <- roc(test$Group, probs,
                    levels = c("Control", "Case"),
                    direction = "<",
                    quiet = TRUE)
      if (is.null(roc_cv)) next
      aucs <- c(aucs, auc(roc_cv))
      all_rocs[[length(all_rocs) + 1]] <- roc_cv
      all_betas[[length(all_betas) + 1]] <- coef(model)
      opt <- coords(roc_cv,
                    x = "best",
                    best.method = "youden",
                    transpose = FALSE)
      sensitivities <- c(sensitivities, opt$sensitivity)
      specificities <- c(specificities, opt$specificity)
      preds_tmp <- tibble(
        PredictedProbability = probs,
        TrueDiagnosis = test$Group,
        Fold = paste0("Rep", i, "_Fold", j)
      )
      all_predictions[[length(all_predictions) + 1]] <- preds_tmp
    }
  }

  if (length(all_rocs) < 2) return(NULL)
  fpr_grid <- seq(0, 1, length.out = 100)
  interp_tprs <- sapply(all_rocs, function(r) {
    x <- 1 - r$specificities
    y <- r$sensitivities
    ord <- order(x)
    approx(x[ord], y[ord],
           xout = fpr_grid,
           rule = 2)$y
  })

  roc_df <- tibble(
    FPR = fpr_grid,
    TPR = rowMeans(interp_tprs),
    TPR_Low = apply(interp_tprs, 1, function(x) quantile(x, 0.025)),
    TPR_High = apply(interp_tprs, 1, function(x) quantile(x, 0.975)),
    Signature = signature_name
  )

  betas_mean <- Reduce("+", all_betas) / length(all_betas)
  return(list(
    roc_df = roc_df,
    auc = mean(aucs),
    auc_ci_low = quantile(aucs, 0.025),
    auc_ci_high = quantile(aucs, 0.975),
    predictions = bind_rows(all_predictions),
    betas = betas_mean[-1],
    intercept = betas_mean[1],
    sensitivity = mean(sensitivities),
    specificity = mean(specificities),
    roc_obj = all_rocs[[1]]
  ))
}

# Compute ROC results 
results <- map2(signatures, names(signatures), calculate_roc_cv)
results <- compact(results)

# Combine ROC curves
roc_df <- bind_rows(map(results, "roc_df"))
auc_vec <- map_dbl(results, "auc")
ordered_signatures <- names(sort(auc_vec, decreasing = TRUE))
roc_df$Signature <- factor(roc_df$Signature, levels = ordered_signatures)

# Define labels for plots
signature_labels <- sapply(ordered_signatures, function(signature) {
  res <- results[[signature]]
  auc_val <- round(res$auc, 3)
  ci_low <- round(res$auc_ci_low, 3)
  ci_high <- round(res$auc_ci_high, 3)
  sens <- round(res$sensitivity, 3)
  spec <- round(res$specificity, 3)
  paste0(signature, " (", auc_val, ", ", ci_low, "-", ci_high, ", ", sens, ", ", spec, ")")
})

# Generate ROC plots
ggplot(roc_df, aes(x = FPR, y = TPR, color = Signature)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = signature_colors[ordered_signatures],
                     labels = signature_labels) +
  scale_fill_manual(values = signature_colors[ordered_signatures], guide = "none") +
  labs(title = "Name",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Signature (AUC, CI 95%, Sensitivity, Specificity)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )
