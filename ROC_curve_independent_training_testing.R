# ROC curve analysis of CSF proteomic biomarker signatures
# Independent training and validation, logistic regression model, optimal threshold determined by Youden index
# Input data described in Data README.md

# Import required libraries
library(pROC)
library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)
library(caret)

# Load and prepare datasets
train_data1 <- read.csv("Path/to/file.csv",
                        header = TRUE, stringsAsFactors = FALSE)
train_data2 <- read.csv("Path/to/file.csv",
                        header = TRUE, stringsAsFactors = FALSE)
validation_data <- read.csv("Path/to/file.csv",
                            header = TRUE, stringsAsFactors = FALSE)
colnames(train_data1) <- trimws(colnames(train_data1))
colnames(train_data2) <- trimws(colnames(train_data2))
colnames(validation_data) <- trimws(colnames(validation_data))
convert_to_numeric <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ suppressWarnings(as.numeric(.))))
}

train_data1 <- convert_to_numeric(train_data1)
train_data2 <- convert_to_numeric(train_data2)
train_data <- bind_rows(train_data1, train_data2)
train_data$Group <- factor(train_data$Group)
validation_data$Group <- factor(validation_data$Group)

# Define biomarker signatures
signatures_list <- list(
  PPAV11 = c("CHI3L1","CYCS","DDAH1","GDA","LRRC4B","NPTX2","PKM","SMOC1","SPON1","YWHAG","YWHAZ"),
  Bader = c("MAPT","PKM","YWHAZ","ALDOC","IMPA1"),
  Wang = c("SMOC1","MAPT","GFAP","SUCLG2","PRDX3","NTN1"),
  Sathe = c("NPTX2","PKM","YWHAG"),
  Tao = c("PLTP","C16orf89","CHRD","FRZB","MGP","ART3","COL18A1","PCDHGC5","CXCL16","B3GALNT1","SHBG","SCRG1","MGAT2","GALNT7","MAN2A2","GM2A","IGHM","MAN1C1","PCSK1N"),
  Shen = c("NPTX2","SMOC1","GFAP","SMOC2","PEA15","TNFRSF1B"),
  Ali = c("TMOD2","SMOC1","YWHAG","TMED2","LRRN1","HHIP","CD248","PPP3R1","PPP1R1A","MIF"),
  VZ = c("ALDOA","PKM","LDHB"),
  Liu = c("AT1B1","SRGN","PRDX3"),
  Campo = c("MMP10","ABL1","SDC4","ITGB2","CLEC5A","TREM1","SPON2","THBD"),
  Guo1 = c("ACHE","YWHAG","PCSK1","MMP10","IRF1"),
  Guo2 = c("YWHAG","SMOC1","PIGR","TMOD2"),
  Hou = c("YWHAZ","LDHA","PKM","CHI3L1","BASP1","SMOC1","ENO2","PRDX1","VSTM2A","F2","SOD1","FN1"),
  Yun  = c("YWHAZ","EPHA5","CABIN1","SST","PPIA","CNTN5","GFAP","POMC","RSPO1","TPT1","MINDY2","PRDX6")
)

# Identify common genes between training and validation datasets 
filter_common_genes <- function(signature, train_data, validation_data) {
  intersect(intersect(signature, colnames(train_data)), colnames(validation_data))
}
signatures_train <- lapply(signatures_list,
                           filter_common_genes,
                           train_data = train_data,
                           validation_data = validation_data)
print(signatures_train)

# Assign colors to signatures
signature_colors <- c(
  PPAV11 = "olivedrab3", Bader = "blue4", Wang = "darkorchid",
  Sathe = "aquamarine3", Tao = "maroon4", Shen = "lightpink3",
  Ali = "cyan3", VZ = "darkseagreen2", Liu = "plum3",
  Campo = "violet", Guo1 = "deeppink2", Guo2 = "gold1",
  Hou = "turquoise4", Yun = "salmon3"
)

# Define function to compute ROC curves
compute_roc_train_validate <- function(genes, signature_name, train_data, validation_data) {
  genes_common <- intersect(genes, colnames(train_data))
  genes_common <- intersect(genes_common, colnames(validation_data))

  if(length(genes_common) == 0) {
    warning(paste("Signature", signature_name, "skipped: no common genes"))
    return(NULL)
  }

  train_subset <- train_data[, c(genes_common, "Group")]
  val_subset <- validation_data[, c(genes_common, "Group")]
  train_genes <- train_subset[, genes_common, drop = FALSE]
  val_genes <- val_subset[, genes_common, drop = FALSE]
  train_subset <- train_subset[rowSums(is.na(train_genes)) == 0, ]
  val_subset <- val_subset[rowSums(is.na(val_genes)) == 0, ]

  if(length(unique(train_subset$Group)) < 2 ||
     length(unique(val_subset$Group)) < 2) {
    warning(paste("Signature", signature_name, "skipped: insufficient classes"))
    return(NULL)
  }
  model <- tryCatch(
    glm(Group ~ ., data = train_subset, family = binomial),
    error = function(e) NULL
  )

  if(is.null(model)) return(NULL)
  probabilities <- predict(model,
                           newdata = val_subset,
                           type = "response")
  roc_obj <- roc(val_subset$Group, probabilities, quiet = TRUE)
  best_threshold <- as.numeric(
    coords(roc_obj, "best", ret = "threshold", best.method = "youden")
  )

  predicted_class <- factor(
    ifelse(probabilities >= best_threshold,
           levels(val_subset$Group)[2],
           levels(val_subset$Group)[1]),
    levels = levels(val_subset$Group)
  )

  predictions <- tibble(
    Probability = probabilities,
    TrueDiagnosis = val_subset$Group,
    PredictedDiagnosis = predicted_class
  )

  conf_mat <- confusionMatrix(predictions$PredictedDiagnosis,
                              predictions$TrueDiagnosis,
                              positive = levels(val_subset$Group)[2])
  sensitivity <- conf_mat$byClass["Sensitivity"]
  specificity <- conf_mat$byClass["Specificity"]
  fpr_grid <- seq(0,1,length.out = 100)
  tpr_interp <- approx(
    1 - roc_obj$specificities,
    roc_obj$sensitivities,
    xout = fpr_grid,
    rule = 2
  )$y
  
  return(list(
    roc_df = tibble(FPR = fpr_grid,
                    TPR = tpr_interp,
                    Signature = signature_name),
    auc = auc(roc_obj),
    roc_obj = roc_obj,
    preds = predictions,
    best_threshold = best_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    betas = coef(model)[-1],
    intercept = coef(model)[1]
  ))
}

# Compute ROC results
results <- map2(signatures_train,
                names(signatures_train),
                ~ compute_roc_train_validate(.x,.y,train_data,validation_data))
names(results) <- names(signatures_train)
results <- compact(results)

# Prepare data for plot
roc_data <- bind_rows(map(results,"roc_df"))
auc_vec <- map_dbl(results,"auc")
signature_order <- names(sort(auc_vec,decreasing = TRUE))
roc_data$Signature <- factor(roc_data$Signature,
                             levels = signature_order)

# Define labels for plots
signature_labels <- sapply(signature_order,function(sig){
  roc_obj <- results[[sig]]$roc_obj
  ci_auc <- ci.auc(roc_obj)
  auc_val <- round(auc(roc_obj),3)
  ci_low <- round(ci_auc[1],3)
  ci_high <- round(ci_auc[3],3)
  sens <- round(results[[sig]]$sensitivity,3)
  spec <- round(results[[sig]]$specificity,3)
  paste0(sig," (",auc_val,", ",ci_low,"-",ci_high,", ",sens,", ",spec,")")
})

# Generate ROC plots
ggplot() +
  geom_line(data = roc_data,
            aes(x = FPR,
                y = TPR,
                color = Signature),
            size = 1) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed",
              color = "gray") +
  scale_color_manual(values = signature_colors[signature_order],
                     labels = signature_labels) +
  labs(title = "Name",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Signature (AUC, CI 95%, Sensitivity, Specificity)") +
  theme_minimal()
preds_all <- map_dfr(results,
                     function(res){
                       if(is.null(res$preds)) return(NULL)
                       res$preds
                     },
                     .id = "Signature")
preds_all$Signature <- factor(preds_all$Signature,
                              levels = signature_order)
