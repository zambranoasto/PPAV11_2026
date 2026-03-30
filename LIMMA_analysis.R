# Limma correlation analysis, proteins abundance adjusted by age and sex

# Load packages
library(tidyverse)
library(limma)
library(pheatmap)
library(viridis)
library(dplyr)

# Load data
expr_matrix <- read.csv("C:\\Users\\Username\\OneDrive\\Name.csv", 
                        row.names = 1, check.names = FALSE)
clinical_data <- read.csv("C:\\Users\\Username\\OneDrive\\Name.csv",
                          stringsAsFactors = FALSE)

# Align subjects between both files
common_ids <- intersect(colnames(expr_matrix), clinical_data$RID)
expr_matrix <- expr_matrix[, common_ids]
clinical_data <- clinical_data %>%
  filter(RID %in% common_ids) %>%
  arrange(match(RID, common_ids))
original_protein_order <- rownames(expr_matrix)

# Define clinical variables to analyze
clinical_vars <- c("FDG", "AV45", "ABETA", "TAU", "PTAU", "CDRSB", "ADAS13", "ADASQ4",
                   "MMSE", "RAVLT_immediate", "RAVLT_learning", "LDELTOTAL", "TRABSCOR",
                   "FAQ", "MOCA", "EcogSPTotal", "Ventricles", "Hippocampus", "WholeBrain",
                   "Entorhinal", "Fusiform", "MidTemp", "ICV")

limma_results <- list()

# Function to run LIMMA per clinical variable 
run_limma_per_trait <- function(trait_name) {
  trait_values <- clinical_data[[trait_name]]
  valid_idx <- which(!is.na(trait_values))
  
  if (length(valid_idx) < 10) return(NULL)  
  expr_sub <- expr_matrix[, valid_idx]
  trait_sub <- trait_values[valid_idx]  
  design <- model.matrix(~ trait_sub + Age + Sex, 
                         data = clinical_data[valid_idx, ])
  colnames(design)[2] <- trait_name
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = trait_name, number = Inf)
  results <- results %>%
    rownames_to_column(var = "Protein") %>%
    mutate(Trait = trait_name)
  return(results)
}

# Run analysis for all clinical variables
for (trait in clinical_vars) {
  res <- run_limma_per_trait(trait)
  if (!is.null(res)) {
    limma_results[[trait]] <- res
  }
}

# Combine results
all_results <- bind_rows(limma_results)

# Create  matrix for heatmap
beta_matrix <- all_results %>%
  select(Protein, Trait, logFC) %>%
  pivot_wider(names_from = Trait, values_from = logFC) %>%
  column_to_rownames("Protein")

# Create label matrix 
label_matrix <- all_results %>%
  select(Protein, Trait, logFC, adj.P.Val) %>%
  mutate(label = ifelse(adj.P.Val < 0.05,
                        paste0(sprintf("%.4f", logFC), "*"),
                        sprintf("%.4f", logFC))) %>%
  select(Protein, Trait, label) %>%
  pivot_wider(names_from = Trait, values_from = label) %>%
  column_to_rownames("Protein")

# Reorder matrices according to original CSV order
beta_matrix <- beta_matrix[original_protein_order, , drop = FALSE]
label_matrix <- label_matrix[original_protein_order, , drop = FALSE]

# 1Define colors and breaks for heatmap
max_abs_beta <- max(abs(beta_matrix), na.rm = TRUE)
breaks <- seq(-max_abs_beta, max_abs_beta, length.out = 101)
colors <- rev(viridis(100, option = "D"))

# Plot heatmap
pheatmap(beta_matrix,
         color = colors,
         breaks = breaks,
         display_numbers = as.matrix(label_matrix),
         number_format = "%.4f",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 9,
         border_color = NA,
         number_color = "white",
         main = "LIMMA per clinical variable (adjusted for age and sex, absolute values)")
