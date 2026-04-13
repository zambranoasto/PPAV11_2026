# Limma correlation analysis
# Proteins abundance adjusted by age and sex
# Input data described in DATA_README.md

# Import required libraries
library(tidyverse)
library(limma)
library(pheatmap)
library(viridis)
library(dplyr)

# Load and prepare datasets
expr_matrix <- read.csv("Path/to/file.csv", 
                        row.names = 1, check.names = FALSE)
clinical_data <- read.csv("Path/to/file.csv",
                          stringsAsFactors = FALSE)
common_ids <- intersect(colnames(expr_matrix), clinical_data$RID)
expr_matrix <- expr_matrix[, common_ids]
clinical_data <- clinical_data %>%
  filter(RID %in% common_ids) %>%
  arrange(match(RID, common_ids))
# To scale abundance values for direct comparison across proteins, add the following line:
# expr_mat_scaled <- t(scale(t(expr_mat))
original_protein_order <- rownames(expr_matrix)

# Define clinical variables
clinical_vars <- c("FDG", "AV45", "ABETA", "TAU", "PTAU", "CDRSB", "ADAS13", "ADASQ4",
                   "MMSE", "RAVLT_immediate", "RAVLT_learning", "LDELTOTAL", "TRABSCOR",
                   "FAQ", "MOCA", "EcogSPTotal", "Ventricles", "Hippocampus", "WholeBrain",
                   "Entorhinal", "Fusiform", "MidTemp", "ICV")

limma_results <- list()

# Define function to run LIMMA 
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

# Run LIMMA analysis
for (trait in clinical_vars) {
  res <- run_limma_per_trait(trait)
  if (!is.null(res)) {
    limma_results[[trait]] <- res
  }
}

# Combine results into a single dataframe
all_results <- bind_rows(limma_results)

# Create  matrix for heatmap visualization
beta_matrix <- all_results %>%
  select(Protein, Trait, logFC) %>%
  pivot_wider(names_from = Trait, values_from = logFC) %>%
  column_to_rownames("Protein")

# Create label matrix for clinical variables
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

# Assign colors and breaks for heatmap
max_abs_beta <- max(abs(beta_matrix), na.rm = TRUE)
breaks <- seq(-max_abs_beta, max_abs_beta, length.out = 101)
colors <- rev(viridis(100, option = "D"))

# Generate heatmap plot
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
         main = "Name")
