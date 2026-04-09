# Principal component analysis (PCA) of PPAV11 
# Input data described in DATA_README.md

# Import required libraries
library(ggplot2)
library(ggfortify)

# Load and prepare dataset
data <- read.csv("Path/to/file.csv")
ID <- data[[1]]               
Group <- as.factor(data[[2]])   
protein_matrix <- data[, -c(1,2)] 
protein_matrix[is.infinite(as.matrix(protein_matrix))] <- NA
rows_complete <- complete.cases(protein_matrix)
protein_matrix_clean <- protein_matrix[rows_complete, ]
Group_clean <- Group[rows_complete]
ID_clean <- ID[rows_complete]

# Run PCA 
pca_result <- prcomp(protein_matrix_clean, scale. = FALSE)

# Create dataframe with PCA scores and metadata
plot_data <- data.frame(protein_matrix_clean)
plot_data$Group <- Group_clean

pca_scores <- as.data.frame(pca_result$x)
pca_scores$Group <- Group_clean
pca_scores$ID <- ID_clean

# Extract PCA loadings and scale 
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$Protein <- rownames(loadings)
scaling_factor <- max(abs(pca_scores$PC1), abs(pca_scores$PC2)) * 0.7
loadings$PC1 <- loadings$PC1 * scaling_factor
loadings$PC2 <- loadings$PC2 * scaling_factor

# Calculate variance 
var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Generate PCA plot
ggplot() +
  stat_ellipse(data = pca_scores, aes(x = PC1, y = PC2, fill = Group), 
               geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = Group), size = 2) +
  #geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
  #             arrow = arrow(length = unit(0.1, "cm")), color = "gray36") +
  #geom_text(data = loadings, aes(x = PC1, y = PC2, label = Protein),
  #          color = "black", size = 2.5, vjust = 1) +
  theme_minimal() +
  labs(title = "PCA PPAV11 Johnson",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  scale_color_manual(values = c("0" = "blue4", "1" = "maroon4")) +
  scale_fill_manual(values = c("0" = "blue4", "1" = "maroon4")) +
  theme(legend.title = element_blank())

# MANOVA statistical analysis to evaluate group separation
manova_result <- manova(cbind(PC1, PC2) ~ Group, data = pca_scores)
summary_manova <- summary(manova_result)
cat("Summary of MANOVA analysis (PC1 and PC2):\n")
print(summary_manova)
