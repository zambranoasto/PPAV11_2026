# Boxplot and statistical analysis of protein abundance og PPAV11
# Comparison across groups stratified by sex and ethnicity
# Input data described in DATA_README.md

# Import required libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# Load and prepare datasets
data <- read.csv("Path/to/file.csv")
data <- data %>%
  mutate(
    GroupSex = paste(Group, Sex, sep = "_"),
    GroupEthnicity = paste(Group, Ethnicity, sep = "_")
  )
data_long <- data %>%
  pivot_longer(
    cols = -c(Sample, Group, Sex, Ethnicity, GroupSex, GroupEthnicity),
    names_to = "Protein",
    values_to = "Abundance"
  )
data_long$Group <- factor(data_long$Group, levels = c("Control", "AD"))
data_long$GroupSex <- factor(data_long$GroupSex)
data_long$GroupEthnicity <- factor(data_long$GroupEthnicity)

# Assign color palette for plots
colors <- c("Control_Male" = "#a6cee3", "Control_Female" = "#1f78b4",
            "AD_Male" = "#b2df8a", "AD_Female" = "#33a02c")

# Perform ANOVA
anova_sex <- data_long %>%
  group_by(Protein) %>%
  anova_test(Abundance ~ GroupSex)

anova_ethnicity <- data_long %>%
  group_by(Protein) %>%
  anova_test(Abundance ~ GroupEthnicity)

# Perform Tukey post hoc comparisons 
stat_sex <- data_long %>%
  group_by(Protein) %>%
  tukey_hsd(Abundance ~ GroupSex) %>%
  add_xy_position(x = "GroupSex")

stat_ethnicity <- data_long %>%
  group_by(Protein) %>%
  tukey_hsd(Abundance ~ GroupEthnicity) %>%
  add_xy_position(x = "GroupEthnicity")

# Generate combined panel plots by Sex
plot_sex <- ggplot(data_long, aes(x = GroupSex, y = Abundance)) +
  geom_boxplot(aes(fill = GroupSex), width = 0.6, outlier.shape = NA, alpha = 0.7, color = "gray30") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
  geom_jitter(aes(color = GroupSex), width = 0.2, size = 0.7, alpha = 0.3, show.legend = FALSE) +
  facet_wrap(~Protein, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  stat_pvalue_manual(stat_sex, label = "p.adj.signif", tip.length = 0.01, size = 2.5) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.border = element_rect(color = "gray80", fill = NA)
  )

# Generate combined panel plots by Ethnicity
plot_ethnicity <- ggplot(data_long, aes(x = GroupEthnicity, y = Abundance)) +
  geom_boxplot(aes(fill = GroupEthnicity), width = 0.6, outlier.shape = NA, alpha = 0.7, color = "gray30") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
  geom_jitter(aes(color = GroupEthnicity), width = 0.2, size = 0.7, alpha = 0.3, show.legend = FALSE) +
  facet_wrap(~Protein, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(stat_ethnicity, label = "p.adj.signif", tip.length = 0.01, size = 2.5) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.border = element_rect(color = "gray80", fill = NA)
  )

# Generate boxplots
plot_sex
plot_ethnicity
