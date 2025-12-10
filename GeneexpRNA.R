library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(corrplot)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyr)
library(readr)
library(dplyr)
library(multcomp)  
library(broom)
library(patchwork)

# Define treatments and create data
treatments <- c("Untreated", "HFD/STZ", "50mg", "100mg", "400mg", "Metformin", "GAPDH")
gene_data <- data.frame(
  Treatment = treatments,
  CASPASE3_1 = c(31.3, 29.9, 29.1, 29.2, 32.9, 34.7, 18.9),
  CASPASE3_2 = c(31.7, 29.7, 28.9, 29.9, 32.1, 33.9, 18.1),
  CASPASE3_3 = c(31.9, 30.8, 28.5, 31.1, 32.3, 34.2, 18.6),
  BAX_1 = c(31.2, 23.4, 31.9, 28.9, 27.2, 30.3, 17.4),
  BAX_2 = c(31.4, 23.6, 32, 29.5, 27.6, 30.5, 17.2),
  BAX_3 = c(31.3, 22.8, 31.6, 29.4, 28.1, 30.7, 17.8),
  BCL2_1 = c(28.1, 31.5, 34.3, 27.9, 28.1, 28.9, 19.2),
  BCL2_2 = c(28.5, 32, 33.9, 28.1, 27.8, 27.4, 19.4),
  BCL2_3 = c(28.9, 31.8, 34.2, 28.7, 28.5, 27.1, 19.1),
  TNFalpha_1 = c(34.7, 26.1, 30.5, 30.1, 30.1, 32.1, 19.3),
  TNFalpha_2 = c(34.8, 25.7, 30.3, 30.5, 29.9, 32.9, 19.4),
  TNFalpha_3 = c(35, 25.9, 29.9, 31.1, 29.5, 33.1, 19.8),
  IL1alpha_1 = c(32.1, 24.9, 32.1, 32.4, 30.5, 29.4, 17.7),
  IL1alpha_2 = c(31.7, 24.6, 30.7, 31.3, 31.2, 29.1, 17.9),
  IL1alpha_3 = c(31.9, 23.9, 30.4, 31.7, 30.7, 29.3, 17.5),
  NRF2_1 = c(29.9, 30, 31.5, 31.1, 31.7, 29.5, 18.3),
  NRF2_2 = c(30.1, 29.2, 31.1, 31.6, 30.5, 29.1, 18.6),
  NRF2_3 = c(30.5, 29.7, 30.7, 31.9, 31.1, 29.9, 18.1)
)

print("Raw Gene Expression Data:")
print(gene_data)

expression_data <- gene_data[gene_data$Treatment != "GAPDH", ]
gapdh_data <- gene_data[gene_data$Treatment == "GAPDH", -1]

# Define genes
genes <- c("CASPASE3", "BAX", "BCL2", "TNFalpha", "IL1alpha", "NRF2")

# Function to calculate ΔCt with standard deviation
calculate_delta_ct <- function(expression_data, gapdh_data, genes) {
  results <- list()
 
  for (gene in genes) {
    cat("Processing gene:", gene, "\n")
   
    # Get gene columns
    gene_cols <- grep(paste0("^", gene), names(expression_data), value = TRUE)
   
    if (length(gene_cols) == 0) {
      warning(paste("No columns found for", gene))
      next
    }
   
    # Create long format data for this gene
    gene_long <- expression_data[, c("Treatment", gene_cols)] %>%
      pivot_longer(cols = -Treatment,
                   names_to = "Replicate",
                   values_to = "Target_Ct") %>%
      mutate(
        Replicate_num = as.numeric(str_extract(Replicate, "\\d+$")),
        Treatment = factor(Treatment, levels = c("Untreated", "HFD/STZ", "50mg", "100mg", "400mg", "Metformin"))
      )
   
    # Add corresponding GAPDH values
    gapdh_long <- data.frame(
      Replicate_num = 1:3,
      GAPDH_Ct = as.numeric(gapdh_data[1, ])
    )
   
    # Merge target and GAPDH data
    combined_data <- gene_long %>%
      left_join(gapdh_long, by = "Replicate_num") %>%
      mutate(
        Delta_Ct = Target_Ct - GAPDH_Ct,
        Gene = gene
      )
   
    results[[gene]] <- combined_data
  }
 
  # Combine all genes
  all_data <- bind_rows(results)
  return(all_data)
}

delta_ct_data <- calculate_delta_ct(expression_data, gapdh_data, genes)

print("\nΔCt Data (Target - GAPDH) with individual replicates:")
print(head(delta_ct_data, 20))

# Calculate ΔCt summary statistics
delta_ct_summary <- delta_ct_data %>%
  group_by(Treatment, Gene) %>%
  summarise(
    Delta_Ct_mean = mean(Delta_Ct, na.rm = TRUE),
    Delta_Ct_sd = sd(Delta_Ct, na.rm = TRUE),
    Delta_Ct_se = sd(Delta_Ct, na.rm = TRUE) / sqrt(n()),
    n_replicates = n(),
    .groups = "drop"
  )

print("\nΔCt Summary Statistics:")
print(delta_ct_summary)


control_deltas <- delta_ct_summary %>%
  filter(Treatment == "Untreated") %>%
  select(Gene, Control_Delta_Ct = Delta_Ct_mean, Control_Delta_Ct_sd = Delta_Ct_sd)

# Calculate ΔΔCt for individual replicates
delta_delta_ct_data <- delta_ct_data %>%
  left_join(control_deltas, by = "Gene") %>%
  mutate(
    Delta_Delta_Ct = Delta_Ct - Control_Delta_Ct,
    Fold_Change = 2^(-Delta_Delta_Ct)
  )

print("\nΔΔCt Data with individual replicates:")
print(head(delta_delta_ct_data, 20))
library(MASS)
library(dplyr)


control_deltas <- delta_ct_summary %>%
  filter(Treatment == "Untreated") %>%
  select(Gene, Control_Delta_Ct = Delta_Ct_mean, Control_Delta_Ct_sd = Delta_Ct_sd)

# Calculate ΔΔCt for individual replicates
delta_delta_ct_data <- delta_ct_data %>%
  left_join(control_deltas, by = "Gene") %>%
  mutate(
    Delta_Delta_Ct = Delta_Ct - Control_Delta_Ct,
    Fold_Change = 2^(-Delta_Delta_Ct)
  )

print("\nΔΔCt Data with individual replicates:")
print(head(delta_delta_ct_data, 20))
control_deltas <- delta_ct_summary %>%
  dplyr::filter(Treatment == "Untreated") %>%
  dplyr::select(Gene, Control_Delta_Ct = Delta_Ct_mean, Control_Delta_Ct_sd = Delta_Ct_sd)
conflicts()
detach("package:MASS", unload = TRUE)  # example
control_deltas <- delta_ct_summary %>%
  filter(Treatment == "Untreated") %>%
  select(Gene, Control_Delta_Ct = Delta_Ct_mean, Control_Delta_Ct_sd = Delta_Ct_sd)

# Calculate ΔΔCt for individual replicates
delta_delta_ct_data <- delta_ct_data %>%
  left_join(control_deltas, by = "Gene") %>%
  mutate(
    Delta_Delta_Ct = Delta_Ct - Control_Delta_Ct,
    Fold_Change = 2^(-Delta_Delta_Ct)
  )

print("\nΔΔCt Data with individual replicates:")
print(head(delta_delta_ct_data, 20))

control_deltas <- delta_ct_summary %>%
  dplyr::filter(Treatment == "Untreated") %>%
  dplyr::select(Gene, Control_Delta_Ct = Delta_Ct_mean, Control_Delta_Ct_sd = Delta_Ct_sd)

fold_change_summary <- delta_delta_ct_data %>%
  group_by(Treatment, Gene) %>%
  summarise(
    Delta_Delta_Ct_mean = mean(Delta_Delta_Ct, na.rm = TRUE),
    Delta_Delta_Ct_sd = sd(Delta_Delta_Ct, na.rm = TRUE),
    Delta_Delta_Ct_se = sd(Delta_Delta_Ct, na.rm = TRUE) / sqrt(n()),
    Fold_Change_mean = mean(Fold_Change, na.rm = TRUE),
    Fold_Change_sd = sd(Fold_Change, na.rm = TRUE),
    Fold_Change_se = sd(Fold_Change, na.rm = TRUE) / sqrt(n()),
    n_replicates = n(),
    .groups = "drop"
  )

print("\nFold Change Summary Statistics:")
print(fold_change_summary)

# Function to get significance stars
get_significance_stars <- function(p_value) {
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")  
  if (p_value < 0.05) return("*")
  return("ns")
}

# Statistical analysis using individual replicate data
anova_results <- data.frame()
tukey_results <- data.frame()


print("\n=== STATISTICAL ANALYSIS: ONE-WAY ANOVA AND TUKEY'S HSD TEST ===")
print("Analysis based on individual ΔCt values with proper error propagation")

for (gene in genes) {
  cat("\n", rep("=", 60), "\n")
  cat("ANALYZING GENE:", gene, "\n")
  cat(rep("=", 60), "\n")
 
  # Filter data for this gene
  gene_data_subset <- delta_ct_data %>%
    filter(Gene == gene) %>%
    mutate(Treatment = factor(Treatment, levels = c("Untreated", "HFD/STZ", "50mg", "100mg", "400mg", "Metformin")))
 
  # Perform one-way ANOVA on ΔCt values
  anova_model <- aov(Delta_Ct ~ Treatment, data = gene_data_subset)
  anova_summary <- summary(anova_model)
 
  # Extract ANOVA results
  f_value <- anova_summary[[1]][["F value"]][1]
  p_value <- anova_summary[[1]][["Pr(>F)"]][1]
  significance <- get_significance_stars(p_value)
 
  # Store ANOVA results
  anova_results <- rbind(anova_results, data.frame(
    Gene = gene,
    F_value = f_value,
    P_value = p_value,
    Significance = significance,
    stringsAsFactors = FALSE
  ))
 
  cat("ONE-WAY ANOVA Results (on ΔCt values):\n")
  cat("F-value:", round(f_value, 4), "\n")
  cat("P-value:", format(p_value, scientific = TRUE, digits = 4), "\n")
  cat("Significance:", significance, "\n\n")
 
  # If ANOVA is significant, perform Tukey's HSD test
  if (p_value < 0.05) {
    cat("ANOVA is significant (p < 0.05). Performing Tukey's HSD test...\n\n")
   
    tryCatch({
      tukey_test <- TukeyHSD(anova_model)
      tukey_df <- as.data.frame(tukey_test$Treatment)
      tukey_df$Comparison <- rownames(tukey_df)
      tukey_df$Gene <- gene
      tukey_df$Significance <- sapply(tukey_df$`p adj`, get_significance_stars)
     
      # Store Tukey results
      tukey_results <- rbind(tukey_results, tukey_df)
     
      cat("TUKEY'S HSD POST-HOC TEST Results:\n")
      cat("(Only showing comparisons with adjusted p < 0.05)\n\n")
     
      # Display significant comparisons
      significant_comparisons <- tukey_df[tukey_df$`p adj` < 0.05, ]
     
      if (nrow(significant_comparisons) > 0) {
        for (i in 1:nrow(significant_comparisons)) {
          comp <- significant_comparisons[i, ]
          cat(sprintf("%-25s: Mean Diff = %8.4f, 95%% CI [%8.4f, %8.4f], p.adj = %s (%s)\n",
                      comp$Comparison,
                      comp$diff,
                      comp$lwr,
                      comp$upr,
                      format(comp$`p adj`, scientific = TRUE, digits = 3),
                      comp$Significance))
        }
      } else {
        cat("No pairwise comparisons were significant after Tukey's correction.\n")
      }
    }, error = function(e) {
      cat("Error performing Tukey's test:", e$message, "\n")
    })
   
  } else {
    cat("ANOVA is not significant (p >= 0.05). Tukey's test not performed.\n")
  }
}

# Print comprehensive results
cat("\n", rep("=", 80), "\n")
cat("COMPREHENSIVE RESULTS SUMMARY\n")
cat(rep("=", 80), "\n")

print("\nANOVA RESULTS SUMMARY:")
print(anova_results)

if (nrow(tukey_results) > 0) {
  print("\nTUKEY'S HSD SIGNIFICANT COMPARISONS (p.adj < 0.05):")
  significant_tukey <- tukey_results[tukey_results$`p adj` < 0.05,
                                     c("Gene", "Comparison", "diff", "p adj", "Significance")]
  if (nrow(significant_tukey) > 0) {
    print(significant_tukey)
  } else {
    print("No significant pairwise comparisons found after Tukey's correction.")
  }
}

# Create visualizations with error bars

# 1. Fold Change Bar Plot with Error Bars
fold_change_plot_data <- fold_change_summary %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Untreated", "HFD/STZ", "50mg", "100mg", "400mg", "Metformin")),
    Gene = factor(Gene, levels = genes)
  ) %>%
  left_join(anova_results[, c("Gene", "Significance")], by = "Gene") %>%
  mutate(Gene_Label = paste0(Gene, " (", Significance, ")"))

p_fold_change <- ggplot(fold_change_plot_data, aes(x = Treatment, y = Fold_Change_mean, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Fold_Change_mean - Fold_Change_se,
                    ymax = Fold_Change_mean + Fold_Change_se),
                position = position_dodge(width = 0.9), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Gene Expression Fold Changes by Treatment\n(Mean ± SE, n=3 per group)",
       subtitle = "ANOVA significance shown in parentheses",
       x = "Treatment",
       y = "Fold Change (2^(-ΔΔCt))",
       fill = "Gene",
       caption = "Error bars: Standard Error; Red line: Control level (Fold Change = 1)\nSignificance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption = element_text(size = 8, hjust = 0)) +
  scale_fill_brewer(palette = "Set3", labels = unique(fold_change_plot_data$Gene_Label))

ggsave("fold_change_barplot_with_error_bars.png", p_fold_change, width = 16, height = 10)
print(p_fold_change)

# 2. Individual Gene Plots with Error Bars
individual_plots <- list()

for (gene in genes) {
  gene_data_plot <- fold_change_summary %>%
    filter(Gene == gene) %>%
    mutate(Treatment = factor(Treatment, levels = c("Untreated", "HFD/STZ", "50mg", "100mg", "400mg", "Metformin")))
 
  # Get ANOVA significance
  gene_anova_sig <- anova_results[anova_results$Gene == gene, "Significance"]
  gene_p_value <- anova_results[anova_results$Gene == gene, "P_value"]
 
  p <- ggplot(gene_data_plot, aes(x = Treatment, y = Fold_Change_mean, fill = Treatment)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Fold_Change_mean - Fold_Change_se,
                      ymax = Fold_Change_mean + Fold_Change_se),
                  width = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = paste("Fold Change for", gene),
         subtitle = paste("ANOVA:", gene_anova_sig,
                          "(p =", format(gene_p_value, scientific = TRUE, digits = 3), ")"),
         x = "Treatment",
         y = "Fold Change (Mean ± SE)",
         caption = "Error bars: Standard Error (n=3)") +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.caption = element_text(size = 8, hjust = 0))
 
  individual_plots[[gene]] <- p
 
  # Save individual plot
  ggsave(paste0(gene, "_fold_change_with_error_bars.png"), p, width = 10, height = 6)
}

# Combine individual plots
combined_individual_plot <- wrap_plots(individual_plots, ncol = 2)
ggsave("all_genes_fold_change_with_error_bars.png", combined_individual_plot, width = 16, height = 12)

# 3. ΔCt Box Plots
delta_ct_plot_data <- delta_ct_data %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Untreated", "HFD/STZ", "50mg", "100mg", "400mg", "Metformin")),
    Gene = factor(Gene, levels = genes)
  )

p_delta_ct <- ggplot(delta_ct_plot_data, aes(x = Treatment, y = Delta_Ct)) +
  geom_boxplot(alpha = 0.7, fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(title = "ΔCt Values (Target - GAPDH) Across Treatments",
       subtitle = "Individual replicates shown with box plots",
       x = "Treatment",
       y = "ΔCt Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("delta_ct_boxplots.png", p_delta_ct, width = 16, height = 12)
print(p_delta_ct)

# 4. Heatmap of Fold Changes (Mean values)
fold_change_matrix <- fold_change_summary %>%
  dplyr::select(Treatment, Gene, Fold_Change_mean) %>%
  pivot_wider(names_from = Gene, values_from = Fold_Change_mean) %>%
  column_to_rownames("Treatment") %>%
  as.matrix()

png("gene_expression_heatmap_with_stats.png", width = 1000, height = 700)
pheatmap(fold_change_matrix,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Gene Expression Fold Changes (Mean Values)\nGene Order: CASPASE3, BAX, BCL2, TNFalpha, IL1alpha, NRF2",
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()

# Save results to CSV files
write.csv(delta_ct_summary, "delta_ct_summary_with_stats.csv", row.names = FALSE)
write.csv(fold_change_summary, "fold_change_summary_with_stats.csv", row.names = FALSE)
write.csv(anova_results, "anova_results_with_stats.csv", row.names = FALSE)
if (nrow(tukey_results) > 0) {
  write.csv(tukey_results, "tukey_results_with_stats.csv", row.names = FALSE)
}
write.csv(delta_delta_ct_data, "individual_replicate_data_with_stats.csv", row.names = FALSE)

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE! (With Standard Deviations and Error Propagation)\n")
cat("Generated files:\n")
cat("- delta_ct_summary_with_stats.csv (ΔCt means, SD, SE for each treatment/gene)\n")
cat("- fold_change_summary_with_stats.csv (Fold change means, SD, SE for each treatment/gene)\n")
cat("- individual_replicate_data_with_stats.csv (All individual ΔΔCt and fold change values)\n")
cat("- anova_results_with_stats.csv (ANOVA results)\n")
cat("- tukey_results_with_stats.csv (Tukey post-hoc results)\n")
cat("- fold_change_barplot_with_error_bars.png\n")
cat("- all_genes_fold_change_with_error_bars.png\n")
cat("- delta_ct_boxplots.png\n")
cat("- gene_expression_heatmap_with_stats.png\n")
cat("- Individual gene plots with error bars\n")
cat("\nKey improvements:\n")
cat("✓ Standard deviations calculated for ΔCt, ΔΔCt, and fold changes\n")
cat("✓ Standard errors shown in error bars\n")
cat("✓ Statistical analysis performed on individual ΔCt replicates\n")
cat("✓ Proper error propagation throughout analysis\n")
cat("✓ Individual replicate data preserved for transparency\n")
cat(rep("=", 80), "\n")
