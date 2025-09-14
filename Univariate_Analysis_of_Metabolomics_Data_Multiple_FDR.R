# Load required libraries
library(dplyr)
library(purrr)

# Install and load additional packages for multiple testing
if (!require("qvalue")) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("qvalue")
  library(qvalue)
}

# Load the data
data7_blend_normalized <- read.csv("discharge.csv", check.names = FALSE)

# Ensure 'Sample' and 'Label' columns are present
if (!all(c("Sample", "Label") %in% colnames(data7_blend_normalized))) {
  stop("The 'Sample' and 'Label' columns must be present in the data.")
}

# Trim whitespace from Label column values
data7_blend_normalized$Label <- trimws(data7_blend_normalized$Label)

# Define the desired order of labels
desired_order <- c('Discharge_Good', 'Discharge_Poor')

# Remove unwanted labels if necessary
mismatched_labels <- setdiff(unique(data7_blend_normalized$Label), desired_order)
if (length(mismatched_labels) > 0) {
  stop(paste("Mismatched labels found:", paste(mismatched_labels, collapse = ", ")))
}

# Variance filtering using MAD/median threshold of 25%
var_filter <- function(df, threshold = 0.25) {
  rsd <- apply(df, 2, function(x) {
    mad(x, na.rm = TRUE) / median(x, na.rm = TRUE)
  })
  df <- df[, rsd >= threshold]
  return(df)
}

# Apply variance filtering
filtered_df <- data7_blend_normalized %>%
  select(-c(Sample, Label)) %>%
  var_filter()

# Add back Sample and Label columns
filtered_df <- data7_blend_normalized %>%
  select(Sample, Label) %>%
  bind_cols(filtered_df)

# Log transformation (log10 + 1)
log_transformed_df <- filtered_df %>%
  mutate(across(-c(Sample, Label), ~ log10(. + 1)))

# Extract group labels
Group <- log_transformed_df$Label

# Perform Wilcoxon rank sum test
wilcoxon_results <- lapply(log_transformed_df[, -c(1, 2)], function(x) {
  if (var(x, na.rm = TRUE) == 0) {
    return(list(p.value = 1))
  }
  result <- tryCatch(
    wilcox.test(x ~ Group, exact = FALSE),
    warning = function(w) {
      message("Warning: ", conditionMessage(w))
      return(list(p.value = 1))
    },
    error = function(e) {
      message("Error: ", conditionMessage(e))
      return(list(p.value = 1))
    }
  )
  return(result)
})

# Extract p-values
p_values <- sapply(wilcoxon_results, function(x) x$p.value)

# Apply multiple testing corrections
bh_values <- p.adjust(p_values, method = "fdr")  # Benjamini-Hochberg
by_values <- p.adjust(p_values, method = "BY")   # Benjamini-Yekutieli

# Calculate Storey's q-values (less conservative)
tryCatch({
  qobj <- qvalue(p_values)
  q_values <- qobj$qvalues
}, error = function(e) {
  message("Q-value calculation failed: ", conditionMessage(e))
  q_values <- rep(NA, length(p_values))
  names(q_values) <- names(p_values)
})

# Get nominally significant metabolites (p <= 0.05)
nominally_significant <- names(p_values)[p_values <= 0.05]

if (length(nominally_significant) == 0) {
  stop("No nominally significant metabolites found (p <= 0.05).")
}

# Function to calculate fold change
calculate_fold_change <- function(df, metabolite, group_col = "Label", group1 = "Discharge_Good", group2 = "Discharge_Poor") {
  group1_values <- df %>% filter(!!sym(group_col) == group1) %>% pull(!!sym(metabolite))
  group2_values <- df %>% filter(!!sym(group_col) == group2) %>% pull(!!sym(metabolite))
  
  mean_group1 <- mean(group1_values, na.rm = TRUE)
  mean_group2 <- mean(group2_values, na.rm = TRUE)
  
  fold_change <- 10^(mean_group2 - mean_group1)  # Convert back from log scale
  return(fold_change)
}

# Generate results for all nominally significant metabolites
results <- map_dfr(nominally_significant, function(metabolite) {
  discharge_good_values <- log_transformed_df %>%
    filter(Label == "Discharge_Good") %>%
    pull(!!sym(metabolite))
  
  discharge_poor_values <- log_transformed_df %>%
    filter(Label == "Discharge_Poor") %>%
    pull(!!sym(metabolite))
  
  mean_good <- mean(discharge_good_values, na.rm = TRUE)
  mean_poor <- mean(discharge_poor_values, na.rm = TRUE)
  fold_change <- calculate_fold_change(log_transformed_df, metabolite)
  
  tibble(
    Metabolite = metabolite,
    P_Value = p_values[metabolite],
    BH_Adjusted_P_Value = bh_values[metabolite],
    BY_Adjusted_P_Value = by_values[metabolite],
    Storey_Q_Value = if(exists("q_values")) q_values[metabolite] else NA,
    Discharge_Good = mean_good,
    Discharge_Poor = mean_poor,
    Fold_Change = fold_change
  )
})

# Save final results
write.csv(results, "significant_metabolites_results.csv", row.names = FALSE)
print("Results saved to significant_metabolites_results.csv")

# Print summary of results
print(paste("Number of nominally significant metabolites:", nrow(results)))
print(paste("Number with BH FDR < 0.05:", sum(results$BH_Adjusted_P_Value < 0.05, na.rm = TRUE)))
print(paste("Number with Storey q-value < 0.05:", sum(results$Storey_Q_Value < 0.05, na.rm = TRUE)))

