# Load required libraries
library(dplyr)
library(purrr)

# Load the data
data7_blend_normalized <- read.csv("discharge.csv", check.names = FALSE)
print("Data loaded:")
print(head(data7_blend_normalized))

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

# Get statistically significant metabolites (p <= 0.05)
significant_metabolites <- names(p_values)[p_values <= 0.05]

# Stop execution if no significant metabolites found
if (length(significant_metabolites) == 0) {
  stop("No significant metabolites found.")
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

# Use log-transformed data for fold change calculation
results <- map_dfr(significant_metabolites, function(metabolite) {
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
    P_Value = p_values[metabolite],  # Keep only raw p-values
    Discharge_Good = mean_good,
    Discharge_Poor = mean_poor,
    Fold_Change = fold_change
  )
})

# Save final results
write.csv(results, "significant_metabolites_results_no_adjustment.csv", row.names = FALSE)
print("Results saved to significant_metabolites_results_no_adjustment.csv")
