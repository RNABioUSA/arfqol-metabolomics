# Load Required Packages
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(rstanarm, ggplot2, bayesplot, viridis, pROC, caret, reshape2)

# Create output directory
output_dir <- "Model_parameters_0618"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Load Data
data <- read.csv('significant_narrowed.csv')
data <- data[, -1]  # Remove first column (Sample ID)
data$Label <- factor(data$Label, levels = c("Discharge_Poor", "Discharge_Good"))
if (any(is.na(data$Label))) stop("Error: Label column contains invalid values")

# Convert metabolite columns to numeric
data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], function(x) as.numeric(as.character(x)))

# Remove NA values
data <- na.omit(data)

# Log-transform & Auto-scale metabolite data
numeric_cols <- sapply(data[, 2:ncol(data)], is.numeric)
data[, 2:ncol(data)][, numeric_cols] <- lapply(data[, 2:ncol(data)][, numeric_cols], function(x) scale(log(x + 1)))

# Set seed for reproducibility
set.seed(123)

# Number of repeats for multiple times k-fold cross-validation
num_repeats <- 20
k_folds <- 5

# Initialize storage for aggregated metrics
auc_values <- c()
mcc_values <- c()
confusion_matrices <- list()
all_true_labels <- c()
all_predicted_probs <- c()

# Perform 20x5-Fold Cross-Validation
for (rep_num in 1:num_repeats) {
  cat("Processing Repeat", rep_num, "...\n")
  folds <- createFolds(data$Label, k = k_folds, list = TRUE)
  
  for (i in 1:k_folds) {
    cat("Processing Fold", i, "in Repeat", rep_num, "...\n")
    
    test_indexes <- folds[[i]]
    train_data <- data[-test_indexes, ]
    test_data <- data[test_indexes, ]
    
    stopifnot(all(!rownames(test_data) %in% rownames(train_data)))
    
    # Train Bayesian logistic regression model
    bayes_model <- stan_glm(Label ~ a + b + d + g + i + j + k,
                            data = train_data,
                            family = binomial(link = "logit"),
                            prior = normal(0, 2.5),
                            prior_intercept = normal(0, 2.5),
                            chains = 5, iter = 4000, warmup = 2000, seed = 12345)
    
    # Posterior Predictions
    posterior_epred_samples <- posterior_epred(bayes_model, newdata = test_data, draws = 2000)
    predicted_probs_mean <- colMeans(posterior_epred_samples)
    
    all_true_labels <- c(all_true_labels, as.numeric(test_data$Label) - 1)
    all_predicted_probs <- c(all_predicted_probs, predicted_probs_mean)
    
    predicted_class <- factor(ifelse(predicted_probs_mean > 0.5, "Discharge_Good", "Discharge_Poor"), 
                              levels = c("Discharge_Poor", "Discharge_Good"))
    
    test_labels <- factor(test_data$Label, levels = c("Discharge_Poor", "Discharge_Good"))
    
    # Compute Confusion Matrix
    confusion <- confusionMatrix(predicted_class, test_labels)
    confusion_matrices <- append(confusion_matrices, list(confusion$table))
    
    # Compute AUC
    roc_curve <- roc(as.numeric(test_data$Label) - 1, predicted_probs_mean, plot = FALSE)
    auc_values <- c(auc_values, roc_curve$auc)
    
    # Compute MCC
    TP <- confusion$table[2, 2]
    TN <- confusion$table[1, 1]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    mcc_numerator <- (TP * TN) - (FP * FN)
    mcc_denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc_values <- c(mcc_values, ifelse(mcc_denominator == 0, 0, mcc_numerator / mcc_denominator))
  }
}

# Compute final aggregated AUC
roc_curve_final <- roc(all_true_labels, all_predicted_probs, plot = FALSE)
auc_value_final <- auc(roc_curve_final)

# … after you've computed roc_curve_final and auc_value_final …

library(pROC)
library(ggplot2)

# Build the ggroc plot with big fonts, red subtitle, no grid, black border
roc_gg <- ggroc(roc_curve_final,
                legacy.axes = TRUE,    # x from 1→0
                size        = 1,       # thick blue line
                colour      = "blue") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "solid", colour = "grey50", size = 1) +
  labs(
    title    = "ROC Curve (20×5-Fold CV)",
    subtitle = sprintf("AUC = %.2f", auc_value_final),
    x        = "Specificity",
    y        = "Sensitivity"
  ) +
  theme_bw(base_size = 24) +
  theme(
    panel.grid    = element_blank(),
    panel.border  = element_rect(colour = "black", fill = NA),
    plot.title    = element_text(face = "bold", size = 28, hjust = 0.5, vjust = -3),
    plot.subtitle = element_text(color = "red",
                                 size  = 24,
                                 hjust = 0.5, vjust = -25,
                                 margin = margin(t = 10)),  # give some space below title
    axis.title    = element_text(face = "bold", size = 24),
    axis.text     = element_text(size = 22)
  )

# Save it
ggsave(
  filename = file.path(output_dir, "roc_curve_20x5fold_improved.png"),
  plot     = roc_gg,
  width    = 8, height = 8, dpi = 600
)



# Compute and Print Performance Metrics
cat("Final Aggregated AUC:", round(auc_value_final, 2), "\n")
cat("Average MCC:", round(mean(mcc_values), 3), "\n")

# Compute Final Confusion Matrix
final_conf_matrix <- Reduce(`+`, confusion_matrices)  # Sum confusion matrices across all runs
conf_matrix_df <- as.data.frame.table(final_conf_matrix)
colnames(conf_matrix_df) <- c("Actual", "Predicted", "Count")

# Ensure factor levels match
conf_matrix_df$Actual <- factor(conf_matrix_df$Actual, levels = c("Discharge_Poor", "Discharge_Good"))
conf_matrix_df$Predicted <- factor(conf_matrix_df$Predicted, levels = c("Discharge_Poor", "Discharge_Good"))

# Define quadrant labels
conf_labels <- data.frame(
  Actual = factor(c("Discharge_Poor", "Discharge_Poor", "Discharge_Good", "Discharge_Good"), 
                  levels = levels(conf_matrix_df$Actual)),
  Predicted = factor(c("Discharge_Poor", "Discharge_Good", "Discharge_Poor", "Discharge_Good"), 
                     levels = levels(conf_matrix_df$Predicted)),
  Label = c("True Negative", "False Positive", "False Negative", "True Positive")
)

# Save the Confusion Matrix Plot
ggsave(file.path(output_dir, "confusion_matrix_20x5fold.png"),
       ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Count)) +
         geom_tile(color = "white") +  
         geom_text(aes(label = Count), 
                   color = "black", size = 12, fontface = "bold", vjust = 1.8) +  
         geom_text(
           data = conf_labels,
           aes(x = Predicted, y = Actual, label = Label),
           inherit.aes = FALSE,            # <<–– prevents looking for Count in conf_labels
           color = "black", size = 9, fontface = "bold", vjust = -0.5
         ) +  
         scale_fill_gradient(low = "deepskyblue", high = "darkorange") +  
         labs(title = "Confusion Matrix (20x5-Fold CV)", x = "Predicted", y = "Actual") +
         scale_x_discrete(labels = c("Poor", "Good")) +  
         scale_y_discrete(labels = c("Poor", "Good")) +  
         theme_classic() +  
         theme(
           text = element_text(size = 24),
           plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
           axis.text.y = element_text(size = 24, angle = 90),
           legend.position = "none"
         ),
       width = 8, height = 8, dpi = 600
)

