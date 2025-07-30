# Full 5-Fold Bayesian Logistic Regression CV Script
# Outputs saved to “Model_parameters4” folder

# 0. Load required packages (install if missing)
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(rstanarm, ggplot2, pROC, caret, reshape2)

# 1. Create output directory
output_dir <- "Model_parameters4"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 2. Load and preprocess data
data <- read.csv("significant_narrowed.csv", stringsAsFactors = FALSE)
data <- data[, -1]  # drop Sample ID
data$Label <- factor(data$Label,
                     levels = c("Discharge_Poor", "Discharge_Good"))
if (any(is.na(data$Label))) stop("Error: Label column contains invalid values")

# ensure metabolite columns are numeric
for (j in 2:ncol(data)) {
  data[[j]] <- as.numeric(as.character(data[[j]]))
}
data <- na.omit(data)

# log-transform + autoscale on columns 2:end
met_cols <- 2:ncol(data)
data[met_cols] <- lapply(data[met_cols], function(x) scale(log(x + 1)))

# 3. Prepare 5‐folds
set.seed(123)
folds <- createFolds(data$Label, k = 5, list = TRUE)

# 4. Initialize storage
auc_values          <- numeric(length(folds))
mcc_values          <- numeric(length(folds))
confusion_matrices  <- vector("list", length(folds))
all_true_labels     <- numeric(0)
all_predicted_probs <- numeric(0)

# 5. Cross‐validation loop
for (i in seq_along(folds)) {
  cat("Processing Fold", i, "of", length(folds), "...\n")
  test_idx  <- folds[[i]]
  train_df  <- data[-test_idx, ]
  test_df   <- data[ test_idx, ]
  
  # fit Bayesian logistic regression
  fit <- stan_glm(
    Label ~ a + b + d + g + i + j + k,
    data            = train_df,
    family          = binomial(link = "logit"),
    prior           = normal(0, 2.5),
    prior_intercept = normal(0, 2.5),
    chains          = 5,
    iter            = 4000,
    warmup          = 2000,
    seed            = 12345
  )
  
  # posterior predictive probabilities
  post_samples      <- posterior_epred(fit, newdata = test_df, draws = 2000)
  probs             <- colMeans(post_samples)
  preds             <- factor(ifelse(probs > 0.5,
                                     "Discharge_Good",
                                     "Discharge_Poor"),
                              levels = c("Discharge_Poor","Discharge_Good"))
  
  # accumulate for aggregated ROC
  all_true_labels     <- c(all_true_labels,
                           as.numeric(test_df$Label) - 1)
  all_predicted_probs <- c(all_predicted_probs, probs)
  
  # confusion matrix & store
  cm <- confusionMatrix(preds, test_df$Label)
  confusion_matrices[[i]] <- cm
  
  # per‐fold AUC
  roc_i         <- roc(as.numeric(test_df$Label) - 1, probs, plot = FALSE)
  auc_values[i] <- as.numeric(auc(roc_i))
  
  # per‐fold MCC
  tbl     <- cm$table
  TP      <- tbl[2,2]; TN <- tbl[1,1]; FP <- tbl[1,2]; FN <- tbl[2,1]
  num     <- (TP * TN) - (FP * FN)
  den     <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc_values[i] <- ifelse(den == 0, 0, num/den)
}

# 6. Aggregate metrics
roc_agg       <- roc(all_true_labels, all_predicted_probs, plot = FALSE)
auc_agg       <- as.numeric(auc(roc_agg))
mean_mcc      <- mean(mcc_values)

cat("Aggregated AUC:", round(auc_agg, 2), "\n")
cat("Average MCC:",   round(mean_mcc, 2), "\n")

# 7. Save aggregated ROC curve (base R)
png(filename = file.path(output_dir, "roc_curve_aggregated.png"),
    width = 4000, height = 4000, res = 600)
plot(roc_agg,
     main  = "ROC Curve (5-Fold CV)",
     col   = "blue",
     lwd   = 2,
     xlab  = "1 - Specificity",
     ylab  = "Sensitivity",
     cex.lab  = 2,
     cex.axis = 2,
     cex.main = 2)
text(0.6, 0.7, paste("AUC =", round(auc_agg, 2)),
     col = "red", cex = 2)
dev.off()

# 8. Build final_conf_matrix by summing the 5 folds
list_tables      <- lapply(confusion_matrices, function(x) x$table)
final_conf_matrix <- Reduce(`+`, list_tables)

# 9. Convert to data.frame for ggplot
conf_matrix_df <- as.data.frame.table(final_conf_matrix)
colnames(conf_matrix_df) <- c("Actual","Predicted","Count")
conf_matrix_df$Actual    <- factor(conf_matrix_df$Actual,
                                   levels = c("Discharge_Poor","Discharge_Good"),
                                   labels = c("Poor","Good"))
conf_matrix_df$Predicted <- factor(conf_matrix_df$Predicted,
                                   levels = c("Discharge_Poor","Discharge_Good"),
                                   labels = c("Poor","Good"))

# text labels
conf_labels <- data.frame(
  Actual    = factor(c("Poor","Poor","Good","Good"), levels = c("Poor","Good")),
  Predicted = factor(c("Poor","Good","Poor","Good"), levels = c("Poor","Good")),
  Label     = c("True Negative","False Positive","False Negative","True Positive")
)

# 10. Save confusion matrix heatmap (ggplot2)
ggsave(
  filename = file.path(output_dir, "confusion_matrix_5fold.png"),
  plot = ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count),
              color = "black", size = 12, fontface = "bold", vjust = 1.8) +
    geom_text(data = conf_labels,
              aes(x = Predicted, y = Actual, label = Label),
              inherit.aes = FALSE,
              color = "black", size = 9, fontface = "bold", vjust = -0.5) +
    scale_fill_gradient(low = "deepskyblue", high = "darkorange") +
    labs(title = sprintf("5-Fold CV Confusion Matrix"),
         x = "Predicted", y = "Actual") +
    scale_x_discrete(labels = c("Poor","Good")) +
    scale_y_discrete(labels = c("Poor","Good")) +
    theme_classic() +
    theme(
      text = element_text(size = 24),  
      plot.title = element_text(hjust = 0.2, size = 24, face = "bold"),
      axis.title.x = element_text(size = 24, face = "bold"),
      axis.title.y = element_text(size = 24, face = "bold"),
      axis.text.x = element_text(size = 24, face = "bold"),
      axis.text.y = element_text(size = 24, face = "bold", angle = 90, vjust = 0.5, hjust = 0.5),  # Rotates y-axis labels
      legend.position = "none"
    ),
  width = 8, height = 8, dpi = 600)

