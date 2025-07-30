# Standard logistic regression with outputs saved to “Standard_Log_Reg” folder

# 0. Define & create output folder
output_dir <- "Standard_Log_Reg"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 1. Load necessary libraries
library(pROC)
library(caret)
library(ggplot2)

# 2. Load the data (make sure your file path is correct)
data <- read.csv("significant_narrowed.csv", check.names = TRUE)

# 3. Convert the 'Label' column to a binary factor with levels "Class0" and "Class1"
if ("Label" %in% colnames(data)) {
  data$Label <- ifelse(data$Label == 'Discharge_Good', 1, 0)
}
data$Label <- factor(data$Label, levels = c(0, 1), labels = c("Class0", "Class1"))

# 4. Remove the 'Sample' column (assuming it's the first column)
data <- data[, -1]

# 5. Apply log transformation to the data (excluding 'Label')
numeric_columns <- setdiff(names(data), "Label")
data[numeric_columns] <- log(data[numeric_columns] + 1)

# 6. Apply autoscaling (standardization)
data[numeric_columns] <- scale(data[numeric_columns])

# 7. Define the formula for logistic regression
formula <- as.formula("Label ~ a + b + d + g + i + j + k")

# 8. Fit the logistic regression model
model <- glm(formula, data = data, family = binomial)

# 9. Predict probabilities (for ROC curve and AUC)
probabilities <- predict(model, type = "response")

# 10. Get predicted class labels (threshold 0.5) and convert to factors
predicted_class <- ifelse(probabilities > 0.5, "Class1", "Class0")
predicted_class <- factor(predicted_class, levels = c("Class0", "Class1"))

# 11. Create a confusion matrix
confusion <- confusionMatrix(predicted_class, data$Label)
print(confusion)

# 12. Cohen's Kappa
kappa_value <- confusion$overall['Kappa']
cat("Cohen's Kappa:", round(kappa_value, 3), "\n")

# 13. Precision & Recall
precision <- confusion$byClass["Precision"]
recall    <- confusion$byClass["Recall"]
cat("Precision: ", round(precision, 3), "\n")
cat("Recall (Sensitivity): ", round(recall, 3), "\n")

# 14. F1 Score
f1_score <- round(confusion$byClass["F1"], 3)
cat("F1 Score: ", f1_score, "\n")

# 15. Matthews Correlation Coefficient (MCC)
TP <- confusion$table[2, 2]
TN <- confusion$table[1, 1]
FP <- confusion$table[1, 2]
FN <- confusion$table[2, 1]
mcc_num <- (TP * TN) - (FP * FN)
mcc_den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
mcc_value <- ifelse(mcc_den == 0, 0, mcc_num / mcc_den)
cat("Matthews Correlation Coefficient (MCC):", round(mcc_value, 3), "\n")

# 16. ROC curve & AUC
roc_curve <- roc(as.numeric(data$Label) - 1, probabilities)
auc_value <- round(auc(roc_curve), 2)
cat("AUC:", auc_value, "\n")

# 17. Prepare data frame for ggplot confusion matrix
cm_table <- confusion$table
conf_matrix_df <- as.data.frame.table(cm_table)
colnames(conf_matrix_df) <- c("Actual", "Predicted", "Count")
conf_matrix_df$Actual    <- factor(conf_matrix_df$Actual,    levels = c("Class0","Class1"), labels = c("Poor","Good"))
conf_matrix_df$Predicted <- factor(conf_matrix_df$Predicted, levels = c("Class0","Class1"), labels = c("Poor","Good"))

conf_labels <- data.frame(
  Actual    = factor(c("Poor","Poor","Good","Good"), levels = c("Poor","Good")),
  Predicted = factor(c("Poor","Good","Poor","Good"), levels = c("Poor","Good")),
  Label     = c("True Negative","False Positive","False Negative","True Positive")
)

# 18. Save the Confusion Matrix plot (PNG)
ggsave(
  filename = file.path(output_dir, "confusion_matrix_standard.png"),
  plot = ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count),
              color = "black", size = 12, fontface = "bold", vjust = 0.6) +
    geom_text(
      data = conf_labels,
      aes(x = Predicted, y = Actual, label = Label),
      inherit.aes = FALSE,
      color = "black", size = 9, fontface = "bold", vjust = -1.2
    ) +
    scale_fill_gradient(low = "deepskyblue", high = "darkorange") +
    labs(
      title = "Confusion Matrix (Standard Logistic)",
      x     = "Predicted",
      y     = "Actual"
    ) +
    scale_x_discrete(labels = c("Poor", "Good")) +
    scale_y_discrete(labels = c("Poor", "Good")) +
    theme_classic() +
    theme(
      text            = element_text(size = 24),
      plot.title      = element_text(hjust = 0.5, size = 24, face = "bold"),
      axis.title.x    = element_text(size = 24, face = "bold"),
      axis.title.y    = element_text(size = 24, face = "bold"),
      axis.text.x     = element_text(size = 24, face = "bold"),
      axis.text.y     = element_text(size = 24, face = "bold",
                                     angle = 90, vjust = 0.5, hjust = 0.5),
      legend.position = "none"
    ),
  width = 8, height = 8, dpi = 600
)

# 19. Save the ROC curve plot (PNG)
# 19. Save an improved, publication-ready ROC curve
library(pROC)
library(ggplot2)

# your roc_curve and auc_value have already been computed above
# now build the ggroc plot exactly as before, but point it at roc_curve
roc_gg <- ggroc(roc_curve,
                legacy.axes = TRUE,  # 1→0 on x-axis
                size        = 2,     # thick blue line
                colour      = "blue") +
  # dashed diagonal
  geom_abline(intercept = 0, slope = 1,
              linetype = "solid",
              colour   = "grey50",
              size     = 1) +
  # title + red subtitle
  labs(
    title    = "ROC Curve (Standard Regression)",
    subtitle = sprintf("AUC = %.2f", auc_value),
    x        = "Specificity",
    y        = "Sensitivity"
  ) +
  # white background + black border
  theme_bw(base_size = 24) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "black", fill = NA),
    plot.title      = element_text(face = "bold", size = 28, hjust = 0.5, vjust = -5),
    plot.subtitle   = element_text(color = "red", size = 24, hjust = 0.5, vjust = -30),
    axis.title      = element_text(face = "bold", size = 24),
    axis.text       = element_text(size = 22)
  )

# save it at high res
ggsave(
  filename = file.path(output_dir, "roc_curve_standard_improved.png"),
  plot     = roc_gg,
  width    = 8,
  height   = 8,
  dpi      = 600
)

