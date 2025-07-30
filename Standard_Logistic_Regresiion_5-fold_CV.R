# 5-Fold Cross-Validation for Standard Logistic Regression
# Outputs saved to “Standard_Log_Reg_CV” folder

# 0. Define & create output folder
cv_output_dir <- "Standard_Log_Reg_CV"
if (!dir.exists(cv_output_dir)) dir.create(cv_output_dir)

# 1. Load necessary libraries
library(caret)
library(pROC)
library(ggplot2)

# 2. Load and preprocess data
data <- read.csv("significant_narrowed.csv", check.names = TRUE)
# Convert Label to binary factor
if ("Label" %in% colnames(data)) {
  data$Label <- ifelse(data$Label == 'Discharge_Good', 1, 0)
}
data$Label <- factor(data$Label, levels = c(0,1), labels = c("Class0","Class1"))
# Drop sample ID column
data <- data[ , -1]
# Log‐transform and autoscale
num_cols <- setdiff(names(data), "Label")
data[ num_cols ] <- log(data[ num_cols ] + 1)
data[ num_cols ] <- scale(data[ num_cols ])

# 3. Set up 5-fold CV
set.seed(123)
k_folds    <- 5
folds_list <- createFolds(data$Label, k = k_folds, list = TRUE)

# 4. Initialize storage
all_true_labels   <- c()
all_pred_probs    <- c()
confusion_tables  <- list()
kappa_vals        <- c()
precision_vals    <- c()
recall_vals       <- c()
f1_vals           <- c()
mcc_vals          <- c()
auc_vals          <- c()

# 5. CV loop
for (i in seq_along(folds_list)) {
  cat("Fold", i, "of", k_folds, "...\n")
  
  test_idx  <- folds_list[[i]]
  train_df  <- data[-test_idx, ]
  test_df   <- data[ test_idx, ]
  
  # Fit standard logistic regression
  fit <- glm(Label ~ a + b + d + g + i + j + k,
             data = train_df, family = binomial)
  
  # Predict on test fold
  probs <- predict(fit, newdata = test_df, type = "response")
  preds <- factor(ifelse(probs > 0.5, "Class1", "Class0"),
                  levels = c("Class0","Class1"))
  
  # Store labels & probs for aggregated ROC
  all_true_labels <- c(all_true_labels, as.numeric(test_df$Label)-1)
  all_pred_probs  <- c(all_pred_probs,  probs)
  
  # Confusion matrix & metrics
  cm     <- confusionMatrix(preds, test_df$Label)
  tbl    <- cm$table
  confusion_tables[[i]] <- tbl
  
  # AUC for this fold
  roc_i  <- roc(as.numeric(test_df$Label)-1, probs, plot = FALSE)
  auc_vals   <- c(auc_vals, as.numeric(auc(roc_i)))
  
  # Kappa, Precision, Recall, F1
  kappa_vals     <- c(kappa_vals,   cm$overall["Kappa"])
  precision_vals <- c(precision_vals, cm$byClass["Precision"])
  recall_vals    <- c(recall_vals,    cm$byClass["Recall"])
  f1_vals        <- c(f1_vals,        cm$byClass["F1"])
  
  # MCC
  TP <- tbl[2,2]; TN <- tbl[1,1]; FP <- tbl[1,2]; FN <- tbl[2,1]
  mcc_num <- (TP*TN) - (FP*FN)
  mcc_den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc_vals <- c(mcc_vals, ifelse(mcc_den==0, 0, mcc_num/mcc_den))
}

# 6. Print average metrics across folds
cat("Average AUC:       ", round(mean(auc_vals), 3), "\n")
cat("Average Kappa:     ", round(mean(kappa_vals),3), "\n")
cat("Average Precision: ", round(mean(precision_vals),3), "\n")
cat("Average Recall:    ", round(mean(recall_vals),3), "\n")
cat("Average F1 Score:  ", round(mean(f1_vals),3), "\n")
cat("Average MCC:       ", round(mean(mcc_vals),3), "\n")

# 7. Aggregated ROC curve & AUC
roc_agg <- roc(all_true_labels, all_pred_probs, plot = FALSE)
auc_agg <- round(auc(roc_agg), 3)

# 8. Aggregated confusion matrix
conf_mat_sum <- Reduce(`+`, confusion_tables)
conf_df <- as.data.frame.table(conf_mat_sum)
colnames(conf_df) <- c("Actual","Predicted","Count")
conf_df$Actual    <- factor(conf_df$Actual,    levels=c("Class0","Class1"), labels=c("Poor","Good"))
conf_df$Predicted <- factor(conf_df$Predicted, levels=c("Class0","Class1"), labels=c("Poor","Good"))
conf_labels <- data.frame(
  Actual    = factor(c("Poor","Poor","Good","Good"), levels=c("Poor","Good")),
  Predicted = factor(c("Poor","Good","Poor","Good"), levels=c("Poor","Good")),
  Label     = c("True Negative","False Positive","False Negative","True Positive")
)

# 9. Save aggregated Confusion Matrix plot
ggsave(
  filename = file.path(cv_output_dir, "confusion_matrix_cv.png"),
  plot = ggplot(conf_df, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white") +
    # BIG counts
    geom_text(aes(label = Count),
              color     = "black",
              size      = 12,
              fontface  = "bold",
              vjust     = 1.8) +
    # BIG quadrant labels
    geom_text(data = conf_labels,
              aes(x = Predicted, y = Actual, label = Label),
              inherit.aes = FALSE,
              color       = "black",
              size        = 9,
              fontface    = "bold",
              vjust       = -0.5) +
    scale_fill_gradient(low  = "deepskyblue",
                        high = "darkorange") +
    labs(title = "5-Fold CV Confusion Matrix",
         x     = "Predicted",
         y     = "Actual") +
    scale_x_discrete(labels = c("Poor","Good")) +
    scale_y_discrete(labels = c("Poor","Good")) +
    theme_classic() +
    theme(
      # BIG everything
      text           = element_text(size = 24),
      plot.title     = element_text(hjust = 0.5,
                                    size  = 24,
                                    face  = "bold"),
      axis.title.x   = element_text(size = 24, face = "bold"),
      axis.title.y   = element_text(size = 24, face = "bold"),
      axis.text.x    = element_text(size = 24, face = "bold"),
      axis.text.y    = element_text(size = 24, face = "bold",
                                    angle = 90,
                                    vjust  = 0.5,
                                    hjust  = 0.5),
      legend.position = "none"
    ),
  width = 8,
  height = 8,
  dpi    = 600
)


library(pROC)
library(ggplot2)

roc_gg <- ggroc(roc_agg,
                legacy.axes = TRUE,
                size        = 2,
                colour      = "blue") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "solid", colour = "grey50") +
  labs(
    title    = "ROC Curve (5-Fold CV)",
    subtitle = sprintf("AUC = %.2f", auc_agg),
    x        = "Specificity",
    y        = "Sensitivity"
  ) +
  theme_bw(base_size = 24) +
  theme(
    panel.grid     = element_blank(),
    panel.border   = element_rect(colour = "black", fill = NA),
    plot.title     = element_text(face = "bold", size = 28, hjust = 0.5, vjust = -5),
    plot.subtitle  = element_text(color = "red", size = 24, hjust = 0.5, vjust = -30),  # ← red subtitle
    axis.title     = element_text(face = "bold", size = 24),
    axis.text      = element_text(size = 22)
  )

ggsave(
  file.path(cv_output_dir, "roc_curve_cv_red_subtitle.png"),
  plot   = roc_gg,
  width  = 8, height = 8, dpi = 600
)

