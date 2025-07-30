# Load required libraries
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(rstanarm, ggplot2, bayesplot, viridis, pROC, caret, reshape2)

# ✅ Create output directory "Bayesian_Model" if it doesn't exist
bayesian_output_dir <- "Bayesian_Model"
if (!dir.exists(bayesian_output_dir)) dir.create(bayesian_output_dir)

# ✅ Load Data
data <- read.csv('significant_narrowed.csv')

# ✅ Data Preprocessing
data <- data[, -1]  # Remove first column (Sample ID)
data$Label <- factor(data$Label, levels = c("Discharge_Poor", "Discharge_Good"))
if (any(is.na(data$Label))) stop("Error: Label column contains invalid values")

# Convert metabolite columns to numeric
data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], function(x) as.numeric(as.character(x)))
data <- na.omit(data)  # Remove NA values

# ✅ Log-transform & Auto-scale metabolite data
numeric_cols <- sapply(data[, 2:ncol(data)], is.numeric)
data[, 2:ncol(data)][, numeric_cols] <- lapply(data[, 2:ncol(data)][, numeric_cols], function(x) scale(log(x + 1)))

# ✅ Bayesian Logistic Regression Model
bayes_model <- stan_glm(Label ~ a + b + d + g + i + j + k,
                        data = data,
                        family = binomial(link = "logit"),
                        prior = normal(0, 2.5),  
                        prior_intercept = normal(0, 2.5),  
                        chains = 5, iter = 4000, warmup = 2000, seed = 12345)

# ✅ Extract Posterior Samples (Last 2000 Iterations)
posterior_epred_samples <- posterior_epred(bayes_model, draws = 2000)
predicted_probs_mean <- colMeans(posterior_epred_samples)

# Convert Labels to binary (0 = Discharge_Poor, 1 = Discharge_Good)
true_labels <- as.numeric(data$Label) - 1

# ✅ Compute Confusion Matrix
predicted_class <- ifelse(predicted_probs_mean > 0.5, "Discharge_Good", "Discharge_Poor")
confusion <- confusionMatrix(factor(predicted_class, levels = c("Discharge_Poor", "Discharge_Good")), data$Label)
print(confusion)

# ✅ Compute Performance Metrics
TP <- confusion$table[2, 2]
TN <- confusion$table[1, 1]
FP <- confusion$table[1, 2]
FN <- confusion$table[2, 1]

mcc_numerator <- (TP * TN) - (FP * FN)
mcc_denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
mcc_value <- ifelse(mcc_denominator == 0, 0, mcc_numerator / mcc_denominator)

precision <- TP / (TP + FP)
recall <- TP / (TP + FN)
f1_score <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))

cat("MCC:", round(mcc_value, 3), "\n")
cat("Precision:", round(precision, 3), "\n")
cat("Recall:", round(recall, 3), "\n")
cat("F1 Score:", round(f1_score, 3), "\n")

# ✅ Compute Final Aggregated AUC
roc_curve_final <- roc(true_labels, predicted_probs_mean, plot = FALSE)
auc_value_final <- auc(roc_curve_final)
cat("Final AUC:", round(auc_value_final, 2), "\n")

# ✅ Save **ROC Curve Plot**
png(file.path(bayesian_output_dir, "roc_curve_bayesian.png"), width = 4000, height = 4000, res = 600)
par(mar = c(5, 5, 4, 2), mgp = c(3, 1.5, 0))
plot(roc_curve_final, 
     main = "ROC Curve (Bayesian Model)", 
     col = "blue", 
     lwd = 2, 
     xlab = "Specificity", 
     ylab = "Sensitivity", 
     cex.lab = 2,  
     cex.axis = 2,  
     cex.main = 2,  
     las = 1)  
text(0.6, 0.7, paste("AUC =", round(auc_value_final, 2)), col = "red", cex = 2, font = 2)
dev.off()

# ✅ Extract Confusion Matrix as Data Frame
final_conf_matrix <- as.table(confusion$table)
conf_matrix_df <- as.data.frame(final_conf_matrix)
colnames(conf_matrix_df) <- c("Actual", "Predicted", "Count")
conf_matrix_df$Count <- as.numeric(conf_matrix_df$Count)

# ✅ Ensure Factor Levels Match
conf_matrix_df$Actual <- factor(conf_matrix_df$Actual, levels = c("Discharge_Poor", "Discharge_Good"))
conf_matrix_df$Predicted <- factor(conf_matrix_df$Predicted, levels = c("Discharge_Poor", "Discharge_Good"))

# ✅ Define Quadrant Labels
conf_labels <- data.frame(
  Actual = factor(c("Discharge_Poor", "Discharge_Poor", "Discharge_Good", "Discharge_Good"), 
                  levels = levels(conf_matrix_df$Actual)),
  Predicted = factor(c("Discharge_Poor", "Discharge_Good", "Discharge_Poor", "Discharge_Good"), 
                     levels = levels(conf_matrix_df$Predicted)),
  Label = c("True Negative", "False Positive", "False Negative", "True Positive")
)

# ✅ Merge Labels into DataFrame
conf_labels <- merge(conf_matrix_df, conf_labels, by = c("Actual", "Predicted"), all.x = TRUE)

# ✅ Save **Confusion Matrix Plot**
ggsave(file.path(bayesian_output_dir, "confusion_matrix_bayesian.png"),
       ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Count)) +
         geom_tile(color = "white") +  
         geom_text(aes(label = Count), color = "black", size = 8, fontface = "bold", vjust = 0.5) +  
         geom_text(data = conf_labels, aes(x = Predicted, y = Actual, label = Label), 
                   color = "black", size = 6, fontface = "bold", vjust = -1.2) +  
         scale_fill_gradient(low = "deepskyblue", high = "darkorange") +  
         labs(title = "Confusion Matrix (Bayesian Model)", x = "Predicted", y = "Actual") +
         scale_x_discrete(labels = c("Poor", "Good")) +  
         scale_y_discrete(labels = c("Poor", "Good")) +  
         theme_classic() +  
         theme(
           text = element_text(size = 16),  
           plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
           axis.title.x = element_text(size = 24, face = "bold"),
           axis.title.y = element_text(size = 24, face = "bold"),
           axis.text.x = element_text(size = 16),
           axis.text.y = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 0.5),  
           legend.position = "none"
         ),
       width = 8, height = 8, dpi = 600)

cat("Plots successfully saved in 'Bayesian_Model' folder! 🚀")
