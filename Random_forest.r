# Load required libraries
library(randomForest)
library(caret)
library(ggplot2)
library(corrplot)

# 1. Load and Explore the Data
# Set check.names = FALSE to prevent R from changing column names
data <- read.csv("Discharge_Significant_met.csv", check.names = FALSE)

# Display basic information about the data
str(data)
head(data)

# 2. Encode the Labels
# Convert 'Label' to a factor (classification labels)
data$Label_encoded <- as.numeric(factor(data$Label, levels = c("Discharge_Good", "Discharge_Poor")))

# 3. Correlation Analysis
# Compute the correlation matrix for all numeric features (without filtering)
cor_matrix <- cor(data[, -c(1, 2)])  # Exclude 'Sample' and 'Label' columns

# Display correlation of features with 'Label_encoded'
cor_with_target <- cor_matrix[, "Label_encoded"]
print(sort(cor_with_target, decreasing = TRUE))

# 4. Prepare Data for Modeling
# Prepare feature matrix X and target vector y
X <- data[, -c(1, 2, ncol(data))]  # Exclude 'Sample', 'Label', and 'Label_encoded'
y <- as.factor(data$Label_encoded)  # Convert 'Label_encoded' to factor for classification

# Split data into training and testing sets (60% train, 40% test)
set.seed(42)
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# 5. Train Random Forest Classifier
rf_model <- randomForest(X_train, y_train, importance = TRUE)

# 6. Evaluate the Model
# Predict on the test set
y_pred <- predict(rf_model, X_test)

# Confusion Matrix and other evaluation metrics
conf_matrix <- confusionMatrix(y_pred, y_test)
print(conf_matrix)

# 7. Feature Importance Visualization
# Extract feature importance
importance_df <- as.data.frame(importance(rf_model))
importance_df$Feature <- rownames(importance_df)

# Sort by importance
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]

# Plot Random Forest feature importance with white background and aligned gridlines
feature_plot <- ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Feature Importance from Random Forest") +
  xlab("Features") +
  ylab("Importance") +
  theme_bw(base_size = 12) +  # Switch to theme_bw() for white background with gridlines
  theme(
    plot.title = element_text(size = 12, face = "bold"),  # Title text size
    axis.title = element_text(size = 9),  # Axis title text size
    axis.text = element_text(size = 6),  # Axis label text size
    axis.text.y = element_text(margin = margin(r = 2)),  # Reduce the margin between bars and labels
    panel.background = element_rect(fill = "white"),  # Explicitly set panel background to white
    plot.background = element_rect(fill = "white")  # Explicitly set plot background to white
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis (Importance) starts at 0
  theme(
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.grid.major.x = element_line(color = "grey80")  # Keep major x-axis gridlines
  )

# Create a folder "Variable_selection" if it doesn't already exist
if (!dir.exists("Variable_selection")) {
  dir.create("Variable_selection")
}

# Save the plot as PNG with high resolution and aligned gridlines
png_file <- "Variable_selection/feature_importance_white_background.png"
ggsave(png_file, plot = feature_plot, width = 10, height = 8, dpi = 600)  # Adjust width, height, and resolution

# Save the plot as PDF with gridlines starting from 0 and a white background
pdf_file <- "Variable_selection/feature_importance_white_background.pdf"
ggsave(pdf_file, plot = feature_plot, width = 10, height = 8)  # Adjust width and height

