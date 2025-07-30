# Load the data with metabolite names preserved
data <- read.csv("significant_clean.csv", check.names = FALSE)

# Remove the 'Sample' and 'Label' columns to retain only metabolite data
metabolite_data <- data[, -c(1, 2)]  # Assuming 'Sample' is column 1 and 'Label' is column 2

# Calculate the correlation matrix using Spearman's method
cor_matrix <- cor(metabolite_data, method = "spearman")

# Set the correlation threshold
threshold <- 0.8

# Get metabolite names
metabolite_names <- colnames(metabolite_data)

# Find pairs with high correlations
high_corr_pairs <- which(abs(cor_matrix) > threshold & cor_matrix != 1, arr.ind = TRUE)
high_corr_list <- data.frame(
  Metabolite1 = metabolite_names[high_corr_pairs[, 1]],
  Metabolite2 = metabolite_names[high_corr_pairs[, 2]],
  Correlation = cor_matrix[high_corr_pairs]
)

# Print the high correlation pairs
print(high_corr_list)

# Save the high correlation pairs to a CSV file (optional)
write.csv(high_corr_list, "high_correlation_metabolites.csv", row.names = FALSE)

# Load the corrplot package
library(corrplot)

# Save the plot as a PNG with smaller size and reduced resolution
png("correlation_matrix_plot_clean_smaller.png", width = 3000, height = 3000, res = 300)  # Smaller plot

# Plot the correlation matrix with larger text labels for better readability
corrplot(cor_matrix, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200), 
         tl.col = "black", tl.cex = 1.0, tl.srt = 90)

# Turn off the PNG device to save the file
dev.off()


