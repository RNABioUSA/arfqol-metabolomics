# Load required libraries
if (!require("rstanarm")) install.packages("rstanarm", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("bayesplot")) install.packages("bayesplot", dependencies = TRUE)
if (!require("viridis")) install.packages("viridis", dependencies = TRUE)

library(rstanarm)
library(ggplot2)
library(bayesplot)
library(viridis)

# Load your data
data <- read.csv('significant_narrowed.csv')

# Exclude the "Sample" column (first column) as it is not relevant for the analysis
data <- data[, -1]

# Convert the Label column to a binary factor with correct levels ("Discharge_Poor" = 0, "Discharge_Good" = 1)
data$Label <- factor(data$Label, levels = c("Discharge_Poor", "Discharge_Good"))
if (any(is.na(data$Label))) {
  stop("Error: The Label column contains values other than 'Discharge_Poor' or 'Discharge_Good'")
}

# Ensure all metabolite columns are numeric
data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], function(x) as.numeric(as.character(x)))

# Check for any NA values in the dataset
if (sum(is.na(data)) > 0) {
  cat("Number of NA values:", sum(is.na(data)), "\n")
  data <- na.omit(data)
}

# --- Log Transformation and Auto-scaling (Standardization) of Metabolite Data ---

# Apply log transformation to metabolite columns (starting from column 2, skipping the label column)
data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], function(x) log(x + 1))  # Log-transform to handle skewness, adding 1 to avoid log(0)

# Apply auto-scaling (standardization) after log transformation
data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], scale)  # Auto-scaling (mean = 0, sd = 1)

# --- Bayesian Logistic Regression Model ---
# Ensure column names in the formula match the actual column names in your dataset

bayes_model <- stan_glm(Label ~ a + b + d + g + i + j + k,
                        data = data,
                        family = binomial(link = "logit"),
                        prior = normal(0, 2.5),  # Weakly informative prior for coefficients
                        prior_intercept = normal(0, 2.5),  # Weakly informative prior for intercept
                        chains = 5, iter = 4000, warmup = 2000, seed = 12345)

# Create a folder to store the results
dir.create("Bayesian_result", showWarnings = FALSE)

# Define a base theme for all plots with axis text size set to 16
base_theme <- theme_minimal(base_size = 24) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(size = 24, color = "#000000"),  # Set x-axis text size to 16 and color to black
        axis.text.y = element_text(size = 24, color = "#000000"))  # Set y-axis text size to 16 and color to black

# Get the number of chains in the model
num_chains <- bayes_model$stanfit@sim$chains

# Use a color-blind friendly palette with enough colors for the number of chains
cb_palette <- viridis(num_chains, option = "D", alpha = 0.4)

# --- Trace Plot (Including the Intercept) ---
traceplot <- mcmc_trace(as.array(bayes_model), pars = c("(Intercept)", "a", "b", "d", "g", "i", "j", "k"), facet_args = list(ncol = 4)) +
  base_theme +
  labs(x = "Iteration", y = "Trace", title = "Trace Plot") +
  theme(axis.text.x = element_text(size = 18, color = "#000000"),  # Set x-axis text size to 12
        axis.text.y = element_text(size = 18, color = "#000000"))  # Set y-axis text size to 12

# Save the trace plot
ggsave("Bayesian_result/traceplot_with_intercept_3000_iterations.png", plot = traceplot, dpi = 600, width = 15, height = 8)

# --- Posterior Plot (Including the Intercept) with 95% Credible Interval ---
posterior_plot_with_intercept <- mcmc_areas(
  as.array(bayes_model), 
  pars = c("(Intercept)", "a", "b", "d", "g", "i", "j", "k"), 
  prob = 0.95  # Show the 95% credible interval instead of the default 50%
) +
  base_theme +
  labs(x = "Posterior Distribution", title = "Posterior Plot")

# Save the posterior plot with the intercept included
ggsave("Bayesian_result/posterior_plot_with_intercept_95CI.png", plot = posterior_plot_with_intercept, dpi = 600, width = 10, height = 8)

# --- Posterior Predictive Plot ---
posterior_predictive_plot <- pp_check(bayes_model, nreps = 64) +  # Use nreps instead of nsamples
  base_theme +
  labs(title = "Posterior Predictive Plot", x = "Observed Outcome", y = "Predicted Density") +
  scale_y_continuous(limits = c(0, 1.5), breaks = c(0.0, 0.5, 1.0, 1.5))  # Set y-axis limits and breaks

# Save the posterior predictive plot
ggsave("Bayesian_result/posterior_predictive_plot_64_samples.png", plot = posterior_predictive_plot, dpi = 600, width = 8, height = 6)

# --- Autocorrelation Plot (Including the Intercept) ---
plot_acf_with_adjustable_scaling <- function(lags = 300, breaks = c(0, 150, 300), axis_text_size = 28, line_color = "#1f77b4") {
  acf_plot <- mcmc_acf(as.array(bayes_model), pars = c("(Intercept)", "a", "b", "d", "g", "i", "j", "k"), lags = lags) +
    base_theme +
    scale_x_continuous(breaks = breaks) +
    labs(x = "Lag", y = "Autocorrelation", title = "Autocorrelation Plot") +
    theme(panel.spacing = unit(2, "lines"), axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
    geom_line(color = line_color, linewidth = 0.5)
  
  # Save the autocorrelation plot
  ggsave(paste0("Bayesian_result/acf_plot_with_intercept_", lags, ".png"), plot = acf_plot, dpi = 600, width = 15, height = 9)
}

# Example autocorrelation plot usage
plot_acf_with_adjustable_scaling(lags = 300, breaks = c(0, 150, 300), axis_text_size = 20, line_color = "#1f77b4")

