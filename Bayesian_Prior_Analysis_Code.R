# Full R Script: Data Prep + Prior Predictive Check for Bayesian Logistic Regression

# 1. Load required libraries
if (!require("rstanarm"))   install.packages("rstanarm",   dependencies = TRUE)
if (!require("bayesplot"))  install.packages("bayesplot",  dependencies = TRUE)
if (!require("ggplot2"))    install.packages("ggplot2",    dependencies = TRUE)

library(rstanarm)
library(bayesplot)
library(ggplot2)

# 2. Read raw data (no prior transformations)
raw <- read.csv("significant_narrowed.csv", check.names = TRUE)

# 3. Extract labels and metabolite matrix
labels    <- raw$Label
metab_raw <- raw[, 3:ncol(raw)]   # assuming cols 3:9 are a, b, d, g, i, j, k

# 4. Check positivity of raw intensities
if (any(metab_raw < 0, na.rm = TRUE)) {
  stop("Negative values found in raw metabolites!")
}

# 5. Log-transform (x → log(x + 1))
metab_log <- log(metab_raw + 1)

# 6. Standardize to mean = 0, sd = 1
metab_scaled <- scale(metab_log)

# 7. Reassemble into a data frame for modeling
data <- data.frame(
  Label = factor(labels, levels = c("Discharge_Poor","Discharge_Good")),
  as.data.frame(metab_scaled)
)

# 8. Fit a prior-only Bayesian logistic model
prior_mod <- stan_glm(
  Label ~ a + b + d + g + i + j + k,
  data            = data,
  family          = binomial(link = "logit"),
  prior           = normal(0, 5),
  prior_intercept = normal(0, 5),
  prior_PD        = TRUE,      # draw from the prior only
  chains          = 4,
  iter            = 2000,
  seed            = 12345,
  refresh         = 0          # suppress Stan progress
)

# 9. Simulate from that prior
y_prior <- posterior_predict(prior_mod)   # simulated 0/1 outcomes
p_prior <- posterior_epred(prior_mod)     # simulated probabilities

# 10. Define a function to compute proportion of zeros
prop_zero <- function(x) mean(x == 0)

# 11. Plot proportion of zeros across 200 prior draws
ppc_stat(
  y    = rep(0, nrow(data)),        # dummy numeric vector, no NAs
  yrep = y_prior[1:200, ],          # subset of prior draws
  stat = prop_zero
) +
  ggtitle("Prior Predictive: Proportion of Zeros per Dataset") +
  ylab("Proportion of Zeros")

# 12. Plot the spread of simulated probabilities (80% intervals)
mcmc_areas(
  p_prior[1:200, ],
  prob = 0.8
) +
  ggtitle("Prior Predictive Probabilities (80% intervals)") +
  xlab("Predicted Probability of 'Good'")

# 13. (Optional) Overlay density of a few simulated label sets
ppc_dens_overlay(
  y    = rep(0, nrow(data)),
  yrep = y_prior[1:20, ]
) +
  ggtitle("Density Overlay: Simulated Labels (first 20 draws)") +
  xlab("Outcome (0 = Poor, 1 = Good)")

