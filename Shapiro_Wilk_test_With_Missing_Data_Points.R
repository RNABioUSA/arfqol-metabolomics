# <<— insert at the very top of run_shapiro.R, before the "if (!exists("df"))" check:
df <- read.csv("Discharge_Compounds_detected.csv",
               check.names      = FALSE,
               stringsAsFactors = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) load your data: here we assume you've already done something like
#    df <- read.csv("your_data.csv", check.names=FALSE, stringsAsFactors=FALSE)
#
#    and that df has:
#      • column 1 = sample IDs
#      • column 2 = labels
#      • columns 3..n = numeric metabolite values (with NAs)
# ─────────────────────────────────────────────────────────────────────────────

# just in case:
if (!exists("df")) stop("Please load your data into a data.frame called `df`")

# ─────────────────────────────────────────────────────────────────────────────
# 2) prepare output directory & results container
# ─────────────────────────────────────────────────────────────────────────────
out_dir <- "shapiro_plots"
if (!dir.exists(out_dir)) dir.create(out_dir)

results <- data.frame(
  metabolite = character(),
  n_obs      = integer(),
  W          = numeric(),
  p.value    = numeric(),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) loop over each metabolite, test & plot
# ─────────────────────────────────────────────────────────────────────────────
met_names <- names(df)[3:ncol(df)]

for (m in met_names) {
  # drop NAs
  vals <- df[[m]]
  vals <- vals[!is.na(vals)]
  n   <- length(vals)
  
  # run Shapiro–Wilk only if 3 ≤ n ≤ 5000
  if (n >= 3 && n <= 5000) {
    t   <- shapiro.test(vals)
    W   <- unname(t$statistic)
    p   <- t$p.value
  } else {
    W <- NA; p <- NA
  }
  
  # record
  results <- rbind(results, data.frame(
    metabolite = m,
    n_obs      = n,
    W          = W,
    p.value    = p,
    stringsAsFactors = FALSE
  ))
  
  # make a filesystem‐safe filename
  fn <- gsub("[^A-Za-z0-9_]", "_", m)
  
  # open PNG device
  png(
    filename = file.path(out_dir, paste0(fn, ".png")),
    width    = 1200, height = 600, res = 150
  )
  par(mfrow = c(1, 2), mar = c(4,4,2,1))
  
  # (a) histogram
  hist(
    vals,
    main = paste("Histogram of", m),
    xlab = m,
    col  = "lightgray",
    border = "white"
  )
  
  # (b) Q‐Q plot
  qqnorm(vals, main = paste("Q-Q Plot of", m))
  qqline(vals)
  
  dev.off()
}

# ─────────────────────────────────────────────────────────────────────────────
# 4) write out your table of Shapiro‐Wilk results
# ─────────────────────────────────────────────────────────────────────────────
write.csv(results,
          file = "shapiro_results.csv",
          row.names = FALSE)

message("Done!  
 • shapiro_results.csv written to your working directory  
 • PNGs in folder: ", normalizePath(out_dir))

