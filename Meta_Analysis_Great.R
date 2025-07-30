library(dplyr)
library(ggplot2)
library(stringr)

# 1. Read sanitized data
data <- read.csv("merged_patient_data_sanitized.csv", stringsAsFactors = FALSE) %>%
  mutate(Label = factor(Label, levels = c("Poor", "Good")))

# 2. Create categorical variables (shortened 'Normal weight' → 'Normal')
data <- data %>%
  mutate(
    BMI = case_when(
      bmi < 18.5 ~ "Underweight",
      bmi >= 18.5 & bmi <= 24.9 ~ "Normal",
      bmi >= 25 & bmi <= 29.9 ~ "Overweight",
      bmi >= 30 ~ "Obese"
    ),
    MAP = case_when(
      APACHE_MAP < 70 ~ "Low",
      APACHE_MAP >= 70 & APACHE_MAP <= 100 ~ "Normal",
      APACHE_MAP > 100 ~ "High"
    ),
    Pulse = case_when(
      APACHE_PULSE < 60 ~ "Low",
      APACHE_PULSE >= 60 & APACHE_PULSE <= 100 ~ "Normal",
      APACHE_PULSE > 100 ~ "High"
    ),
    Respiratory_Rate = case_when(
      APACHE_RESP < 12 ~ "Low",
      APACHE_RESP >= 12 & APACHE_RESP <= 20 ~ "Normal",
      APACHE_RESP > 20 ~ "High"
    ),
    APACHE_score = case_when(
      apache_score <= 50 ~ "Mild",
      apache_score > 50 & apache_score <= 75 ~ "Moderate",
      apache_score > 75 ~ "Severe"
    ),
    PaO2_status = case_when(
      APACHE_PAO2 < 80 ~ "Hypoxia",
      APACHE_PAO2 <= 100 ~ "Normoxia",
      APACHE_PAO2 > 100 ~ "Hyperoxia"
    ),
    FiO2_status = case_when(
      APACHE_FIO2 <= 0.21 ~ "Room air",
      APACHE_FIO2 <= 0.50 ~ "Moderate FiO₂",
      APACHE_FIO2 > 0.50 ~ "High FiO₂"
    ),
    PaCO2_status = case_when(
      APACHE_PACO2 < 35 ~ "Hypocapnia",
      APACHE_PACO2 <= 45 ~ "Normocapnia",
      APACHE_PACO2 > 45 ~ "Hypercapnia"
    )
  )

# 3. Totals for Fisher’s tests
total_good <- sum(data$Label == "Good")
total_poor <- sum(data$Label == "Poor")

# 4. Range lookup (updated level to match shortened 'Normal')
range_map <- list(
  BMI = tibble(level = c("Underweight","Normal","Overweight","Obese"), range = c("< 18.5","18.5–24.9","25–29.9",">= 30")),
  MAP = tibble(level = c("Low","Normal","High"), range = c("< 70","70–100","> 100")),
  Pulse = tibble(level = c("Low","Normal","High"), range = c("< 60","60–100","> 100")),
  Respiratory_Rate = tibble(level = c("Low","Normal","High"), range = c("< 12","12–20","> 20")),
  APACHE_score = tibble(level = c("Mild","Moderate","Severe"), range = c("0–50","51–75","> 75")),
  PaO2_status = tibble(level = c("Hypoxia","Normoxia","Hyperoxia"), range = c("< 80","80–100","> 100")),
  FiO2_status = tibble(level = c("Room air","Moderate FiO₂","High FiO₂"), range = c("<= 0.21","0.22–0.50","> 0.50")),
  PaCO2_status = tibble(level = c("Hypocapnia","Normocapnia","Hypercapnia"), range = c("< 35","35–45","> 45"))
)

# 5. Plot setup
output_dir <- "Meta_results4"
if (!dir.exists(output_dir)) dir.create(output_dir)
custom_colors <- c("Poor" = "#666666", "Good" = "#CCCCCC")
cats <- names(range_map)
bar_width <- 0.5

# 6. Generate plots
for (cat_var in cats) {
  pvals <- data %>%
    group_by(.data[[cat_var]]) %>%
    summarise(in_good = sum(Label == "Good"), in_poor = sum(Label == "Poor"), .groups = "drop") %>%
    rowwise() %>%
    mutate(p.value = fisher.test(matrix(c(in_good, in_poor, total_good - in_good, total_poor - in_poor), nrow = 2))$p.value) %>%
    ungroup()
  
  df_perc <- data %>%
    count(.data[[cat_var]], Label) %>%
    group_by(Label) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup() %>%
    rename(level = .data[[cat_var]]) %>%
    left_join(range_map[[cat_var]], by = "level") %>%
    mutate(level = factor(level, levels = range_map[[cat_var]]$level))
  
  y_p <- df_perc %>%
    group_by(level) %>%
    summarise(y = max(percentage) + 1, .groups = "drop")
  
  ann <- pvals %>%
    rename(level = .data[[cat_var]]) %>%
    left_join(y_p, by = "level") %>%
    mutate(label = paste0("p = ", signif(p.value, 2)))
  
  unit_label <- switch(cat_var,
                       BMI = " (kg/m²)",
                       MAP = " (mmHg)",
                       Pulse = " (bpm)",
                       Respiratory_Rate = " (b/m)",
                       APACHE_score = "",
                       PaO2_status = " (mmHg)",
                       FiO2_status = "",
                       PaCO2_status = " (mmHg)"
  )
  
  x_lab_expr <- switch(cat_var,
                       BMI = bquote(bold(BMI) * .(unit_label)),
                       MAP = bquote(bold(MAP) * .(unit_label)),
                       Pulse = bquote(bold(Pulse) * .(unit_label)),
                       Respiratory_Rate = bquote(bold("Respiratory Rate") * .(unit_label)),
                       APACHE_score = bquote(bold("APACHE Score") * .(unit_label)),
                       PaO2_status = bquote(bold(Pa*O[2]) * .(unit_label)),
                       FiO2_status = bquote(bold(Fi*O[2]) * .(unit_label)),
                       PaCO2_status = bquote(bold(Pa*CO[2]) * .(unit_label))
  )
  
  plot_title <- paste0(gsub("_", " ", cat_var), "")
  
  label_df <- df_perc %>% distinct(level, range) %>% mutate(y = 3)
  
  p <- ggplot(df_perc, aes(x = level, y = percentage, fill = Label)) +
    geom_col(width = bar_width, position = position_dodge(width = bar_width)) +
    scale_x_discrete(expand = expansion(add = c(bar_width * 1.5, bar_width / 2))) +
    scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = custom_colors) +
    labs(title = plot_title, x = x_lab_expr, y = expression(bold(Percentage)), fill = "Label") +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA)
    ) +
    geom_text(data = label_df, aes(x = level, y = y, label = range), inherit.aes = FALSE, size = 4, color = "black") +
    geom_text(data = ann, aes(x = level, y = y, label = label), inherit.aes = FALSE, size = 5, vjust = -0.3)
  
  ggsave(file.path(output_dir, paste0(cat_var, "_with_ranges_and_p.png")), plot = p, dpi = 600, width = 6, height = 4, bg = "transparent")
}

message("Done: eight plots saved with centered range labels at bar base and p-values above.")

