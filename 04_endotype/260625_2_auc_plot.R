rm(list = ls())

library(tidyverse)

base_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_auc_sel"
plot_dir <- file.path(base_dir, "plot")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

auc_feature_file <- file.path(base_dir, "auc_feature_space_summary.csv")
auc_model_file <- file.path(base_dir, "auc_validation_summary.csv")
auc_long_file <- file.path(base_dir, "auc_validation_long.csv")

auc_feature <- read.csv(auc_feature_file, check.names = FALSE)
auc_model <- read.csv(auc_model_file, check.names = FALSE)
auc_long <- read.csv(auc_long_file, check.names = FALSE)

task_levels <- c("D1_M3", "D1_M7", "D1_M28", "D4_M28", "D7_M28")
task_labels <- c("D1 3-day", "D1 7-day", "D1 28-day", "D4 28-day", "D7 28-day")
feature_levels <- c("CTS + microbe", "CTS", "microbe", "APACHEII")
feature_labels <- c(
  "CTS + microbe" = "CTS + Microbes (Integrated)",
  "CTS" = "CTS (Host endotype)",
  "microbe" = "Microbes",
  "APACHEII" = "APACHE II (Severity score)"
)
feature_cols <- c(
  "CTS + microbe" = "#E41A1C",
  "CTS" = "#1F78B4",
  "microbe" = "#FDAE61",
  "APACHEII" = "#8C8C8C"
)
feature_shapes <- c(
  "CTS + microbe" = 8,
  "CTS" = 16,
  "microbe" = 17,
  "APACHEII" = 15
)
feature_linetypes <- c(
  "CTS + microbe" = "solid",
  "CTS" = "solid",
  "microbe" = "solid",
  "APACHEII" = "dotted"
)

#### auc_plot ####

prepare_auc_plot_df <- function(df) {
  df %>%
    mutate(
      task_id = factor(task_id, levels = task_levels),
      task_label = factor(task_label, levels = task_labels),
      feature_space = factor(feature_space, levels = feature_levels),
      feature_label = factor(
        feature_labels[as.character(feature_space)],
        levels = feature_labels[feature_levels]
      )
    ) %>%
    arrange(task_id, feature_space)
}

make_auc_trend_plot <- function(plot_df,
                                title = NULL,
                                legend_position = c(0.25, 0.16),
                                y_limits = c(0.48, 0.80),
                                y_breaks = seq(0.50, 0.80, by = 0.05)) {
  ggplot(
    plot_df,
    aes(
      x = task_label,
      y = median_auc,
      group = feature_space,
      color = feature_space,
      shape = feature_space,
      linetype = feature_space
    )
  ) +
    geom_ribbon(
      data = filter(plot_df, feature_space == "CTS + microbe"),
      aes(x = task_label, ymin = q25_auc, ymax = q75_auc, fill = feature_space, group = 1),
      inherit.aes = FALSE,
      alpha = 0.10,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.9, stroke = 0.7) +
    scale_color_manual(values = feature_cols, labels = feature_labels) +
    scale_fill_manual(values = feature_cols, guide = "none") +
    scale_shape_manual(values = feature_shapes, labels = feature_labels) +
    scale_linetype_manual(values = feature_linetypes, labels = feature_labels) +
    scale_x_discrete(
      expand = expansion(add = c(0.25, 0.25))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      title = title,
      x = "Sequential Prediction Task",
      y = "Model Performance (Median AUC)",
      color = NULL,
      shape = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid.minor = element_blank(),
      legend.position = legend_position,
      legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "grey80"),
      legend.key = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.75, "cm"),
      legend.spacing.y  = unit(0.02, "cm"),
      legend.margin     = margin(2, 4, 2, 4),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 5.5, r = 18, b = 5.5, l = 5.5)
    )
}

auc_plot_df <- prepare_auc_plot_df(auc_feature)

auc_total_no_d1_svm_df <- auc_long %>%
  filter(!(day == "D1" & algorithm == "svm")) %>%
  group_by(
    task_id, task_label, day, horizon_days,
    feature_space, feature_space_label
  ) %>%
  summarise(
    n_values = sum(!is.na(auc)),
    median_auc = median(auc, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    q25_auc = quantile(auc, 0.25, na.rm = TRUE),
    q75_auc = quantile(auc, 0.75, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE),
    n_total = first(n_total),
    n_errors = sum(!is.na(error)),
    .groups = "drop"
  ) %>%
  mutate(
    auc_label = sprintf("%.3f (%.3f-%.3f)", median_auc, q25_auc, q75_auc)
  ) %>%
  prepare_auc_plot_df()

auc_total_no_svm_df <- auc_long %>%
  filter(algorithm != "svm") %>%
  group_by(
    task_id, task_label, day, horizon_days,
    feature_space, feature_space_label
  ) %>%
  summarise(
    n_values = sum(!is.na(auc)),
    median_auc = median(auc, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    q25_auc = quantile(auc, 0.25, na.rm = TRUE),
    q75_auc = quantile(auc, 0.75, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE),
    n_total = first(n_total),
    n_errors = sum(!is.na(error)),
    .groups = "drop"
  ) %>%
  mutate(
    auc_label = sprintf("%.3f (%.3f-%.3f)", median_auc, q25_auc, q75_auc)
  ) %>%
  prepare_auc_plot_df()

auc_model_plot_df <- auc_model %>%
  prepare_auc_plot_df() %>%
  mutate(
    algorithm_label = factor(
      algorithm_label,
      levels = c("XGBoost", "Random Forest", "SVM-RBF", "Elastic Net")
    )
  )

write.csv(
  auc_plot_df,
  file.path(plot_dir, "plotting_data_auc_feature_space_summary.csv"),
  row.names = FALSE
)

write.csv(
  auc_model_plot_df,
  file.path(plot_dir, "plotting_data_auc_model_feature_space_summary.csv"),
  row.names = FALSE
)

write.csv(
  auc_total_no_d1_svm_df,
  file.path(plot_dir, "plotting_data_auc_total_excluding_d1_svm.csv"),
  row.names = FALSE
)

write.csv(
  auc_total_no_svm_df,
  file.path(plot_dir, "plotting_data_auc_total_excluding_svm.csv"),
  row.names = FALSE
)

p_auc <- make_auc_trend_plot(auc_plot_df)
p_auc
ggsave(
  file.path(plot_dir, "260625_auc_trend.pdf"),
  p_auc,
  width = 4.8,
  height = 3
)

p_auc_total_no_d1_svm <- make_auc_trend_plot(
  auc_total_no_d1_svm_df,
  title = "Overall AUC Trend (D1 Excluding SVM)"
)

ggsave(
  file.path(plot_dir, "260625_auc_trend_overall_excluding_d1_svm.pdf"),
  p_auc_total_no_d1_svm,
  width = 4.8,
  height = 3
)

p_auc_total_no_svm <- make_auc_trend_plot(
  auc_total_no_svm_df,
  title = "Overall AUC Trend (Excluding SVM)"
)

ggsave(
  file.path(plot_dir, "260625_auc_trend_overall_excluding_svm.pdf"),
  p_auc_total_no_svm,
  width = 4.8,
  height = 3
)

for (algorithm_name in levels(auc_model_plot_df$algorithm_label)) {
  model_plot_df <- auc_model_plot_df %>%
    filter(algorithm_label == algorithm_name)

  p_auc_model <- make_auc_trend_plot(
    model_plot_df,
    title = paste0(algorithm_name, " AUC Trend"),
    legend_position = c(0.25, 0.16),
    y_limits = c(0.40, 0.80),
    y_breaks = seq(0.40, 0.80, by = 0.05)
  )

  model_file <- algorithm_name %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")

  ggsave(
    file.path(plot_dir, paste0("260625_auc_trend_", model_file, ".pdf")),
    p_auc_model,
    width = 4.8,
    height = 3
  )
}

message("AUC plots written to: ", plot_dir)
