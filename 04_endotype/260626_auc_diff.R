rm(list = ls())

library(tidyverse)

out_dir <- file.path(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260626_auc_diff"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

auc_df <- read.csv(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_auc_sel/plot/plotting_data_auc_feature_space_summary.csv",
  check.names = FALSE
)

#### write auc table ####

auc_label_wide <- auc_df %>%
  dplyr::select(task_label, feature_space_label, auc_label) %>%
  pivot_wider(
    id_cols = task_label,
    names_from = feature_space_label,
    values_from = auc_label
  )

auc_gain_wide <- auc_df %>%
  dplyr::select(task_label, feature_space_label, median_auc) %>%
  pivot_wider(
    id_cols = task_label,
    names_from = feature_space_label,
    values_from = median_auc
  ) %>%
  transmute(
    task_label,
    `Incremental AUC Gain` = `CTS + Microbes` - `APACHE II`
  )

auc_df2 <- auc_label_wide %>%
  left_join(auc_gain_wide, by = "task_label") %>%
  mutate(`Incremental AUC Gain` = sprintf("%.3f", `Incremental AUC Gain`))

average_horizon_row <- auc_df %>%
  group_by(feature_space_label) %>%
  summarise(median_auc = median(median_auc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = feature_space_label,
    values_from = median_auc
  ) %>%
  mutate(
    task_label = "Average Performance across Horizons",
    `Incremental AUC Gain` = sprintf("%.3f", `CTS + Microbes` - `APACHE II`),
    across(
      -c(task_label, `Incremental AUC Gain`),
      ~sprintf("%.3f", .x)
    )
  ) %>%
  dplyr::select(all_of(names(auc_df2)))

auc_df2 <- bind_rows(auc_df2, average_horizon_row)

write.csv(
  auc_df2,
  file.path(out_dir, "auc_feature_space_wide.csv"),
  row.names = FALSE
)

#### draw auc gain plot ####

task_levels <- c("D1 3-day", "D1 7-day", "D1 28-day", "D4 28-day", "D7 28-day")
task_axis_labels <- c(
  "D1 3-day" = "D1\n3-day",
  "D1 7-day" = "D1\n7-day",
  "D1 28-day" = "D1\n28-day",
  "D4 28-day" = "D4\n28-day",
  "D7 28-day" = "D7\n28-day"
)

auc_gain_plot_df <- auc_gain_wide %>%
  mutate(
    task_label = factor(task_label, levels = task_levels),
    gain_label = sprintf("%+.3f", `Incremental AUC Gain`)
  )

p_auc_gain <- ggplot(
  auc_gain_plot_df,
  aes(x = task_label, y = `Incremental AUC Gain`)
) +
  geom_col(width = 0.62, fill = "#DF6B21") +
  geom_text(
    aes(label = gain_label),
    vjust = -0.25,
    size = 4.2
  ) +
  scale_x_discrete(labels = task_axis_labels) +
  scale_y_continuous(
    limits = c(0, max(auc_gain_plot_df$`Incremental AUC Gain`, na.rm = TRUE) * 1.25),
    breaks = seq(0, 0.25, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    # title = "Integrated model\nAUC difference",
    x = NULL,
    y = "Gain over APACHE II"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, lineheight = 0.95),
    # axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "grey45", linewidth = 0.6),
    axis.ticks = element_line(color = "grey45", linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(
  file.path(out_dir, "auc_incremental_gain_over_apacheii.pdf"),
  p_auc_gain,
  width = 3.8,
  height = 4.0
)

#### rank plot ####

auc_model_df <- read.csv(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_auc_sel/plot/plotting_data_auc_model_feature_space_summary.csv",
  check.names = FALSE
)

rank_feature_levels <- c("APACHE II", "Microbes", "CTS", "CTS + Microbes")
rank_algorithm_levels <- c("glmnet", "svm", "randomForest", "xgboost")
rank_algorithm_labels <- c(
  "glmnet" = "GLMNET",
  "svm" = "SVM",
  "randomForest" = "Random Forest",
  "xgboost" = "XGBoost"
)

auc_rank_long <- auc_model_df %>%
  mutate(
    feature_space_label = factor(feature_space_label, levels = rank_feature_levels),
    algorithm = factor(algorithm, levels = rank_algorithm_levels)
  )

auc_rank_summary <- auc_rank_long %>%
  group_by(algorithm, feature_space_label) %>%
  summarise(
    median_AUC = median(median_auc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(algorithm) %>%
  arrange(desc(median_AUC), feature_space_label, .by_group = TRUE) %>%
  mutate(median_AUC_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    algorithm_label = factor(
      rank_algorithm_labels[as.character(algorithm)],
      levels = rev(rank_algorithm_labels[rank_algorithm_levels])
    ),
    feature_space_label = factor(feature_space_label, levels = rank_feature_levels),
    rank_label = paste0(
      sprintf("%.3f", median_AUC), "; ",
      "rank #", median_AUC_rank
    )
  )

auc_rank_summary_wide = auc_rank_summary %>% 
  pivot_wider(id_cols = algorithm_label,
              names_from = feature_space_label,
              values_from = rank_label
              )
auc_rank_summary_wide = auc_rank_summary_wide[,c(1:3,5,4)]

write.csv(
  auc_rank_summary_wide,
  file.path(out_dir, "auc_median_rank_by_model_feature.csv"),
  row.names = FALSE
)

p_auc_rank <- ggplot(
  auc_rank_summary,
  aes(x = feature_space_label, y = algorithm_label, fill = median_AUC)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = rank_label), size = 3.0, lineheight = 0.9) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#DF6B21",
    name = "Median\nAUC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(
  file.path(out_dir, "auc_median_rank_heatmap_by_model_feature.pdf"),
  p_auc_rank,
  width = 5.4,
  height = 2.8
)

#### results excluding SVM ####

auc_no_svm_df <- read.csv(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_auc_sel/plot/plotting_data_auc_total_excluding_svm.csv",
  check.names = FALSE
)

auc_no_svm_label_wide <- auc_no_svm_df %>%
  dplyr::select(task_label, feature_space_label, auc_label) %>%
  pivot_wider(
    id_cols = task_label,
    names_from = feature_space_label,
    values_from = auc_label
  )

auc_no_svm_gain_wide <- auc_no_svm_df %>%
  dplyr::select(task_label, feature_space_label, median_auc) %>%
  pivot_wider(
    id_cols = task_label,
    names_from = feature_space_label,
    values_from = median_auc
  ) %>%
  transmute(
    task_label,
    `Incremental AUC Gain` = `CTS + Microbes` - `APACHE II`
  )

auc_no_svm_table <- auc_no_svm_label_wide %>%
  left_join(auc_no_svm_gain_wide, by = "task_label") %>%
  mutate(`Incremental AUC Gain` = sprintf("%.3f", `Incremental AUC Gain`))

average_no_svm_horizon_row <- auc_no_svm_df %>%
  group_by(feature_space_label) %>%
  summarise(median_auc = median(median_auc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = feature_space_label,
    values_from = median_auc
  ) %>%
  mutate(
    task_label = "Average Performance across Horizons",
    `Incremental AUC Gain` = sprintf("%.3f", `CTS + Microbes` - `APACHE II`),
    across(
      -c(task_label, `Incremental AUC Gain`),
      ~sprintf("%.3f", .x)
    )
  ) %>%
  dplyr::select(all_of(names(auc_no_svm_table)))

auc_no_svm_table <- bind_rows(auc_no_svm_table, average_no_svm_horizon_row)

write.csv(
  auc_no_svm_table,
  file.path(out_dir, "auc_feature_space_wide_excluding_svm.csv"),
  row.names = FALSE
)

auc_no_svm_gain_plot_df <- auc_no_svm_gain_wide %>%
  mutate(
    task_label = factor(task_label, levels = task_levels),
    gain_label = sprintf("%+.3f", `Incremental AUC Gain`),
    label_vjust = if_else(`Incremental AUC Gain` >= 0, -0.25, 1.25)
  )

auc_no_svm_gain_range <- range(
  c(0, auc_no_svm_gain_plot_df$`Incremental AUC Gain`),
  na.rm = TRUE
)
auc_no_svm_gain_padding <- diff(auc_no_svm_gain_range) * 0.12
auc_no_svm_gain_limits <- auc_no_svm_gain_range +
  c(-auc_no_svm_gain_padding, auc_no_svm_gain_padding)

p_auc_no_svm_gain <- ggplot(
  auc_no_svm_gain_plot_df,
  aes(x = task_label, y = `Incremental AUC Gain`)
) +
  geom_col(width = 0.62, fill = "#DF6B21") +
  geom_text(aes(label = gain_label, vjust = label_vjust), size = 4.2) +
  scale_x_discrete(labels = task_axis_labels) +
  scale_y_continuous(
    limits = auc_no_svm_gain_limits,
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Gain over APACHE II") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, color = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "grey45", linewidth = 0.6),
    axis.ticks = element_line(color = "grey45", linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(
  file.path(out_dir, "auc_incremental_gain_over_apacheii_excluding_svm.pdf"),
  p_auc_no_svm_gain,
  width = 3.8,
  height = 4.0,
  device = cairo_pdf
)

auc_rank_summary_no_svm <- auc_rank_summary %>%
  filter(as.character(algorithm) != "svm") %>%
  droplevels()

auc_rank_summary_no_svm_wide <- auc_rank_summary_no_svm %>%
  pivot_wider(
    id_cols = algorithm_label,
    names_from = feature_space_label,
    values_from = rank_label
  )
auc_rank_summary_no_svm_wide <- auc_rank_summary_no_svm_wide[, c(1:3, 5, 4)]

write.csv(
  auc_rank_summary_no_svm_wide,
  file.path(out_dir, "auc_median_rank_by_model_feature_excluding_svm.csv"),
  row.names = FALSE
)

p_auc_rank_no_svm <- ggplot(
  auc_rank_summary_no_svm,
  aes(x = feature_space_label, y = algorithm_label, fill = median_AUC)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = rank_label), size = 3.0, lineheight = 0.9) +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#DF6B21",
    name = "Median\nAUC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(
  file.path(out_dir, "auc_median_rank_heatmap_by_model_feature_excluding_svm.pdf"),
  p_auc_rank_no_svm,
  width = 5.4,
  height = 2.4,
  device = cairo_pdf
)
