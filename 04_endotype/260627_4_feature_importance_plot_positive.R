rm(list = ls())

library(tidyverse)

out_dir <- file.path(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260627_feature_positive/"
)
plot_dir <- file.path(out_dir, "plot")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path("/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_feature_sel/feature_importance_permutation_summary.csv")
if (!file.exists(summary_file)) {
  stop("Missing feature importance summary: ", summary_file)
}

importance_summary <- read.csv(summary_file, check.names = FALSE)

parse_env_list <- function(var, default) {
  value <- Sys.getenv(var, default)
  strsplit(value, ",", fixed = TRUE)[[1]] %>%
    trimws() %>%
    discard(~.x == "")
}

forecast_tasks <- tibble::tribble(
  ~task_id,  ~task_label,  ~day,  ~horizon_days,
  "D1_M3",  "D1 3-day",   "D1",  3,
  "D1_M7",  "D1 7-day",   "D1",  7,
  "D1_M28", "D1 28-day",  "D1",  28,
  "D4_M28", "D4 28-day",  "D4",  28,
  "D7_M28", "D7 28-day",  "D7",  28
)

importance_algorithms <- tibble::tribble(
  ~algorithm,      ~algorithm_label,
  "randomForest", "Random Forest",
  "xgboost",      "XGBoost",
  "svm",          "SVM-RBF",
  "glmnet",       "Elastic Net"
)

task_levels <- forecast_tasks$task_id
task_labels <- setNames(forecast_tasks$task_label, forecast_tasks$task_id)
algorithm_levels <- importance_algorithms$algorithm
algorithm_labels <- setNames(
  importance_algorithms$algorithm_label,
  importance_algorithms$algorithm
)

plot_top_n <- as.integer(Sys.getenv("FI_PLOT_TOP_N", "8"))
main_plot_top_n <- as.integer(Sys.getenv("FI_MAIN_PLOT_TOP_N", "10"))
plot_excluded_features <- parse_env_list(
  "FI_PLOT_EXCLUDE_FEATURES",
  "Total,Viruses,Bacteria/Fungi"
)

integrated_frequency_by_task <- importance_summary %>%
  filter(!feature_label %in% plot_excluded_features) %>%
  group_by(task_id, task_label, day, horizon_days, feature_label, feature_type) %>%
  summarise(
    integrated_positive_importance_frequency = mean(positive_importance_frequency, na.rm = TRUE),
    max_model_positive_importance_frequency = max(positive_importance_frequency, na.rm = TRUE),
    median_importance_across_models = median(median_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(task_id) %>%
  arrange(
    desc(integrated_positive_importance_frequency),
    desc(max_model_positive_importance_frequency),
    desc(median_importance_across_models),
    feature_label,
    .by_group = TRUE
  ) %>%
  mutate(integrated_frequency_rank = row_number()) %>%
  ungroup()

main_plot_features <- integrated_frequency_by_task %>%
  group_by(feature_label, feature_type) %>%
  summarise(
    overall_integrated_frequency = mean(integrated_positive_importance_frequency, na.rm = TRUE),
    max_integrated_frequency = max(integrated_positive_importance_frequency, na.rm = TRUE),
    max_model_frequency = max(max_model_positive_importance_frequency, na.rm = TRUE),
    max_importance = max(median_importance_across_models, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(overall_integrated_frequency),
    desc(max_integrated_frequency),
    desc(max_model_frequency),
    desc(max_importance),
    feature_label
  ) %>%
  slice_head(n = main_plot_top_n) %>%
  mutate(
    overall_rank = row_number(),
    feature_type_label = recode(feature_type, "CTSg" = "CTS Gene"),
    feature_display = paste0(overall_rank, ". ", feature_label)
  )

main_frequency_plot_df <- tidyr::expand_grid(
  main_plot_features %>%
    select(feature_label, feature_type, feature_type_label, feature_display, overall_rank),
  forecast_tasks %>%
    select(task_id, task_label, day, horizon_days)
) %>%
  left_join(
    integrated_frequency_by_task,
    by = c("task_id", "task_label", "day", "horizon_days",
           "feature_label", "feature_type")
  ) %>%
  mutate(
    task_id = factor(task_id, levels = task_levels),
    task_label = factor(task_labels[as.character(task_id)], levels = task_labels),
    task_order = as.integer(task_id),
    feature_display = factor(
      feature_display,
      levels = rev(main_plot_features$feature_display)
    )
  ) %>%
  arrange(overall_rank, task_id)

main_freq_wide <- main_frequency_plot_df %>% 
  dplyr::select(
    feature_label,
    feature_type_label,
    task_label,
    integrated_positive_importance_frequency
  ) %>%
  pivot_wider(
    id_cols = c(feature_label, feature_type_label),
    names_from = task_label,
    values_from = integrated_positive_importance_frequency
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write.csv(main_freq_wide, file = file.path(out_dir, "freq_wide.csv"), row.names = F)

feature_type_annotation_df <- main_frequency_plot_df %>%
  distinct(feature_display, feature_type_label)

p_integrated_top10_rank <- ggplot(
  main_frequency_plot_df,
  aes(x = task_order, y = feature_display, fill = integrated_frequency_rank)
) +
  geom_point(
    data = feature_type_annotation_df,
    aes(x = 0.55, y = feature_display, color = feature_type_label),
    inherit.aes = FALSE,
    shape = 15,
    size = 5.2
  ) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = if_else(
      is.na(integrated_frequency_rank),
      "",
      as.character(integrated_frequency_rank)
    )),
    size = 2.8
  ) +
  scale_fill_gradientn(
    colours = c("#9E0142", "#D7191C", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF"),
    limits = c(1, max(integrated_frequency_by_task$integrated_frequency_rank, na.rm = TRUE)),
    na.value = "grey95",
    name = "Integrated\nRank"
  ) +
  scale_color_manual(
    values = c("CTS Gene" = "#1F78B4", "Microbe" = "#FDAE61"),
    name = "Feature type"
  ) +
  scale_x_continuous(
    breaks = seq_along(task_labels),
    labels = task_labels,
    limits = c(0.35, length(task_labels) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    # title = "Overall Top 10 Integrated Feature Importance Rank",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

top_features_for_plot <- importance_summary %>%
  filter(!feature_label %in% plot_excluded_features) %>%
  group_by(feature_space, feature_space_label, task_id, task_label, feature, feature_label) %>%
  summarise(
    max_median_importance = max(median_importance, na.rm = TRUE),
    max_positive_importance_frequency = max(positive_importance_frequency, na.rm = TRUE),
    mean_median_importance = mean(median_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(feature_space, task_id) %>%
  arrange(
    desc(max_median_importance),
    desc(max_positive_importance_frequency),
    feature_label,
    .by_group = TRUE
  ) %>%
  slice_head(n = plot_top_n) %>%
  mutate(
    feature_order = row_number(),
    feature_panel_id = paste(feature_space, task_id, feature, sep = "___")
  ) %>%
  ungroup()

importance_plot_df <- importance_summary %>%
  inner_join(
    top_features_for_plot %>%
      select(feature_space, task_id, feature, feature_order, feature_panel_id),
    by = c("feature_space", "task_id", "feature")
  ) %>%
  mutate(
    task_id = factor(task_id, levels = task_levels),
    task_label = factor(task_labels[as.character(task_id)], levels = task_labels),
    algorithm = factor(algorithm, levels = algorithm_levels),
    algorithm_label = factor(
      algorithm_labels[as.character(algorithm)],
      levels = algorithm_labels
    ),
    feature_panel_id = factor(
      feature_panel_id,
      levels = top_features_for_plot %>%
        arrange(feature_space, task_id, desc(feature_order)) %>%
        pull(feature_panel_id)
    )
  ) %>%
  arrange(feature_space, task_id, feature_order, algorithm)

feature_panel_labels <- importance_plot_df %>%
  distinct(feature_panel_id, feature_label) %>%
  deframe()

p_importance_positive_freq <- ggplot(
  importance_plot_df,
  aes(x = algorithm_label, y = feature_panel_id, fill = positive_importance_frequency)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", positive_importance_frequency)), size = 2.3) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    name = "Positive Importance\nFrequency"
  ) +
  scale_y_discrete(labels = feature_panel_labels) +
  facet_wrap(vars(feature_space_label, task_label), nrow = 2, scales = "free_y") +
  labs(
    title = "Positive Permutation Importance Frequency",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

#### results excluding SVM ####

importance_summary_no_svm <- importance_summary %>%
  filter(algorithm != "svm")

integrated_positive_by_task_no_svm <- importance_summary_no_svm %>%
  filter(!feature_label %in% plot_excluded_features) %>%
  group_by(task_id, task_label, day, horizon_days, feature_label, feature_type) %>%
  summarise(
    integrated_positive_importance_frequency = mean(positive_importance_frequency, na.rm = TRUE),
    max_model_positive_importance_frequency = max(positive_importance_frequency, na.rm = TRUE),
    median_importance_across_models = median(median_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(task_id) %>%
  arrange(
    desc(integrated_positive_importance_frequency),
    desc(max_model_positive_importance_frequency),
    desc(median_importance_across_models),
    feature_label,
    .by_group = TRUE
  ) %>%
  mutate(integrated_frequency_rank = row_number()) %>%
  ungroup()

main_positive_features_no_svm <- integrated_positive_by_task_no_svm %>%
  group_by(feature_label, feature_type) %>%
  summarise(
    overall_integrated_frequency = mean(integrated_positive_importance_frequency, na.rm = TRUE),
    max_integrated_frequency = max(integrated_positive_importance_frequency, na.rm = TRUE),
    max_model_frequency = max(max_model_positive_importance_frequency, na.rm = TRUE),
    max_importance = max(median_importance_across_models, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(overall_integrated_frequency),
    desc(max_integrated_frequency),
    desc(max_model_frequency),
    desc(max_importance),
    feature_label
  ) %>%
  slice_head(n = main_plot_top_n) %>%
  mutate(
    overall_rank = row_number(),
    feature_type_label = recode(feature_type, "CTSg" = "CTS Gene"),
    feature_display = paste0(overall_rank, ". ", feature_label)
  )

main_positive_plot_no_svm_df <- tidyr::expand_grid(
  main_positive_features_no_svm %>%
    select(feature_label, feature_type, feature_type_label, feature_display, overall_rank),
  forecast_tasks %>%
    select(task_id, task_label, day, horizon_days)
) %>%
  left_join(
    integrated_positive_by_task_no_svm,
    by = c("task_id", "task_label", "day", "horizon_days",
           "feature_label", "feature_type")
  ) %>%
  mutate(
    task_id = factor(task_id, levels = task_levels),
    task_label = factor(task_labels[as.character(task_id)], levels = task_labels),
    task_order = as.integer(task_id),
    feature_display = factor(
      feature_display,
      levels = rev(main_positive_features_no_svm$feature_display)
    )
  ) %>%
  arrange(overall_rank, task_id)

main_positive_wide_no_svm <- main_positive_plot_no_svm_df %>%
  dplyr::select(
    feature_label,
    feature_type_label,
    task_label,
    integrated_positive_importance_frequency
  ) %>%
  pivot_wider(
    id_cols = c(feature_label, feature_type_label),
    names_from = task_label,
    values_from = integrated_positive_importance_frequency
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

write.csv(
  main_positive_wide_no_svm,
  file.path(out_dir, "freq_wide_excluding_svm.csv"),
  row.names = FALSE
)

positive_type_annotation_no_svm_df <- main_positive_plot_no_svm_df %>%
  distinct(feature_display, feature_type_label)

p_integrated_positive_rank_no_svm <- ggplot(
  main_positive_plot_no_svm_df,
  aes(x = task_order, y = feature_display, fill = integrated_frequency_rank)
) +
  geom_point(
    data = positive_type_annotation_no_svm_df,
    aes(x = 0.55, y = feature_display, color = feature_type_label),
    inherit.aes = FALSE,
    shape = 15,
    size = 5.2
  ) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = if_else(
      is.na(integrated_frequency_rank),
      "",
      as.character(integrated_frequency_rank)
    )),
    size = 2.8
  ) +
  scale_fill_gradientn(
    colours = c("#9E0142", "#D7191C", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF"),
    limits = c(1, max(integrated_positive_by_task_no_svm$integrated_frequency_rank, na.rm = TRUE)),
    na.value = "grey95",
    name = "Integrated\nRank"
  ) +
  scale_color_manual(
    values = c("CTS Gene" = "#1F78B4", "Microbe" = "#FDAE61"),
    name = "Feature type"
  ) +
  scale_x_continuous(
    breaks = seq_along(task_labels),
    labels = task_labels,
    limits = c(0.35, length(task_labels) + 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

top_positive_features_no_svm <- importance_summary_no_svm %>%
  filter(!feature_label %in% plot_excluded_features) %>%
  group_by(feature_space, feature_space_label, task_id, task_label, feature, feature_label) %>%
  summarise(
    max_median_importance = max(median_importance, na.rm = TRUE),
    max_positive_importance_frequency = max(positive_importance_frequency, na.rm = TRUE),
    mean_median_importance = mean(median_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(feature_space, task_id) %>%
  arrange(
    desc(max_median_importance),
    desc(max_positive_importance_frequency),
    feature_label,
    .by_group = TRUE
  ) %>%
  slice_head(n = plot_top_n) %>%
  mutate(
    feature_order = row_number(),
    feature_panel_id = paste(feature_space, task_id, feature, sep = "___")
  ) %>%
  ungroup()

algorithm_levels_no_svm <- setdiff(algorithm_levels, "svm")

positive_importance_plot_no_svm_df <- importance_summary_no_svm %>%
  inner_join(
    top_positive_features_no_svm %>%
      select(feature_space, task_id, feature, feature_order, feature_panel_id),
    by = c("feature_space", "task_id", "feature")
  ) %>%
  mutate(
    task_id = factor(task_id, levels = task_levels),
    task_label = factor(task_labels[as.character(task_id)], levels = task_labels),
    algorithm = factor(algorithm, levels = algorithm_levels_no_svm),
    algorithm_label = factor(
      algorithm_labels[as.character(algorithm)],
      levels = algorithm_labels[algorithm_levels_no_svm]
    ),
    feature_panel_id = factor(
      feature_panel_id,
      levels = top_positive_features_no_svm %>%
        arrange(feature_space, task_id, desc(feature_order)) %>%
        pull(feature_panel_id)
    )
  ) %>%
  arrange(feature_space, task_id, feature_order, algorithm)

positive_panel_no_svm_labels <- positive_importance_plot_no_svm_df %>%
  distinct(feature_panel_id, feature_label) %>%
  deframe()

p_importance_positive_freq_no_svm <- ggplot(
  positive_importance_plot_no_svm_df,
  aes(x = algorithm_label, y = feature_panel_id, fill = positive_importance_frequency)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", positive_importance_frequency)), size = 2.3) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    name = "Positive Importance\nFrequency"
  ) +
  scale_y_discrete(labels = positive_panel_no_svm_labels) +
  facet_wrap(vars(feature_space_label, task_label), nrow = 2, scales = "free_y") +
  labs(
    title = "Positive Permutation Importance Frequency (Excluding SVM)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

write.csv(
  main_positive_plot_no_svm_df,
  file.path(plot_dir, "plotting_data_positive_top10_excluding_svm.csv"),
  row.names = FALSE
)

write.csv(
  positive_importance_plot_no_svm_df,
  file.path(plot_dir, "plotting_data_positive_importance_excluding_svm.csv"),
  row.names = FALSE
)

ggsave(
  file.path(plot_dir, "260627_integrated_feature_top10_rank_heatmap.pdf"),
  p_integrated_top10_rank,
  width = 4,
  height = 3
)
ggsave(
  file.path(plot_dir, "260627_feature_importance_positive_frequency_heatmap.pdf"),
  p_importance_positive_freq,
  width = 11.5,
  height = 5.6
)
ggsave(
  file.path(plot_dir, "260627_integrated_feature_top10_rank_heatmap_excluding_svm.pdf"),
  p_integrated_positive_rank_no_svm,
  width = 4,
  height = 3
)
ggsave(
  file.path(plot_dir, "260627_feature_importance_positive_frequency_heatmap_excluding_svm.pdf"),
  p_importance_positive_freq_no_svm,
  width = 11.5,
  height = 5.6
)

message("Feature importance plots written to: ", plot_dir)
