rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(ggpubr)
  library(rstatix)
})

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)
load("Inputs/1211_metadata.rdata")
gsva_mat <- read.csv(
  file = "Outputs/Supplementary_data_6_module_gsva.csv",
  row.names = 1,
  check.names = FALSE
)


out_dir <- file.path("/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260620_module_inflammatory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_labels <- c(
  up1 = "Inflammatory signaling",
  up2 = "Phagolysosome function",
  up3 = "IFN-1 signaling",
  dw1 = "Ribosome biogenesis"
)

timepoint_levels <- c("D1", "D4", "D7", "D14", "D21")

meta_gsva <- df_long %>%
  dplyr::select(HumanID, Timepoint, SampleID, Mortality28d) %>%
  distinct() %>%
  filter(SampleID %in% colnames(gsva_mat)) %>%
  arrange(match(SampleID, colnames(gsva_mat))) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoint_levels),
    Mortality28d = factor(Mortality28d, levels = c("0", "1"))
  )

missing_meta <- setdiff(colnames(gsva_mat), meta_gsva$SampleID)
if (length(missing_meta) > 0) {
  warning(
    "The following GSVA samples are missing from df_long and were dropped: ",
    paste(missing_meta, collapse = ", ")
  )
}

gsva_plot_mat <- gsva_mat[, meta_gsva$SampleID, drop = FALSE]
stopifnot(identical(meta_gsva$SampleID, colnames(gsva_plot_mat)))

plot_df <- bind_cols(meta_gsva, as_tibble(t(gsva_plot_mat), .name_repair = "minimal"))

plot_long <- plot_df %>%
  pivot_longer(
    cols = all_of(names(module_labels)),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = factor(Module, levels = names(module_labels), labels = module_labels),
    MortalityGroup = factor(
      Mortality28d,
      levels = c("0", "1"),
      labels = c("Survival", "Mortality")
    )
  )

write.csv(
  plot_long,
  file = file.path(out_dir, "module_gsva_plot_long.csv"),
  row.names = FALSE
)

group_summary <- plot_long %>%
  group_by(Module, Timepoint, Mortality28d, MortalityGroup) %>%
  summarise(
    n = sum(!is.na(Score)),
    mean = mean(Score, na.rm = TRUE),
    sd = sd(Score, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

wilcox_df <- plot_long %>%
  filter(!is.na(Score), !is.na(Mortality28d)) %>%
  group_by(Module, Timepoint) %>%
  filter(n_distinct(Mortality28d) == 2) %>%
  wilcox_test(Score ~ Mortality28d) %>%
  ungroup() %>%
  mutate(
    p.label = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

y_pos_df <- group_summary %>%
  group_by(Module, Timepoint) %>%
  summarise(
    y.position = max(mean + se, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    plot_long %>%
      group_by(Module, Timepoint) %>%
      summarise(
        rng = diff(range(Score, na.rm = TRUE)),
        .groups = "drop"
      ),
    by = c("Module", "Timepoint")
  ) %>%
  mutate(
    rng = if_else(is.finite(rng) & rng > 0, rng, 0.05),
    y.position = y.position + 0.035 * rng
  ) %>%
  dplyr::select(Module, Timepoint, y.position)

wilcox_df <- wilcox_df %>%
  left_join(y_pos_df, by = c("Module", "Timepoint"))

write.csv(
  group_summary,
  file = file.path(out_dir, "module_gsva_group_summary.csv"),
  row.names = FALSE
)
write.csv(
  wilcox_df,
  file = file.path(out_dir, "module_gsva_wilcox_by_timepoint.csv"),
  row.names = FALSE
)

plot_one_module <- function(module_name) {
  plot_data <- plot_long %>% filter(Module == module_name)
  stat_data <- wilcox_df %>% filter(Module == module_name)

  ggplot(
    plot_data,
    aes(
      x = Timepoint,
      y = Score,
      color = Mortality28d,
      group = interaction(HumanID, Mortality28d)
    )
  ) +
    annotate(
      "rect",
      xmin = 1,
      xmax = 2,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey",
      alpha = 0.28
    ) +
    stat_summary(
      aes(group = Mortality28d),
      fun = mean,
      geom = "line",
      linewidth = 1.3
    ) +
    stat_summary(
      aes(group = Mortality28d),
      fun = mean,
      geom = "point",
      size = 3
    ) +
    stat_summary(
      aes(group = Mortality28d),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15,
      linewidth = 0.6
    ) +
    stat_pvalue_manual(
      stat_data,
      label = "p.label",
      x = "Timepoint",
      y.position = "y.position",
      tip.length = 0,
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c("0" = "#4575B4", "1" = "#D73027"),
      name = NULL,
      labels = c("0" = "Survival", "1" = "Mortality")
    ) +
    scale_x_discrete(expand = expansion(add = c(0.12, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    labs(x = "Timepoint", y = "Module GSVA score", title = module_name) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = "top"
    )
}

module_plot_names <- set_names(names(module_labels), module_labels)

module_plots <- imap(
  module_plot_names,
  function(module_code, module_name) {
    p <- plot_one_module(module_name)
    ggsave(
      filename = file.path(out_dir, paste0("module_gsva_", module_code, "_mortality_timepoint_lineplot.pdf")),
      plot = p,
      width = 3.3,
      height = 3.1
    )
    p
  }
)

d1_plot_long <- plot_long %>%
  filter(Timepoint == "D1")

d1_p_df <- wilcox_df %>%
  filter(Timepoint == "D1") %>%
  left_join(
    d1_plot_long %>%
      group_by(Module) %>%
      summarise(
        y.position = max(Score, na.rm = TRUE) +
          0.06 * diff(range(Score, na.rm = TRUE)),
        .groups = "drop"
      ),
    by = "Module",
    suffix = c("", ".d1")
  ) %>%
  mutate(y.position = coalesce(y.position.d1, y.position)) %>%
  dplyr::select(Module, p, p.label, y.position)

p_d1_box <- ggplot(
  d1_plot_long,
  aes(x = Module, y = Score, fill = Mortality28d, color = Mortality28d)
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.28,
    linewidth = 0.55,
    position = position_dodge(width = 0.72)
  ) +
  geom_boxplot(
    width = 0.65,
    linewidth = 0.55,
    outlier.shape = 21,
    outlier.size = 1.8,
    position = position_dodge(width = 0.72),
    alpha = 0.9
  ) +
  geom_text(
    data = d1_p_df,
    aes(x = Module, y = y.position, label = p.label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.2
  ) +
  scale_fill_manual(
    values = c("0" = "#4575B4", "1" = "#D73027"),
    name = NULL,
    labels = c("0" = "Survival", "1" = "Mortality")
  ) +
  scale_color_manual(
    values = c("0" = "#4575B4", "1" = "#D73027"),
    name = NULL,
    labels = c("0" = "Survival", "1" = "Mortality")
  ) +
  scale_x_discrete(
    labels = c(
      "Inflammatory signaling" = "Inflammatory\nsignaling",
      "Phagolysosome function" = "Phagolysosome\nfunction",
      "IFN-1 signaling" = "IFN-1\nsignaling",
      "Ribosome biogenesis" = "Ribosome\nbiogenesis"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    x = NULL,
    y = "GSVA Score",
    title = "D1 Host Transcriptomic GSVA Scores"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 12),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "top"
  )

p_d1_box

ggsave(
  filename = file.path(out_dir, "module_gsva_D1_mortality_boxplot.pdf"),
  plot = p_d1_box,
  width = 5.2,
  height = 3.6
)
ggsave(
  filename = file.path(out_dir, "module_gsva_D1_mortality_boxplot.png"),
  plot = p_d1_box,
  width = 5.2,
  height = 3.6,
  dpi = 300
)
