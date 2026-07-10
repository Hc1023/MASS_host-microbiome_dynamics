rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

input_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
output_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs"
out_dir <- file.path(output_dir, "260618_metadata_overview")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(input_dir, "Inputs/1211_metadata.rdata"))

df_right <- df_long %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7", "D14", "D21")),
    Mortality28d = factor(Mortality28d, levels = c("0", "1"))
  ) %>%
  count(Timepoint, Mortality28d)

df_total <- df_right %>%
  summarise(n_total = sum(n), .by = Timepoint)

p_right <- ggplot(df_right, aes(x = Timepoint, y = n, fill = Mortality28d)) +
  geom_col(
    alpha = 0.9,
    color = "white",
    linewidth = 0.4,
    width = 0.6
  ) +
  geom_text(
    data = df_total,
    aes(x = Timepoint, y = n_total, label = paste0("n=", n_total)),
    inherit.aes = FALSE,
    vjust = -0.8,
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_manual(
    values = c("0" = "#4575B4", "1" = "#D73027"),
    labels = c("0" = "Survival", "1" = "Mortality")
  ) +
  labs(x = "Study Timepoint", y = "Number of Samples", fill = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.78, 0.9),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size = unit(0.45, "cm"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  scale_y_continuous(
    limits = c(0, 450),
    breaks = seq(0, 450, 100),
    expand = expansion(mult = c(0, 0))
  )

p_right

# write_csv(df_right, file.path(out_dir, "p_right_sample_counts.csv"))
ggsave(
  filename = file.path(out_dir, "p_right_sample_counts.pdf"),
  plot = p_right,
  width = 3,
  height = 3.7
)
