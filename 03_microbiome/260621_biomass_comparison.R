rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)

out_dir <- file.path("/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260622_microbe_comparison_D1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


load('Inputs/1211_metadata.rdata')
load('Inputs/1616_microbe.rdata')

# D1 mortality comparison: Total, Viruses, EBV and HCMV.
plot_d1_biomass = function(df_long, data, mapid_vec){
  df_d1 = df_long %>%
    filter(Timepoint == "D1") %>%
    mutate(
      Mortality28d = factor(
        Mortality28d,
        levels = c("0", "1"),
        labels = c("Survival", "Mortality")
      )
    )

  data_d1 = data[, df_d1$SampleID, drop = FALSE]
  virus_rows = mapid_vec[rownames(data_d1)] == "Viruses"

  plot_df = df_d1 %>%
    mutate(
      Total = colSums(data_d1),
      Viruses = colSums(data_d1[virus_rows, , drop = FALSE]),
      EBV = as.numeric(data_d1["HHV-4", ]),
      HCMV = as.numeric(data_d1["HCMV", ])
    ) %>%
    select(SampleID, Mortality28d, Total, Viruses, EBV, HCMV) %>%
    pivot_longer(
      cols = c(Total, Viruses, EBV, HCMV),
      names_to = "Feature",
      values_to = "Mass"
    ) %>%
    mutate(
      Feature = factor(Feature, levels = c("Total", "Viruses", "EBV", "HCMV")),
      Mass_log2 = log2(Mass + 1)
    )

  ggplot(plot_df, aes(x = Feature, y = Mass_log2, fill = Mortality28d)) +
    geom_boxplot(
      aes(color = Mortality28d),
      outlier.shape = NA,
      alpha = 0.75,
      width = 0.55,
      position = position_dodge(width = 0.72)
    ) +
    geom_point(
      aes(color = Mortality28d),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72),
      alpha = 0.65,
      size = 1.4
    ) +
    stat_compare_means(
      aes(group = Mortality28d),
      method = "wilcox.test",
      label = "p",
      hide.ns = TRUE,
      color = "#c43932",
      size = 4
    ) +
    scale_color_manual(values = c("Survival" = "#4575B4", "Mortality" = "#D73027")) +
    scale_fill_manual(values = c("Survival" = "#4575B4", "Mortality" = "#D73027")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12),
      plot.title = element_text(face = "bold"),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) +
    labs(
      # title = "Day 1 Pathogen Abundance and Mass",
      x = NULL,
      y = "log2(mass + 1)",
      fill = "28-day Mortality",
      color = "28-day Mortality"
    )
}

pd1 = plot_d1_biomass(df_long = df_long, data = data, mapid_vec = mapid_vec)

ggsave(file.path(out_dir, "D1_Total_Viruses_EBV_HCMV_boxplot.pdf"),
       pd1, width = 5, height = 3.3)

