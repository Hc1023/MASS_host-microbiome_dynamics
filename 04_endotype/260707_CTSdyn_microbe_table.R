rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(purrr)
library(ggplot2)
library(stringr)
library(forcats)
library(grid)
library(survival)

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)

out_dir <- file.path("/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260707_CTSdyn_microbe_tab")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


load('Inputs/1616_meta_model.rdata')

df1 = meta_model %>% dplyr::select(HumanID, Mortality28d, D1_CTS, D4_CTS, D7_CTS)
df1 %<>% na.omit()

# CTS3 is the low-risk endotype; CTS1/2 are the high-risk endotypes.
# Therefore CTS3 -> CTS1/2 is deterioration and CTS1/2 -> CTS3 is improvement.
# Stable groups additionally require the same risk class at Day 4.
df1$grp = "Mixed"
idx = df1$D1_CTS == 3 & df1$D7_CTS != 3
df1$grp[idx] = 'Deteriorating'
idx = df1$D1_CTS != 3 & df1$D7_CTS == 3
df1$grp[idx] = 'Improving'
idx = df1$D1_CTS == 3 & df1$D4_CTS == 3 & df1$D7_CTS == 3
df1$grp[idx] = 'Stable low-risk'
idx = df1$D1_CTS != 3 & df1$D4_CTS != 3 & df1$D7_CTS != 3
df1$grp[idx] = 'Stable high-risk'



df_sum <- df1 %>%
  filter(!is.na(grp), !is.na(Mortality28d)) %>%
  group_by(grp) %>%
  summarise(
    n = n(),
    death_rate = mean(Mortality28d == 1),
    .groups = "drop"
  )

df_sum <- df_sum %>%
  mutate(
    death_pct = death_rate * 100,
    n_label = paste0("n = ", n)
  )
df_sum$grp = factor(df_sum$grp, 
                    levels = c("Stable high-risk", "Stable low-risk",
                               "Deteriorating", "Improving", "Mixed"))
p = ggplot(df_sum, aes(x = grp, y = death_pct)) +
  geom_col(
    fill = "#D55E00",
    color = 'black',
    width = 0.6,
    alpha = 0.85
  ) +
  
  geom_text(
    aes(label = n_label),
    # vjust = -0.5,
    hjust = -0.1,
    size = 3,
    color = "black"
  ) +
  
  scale_y_continuous(
    name = "Mortality (%)",
    limits = c(0, max(df_sum$death_pct) * 1.3),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 11),
    # plot.margin = margin(5.5, 5.5, 10, 5.5),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL) +
  coord_flip()

# print(p)


#### microbe changes ####

tmp = meta_model %>% dplyr::select(HumanID, matches('D[147]_mi'))
df1_mi = df1 %>% left_join(tmp, by = 'HumanID')

delta_long <- df1_mi %>%
  pivot_longer(
    cols = matches("^D(1|4|7)_(mi2?|mi)_"),
    names_to = c("Day", "feature"),
    names_pattern = "^D(1|4|7)_(.+)$",
    values_to = "value"
  ) %>%
  mutate(Day = paste0("D", Day)) %>%
  pivot_wider(names_from = Day, values_from = value) %>%
  mutate(delta = D7 - D1)

delta_long$feature <- gsub("^(mi2|mi)_", "", delta_long$feature)
delta_long %<>% filter(grp != 'Mixed')

#### log2FC ####


library(rstatix)

res_delta <- delta_long %>%
  filter(!is.na(D1), !is.na(D7)) %>%
  group_by(grp, feature) %>%
  summarise(
    n = n(),
    median_delta = median(delta, na.rm = TRUE),
    mean_delta   = mean(delta, na.rm = TRUE),
    se_delta     = sd(delta, na.rm = TRUE) / sqrt(n()),
    p_value = tryCatch(
      wilcox.test(D7, D1, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

mean_delta_mat <- res_delta %>%
  dplyr::select(feature, grp, mean_delta) %>%
  tidyr::pivot_wider(
    names_from  = grp,
    values_from = mean_delta
  )

pval_mat <- res_delta %>%
  dplyr::select(feature, grp, p_value) %>%
  tidyr::pivot_wider(
    names_from  = grp,
    values_from = p_value
  )

pval_mat_2 = pval_mat %>% filter(!is.na(Deteriorating)) %>% filter(Deteriorating<1)
pval_mat_2 %<>% arrange(Deteriorating)


bubble_df <- pval_mat_2 %>%
  pivot_longer(
    cols = -feature,
    names_to = "grp",
    values_to = "p"
  ) %>%
  left_join(
    mean_delta_mat %>%
      pivot_longer(
        cols = -feature,
        names_to = "grp",
        values_to = "delta"
      ),
    by = c("feature", "grp")
  ) %>%
  # filter(!is.na(p)) %>%
  mutate(
    neg_log10_p = -log10(p)
  )

bubble_df$feature = factor(bubble_df$feature, levels = rev(pval_mat_2$feature))

bubble_df <- bubble_df %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1  ~ ".",
      TRUE ~ ""
    )
  )

levels(bubble_df$feature)[levels(bubble_df$feature) == "HHV.4"] <- "EBV"
levels(bubble_df$feature)[levels(bubble_df$feature) == "Influenza.A"] <- "Influenza A"

levels(bubble_df$feature)[levels(bubble_df$feature) == "Bacf"] <- "Bacteria/Fungi"


#### Heatmap and source tables ####

# Keep the group order and labels explicit so that the clinical direction cannot
# be reversed accidentally in either the plot or the exported tables.
group_order <- c(
  "Deteriorating",
  "Improving",
  "Stable high-risk",
  "Stable low-risk"
)

group_labels <- c(
  "Deteriorating" = "CTS3 -> 1/2",
  "Improving" = "CTS1/2 -> 3",
  "Stable high-risk" = "Stable CTS1/2",
  "Stable low-risk" = "Stable CTS3"
)

feature_labels <- c(
  "HHV.4" = "EBV",
  "Influenza.A" = "Influenza A",
  "Bacf" = "Bacteria/Fungi"
)

heatmap_features <- pval_mat_2$feature

heatmap_stats <- res_delta %>%
  filter(feature %in% heatmap_features, grp %in% group_order) %>%
  mutate(
    feature_raw = feature,
    feature = recode(feature, !!!feature_labels),
    trajectory_group = grp,
    cts_transition = unname(group_labels[grp]),
    risk_direction = case_when(
      grp == "Deteriorating" ~ "Low-risk CTS3 to high-risk CTS1/2",
      grp == "Improving" ~ "High-risk CTS1/2 to low-risk CTS3",
      grp == "Stable high-risk" ~ "Persistently high-risk CTS1/2",
      grp == "Stable low-risk" ~ "Persistently low-risk CTS3"
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  arrange(
    match(feature_raw, heatmap_features),
    match(trajectory_group, group_order)
  ) %>%
  select(
    feature, feature_raw, trajectory_group, cts_transition, risk_direction,
    n, mean_delta, se_delta, median_delta, p_value, p_adj, significance
  )

heatmap_table <- heatmap_stats %>%
  mutate(
    feature = factor(feature, levels = unique(feature)),
    cts_transition = factor(cts_transition, levels = unname(group_labels)),
    value = paste0(sprintf("%.2f", mean_delta), significance)
  ) %>%
  select(feature, cts_transition, value) %>%
  pivot_wider(names_from = cts_transition, values_from = value) %>%
  arrange(feature)

transition_summary <- df_sum %>%
  filter(as.character(grp) %in% group_order) %>%
  mutate(
    trajectory_group = as.character(grp),
    cts_transition = unname(group_labels[trajectory_group]),
    risk_direction = case_when(
      trajectory_group == "Deteriorating" ~ "Low-risk CTS3 to high-risk CTS1/2",
      trajectory_group == "Improving" ~ "High-risk CTS1/2 to low-risk CTS3",
      trajectory_group == "Stable high-risk" ~ "Persistently high-risk CTS1/2",
      trajectory_group == "Stable low-risk" ~ "Persistently low-risk CTS3"
    )
  ) %>%
  arrange(match(trajectory_group, group_order)) %>%
  transmute(
    trajectory_group, cts_transition, risk_direction,
    sample_size = n, deaths = round(death_rate * n), mortality_percent = death_pct
  )

write_csv(
  heatmap_stats,
  file.path(out_dir, "CTS_transition_microbe_heatmap_statistics.csv")
)
write_csv(
  heatmap_table,
  file.path(out_dir, "CTS_transition_microbe_heatmap_display_table.csv")
)
write_csv(
  transition_summary,
  file.path(out_dir, "CTS_transition_group_summary.csv")
)

heatmap_plot <- ggplot(
  bubble_df %>% mutate(grp = factor(grp, levels = group_order)),
  aes(x = grp, y = feature)
) +
  geom_tile(
    aes(fill = delta),
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", delta), sig)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = expression(log[2]((D7 + 1) / (D1 + 1)))
  ) +
  scale_x_discrete(labels = group_labels, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "CTS_transition_microbe_heatmap.pdf"),
  heatmap_plot,
  width = 5,
  height = 3
)
ggsave(
  file.path(out_dir, "CTS_transition_microbe_heatmap.png"),
  heatmap_plot,
  width = 5,
  height = 3,
  dpi = 300
)


#### Styled Excel report ####

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("The 'openxlsx' package is required to create endotype_transition.xlsx")
}

# Section 1 uses the same features and order as the heatmap. Stars are based on
# the raw paired Wilcoxon P value; the BH-adjusted P value is also reported.
deteriorating_table <- heatmap_stats %>%
  filter(trajectory_group == "Deteriorating") %>%
  mutate(
    interpretation = case_when(
      p_value < 0.05 & mean_delta > 0 ~ "Significant increase",
      p_value < 0.05 & mean_delta < 0 ~ "Significant decrease",
      p_value < 0.1 & mean_delta > 0 ~ "Increasing trend",
      p_value < 0.1 & mean_delta < 0 ~ "Decreasing trend",
      TRUE ~ "No significant paired change"
    )
  ) %>%
  transmute(
    `Microbial feature` = feature,
    `Mean Delta log2` = mean_delta,
    `SE of Delta` = se_delta,
    `P value` = p_value,
    `BH-adjusted P` = p_adj,
    `Significance` = significance,
    `Statistical interpretation` = interpretation
  )

trajectory_deltas <- res_delta %>%
  filter(
    grp %in% group_order,
    feature %in% c("Total", "Bacf", "Viruses")
  ) %>%
  mutate(feature = recode(feature, "Bacf" = "Bacteria/Fungi")) %>%
  select(grp, feature, mean_delta) %>%
  pivot_wider(names_from = feature, values_from = mean_delta)

trajectory_table <- transition_summary %>%
  left_join(trajectory_deltas, by = c("trajectory_group" = "grp")) %>%
  transmute(
    `Dynamic host-state trajectory` = paste0(cts_transition, " (", trajectory_group, ")"),
    `Sample size (n)` = sample_size,
    `Deaths (n)` = deaths,
    `28-day mortality` = mortality_percent / 100,
    `Total Delta log2` = Total,
    `Bacteria/Fungi Delta log2` = `Bacteria/Fungi`,
    `Viral Delta log2` = Viruses
  )

wb <- openxlsx::createWorkbook(creator = "MASS sepsis microbiome analysis")
openxlsx::addWorksheet(wb, "Endotype Transitions", gridLines = FALSE)
openxlsx::addWorksheet(wb, "Full Statistics", gridLines = FALSE)
openxlsx::addWorksheet(wb, "Heatmap", gridLines = FALSE)

report_sheet <- "Endotype Transitions"
dark_blue <- "#1F4E78"
medium_blue <- "#4472C4"
header_blue <- "#2F5597"
light_blue <- "#D9EAF7"
light_grey <- "#F2F2F2"
border_grey <- "#D9D9D9"

title_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 20, fontColour = dark_blue,
  textDecoration = "bold", halign = "left", valign = "center"
)
subtitle_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 11, fontColour = "#595959",
  textDecoration = "italic", halign = "left", valign = "center"
)
section_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 12, fontColour = "#FFFFFF",
  textDecoration = "bold", fgFill = medium_blue,
  halign = "left", valign = "center"
)
header_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "#FFFFFF",
  textDecoration = "bold", fgFill = header_blue,
  halign = "center", valign = "center", wrapText = TRUE,
  border = "Bottom", borderStyle = "medium", borderColour = "#FFFFFF"
)
body_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 10, valign = "center",
  border = "Bottom", borderStyle = "thin", borderColour = border_grey
)
band_style <- openxlsx::createStyle(fgFill = light_grey)
numeric_style <- openxlsx::createStyle(numFmt = "0.00", halign = "right")
pvalue_style <- openxlsx::createStyle(numFmt = "0.0000", halign = "right")
percent_style <- openxlsx::createStyle(numFmt = "0.0%", halign = "right")
integer_style <- openxlsx::createStyle(numFmt = "0", halign = "right")
definition_label_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 10, textDecoration = "bold",
  fgFill = light_blue, valign = "top"
)
definition_text_style <- openxlsx::createStyle(
  fontName = "Arial", fontSize = 10, wrapText = TRUE, valign = "top"
)

openxlsx::mergeCells(wb, report_sheet, cols = 1:7, rows = 1)
openxlsx::writeData(wb, report_sheet,
  "Host Endotype Transitions & Microbial Expansion", startRow = 1, startCol = 1
)
openxlsx::addStyle(wb, report_sheet, title_style, rows = 1, cols = 1:7, gridExpand = TRUE)
openxlsx::setRowHeights(wb, report_sheet, rows = 1, heights = 30)

openxlsx::mergeCells(wb, report_sheet, cols = 1:7, rows = 2)
openxlsx::writeData(wb, report_sheet,
  "Paired longitudinal change in spike-in-normalized absolute microbial burden (Day 7 - Day 1)",
  startRow = 2, startCol = 1
)
openxlsx::addStyle(wb, report_sheet, subtitle_style, rows = 2, cols = 1:7, gridExpand = TRUE)
openxlsx::setRowHeights(wb, report_sheet, rows = 2, heights = 22)

openxlsx::mergeCells(wb, report_sheet, cols = 1:7, rows = 4)
openxlsx::writeData(wb, report_sheet,
  "1. MICROBIAL EXPANSION IN THE CTS3 -> CTS1/2 DETERIORATING GROUP",
  startRow = 4, startCol = 1
)
openxlsx::addStyle(wb, report_sheet, section_style, rows = 4, cols = 1:7, gridExpand = TRUE)
openxlsx::setRowHeights(wb, report_sheet, rows = 4, heights = 24)

openxlsx::writeData(
  wb, report_sheet, deteriorating_table,
  startRow = 5, startCol = 1, headerStyle = header_style, withFilter = FALSE
)
section1_rows <- 6:(5 + nrow(deteriorating_table))
openxlsx::addStyle(wb, report_sheet, body_style, rows = section1_rows, cols = 1:7, gridExpand = TRUE)
openxlsx::addStyle(wb, report_sheet, numeric_style, rows = section1_rows, cols = 2:3, gridExpand = TRUE, stack = TRUE)
openxlsx::addStyle(wb, report_sheet, pvalue_style, rows = section1_rows, cols = 4:5, gridExpand = TRUE, stack = TRUE)
for (r in section1_rows[seq_along(section1_rows) %% 2 == 0]) {
  openxlsx::addStyle(wb, report_sheet, band_style, rows = r, cols = 1:7, gridExpand = TRUE, stack = TRUE)
}

section2_row <- 5 + nrow(deteriorating_table) + 2
openxlsx::mergeCells(wb, report_sheet, cols = 1:7, rows = section2_row)
openxlsx::writeData(wb, report_sheet,
  "2. MICROBIAL TRAJECTORIES AND 28-DAY OUTCOMES BY HOST ENDOTYPE TRANSITION",
  startRow = section2_row, startCol = 1
)
openxlsx::addStyle(wb, report_sheet, section_style, rows = section2_row, cols = 1:7, gridExpand = TRUE)
openxlsx::setRowHeights(wb, report_sheet, rows = section2_row, heights = 24)

openxlsx::writeData(
  wb, report_sheet, trajectory_table,
  startRow = section2_row + 1, startCol = 1, headerStyle = header_style, withFilter = FALSE
)
section2_rows <- (section2_row + 2):(section2_row + 1 + nrow(trajectory_table))
openxlsx::addStyle(wb, report_sheet, body_style, rows = section2_rows, cols = 1:7, gridExpand = TRUE)
openxlsx::addStyle(wb, report_sheet, integer_style, rows = section2_rows, cols = 2:3, gridExpand = TRUE, stack = TRUE)
openxlsx::addStyle(wb, report_sheet, percent_style, rows = section2_rows, cols = 4, gridExpand = TRUE, stack = TRUE)
openxlsx::addStyle(wb, report_sheet, numeric_style, rows = section2_rows, cols = 5:7, gridExpand = TRUE, stack = TRUE)
for (r in section2_rows[seq_along(section2_rows) %% 2 == 0]) {
  openxlsx::addStyle(wb, report_sheet, band_style, rows = r, cols = 1:7, gridExpand = TRUE, stack = TRUE)
}

definitions_row <- section2_row + nrow(trajectory_table) + 4
openxlsx::mergeCells(wb, report_sheet, cols = 1:7, rows = definitions_row)
openxlsx::writeData(wb, report_sheet, "DEFINITIONS & CONVENTIONS", startRow = definitions_row, startCol = 1)
openxlsx::addStyle(wb, report_sheet, section_style, rows = definitions_row, cols = 1:7, gridExpand = TRUE)

definitions <- tibble::tribble(
  ~Term, ~Definition,
  "CTS3", "Low-risk Clinical Transcriptomic State.",
  "CTS1 & CTS2", "High-risk Clinical Transcriptomic States.",
  "Deteriorating", "CTS3 at Day 1 changing to CTS1/2 at Day 7.",
  "Improving", "CTS1/2 at Day 1 changing to CTS3 at Day 7.",
  "Delta (Day 7 - Day 1)", "Change in log2(mass + 1); positive values indicate increasing microbial burden.",
  "Statistical test", "Paired two-sided Wilcoxon signed-rank test. Stars use raw P values: *** <0.001, ** <0.01, * <0.05, . <0.10."
)
definition_data_rows <- (definitions_row + 1):(definitions_row + nrow(definitions))
openxlsx::writeData(wb, report_sheet, definitions, startRow = definitions_row + 1, startCol = 1, colNames = FALSE)
for (r in definition_data_rows) {
  openxlsx::mergeCells(wb, report_sheet, cols = 2:7, rows = r)
}
openxlsx::addStyle(wb, report_sheet, definition_label_style, rows = definition_data_rows, cols = 1, gridExpand = TRUE)
openxlsx::addStyle(wb, report_sheet, definition_text_style, rows = definition_data_rows, cols = 2:7, gridExpand = TRUE)
openxlsx::setRowHeights(wb, report_sheet, rows = definition_data_rows, heights = 26)

openxlsx::setColWidths(wb, report_sheet, cols = 1, widths = 34)
openxlsx::setColWidths(wb, report_sheet, cols = 2:6, widths = c(15, 13, 13, 16, 12))
openxlsx::setColWidths(wb, report_sheet, cols = 7, widths = 32)
openxlsx::setRowHeights(wb, report_sheet, rows = 5, heights = 36)
openxlsx::setRowHeights(wb, report_sheet, rows = section2_row + 1, heights = 42)
openxlsx::freezePane(wb, report_sheet, firstActiveRow = 5)
openxlsx::pageSetup(
  wb, report_sheet, orientation = "landscape", fitToWidth = 1, fitToHeight = 1,
  paperSize = 9
)

# Machine-readable numeric source sheet.
openxlsx::writeDataTable(
  wb, "Full Statistics", heatmap_stats,
  startRow = 1, startCol = 1, tableStyle = "TableStyleMedium2",
  tableName = "HeatmapStatistics", withFilter = TRUE
)
openxlsx::freezePane(wb, "Full Statistics", firstRow = TRUE, firstCol = TRUE)
openxlsx::setColWidths(wb, "Full Statistics", cols = 1:ncol(heatmap_stats), widths = "auto")
openxlsx::setColWidths(wb, "Full Statistics", cols = 5, widths = 38)
openxlsx::addStyle(
  wb, "Full Statistics", numeric_style,
  rows = 2:(nrow(heatmap_stats) + 1), cols = c(7:9), gridExpand = TRUE, stack = TRUE
)
openxlsx::addStyle(
  wb, "Full Statistics", pvalue_style,
  rows = 2:(nrow(heatmap_stats) + 1), cols = c(10:11), gridExpand = TRUE, stack = TRUE
)

# Keep the publication heatmap in the workbook as a separate, uncluttered sheet.
openxlsx::insertImage(
  wb, "Heatmap", file.path(out_dir, "CTS_transition_microbe_heatmap.png"),
  startRow = 1, startCol = 1, width = 10, height = 6, units = "in"
)
openxlsx::setColWidths(wb, "Heatmap", cols = 1:12, widths = 12)
openxlsx::setRowHeights(wb, "Heatmap", rows = 1:30, heights = 18)

openxlsx::saveWorkbook(
  wb,
  file.path(out_dir, "endotype_transition.xlsx"),
  overwrite = TRUE
)



#### 折线图  ####

feature_x = 'Klebsiella'
feature_x = 'Total'

delta_long$feature[delta_long$feature == 'HHV.4'] = 'EBV'
delta_long$feature = gsub('\\.',' ',delta_long$feature)
delta_long$feature[delta_long$feature == 'SARS CoV 2'] = 'SARS-CoV-2'
delta_long$feature[delta_long$feature == 'HAdV B'] = 'HAdV-B'
delta_long$feature[delta_long$feature == 'HHV 6B'] = 'HHV-6B'
delta_long$feature[delta_long$feature == 'Bacf'] = 'Bacteria/Fungi'

draw_feature_x = function(feature_x){
  group_cols <- c(
    "Deteriorating" = "#8DA0CB",
    "Improving"     = "#D95F02"
  )

  dfx = delta_long %>%
    filter(feature == feature_x, grp %in% names(group_cols)) %>%
    mutate(grp = factor(grp, levels = names(group_cols)))
  
  df_long <- dfx %>%
    pivot_longer(
      cols = c(D1, D4, D7),
      names_to = "Day",
      values_to = "value"
    ) %>%
    mutate(
      Day = factor(Day, levels = c("D1", "D4", "D7")),
      Day_num = as.numeric(Day)
    )

  stat_line <- df_long %>%
    group_by(Day, grp) %>%
    summarise(
      m = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    ) %>%
    mutate(
      lower = m - se,
      upper = m + se
    )

  y_range <- range(c(stat_line$lower, stat_line$upper), na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  
  stat_grp <- dfx %>%
    filter(!is.na(D1), !is.na(D7)) %>%
    group_by(grp) %>%
    summarise(
      n = n(),
      log2FC = mean(D7 - D1, na.rm = TRUE),
      p = tryCatch(
        wilcox.test(D7, D1, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      sig = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ ".",
        TRUE ~ "ns"
      ),
      label = paste0("log2FC = ", sprintf("%.2f", log2FC), " ", sig),
      x.position = 2,
      y.position = case_when(
        grp == "Deteriorating" ~ y_range[2] - y_span * 0.15,
        grp == "Improving" ~ y_range[1] + y_span * 0.05,
        TRUE ~ NA_real_
      )
    )
  
  p = ggplot(df_long,
             aes(x = Day_num, y = value,
                 color = grp, group = grp)) +
    
    ## mean ± SE
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15,
      linewidth = 0.6
    ) +
    
    ## 组均值线（加 linetype）
    stat_summary(
      fun = mean,
      geom = "line",
      linewidth = 0.8
    ) +
    
    ## 组均值点
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2.5
    ) +
    
    theme_bw() +
    labs(
      x = NULL,
      y = "log2(mass+1)",
      color = "Trajectory group"
    ) +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'none') +
    ggtitle(feature_x) +
    geom_text(
      data = stat_grp,
      aes(x = x.position, y = y.position, label = label, color = grp),
      inherit.aes = FALSE,
      size = 3.2,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = group_cols
    ) +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("D1", "D4", "D7"),
      limits = c(0.9, 3.1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
  p
  return(p)
}

features_to_plot <- c("Total", "Bacteria/Fungi", "Viruses")
plots <- setNames(map(features_to_plot, draw_feature_x), features_to_plot)
draw_feature_x("Total")
draw_feature_x("Bacteria/Fungi")
draw_feature_x("Viruses")

group_cols <- c(
  "Deteriorating" = "#8DA0CB",
  "Improving"     = "#D95F02"
)

legend_df <- tibble(
  Day = rep(c(1, 2), 2),
  value = rep(c(1, 1), 2),
  grp = factor(
    rep(names(group_cols), each = 2),
    levels = names(group_cols)
  )
)

legend_plot <- ggplot(
  legend_df,
  aes(x = Day, y = value, color = grp, group = grp)
) +
  geom_line(linewidth = 0.8, show.legend = TRUE) +
  geom_point(size = 2.5, show.legend = TRUE) +
  scale_color_manual(
    values = group_cols,
    labels = c(
      "Deteriorating (CTS3 -> 1/2)",
      "Improving (CTS1/2 -> 3)"
    )
  ) +
  labs(color = NULL) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.1, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 0.8, size = 2.5)
    )
  )

legend_grob <- ggplotGrob(legend_plot)
legend_grob <- legend_grob$grobs[[which(
  vapply(legend_grob$grobs, function(x) x$name, character(1)) == "guide-box"
)]]

ggsave(
  filename = file.path(out_dir, "trajectory_two_groups_legend.pdf"),
  plot = legend_grob,
  width = 4.8,
  height = 0.35
)

walk2(
  plots,
  names(plots),
  ~ ggsave(
    filename = file.path(
      out_dir,
      paste0("trajectory_two_groups_", str_replace_all(.y, "[/ ]+", "_"), ".pdf")
    ),
    plot = .x,
    width = 2.3,
    height = 2.5
  )
)
