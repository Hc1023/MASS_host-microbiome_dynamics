# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

project_dir <- "."

analysis_dir <- repo_path(project_dir, "Outputs")

plot_dir <- repo_path(analysis_dir, "plot")

repo_dir(plot_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "260901_Step4_2_CMAISE_subgroups"

performance <- read.csv("Outputs/FigureS11C_analysis_Step4_2_CMAISE_subgroups_performance.csv", check.names = FALSE)

comparisons <- read.csv("Outputs/FigureS11C_analysis_Step4_2_CMAISE_subgroups_comparisons.csv", check.names = FALSE)

stopifnot(nrow(performance) == 16L, nrow(comparisons) == 24L)

model_colors <- c(M0 = "#777777", M1 = "#4575B4", M2 = "#7B3294", M3 = "#D95F02")

main_labels <- c(M0 = "SOFA", M1 = "SOFA +\nCTS genes", M2 = "SOFA +\ntrajectory genes", M3 = "SOFA +\ncombined")

main <- performance %>% filter(Model %in% names(main_labels)) %>% mutate(Timepoint = factor(Timepoint, 
    c("D3", "D5")), Subgroup = factor(Subgroup, c("Lung infection", "Non-lung infection")), Model = factor(Model, 
    names(main_labels)), Label = factor(main_labels[as.character(Model)], levels = main_labels), AUC_label = sprintf("%.3f", 
    AUC))

panel_tops <- main %>% group_by(Timepoint, Subgroup) %>% summarise(Panel_top = max(AUC_CI_upper + 0.025), 
    .groups = "drop")

brackets <- comparisons %>% filter(Contrast %in% c("M2 - M1", "M3 - M2")) %>% mutate(Timepoint = factor(Timepoint, 
    c("D3", "D5")), Subgroup = factor(Subgroup, c("Lung infection", "Non-lung infection")), x = ifelse(Contrast == 
    "M2 - M1", 2, 3), xend = x + 1, label = sprintf("AUC %+.3f", DeltaAUC), color = ifelse(Contrast == 
    "M2 - M1", "#285A84", "#A64B00")) %>% left_join(panel_tops, by = c("Timepoint", "Subgroup")) %>% 
    mutate(y = Panel_top + ifelse(Contrast == "M2 - M1", 0.05, 0.11))

plot_upper <- max(brackets$y) + 0.055

p_main <- ggplot(main, aes(Label, AUC, color = Model)) + geom_hline(yintercept = 0.5, linetype = "dashed", 
    color = "grey70") + geom_errorbar(aes(ymin = AUC_CI_lower, ymax = AUC_CI_upper), width = 0.08, linewidth = 0.65) + 
    geom_point(aes(fill = Model), shape = 23, size = 3.8, stroke = 0.75, color = "black") + geom_text(aes(y = AUC_CI_upper + 
    0.03, label = AUC_label), color = "black", size = 3, fontface = "bold") + geom_segment(data = brackets, 
    aes(x = x, xend = xend, y = y, yend = y, color = color), inherit.aes = FALSE, linewidth = 0.6) + 
    geom_segment(data = brackets, aes(x = x, xend = x, y = y - 0.012, yend = y, color = color), inherit.aes = FALSE, 
        linewidth = 0.6) + geom_segment(data = brackets, aes(x = xend, xend = xend, y = y - 0.012, yend = y, 
    color = color), inherit.aes = FALSE, linewidth = 0.6) + geom_text(data = brackets, aes(x = (x + xend)/2, 
    y = y + 0.009, label = label, color = color), inherit.aes = FALSE, size = 2.8, vjust = 0) + facet_grid(Subgroup ~ 
    Timepoint, scales = "free_y", labeller = labeller(Timepoint = as_labeller(c(D3 = "CMAISE Day 3", 
    D5 = "CMAISE Day 5")), Subgroup = label_value)) + scale_color_manual(values = c(model_colors, `#285A84` = "#285A84", 
    `#A64B00` = "#A64B00"), guide = "none") + scale_fill_manual(values = model_colors, guide = "none") + 
    coord_cartesian(ylim = c(0.45, NA), clip = "off") + scale_y_continuous(breaks = c(0.5, 0.75, 1)) + 
    labs(x = NULL, y = "Subgroup AUC", title = "CMAISE validation by infection source", subtitle = "Frozen MASS models; error bars show 95% CI") + 
    theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), 
    plot.subtitle = element_text(hjust = 0.5), strip.text = element_text(face = "bold"), axis.text.x = element_text(color = "black", 
        size = 8, lineheight = 0.88, angle = 45, hjust = 1, vjust = 1))

ggsave(repo_output(repo_path(plot_dir, paste0(prefix, "_FigS_subgroup_validation.pdf")), "FigureS11C_plot.R"), 
    p_main, width = 6, height = 5.5, device = "pdf")

