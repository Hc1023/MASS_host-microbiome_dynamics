# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

project_dir <- "."

analysis_dir <- repo_path(project_dir, "Outputs")

plot_dir <- repo_path(analysis_dir, "plot")

repo_dir(plot_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "260901_Step4_CMAISE_validation"

performance <- read.csv("Outputs/Figure6D_analysis_Step4_CMAISE_validation_performance.csv", check.names = FALSE)

comparisons <- read.csv("Outputs/Figure6D_analysis_Step4_CMAISE_validation_comparisons.csv", check.names = FALSE)

stopifnot(nrow(performance) == 8L, nrow(comparisons) == 12L)

model_levels <- c("M0", "M1", "M2", "M3", "M4")

main_model_levels <- c("M0", "M1", "M2", "M3")

model_labels <- c(M0 = "SOFA", M1 = "SOFA +\nCTS genes", M2 = "SOFA +\ntrajectory genes", M3 = "SOFA +\ncombined")

model_colors <- c(M0 = "#777777", M1 = "#4575B4", M2 = "#7B3294", M3 = "#D95F02")

plot_data <- performance %>% filter(Model %in% main_model_levels) %>% mutate(Timepoint = factor(Timepoint, 
    c("D3", "D5")), Model = factor(Model, main_model_levels), Label = factor(model_labels[as.character(Model)], 
    levels = model_labels), AUC_label = sprintf("%.3f", AUC))

brackets <- comparisons %>% filter(Contrast %in% c("M2 - M1", "M3 - M2")) %>% mutate(Timepoint = factor(Timepoint, 
    c("D3", "D5")), x = ifelse(Contrast == "M2 - M1", 2, 3), xend = x + 1, y = ifelse(Contrast == "M2 - M1", 
    0.935, 0.975) + ifelse(Timepoint == "D3", -0.02, 0.01), label = sprintf("AUC %+.3f", DeltaAUC), color = ifelse(Contrast == 
    "M2 - M1", "#285A84", "#A64B00"))

p_main <- ggplot(plot_data, aes(Label, AUC, color = Model)) + geom_hline(yintercept = 0.5, linetype = "dashed", 
    color = "grey70") + geom_errorbar(aes(ymin = AUC_CI_lower, ymax = AUC_CI_upper), width = 0.09, linewidth = 0.7) + 
    geom_point(aes(fill = Model), shape = 23, size = 4.2, stroke = 0.8, color = "black") + geom_text(aes(y = AUC_CI_upper + 
    0.025, label = AUC_label), color = "black", size = 3.5, fontface = "bold") + geom_segment(data = brackets, 
    aes(x = x, xend = xend, y = y, yend = y, color = color), inherit.aes = FALSE, linewidth = 0.6) + 
    geom_segment(data = brackets, aes(x = x, xend = x, y = y - 0.012, yend = y, color = color), inherit.aes = FALSE, 
        linewidth = 0.6) + geom_segment(data = brackets, aes(x = xend, xend = xend, y = y - 0.012, yend = y, 
    color = color), inherit.aes = FALSE, linewidth = 0.6) + geom_text(data = brackets, aes(x = (x + xend)/2, 
    y = y + 0.009, label = label, color = color), inherit.aes = FALSE, size = 3, vjust = 0) + facet_wrap(~Timepoint, 
    ncol = 1, labeller = as_labeller(c(D3 = "CMAISE Day 3", D5 = "CMAISE Day 5"))) + scale_color_manual(values = c(model_colors, 
    `#285A84` = "#285A84", `#A64B00` = "#A64B00"), guide = "none") + scale_fill_manual(values = model_colors, 
    guide = "none") + coord_cartesian(ylim = c(0.48, 1.01), clip = "off") + labs(x = NULL, y = "External-validation AUC", 
    title = "External validation of\nDay-4 host models", subtitle = "Direct application in CMAISE;\npoints and error bars show AUC and 95% CI") + 
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), 
    plot.subtitle = element_text(hjust = 0.5, size = 10), strip.text = element_text(face = "bold"), axis.text.x = element_text(color = "black", 
        size = 8.8, lineheight = 0.88, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black"), 
    plot.margin = margin(8, 12, 8, 10))

ggsave(repo_output(repo_path(plot_dir, paste0(prefix, "_Figure6D.pdf")), "Figure6D_plot.R"), 
    p_main, width = 3.7, height = 6.5, device = "pdf")

