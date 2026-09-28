# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

project_dir <- "."

analysis_dir <- repo_path(project_dir, "Outputs")

plot_dir <- repo_path(analysis_dir, "plot")

repo_dir(plot_dir, recursive = TRUE, showWarnings = FALSE)

performance_file <- "Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_performance.csv"

contrasts_file <- "Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_contrasts.csv"

repeat_file <- "Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_repeat_performance.csv"

distribution_file <- "Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_contrast_distributions.csv"

stopifnot(repo_exists(performance_file), repo_exists(contrasts_file), repo_exists(repeat_file), repo_exists(distribution_file))

performance <- read.csv(repo_input(performance_file), check.names = FALSE)

contrasts <- read.csv(repo_input(contrasts_file), check.names = FALSE)

repeat_performance <- read.csv(repo_input(repeat_file), check.names = FALSE)

contrast_distributions <- read.csv(repo_input(distribution_file), check.names = FALSE)

all_model_levels <- paste0("M", 0:5)

main_model_levels <- paste0("M", 0:3)

contrast_levels <- c("M1 - M0", "M2 - M1", "M3 - M2", "M3 - M0", "M4 - M0", "M5 - M4", "M5 - M2", "M5 - M3")

incremental_plot_levels <- c("M1 - M0", "M2 - M1", "M3 - M2", "M3 - M0")

stopifnot(identical(performance$Model, all_model_levels), identical(contrasts$Contrast, contrast_levels), 
    setequal(unique(repeat_performance$Model), all_model_levels), setequal(unique(contrast_distributions$Contrast), 
        contrast_levels), length(unique(repeat_performance$Repeat)) == performance$Repeats[1])

performance_qc <- repeat_performance %>% group_by(Model) %>% summarise(Mean_CV_AUC = mean(AUC), AUC_Q25 = quantile(AUC, 
    0.25), AUC_Q75 = quantile(AUC, 0.75), Mean_CV_Brier = mean(Brier), Brier_Q25 = quantile(Brier, 0.25), 
    Brier_Q75 = quantile(Brier, 0.75), .groups = "drop")

contrast_qc <- contrast_distributions %>% group_by(Contrast) %>% summarise(Mean_paired_DeltaAUC = mean(DeltaAUC), 
    Mean_paired_DeltaBrier = mean(DeltaBrier), .groups = "drop")

merged_performance_qc <- performance %>% select(Model, Mean_CV_AUC, AUC_Q25, AUC_Q75, Mean_CV_Brier, 
    Brier_Q25, Brier_Q75) %>% inner_join(performance_qc, by = "Model", suffix = c("_saved", "_recalculated"))

merged_contrast_qc <- contrasts %>% select(Contrast, Mean_paired_DeltaAUC, Mean_paired_DeltaBrier) %>% 
    inner_join(contrast_qc, by = "Contrast", suffix = c("_saved", "_recalculated"))

if (max(abs(as.matrix(merged_performance_qc[, grep("_saved$", names(merged_performance_qc))]) - as.matrix(merged_performance_qc[, 
    grep("_recalculated$", names(merged_performance_qc))]))) > 1e-12 || max(abs(as.matrix(merged_contrast_qc[, 
    grep("_saved$", names(merged_contrast_qc))]) - as.matrix(merged_contrast_qc[, grep("_recalculated$", 
    names(merged_contrast_qc))]))) > 1e-12) {
    stop("Saved detailed outputs do not match saved summaries")
}

viral_model_levels <- c("M0", "M4", "M2", "M5")

viral_model_labels <- c(M0 = "SOFA", M4 = "SOFA +\nEBV/HCMV", M2 = "SOFA +\nCTS genes", M5 = "SOFA +\nCTS genes +\nEBV/HCMV")

viral_model_colors <- c(M0 = "#7F7F7F", M4 = "#D95F02", M2 = "#4575B4", M5 = "#B2182B")

viral_repeat_data <- repeat_performance %>% filter(Model %in% viral_model_levels) %>% mutate(Model = factor(Model, 
    levels = viral_model_levels), Label = factor(viral_model_labels[as.character(Model)], levels = viral_model_labels))

viral_data <- performance %>% filter(Model %in% viral_model_levels) %>% left_join(repeat_performance %>% 
    filter(Model %in% viral_model_levels) %>% group_by(Model) %>% summarise(Max_repeat_AUC = max(AUC), 
    .groups = "drop"), by = "Model") %>% mutate(Model = factor(Model, levels = viral_model_levels), Label = factor(viral_model_labels[as.character(Model)], 
    levels = viral_model_labels), Mean_label = sprintf("%.3f", Mean_CV_AUC))

viral_min <- min(viral_repeat_data$AUC)

viral_max <- max(viral_repeat_data$AUC)

viral_span <- viral_max - viral_min

if (!is.finite(viral_span) || viral_span <= 0) viral_span <- 0.1

viral_label_offset <- 0.06 * viral_span

viral_bracket_step <- 0.13 * viral_span

viral_bracket_tick <- 0.028 * viral_span

viral_first_y <- viral_max + 0.18 * viral_span

viral_second_y <- viral_first_y + viral_bracket_step

viral_data$Mean_label_y <- viral_data$Max_repeat_AUC + viral_label_offset

viral_lower <- viral_min - 0.08 * viral_span

viral_upper <- viral_second_y + 0.17 * viral_span

viral_delta <- contrasts$Mean_paired_DeltaAUC[contrasts$Contrast == "M4 - M0"]

viral_after_cts_delta <- contrasts$Mean_paired_DeltaAUC[contrasts$Contrast == "M5 - M2"]

viral_delta_label <- sprintf("AUC %+.3f", viral_delta)

viral_after_cts_delta_label <- sprintf("AUC %+.3f", viral_after_cts_delta)

viral_annotation_gap <- 0.09 * viral_span

p_viral_cts <- ggplot(viral_repeat_data, aes(x = Label, y = AUC, color = Model)) + geom_boxplot(aes(fill = Model), 
    width = 0.32, outlier.shape = NA, alpha = 0.22, linewidth = 0.65, show.legend = FALSE) + geom_jitter(size = 1.7, 
    alpha = 0.55, show.legend = FALSE, position = position_jitter(width = 0.085, height = 0, seed = 260831)) + 
    geom_point(data = viral_data, aes(x = Label, y = Mean_CV_AUC, fill = Model), inherit.aes = FALSE, 
        shape = 23, size = 4.3, stroke = 0.8, color = "black", show.legend = FALSE) + geom_text(data = viral_data, 
    aes(x = Label, y = Mean_label_y, label = Mean_label), inherit.aes = FALSE, size = 3.5, fontface = "bold") + 
    annotate("segment", x = 1, xend = 2, y = viral_first_y, yend = viral_first_y, linewidth = 0.6, color = "#A64200") + 
    annotate("segment", x = c(1, 2), xend = c(1, 2), y = viral_first_y - viral_bracket_tick, yend = viral_first_y, 
        linewidth = 0.6, color = "#A64200") + annotate("text", x = 1.5, y = viral_first_y + viral_bracket_tick, 
    label = viral_delta_label, vjust = 0, size = 3.1, color = "#A64200") + annotate("text", x = 1.5, 
    y = viral_first_y + viral_bracket_tick + viral_annotation_gap, label = "D4 EBV/HCMV", vjust = -0.3, 
    size = 3.1, color = "#A64200") + annotate("segment", x = 3, xend = 4, y = viral_second_y, yend = viral_second_y, 
    linewidth = 0.6, color = "#8E1421") + annotate("segment", x = c(3, 4), xend = c(3, 4), y = viral_second_y - 
    viral_bracket_tick, yend = viral_second_y, linewidth = 0.6, color = "#8E1421") + annotate("text", 
    x = 3.5, y = viral_second_y + viral_bracket_tick, label = viral_after_cts_delta_label, vjust = 0, 
    size = 3.1, color = "#8E1421") + annotate("text", x = 3.5, y = viral_second_y + viral_bracket_tick + 
    viral_annotation_gap, label = "D4 EBV/HCMV beyond CTS", vjust = -0.3, size = 3.1, color = "#8E1421") + 
    scale_color_manual(values = viral_model_colors) + scale_fill_manual(values = viral_model_colors) + 
    coord_cartesian(ylim = c(viral_lower, viral_upper), clip = "off") + labs(title = "Day-4 viral and CTS reassessment", 
    subtitle = sprintf("%d repeats x stratified %d-fold outer CV; n = %d (%d deaths, %d survivors)", 
        performance$Repeats[1], performance$Outer_folds_per_repeat[1], performance$N[1], performance$Deaths[1], 
        performance$Survivors[1]), x = NULL, y = "Outer-CV AUC", caption = "Small points are repeat-level AUCs; boxes show median and IQR; diamonds show means.") + 
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), 
    plot.subtitle = element_text(hjust = 0.5, size = 10.5), axis.text.x = element_text(color = "black", 
        size = 8.5, lineheight = 0.92, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black"), 
    plot.caption = element_text(size = 8.5), plot.margin = margin(10, 12, 8, 10))

ggsave(repo_output(repo_path(plot_dir, "260831_Step2_D4_D1_FigS_viral_CTS_AUC.pdf"), "FigureS11A_plot.R"), 
    p_viral_cts, width = 4.5, height = 4.2, device = "pdf")

message("Saved Figure 6B and both supplementary figures to: ", plot_dir)

