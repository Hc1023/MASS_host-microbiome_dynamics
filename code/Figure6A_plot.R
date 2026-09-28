# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(pROC)
    library(tidyr)
})

project_dir <- "."

analysis_dir <- repo_path(project_dir, "Outputs")

plot_dir <- repo_path(analysis_dir, "plot")

repo_dir(plot_dir, recursive = TRUE, showWarnings = FALSE)

performance_file <- repo_path(analysis_dir, "Inputs/260829_Step1_D1_performance.csv")

contrasts_file <- repo_path(analysis_dir, "Inputs/260829_Step1_D1_contrasts.csv")

oof_file <- repo_path(analysis_dir, "Inputs/260829_Step1_D1_oof_predictions.csv")

stopifnot(repo_exists(performance_file), repo_exists(contrasts_file), repo_exists(oof_file))

performance <- read.csv(repo_input(performance_file), check.names = FALSE)

contrasts <- read.csv(repo_input(contrasts_file), check.names = FALSE)

oof <- read.csv(repo_input(oof_file), check.names = FALSE)

model_levels <- c("M0", "M1", "M2", "M3")

contrast_levels <- c("M1 - M0", "M2 - M0", "M3 - M1", "M3 - M2")

required_probability_cols <- paste0(model_levels, "_probability")

stopifnot(identical(performance$Model, model_levels), identical(contrasts$Contrast, contrast_levels), 
    all(c("HumanID", "Mortality28d", "Repeat", "Outer_fold", required_probability_cols) %in% names(oof)), 
    length(unique(oof$Repeat)) == 20L, all(count(oof, Repeat, HumanID)$n == 1L))

calc_auc <- function(y, probability) {
    as.numeric(pROC::auc(pROC::roc(response = y, predictor = probability, levels = c(0, 1), direction = "<", 
        quiet = TRUE)))
}

repeat_performance <- lapply(sort(unique(oof$Repeat)), function(repeat_id) {
    repeat_data <- oof %>% filter(Repeat == repeat_id)
    bind_rows(lapply(model_levels, function(model) {
        probability <- repeat_data[[paste0(model, "_probability")]]
        tibble(Repeat = repeat_id, Model = model, AUC = calc_auc(repeat_data$Mortality28d, probability), 
            Brier = mean((repeat_data$Mortality28d - probability)^2))
    }))
}) %>% bind_rows()

contrast_definitions <- tibble(Contrast = contrast_levels, New_model = c("M1", "M2", "M3", "M3"), Reference_model = c("M0", 
    "M0", "M1", "M2"))

contrast_distributions <- lapply(seq_len(nrow(contrast_definitions)), function(i) {
    definition <- contrast_definitions[i, ]
    wide <- repeat_performance %>% filter(Model %in% c(definition$New_model, definition$Reference_model)) %>% 
        pivot_wider(names_from = Model, values_from = c(AUC, Brier))
    tibble(Repeat = wide$Repeat, Contrast = definition$Contrast, DeltaAUC = wide[[paste0("AUC_", definition$New_model)]] - 
        wide[[paste0("AUC_", definition$Reference_model)]], DeltaBrier = wide[[paste0("Brier_", definition$New_model)]] - 
        wide[[paste0("Brier_", definition$Reference_model)]])
}) %>% bind_rows() %>% mutate(Contrast = factor(Contrast, levels = contrast_levels))

plot_performance_qc <- repeat_performance %>% group_by(Model) %>% summarise(Mean_CV_AUC = mean(AUC), 
    AUC_Q25 = quantile(AUC, 0.25), AUC_Q75 = quantile(AUC, 0.75), Mean_CV_Brier = mean(Brier), Brier_Q25 = quantile(Brier, 
        0.25), Brier_Q75 = quantile(Brier, 0.75), .groups = "drop")

plot_contrast_qc <- contrast_distributions %>% group_by(Contrast) %>% summarise(Mean_paired_DeltaAUC = mean(DeltaAUC), 
    Mean_paired_DeltaBrier = mean(DeltaBrier), .groups = "drop") %>% mutate(Contrast = as.character(Contrast))

merged_performance_qc <- performance %>% select(Model, Mean_CV_AUC, AUC_Q25, AUC_Q75, Mean_CV_Brier, 
    Brier_Q25, Brier_Q75) %>% inner_join(plot_performance_qc, by = "Model", suffix = c("_saved", "_reconstructed"))

merged_contrast_qc <- contrasts %>% select(Contrast, Mean_paired_DeltaAUC, Mean_paired_DeltaBrier) %>% 
    inner_join(plot_contrast_qc, by = "Contrast", suffix = c("_saved", "_reconstructed"))

if (max(abs(as.matrix(merged_performance_qc[, grep("_saved$", names(merged_performance_qc))]) - as.matrix(merged_performance_qc[, 
    grep("_reconstructed$", names(merged_performance_qc))]))) > 1e-12 || max(abs(as.matrix(merged_contrast_qc[, 
    grep("_saved$", names(merged_contrast_qc))]) - as.matrix(merged_contrast_qc[, grep("_reconstructed$", 
    names(merged_contrast_qc))]))) > 1e-12) {
    stop("Plotting reconstruction does not match saved analysis summaries")
}

model_labels <- c(M0 = "SOFA", M1 = "SOFA +\nEBV/HCMV", M2 = "SOFA +\nCTS genes", M3 = "SOFA + EBV/HCMV\n+ CTS genes")

model_colors <- c(M0 = "#7F7F7F", M1 = "#D95F02", M2 = "#1B9E77", M3 = "#7B3294")

repeat_maxima <- repeat_performance %>% group_by(Model) %>% summarise(Max_repeat_AUC = max(AUC), .groups = "drop")

main_data <- performance %>% left_join(repeat_maxima, by = "Model") %>% mutate(Model = factor(Model, 
    levels = model_levels), Label = factor(model_labels[as.character(Model)], levels = model_labels), 
    Mean_label = sprintf("%.3f", Mean_CV_AUC), Mean_label_y = Max_repeat_AUC + 0.0025)

main_repeat_data <- repeat_performance %>% mutate(Model = factor(Model, levels = model_levels), Label = factor(model_labels[as.character(Model)], 
    levels = model_labels))

viral_delta <- contrasts$Mean_paired_DeltaAUC[contrasts$Contrast == "M1 - M0"]

cts_delta <- contrasts$Mean_paired_DeltaAUC[contrasts$Contrast == "M3 - M1"]

viral_label <- sprintf("EBV/HCMV\nAUC %+.3f", viral_delta)

cts_label <- sprintf("CTS genes\nAUC %+.3f", cts_delta)

p_main <- ggplot(main_repeat_data, aes(x = Label, y = AUC, color = Model)) + geom_boxplot(aes(fill = Model), 
    width = 0.32, outlier.shape = NA, alpha = 0.22, linewidth = 0.65, show.legend = FALSE) + geom_jitter(size = 1.7, 
    alpha = 0.55, show.legend = FALSE, position = position_jitter(width = 0.085, height = 0, seed = 260829)) + 
    geom_point(data = main_data, aes(x = Label, y = Mean_CV_AUC, fill = Model), inherit.aes = FALSE, 
        shape = 23, size = 4.3, stroke = 0.8, color = "black", show.legend = FALSE) + geom_text(data = main_data, 
    aes(x = Label, y = Mean_label_y, label = Mean_label), inherit.aes = FALSE, size = 3.5, fontface = "bold") + 
    annotate("segment", x = 1, xend = 2, y = 0.703, yend = 0.703, linewidth = 0.6, color = "#A64200") + 
    annotate("segment", x = 1, xend = 1, y = 0.7005, yend = 0.703, linewidth = 0.6, color = "#A64200") + 
    annotate("segment", x = 2, xend = 2, y = 0.7005, yend = 0.703, linewidth = 0.6, color = "#A64200") + 
    annotate("text", x = 1.5, y = 0.705, label = viral_label, vjust = 0, size = 3.2, color = "#A64200", 
        lineheight = 0.95) + annotate("segment", x = 2, xend = 4, y = 0.716, yend = 0.716, linewidth = 0.6, 
    color = "#5A226C") + annotate("segment", x = 2, xend = 2, y = 0.7135, yend = 0.716, linewidth = 0.6, 
    color = "#5A226C") + annotate("segment", x = 4, xend = 4, y = 0.7135, yend = 0.716, linewidth = 0.6, 
    color = "#5A226C") + annotate("text", x = 3, y = 0.718, label = cts_label, vjust = 0, size = 3.2, 
    color = "#5A226C", lineheight = 0.95) + scale_color_manual(values = model_colors) + scale_fill_manual(values = model_colors) + 
    coord_cartesian(ylim = c(0.62, 0.728), clip = "off") + scale_y_continuous(breaks = seq(0.62, 0.72, 
    0.02)) + labs(title = "Day-1 admission prediction of 28-day mortality", subtitle = sprintf("20 repeats x stratified 5-fold outer CV n = %d (%d deaths, %d survivors)", 
    performance$N[1], performance$Deaths[1], performance$Survivors[1]), x = NULL, y = "Outer-CV AUC", 
    caption = "Small points are repeat-level AUCs; boxes show median and IQR;\ndiamonds show means.") + 
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), 
    plot.subtitle = element_text(hjust = 0.5, size = 10.5), axis.text.x = element_text(color = "black", 
        size = 9, angle = 45, hjust = 1, vjust = 1, lineheight = 0.9), axis.text.y = element_text(color = "black"), 
    plot.caption = element_text(size = 8.5), plot.margin = margin(10, 12, 8, 10))

ggsave(repo_output(repo_path(plot_dir, "260829_Step1_D1_Figure6A.pdf"), "Figure6A_plot.R"), p_main, 
    width = 3.7, height = 4.2, device = "pdf")

write.csv(repeat_performance,"Outputs/Figure6A_repeat_performance.csv",row.names=FALSE)
write.csv(contrast_distributions,"Outputs/Figure6A_paired_contrasts.csv",row.names=FALSE)
