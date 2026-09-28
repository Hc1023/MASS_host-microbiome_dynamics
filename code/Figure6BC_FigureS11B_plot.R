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

prefix <- "260901_Step3_D4_signature"

read_output <- function(suffix) {
    file <- repo_path(analysis_dir, paste0(prefix, "_", suffix, ".csv"))
    if (!repo_exists(file)) 
        stop("Missing analysis output: ", file)
    read.csv(repo_input(file), check.names = FALSE)
}

performance <- read_output("performance")

contrasts <- read_output("contrasts")

repeat_performance <- read_output("repeat_performance")

contrast_distributions <- read_output("contrast_distributions")

outer_models <- read_output("outer_model_summary")

final_coefficients <- read_output("final_coefficients")

final_genes <- read_output("final_genes")

step2_dir <- repo_path(project_dir, "Outputs")

step2_performance <- read.csv(repo_input("Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_performance.csv"), 
    check.names = FALSE)

step2_repeat_performance <- read.csv(repo_input("Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_repeat_performance.csv"), 
    check.names = FALSE)

step2_contrasts <- read.csv(repo_input("Outputs/Figure6B_FigureS11A_analysis_Step2_D4_D1_contrasts.csv"), 
    check.names = FALSE)

model_levels <- c("M0", "M1", "M2", "M3")

main_model_levels <- c("M0", "M1", "M2", "M3")

contrast_levels <- c("M1 - M0", "M2 - M0", "M3 - M0", "M2 - M1", "M3 - M2", "M3 - M1")

stopifnot(identical(performance$Model, model_levels), identical(contrasts$Contrast, contrast_levels), 
    setequal(unique(repeat_performance$Model), model_levels), setequal(unique(contrast_distributions$Contrast), 
        contrast_levels), n_distinct(repeat_performance$Repeat) == 20L, nrow(outer_models) == 300L)

performance_qc <- repeat_performance %>% group_by(Model) %>% summarise(Mean_CV_AUC = mean(AUC), AUC_Q25 = quantile(AUC, 
    0.25), AUC_Q75 = quantile(AUC, 0.75), Mean_CV_Brier = mean(Brier), Brier_Q25 = quantile(Brier, 0.25), 
    Brier_Q75 = quantile(Brier, 0.75), .groups = "drop")

contrast_qc <- contrast_distributions %>% group_by(Contrast) %>% summarise(Mean_DeltaAUC = mean(DeltaAUC), 
    Mean_DeltaBrier = mean(DeltaBrier), .groups = "drop")

if (max(abs(performance$Mean_CV_AUC - performance_qc$Mean_CV_AUC[match(performance$Model, performance_qc$Model)])) > 
    1e-12 || max(abs(contrasts$Mean_DeltaAUC - contrast_qc$Mean_DeltaAUC[match(contrasts$Contrast, contrast_qc$Contrast)])) > 
    1e-12) {
    stop("Detailed plot inputs do not reproduce the saved summaries")
}

figure_model_levels <- paste0("M", 0:5)

model_labels <- c(M0 = "SOFA", M1 = "SOFA +\n(D1) CTS genes", M2 = "SOFA +\nCTS genes", M3 = "SOFA +\nCTS genes +\ntrajectory modules", 
    M4 = "SOFA +\ntrajectory genes", M5 = "SOFA + \ntrajectory genes\n+ CTS genes")

model_colors <- c(M0 = "#777777", M1 = "#1B9E77", M2 = "#4575B4", M3 = "#8C6BB1", M4 = "#7B3294", M5 = "#D95F02")

auc_data <- bind_rows(step2_repeat_performance %>% filter(Model %in% paste0("M", 0:3)) %>% select(Repeat, 
    Model, AUC, Brier), repeat_performance %>% filter(Model %in% c("M2", "M3")) %>% mutate(Model = recode(Model, 
    M2 = "M4", M3 = "M5")) %>% select(Repeat, Model, AUC, Brier)) %>% mutate(Model = factor(Model, levels = figure_model_levels), 
    Label = factor(model_labels[as.character(Model)], levels = model_labels))

auc_summary <- bind_rows(step2_performance %>% filter(Model %in% paste0("M", 0:3)) %>% select(Model, 
    Mean_CV_AUC), performance %>% filter(Model %in% c("M2", "M3")) %>% mutate(Model = recode(Model, M2 = "M4", 
    M3 = "M5")) %>% select(Model, Mean_CV_AUC)) %>% left_join(auc_data %>% group_by(Model) %>% summarise(Max_AUC = max(AUC), 
    .groups = "drop"), by = "Model") %>% mutate(Model = factor(Model, levels = figure_model_levels), 
    Label = factor(model_labels[as.character(Model)], levels = model_labels), Mean_label = sprintf("%.3f", 
        Mean_CV_AUC))

if (n_distinct(auc_data$Repeat) != 20L || any(count(auc_data, Repeat, Model)$n != 1L) || step2_performance$N[1] != 
    performance$N[1] || abs(step2_performance$Mean_CV_AUC[step2_performance$Model == "M2"] - performance$Mean_CV_AUC[performance$Model == 
    "M1"]) > 1e-12) {
    stop("Step2 and Step3 performance outputs are not aligned")
}

auc_min <- min(auc_data$AUC)

auc_max <- max(auc_data$AUC)

auc_span <- auc_max - auc_min

if (!is.finite(auc_span) || auc_span <= 0) auc_span <- 0.1

auc_summary$Mean_label_y <- auc_summary$Max_AUC + 0.055 * auc_span

bracket_tick <- 0.022 * auc_span

bracket_cts_y <- auc_max + 0.17 * auc_span

bracket_1_y <- auc_max + 0.34 * auc_span

bracket_2_y <- bracket_1_y

plot_upper <- auc_max + 0.54 * auc_span

contrast_value <- function(name) contrasts$Mean_DeltaAUC[match(name, contrasts$Contrast)]

step2_contrast_value <- function(name) step2_contrasts$Mean_paired_DeltaAUC[match(name, step2_contrasts$Contrast)]

p_main <- ggplot(auc_data, aes(Label, AUC, color = Model)) + geom_boxplot(aes(fill = Model), width = 0.3, 
    outlier.shape = NA, alpha = 0.22, linewidth = 0.65, show.legend = FALSE) + geom_jitter(size = 1.65, 
    alpha = 0.55, position = position_jitter(width = 0.08, height = 0, seed = 260901), show.legend = FALSE) + 
    geom_point(data = auc_summary, aes(Label, Mean_CV_AUC, fill = Model), inherit.aes = FALSE, shape = 23, 
        size = 4.2, stroke = 0.8, color = "black", show.legend = FALSE) + geom_text(data = auc_summary, 
    aes(Label, Mean_label_y, label = Mean_label), inherit.aes = FALSE, size = 3.4, fontface = "bold") + 
    annotate("segment", x = 2, xend = 3, y = bracket_cts_y, yend = bracket_cts_y, color = "#285A84", 
        linewidth = 0.6) + annotate("segment", x = c(2, 3), xend = c(2, 3), y = bracket_cts_y - bracket_tick, 
    yend = bracket_cts_y, color = "#285A84", linewidth = 0.6) + annotate("text", x = 2.5, y = bracket_cts_y + 
    bracket_tick, label = sprintf("AUC %+.3f", step2_contrast_value("M2 - M1")), vjust = 0, size = 3.1, 
    color = "#285A84") + annotate("segment", x = 3, xend = 4, y = bracket_1_y, yend = bracket_1_y, color = "#8C6BB1", 
    linewidth = 0.6) + annotate("segment", x = c(3, 4), xend = c(3, 4), y = bracket_1_y - bracket_tick, 
    yend = bracket_1_y, color = "#8C6BB1", linewidth = 0.6) + annotate("text", x = 3.5, y = bracket_1_y + 
    bracket_tick, label = sprintf("AUC %+.3f", step2_contrast_value("M3 - M2")), vjust = 0, size = 3.1, 
    color = "#8C6BB1") + annotate("segment", x = 5, xend = 6, y = bracket_2_y, yend = bracket_2_y, color = "#A64B00", 
    linewidth = 0.6) + annotate("segment", x = c(5, 6), xend = c(5, 6), y = bracket_2_y - bracket_tick, 
    yend = bracket_2_y, color = "#A64B00", linewidth = 0.6) + annotate("text", x = 5.5, y = bracket_2_y + 
    bracket_tick, label = sprintf("AUC %+.3f", contrast_value("M3 - M2")), vjust = 0, size = 3.1, color = "#A64B00") + 
    scale_color_manual(values = model_colors) + scale_fill_manual(values = model_colors) + coord_cartesian(ylim = c(auc_min - 
    0.08 * auc_span, plot_upper), clip = "off") + labs(x = NULL, y = "Outer-CV AUC", title = "Day-4 reassessment and gene compression", 
    subtitle = sprintf("%d repeats x stratified %d-fold outer CV; n = %d (%d deaths, %d survivors)", 
        performance$Repeats[1], performance$Outer_folds_per_repeat[1], performance$N[1], performance$Deaths[1], 
        performance$Survivors[1])) + theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", 
    hjust = 0.5, size = 14), plot.subtitle = element_text(hjust = 0.5, size = 10.2), axis.text.x = element_text(color = "black", 
    size = 9.2, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black"), plot.margin = margin(9, 
    12, 8, 9))

p_main <- p_main + labs(title=NULL,subtitle=NULL)

ggsave(repo_output(repo_path(plot_dir, paste0(prefix, "_Figure6B.pdf")), "Figure6BC_FigureS11B_plot.R"), 
    p_main, width = 4.8, height = 4.1, device = "pdf")

compression_levels <- c("M1", "M2", "M3")

compression_labels <- c(M1 = "CTS", M2 = "Trajectory", M3 = "Combined")

compression_colors <- c(M1 = "#4575B4", M2 = "#7B3294", M3 = "#D95F02")

count_data <- outer_models %>% filter(Model %in% compression_levels) %>% mutate(Model = factor(Model, 
    levels = compression_levels), x = as.numeric(Model))

stopifnot(all(count(count_data, Model)$n == 100L))

count_summary <- count_data %>% group_by(Model, x) %>% summarise(Median = median(Selected_gene_n), Q25 = quantile(Selected_gene_n, 
    0.25), Q75 = quantile(Selected_gene_n, 0.75), .groups = "drop")

compression_axis_labels <- setNames(sprintf("%s\nOuter median = %.0f", unname(compression_labels[as.character(count_summary$Model)]), 
    count_summary$Median), as.character(count_summary$Model))

final_fit_summary <- final_genes %>% filter(Model %in% compression_levels) %>% count(Model, name = "Final_gene_n") %>% 
    mutate(Model = factor(Model, levels = compression_levels), x = as.numeric(Model), Final_x = x + 0.2, 
        Detail = ifelse(Model == "M3", "14 trajectory\n+ 2 CTS", NA_character_), Final_label = sprintf("Final = %d", 
            Final_gene_n))

stopifnot(setequal(final_fit_summary$Final_gene_n, c(16L, 17L)))

p_compression <- ggplot(count_data, aes(x, Selected_gene_n, group = Model)) + geom_boxplot(aes(fill = Model), 
    width = 0.3, outlier.shape = NA, alpha = 0.24, linewidth = 0.65, show.legend = FALSE) + geom_jitter(aes(color = Model), 
    width = 0.075, height = 0, size = 1.35, alpha = 0.5, show.legend = FALSE) + geom_point(data = final_fit_summary, 
    aes(x = Final_x, y = Final_gene_n, fill = Model), inherit.aes = FALSE, shape = 23, size = 5.2, stroke = 0.9, 
    color = "black", show.legend = FALSE) + scale_x_continuous(breaks = seq_along(compression_levels), 
    labels = unname(compression_axis_labels[compression_levels]), limits = c(0.55, 3.72)) + scale_y_continuous(expand = expansion(mult = c(0.02, 
    0.08))) + scale_fill_manual(values = compression_colors) + scale_color_manual(values = compression_colors) + 
    labs(x = NULL, y = "Selected genes per outer model", title = "Gene-compression stability", subtitle = "100 outer models; boxes show median and IQR") + 
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), 
    axis.text.x = element_text(color = "black", lineheight = 0.95), axis.text.y = element_text(color = "black"), 
    plot.margin = margin(8, 16, 8, 8))

ggsave(repo_output(repo_path(plot_dir, paste0(prefix, "_FigS2_gene_compression.pdf")), "Figure6BC_FigureS11B_plot.R"), 
    p_compression, width = 4, height = 4, device = "pdf")

coefficient_model_labels <- c(M1 = "CTS", M2 = "trajectory", M3 = "combined")

make_coefficient_plot <- function(models, title) {
    dat <- final_coefficients %>% filter(Model %in% models, Gene_source %in% c("CTS", "trajectory")) %>% 
        mutate(Model = factor(Model, levels = models), Model_label = factor(coefficient_model_labels[as.character(Model)], 
            levels = unname(coefficient_model_labels[models])), Display_gene = paste(Model, Gene, sep = "__")) %>% 
        group_by(Model) %>% arrange(Coefficient, .by_group = TRUE) %>% mutate(Display_gene = factor(Display_gene, 
        levels = unique(Display_gene))) %>% ungroup()
    gene_labels <- setNames(dat$Gene, as.character(dat$Display_gene))
    ggplot(dat, aes(Coefficient, Display_gene, color = Gene_source)) + geom_vline(xintercept = 0, color = "grey65", 
        linetype = "dashed") + geom_segment(aes(x = 0, xend = Coefficient, yend = Display_gene), linewidth = 0.65) + 
        geom_point(size = 2.3) + facet_wrap(~Model_label, scales = "free_y", nrow = 1) + scale_y_discrete(labels = gene_labels) + 
        scale_color_manual(values = c(CTS = "#4575B4", trajectory = "#D95F02"), name = "Gene source") + 
        labs(x = "Full-data elastic-net coefficient", y = NULL, title = title) + theme_classic(base_size = 9.5) + 
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), strip.text = element_text(face = "bold"), 
            legend.position = "bottom")
}

p_coef_cts <- make_coefficient_plot("M1", "Final full-data selected CTS genes")

p_coef_trajectory_combined <- make_coefficient_plot(c("M2", "M3"), "Final full-data selected trajectory and combined genes")

ggsave(repo_output(repo_path(plot_dir, paste0(prefix, "_FigS3A_CTS_final_coefficients.pdf")), "Figure6BC_FigureS11B_plot.R"), 
    p_coef_cts, width = 4.8, height = 6.5, device = "pdf")

ggsave(repo_output(repo_path(plot_dir, paste0(prefix, "_FigS3B_trajectory_combined_final_coefficients.pdf")), 
    "Figure6BC_FigureS11B_plot.R"), p_coef_trajectory_combined, width = 4, height = 5, device = "pdf")

