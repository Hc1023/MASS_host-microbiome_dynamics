# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
    library(openxlsx)
    library(ggrepel)
    library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)

script_dir <- "."

out_dir <- repo_path(script_dir, "Outputs")

repo_dir(out_dir, showWarnings = FALSE, recursive = TRUE)

nc_file <- "Inputs/NC_RNA_species_summary.csv"

raw_count_file <- "Inputs/PBMC_species_count_all.anno.csv"

expected_nc_controls <- 50

freq_cutoff <- 0.5

mean_ra_pct_cutoff <- 0.1

blrb_mean_ra_cutoff <- mean_ra_pct_cutoff/100

theme_pub <- function(base_size = 10) {
    theme_bw(base_size = base_size) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.25, 
        color = "grey88"), plot.title = element_text(face = "bold", size = base_size + 2), plot.subtitle = element_text(color = "grey30"), 
        axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"), strip.background = element_rect(fill = "grey95", 
            color = "grey70"), strip.text = element_text(face = "bold"))
}

save_plot <- function(plot, name, width, height) {
    ggsave(repo_output(repo_path(out_dir, paste0(name, ".pdf")), "FigureS1B_analysis.R"), plot, width = width, 
        height = height, bg = "white")
    ggsave(repo_output(repo_path(out_dir, paste0(name, ".png")), "FigureS1B_analysis.R"), plot, width = width, 
        height = height, dpi = 300, bg = "white")
}

message("Reading NC species summary...")

nc_raw <- read.csv(repo_input(nc_file), check.names = FALSE)

required_nc_cols <- c("Species", "Type", "Genus", "Detected_Samples_n", "Total_NC_Samples_N", "Detection_Frequency", 
    "Median_RelAbund_pct", "Mean_RelAbund_pct")

stopifnot(all(required_nc_cols %in% colnames(nc_raw)))

nc <- nc_raw %>% mutate(across(c(Detected_Samples_n, Total_NC_Samples_N, Detection_Frequency, Median_RelAbund_pct, 
    Q1_RelAbund_pct, Q3_RelAbund_pct, IQR_RelAbund_pct, Mean_RelAbund_pct, Std_RelAbund_pct, Mean_RA_DetectedOnly_pct, 
    Std_RA_DetectedOnly_pct), as.numeric), Type = replace_na(Type, "Unknown"), Genus = replace_na(Genus, 
    "Unknown"), Species = str_replace_all(Species, "\\s+", "_"), is_recurrent_background = Detection_Frequency > 
    freq_cutoff & Mean_RelAbund_pct > mean_ra_pct_cutoff, background_status = if_else(is_recurrent_background, 
    "Recurrent background: frequency >50% and mean RA >0.1%", "Below recurrent-background threshold"), 
    neglog10_mean_ra = -log10(pmax(Mean_RelAbund_pct, .Machine$double.xmin)), Species_label = str_replace_all(Species, 
        "_", " "))

nc_control_n <- sort(unique(na.omit(nc$Total_NC_Samples_N)))

if (length(nc_control_n) != 1 || nc_control_n != expected_nc_controls) {
    warning("The NC summary file reports Total_NC_Samples_N = ", paste(nc_control_n, collapse = ", "), 
        ", not the expected n = ", expected_nc_controls, ". Plots and summaries use the values present in the file.")
}

nc_recurrent <- nc %>% filter(is_recurrent_background) %>% arrange(desc(Detection_Frequency), desc(Mean_RelAbund_pct), 
    desc(Detected_Samples_n))

nc_retained <- nc %>% filter(!is_recurrent_background) %>% arrange(desc(Detection_Frequency), desc(Mean_RelAbund_pct), 
    desc(Detected_Samples_n))

write.csv(nc, repo_output(repo_path(out_dir, "NC_species_summary_with_background_flags.csv"), "FigureS1B_analysis.R"), 
    row.names = FALSE)

write.csv(nc_recurrent, repo_output(repo_path(out_dir, "NC_recurrent_background_taxa.csv"), "FigureS1B_analysis.R"), 
    row.names = FALSE)

write.csv(nc_retained, repo_output(repo_path(out_dir, "NC_taxa_below_recurrent_background_threshold.csv"), 
    "FigureS1B_analysis.R"), row.names = FALSE)

summary_tbl <- tibble(metric = c("expected_NC_controls_from_methods_text", "NC_controls_reported_in_input_file", 
    "NC_species_annotations", "NC_taxa_frequency_gt_50pct", "NC_taxa_mean_RA_gt_0.1pct", "NC_recurrent_background_taxa_frequency_gt_50pct_and_mean_RA_gt_0.1pct", 
    "NC_taxa_below_recurrent_background_threshold"), value = c(expected_nc_controls, paste(nc_control_n, 
    collapse = ";"), nrow(nc), sum(nc$Detection_Frequency > freq_cutoff, na.rm = TRUE), sum(nc$Mean_RelAbund_pct > 
    mean_ra_pct_cutoff, na.rm = TRUE), nrow(nc_recurrent), nrow(nc_retained)))

write.csv(summary_tbl, repo_output(repo_path(out_dir, "NC_background_filter_summary.csv"), "FigureS1B_analysis.R"), 
    row.names = FALSE)

top_labels <- nc %>% filter(is_recurrent_background) %>% slice_max(order_by = Mean_RelAbund_pct * Detection_Frequency, 
    n = 5, with_ties = FALSE)

p_threshold <- ggplot(nc, aes(Detection_Frequency, Mean_RelAbund_pct)) + geom_hline(yintercept = mean_ra_pct_cutoff, 
    linetype = "dashed", linewidth = 0.45, color = "grey35") + geom_vline(xintercept = freq_cutoff, linetype = "dashed", 
    linewidth = 0.45, color = "grey35") + geom_point(aes(color = background_status, size = Detected_Samples_n), 
    alpha = 0.72) + ggrepel::geom_text_repel(data = top_labels, aes(label = Species_label), size = 2.5, 
    max.overlaps = Inf, min.segment.length = 0, box.padding = 0.25, seed = 1718) + scale_x_continuous(labels = scales::percent_format(accuracy = 1), 
    limits = c(0, 1.02)) + scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01), breaks = c(0.001, 
    0.01, 0.03, 0.1, 0.3, 1, 3), labels = function(x) paste0(x, "%")) + scale_color_manual(values = c(`Below recurrent-background threshold` = "#4C78A8", 
    `Recurrent background: frequency >50% and mean RA >0.1%` = "#D55E00")) + scale_size_continuous(range = c(1, 
    4.5), guide = "none") + labs(x = "Detection frequency in NC plasma controls", y = "Mean relative abundance in NC controls", 
    title = "Annotation-Level Background Taxa Defined from NC Controls", subtitle = paste0(nrow(nc_recurrent), 
        " recurrent taxa exceed both thresholds; input file reports NC n = ", paste(nc_control_n, collapse = "/")), 
    color = NULL) + theme_pub(10) + theme(legend.position = "bottom")

save_plot(p_threshold, "plot_NC_frequency_vs_meanRA_thresholds", 8.2, 6.2)

top_n <- 50

top_background <- nc_recurrent %>% arrange(desc(Detection_Frequency), desc(Mean_RelAbund_pct), desc(Detected_Samples_n)) %>% 
    slice_head(n = top_n) %>% mutate(Species_label = factor(Species_label, levels = rev(Species_label)), 
    Mean_RA_bin = cut(Mean_RelAbund_pct, breaks = c(-Inf, 0.1, 0.3, 1, Inf), labels = c("0.1-0.3%", "0.3-1%", 
        "1-3%", ">=3%"), right = FALSE))

p_top_prev <- ggplot(top_background, aes(Detection_Frequency, Species_label, fill = Mean_RA_bin)) + geom_col(width = 0.72) + 
    geom_errorbar(aes(xmin = pmax(CI95_Lower, 0), xmax = pmin(CI95_Upper, 1)), orientation = "y", width = 0.18, 
        linewidth = 0.3, color = "grey25") + scale_x_continuous(labels = scales::percent_format(accuracy = 1), 
    limits = c(0, 1.02)) + scale_fill_manual(values = c(`0.1-0.3%` = "#8DD3C7", `0.3-1%` = "#80B1D3", 
    `1-3%` = "#FDB462", `>=3%` = "#B2182B"), drop = FALSE) + labs(x = "Detection frequency in NC controls with 95% CI", 
    y = NULL, title = "Top Recurrent Background Species in NC Controls", subtitle = "Species shown passed frequency >50% and mean relative abundance >0.1%", 
    fill = "Mean RA") + theme_pub(9) + theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

save_plot(p_top_prev, "plot_NC_top50_recurrent_background_species", 8.2, 9.5)

p_dist_freq <- ggplot(nc, aes(Detection_Frequency, fill = background_status)) + geom_histogram(binwidth = 1/60, 
    boundary = 0, color = "white", linewidth = 0.15) + geom_vline(xintercept = freq_cutoff, linetype = "dashed", 
    linewidth = 0.45, color = "grey25") + scale_x_continuous(labels = function(x) paste0(100 * x)) + 
    scale_fill_manual(values = c(`Below recurrent-background threshold` = "#4C78A8", `Recurrent background: frequency >50% and mean RA >0.1%` = "#D55E00")) + 
    labs(x = "Detection frequency (%)", y = "Number of species", title = "NC Prevalence Distribution", 
        fill = NULL) + theme_pub(10) + theme(legend.position = "none")

p_dist_ra <- ggplot(nc, aes(Mean_RelAbund_pct, fill = background_status)) + geom_histogram(bins = 45, 
    color = "white", linewidth = 0.15) + geom_vline(xintercept = mean_ra_pct_cutoff, linetype = "dashed", 
    linewidth = 0.45, color = "grey25") + scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.01), 
    breaks = c(0, 0.01, 0.03, 0.1, 0.3, 1, 3), labels = function(x) paste0(x)) + scale_fill_manual(values = c(`Below recurrent-background threshold` = "#4C78A8", 
    `Recurrent background: frequency >50% and mean RA >0.1%` = "#D55E00")) + labs(x = "Mean relative abundance (%)", 
    y = "Number of species", title = "NC Abundance Distribution", fill = NULL) + theme_pub(10) + theme(legend.position = "bottom")

p_dist <- p_dist_freq + p_dist_ra + plot_layout(widths = c(1, 1), guides = "collect") & theme(legend.position = "bottom")

save_plot(p_dist, "plot_NC_background_taxa_distributions", 5.5, 2.5)

message("Applying clinical sample-level relative-abundance filter...")

raw_counts <- read.csv(repo_input(raw_count_file), check.names = FALSE, stringsAsFactors = FALSE, colClasses = c("character", 
    "character", "character", "NULL", "character", "NULL", "numeric")) %>% rename(Sample = `#Sample`) %>% 
    mutate(Species = str_replace_all(Species, "\\s+", "_"), SMRN = as.numeric(SMRN), SMRN = replace_na(SMRN, 
        0), Type = replace_na(Type, "Unknown"), Genus = replace_na(Genus, "Unknown"), is_BLRB_whole_blood_HC = str_detect(Sample, 
        "BLRB")) %>% filter(!is.na(Sample), Sample != "", !is.na(Species), Species != "")

clinical_calls <- raw_counts %>% filter(!is_BLRB_whole_blood_HC) %>% group_by(Sample, Species, Type, 
    Genus) %>% summarise(SMRN = sum(SMRN), .groups = "drop") %>% group_by(Sample) %>% mutate(sample_total_SMRN = sum(SMRN), 
    clinical_RA_pct = if_else(sample_total_SMRN > 0, 100 * SMRN/sample_total_SMRN, 0)) %>% ungroup() %>% 
    left_join(nc %>% select(Species, NC_detection_frequency = Detection_Frequency, NC_mean_RA_pct = Mean_RelAbund_pct, 
        NC_median_RA_pct = Median_RelAbund_pct, NC_recurrent_background = is_recurrent_background), by = "Species") %>% 
    mutate(NC_recurrent_background = replace_na(NC_recurrent_background, FALSE), NC_median_RA_pct = replace_na(NC_median_RA_pct, 
        0), sample_level_pass = clinical_RA_pct >= NC_median_RA_pct, filter_step = case_when(NC_recurrent_background ~ 
        "Removed by global recurrent-background filter", !sample_level_pass ~ "Removed by sample-level median-RA filter", 
        TRUE ~ "Retained after both NC filters"), Species_label = str_replace_all(Species, "_", " "))

clinical_after_global <- clinical_calls %>% filter(!NC_recurrent_background)

clinical_after_sample <- clinical_after_global %>% filter(sample_level_pass)

clinical_filter_summary <- tibble(step = c("Raw clinical sample-species calls", "After global recurrent-background removal", 
    "After sample-level median-RA filter"), sample_species_calls = c(nrow(clinical_calls), nrow(clinical_after_global), 
    nrow(clinical_after_sample)), species_annotations = c(n_distinct(clinical_calls$Species), n_distinct(clinical_after_global$Species), 
    n_distinct(clinical_after_sample$Species)), clinical_samples = c(n_distinct(clinical_calls$Sample), 
    n_distinct(clinical_after_global$Sample), n_distinct(clinical_after_sample$Sample)))

write.csv(clinical_calls, repo_output(repo_path(out_dir, "clinical_sample_species_calls_with_NC_filter_flags.csv"), 
    "FigureS1B_analysis.R"), row.names = FALSE)

write.csv(clinical_after_sample, repo_output(repo_path(out_dir, "clinical_sample_species_calls_retained_after_NC_filters.csv"), 
    "FigureS1B_analysis.R"), row.names = FALSE)

write.csv(clinical_filter_summary, repo_output(repo_path(out_dir, "clinical_sample_level_filter_summary.csv"), 
    "FigureS1B_analysis.R"), row.names = FALSE)

summary_tbl <- bind_rows(summary_tbl, tibble(metric = c("clinical_raw_sample_species_calls", "clinical_raw_species_annotations_in_count_file", 
    "clinical_calls_after_global_recurrent_background_removal", "clinical_species_after_global_recurrent_background_removal", 
    "clinical_calls_after_sample_level_median_RA_filter", "clinical_species_after_sample_level_median_RA_filter"), 
    value = as.character(c(nrow(clinical_calls), n_distinct(clinical_calls$Species), nrow(clinical_after_global), 
        n_distinct(clinical_after_global$Species), nrow(clinical_after_sample), n_distinct(clinical_after_sample$Species)))))

write.csv(summary_tbl, repo_output(repo_path(out_dir, "NC_background_filter_summary.csv"), "FigureS1B_analysis.R"), 
    row.names = FALSE)

clinical_step_long <- clinical_filter_summary %>% mutate(step = factor(step, levels = step)) %>% pivot_longer(c(sample_species_calls, 
    species_annotations), names_to = "metric", values_to = "n") %>% mutate(metric = recode(metric, sample_species_calls = "Sample-species calls", 
    species_annotations = "Species annotations"))

p_clinical_filter_counts <- ggplot(clinical_step_long, aes(step, n, fill = metric)) + geom_col(position = position_dodge(width = 0.72), 
    width = 0.62) + geom_text(aes(label = scales::comma(n)), position = position_dodge(width = 0.72), 
    vjust = -0.22, size = 3) + scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 
    0.14))) + scale_fill_manual(values = c(`Sample-species calls` = "#4C78A8", `Species annotations` = "#D55E00")) + 
    labs(x = NULL, y = "Count", title = "Clinical Calls Are Reduced by the Two NC-Based Annotation Filters", 
        subtitle = "The sample-level filter removes calls with clinical RA below the same species median RA in NC plasma controls", 
        fill = NULL) + theme_pub(10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 18, 
    hjust = 1))

save_plot(p_clinical_filter_counts, "plot_clinical_call_counts_across_NC_filters", 8.4, 5.2)

sample_filter_plot_data <- clinical_after_global %>% mutate(sample_level_status = if_else(sample_level_pass, 
    "Retained: clinical RA >= NC median RA", "Removed: clinical RA < NC median RA"))

sample_filter_labels <- sample_filter_plot_data %>% filter(sample_level_pass) %>% group_by(Species, Species_label) %>% 
    summarise(max_clinical_RA_pct = max(clinical_RA_pct, na.rm = TRUE), NC_median_RA_pct = max(NC_median_RA_pct, 
        na.rm = TRUE), retained_calls = n(), .groups = "drop") %>% slice_max(order_by = max_clinical_RA_pct, 
    n = 12, with_ties = FALSE)

p_sample_filter <- ggplot(sample_filter_plot_data, aes(NC_median_RA_pct, clinical_RA_pct, color = sample_level_status)) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.45, color = "grey35") + 
    geom_point(alpha = 0.42, size = 1.1) + ggrepel::geom_text_repel(data = sample_filter_plot_data %>% 
    semi_join(sample_filter_labels, by = "Species") %>% group_by(Species, Species_label) %>% slice_max(clinical_RA_pct, 
    n = 1, with_ties = FALSE) %>% ungroup(), aes(label = Species_label), size = 2.35, max.overlaps = Inf, 
    min.segment.length = 0, seed = 1720, show.legend = FALSE) + scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.001), 
    breaks = c(0, 0.001, 0.01, 0.03, 0.1, 0.3, 1, 3), labels = function(x) paste0(x)) + scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.001), 
    breaks = c(0, 0.001, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100), labels = function(x) paste0(x)) + 
    scale_color_manual(values = c(`Retained: clinical RA >= NC median RA` = "#1B9E77", `Removed: clinical RA < NC median RA` = "#D55E00")) + 
    labs(x = "Median relative abundance in NC plasma controls (%)", y = "Relative abundance in \neach clinical sample (%)", 
        title = "Sample-Level Background Filter Uses Species-Specific NC Median RA", subtitle = paste0(scales::comma(nrow(clinical_after_sample)), 
            " of ", scales::comma(nrow(clinical_after_global)), " post-global-filter sample-species calls retained"), 
        color = NULL) + theme_pub(10) + theme(legend.position = "bottom")

save_plot(p_sample_filter, "plot_sample_level_RA_filter_against_NC_median", 5.6, 4)

