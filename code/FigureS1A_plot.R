# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(vegan)
})

out_dir <- "Outputs"

repo_dir(out_dir, showWarnings = FALSE, recursive = TRUE)

counts_file <- "Inputs/1708_spcount_all_long.csv"

timepoint_levels <- c("HC", "D1", "D4", "D7", "D14", "D21")

timepoint_palette <- c(HC = "#4D4D4D", D1 = "#1B9E77", D4 = "#D95F02", D7 = "#7570B3", D14 = "#E7298A", 
    D21 = "#66A61E")

message("Reading counts...")

counts_long <- read.csv(repo_input(counts_file), check.names = FALSE, stringsAsFactors = FALSE) %>% mutate(SMRN = as.numeric(SMRN), 
    SMRN = replace_na(SMRN, 0)) %>% filter(!is.na(SampleID), SampleID != "", !is.na(Species), Species != 
    "", SMRN > 0) %>% group_by(SampleID, Species) %>% summarise(SMRN = sum(SMRN), .groups = "drop")

sample_meta <- counts_long %>% group_by(SampleID) %>% summarise(Total_reads = sum(SMRN), Observed_species = n_distinct(Species), 
    .groups = "drop") %>% mutate(Is_HC = str_detect(SampleID, "^(BLRB|RB)"), Timepoint = if_else(Is_HC, 
    "HC", str_extract(SampleID, "D[0-9]+$")), Timepoint = replace_na(Timepoint, "Other"), Timepoint = factor(Timepoint, 
    levels = c(timepoint_levels, "Other")), PatientID = if_else(Is_HC, SampleID, str_remove(SampleID, 
    "_D[0-9]+$"))) %>% arrange(Timepoint, SampleID)

write.csv(sample_meta, repo_output(repo_path(out_dir, "sample_depth_richness_summary.csv"), "FigureS1A_plot.R"), 
    row.names = FALSE)

reference_depth <- 1500

shannon_rarefaction_reps <- 50

reference_coverage <- mean(sample_meta$Total_reads >= reference_depth)

reference_meta <- sample_meta %>% filter(Total_reads >= reference_depth)

cutoff_summary <- tibble(reference_depth = reference_depth, coverage = reference_coverage, total_samples = nrow(sample_meta), 
    samples_at_or_above_reference = nrow(reference_meta), fraction_at_or_above_reference = nrow(reference_meta)/nrow(sample_meta))

write.csv(cutoff_summary, repo_output(repo_path(out_dir, "rarefaction_depth_cutoff_summary.csv"), "FigureS1A_plot.R"), 
    row.names = FALSE)

candidate_cutoffs <- tibble(depth_cutoff = c(floor(quantile(sample_meta$Total_reads, 0.2, names = FALSE)), 
    floor(quantile(sample_meta$Total_reads, 0.5, names = FALSE)), 1500, 1800, 2000)) %>% distinct() %>% 
    mutate(retained_samples = map_int(depth_cutoff, ~sum(sample_meta$Total_reads >= .x)), retained_fraction = retained_samples/nrow(sample_meta))

write.csv(candidate_cutoffs, repo_output(repo_path(out_dir, "rarefaction_candidate_cutoffs.csv"), "FigureS1A_plot.R"), 
    row.names = FALSE)

counts_wide <- counts_long %>% filter(SampleID %in% sample_meta$SampleID) %>% pivot_wider(names_from = Species, 
    values_from = SMRN, values_fill = 0) %>% arrange(match(SampleID, sample_meta$SampleID))

count_mat <- counts_wide %>% select(-SampleID) %>% as.matrix()

rownames(count_mat) <- counts_wide$SampleID

max_depth <- max(sample_meta$Total_reads)

global_depths <- sort(unique(c(round(seq(1, reference_depth, length.out = 45)), reference_depth)))

if (FALSE) {
    message("Calculating Shannon rarefaction curves...")
    set.seed(1708)
    curve_list <- vector("list", length(global_depths))
    for (j in seq_along(global_depths)) {
        depth <- global_depths[j]
        sample_ids <- rownames(count_mat)[rowSums(count_mat) >= depth]
        if (length(sample_ids) == 0) {
            next
        }
        shannon_mat <- matrix(NA_real_, nrow = length(sample_ids), ncol = shannon_rarefaction_reps)
        rownames(shannon_mat) <- sample_ids
        for (rep_i in seq_len(shannon_rarefaction_reps)) {
            rare_counts <- vegan::rrarefy(count_mat[sample_ids, , drop = FALSE], sample = depth)
            shannon_mat[, rep_i] <- vegan::diversity(rare_counts, index = "shannon")
        }
        curve_list[[j]] <- tibble(SampleID = sample_ids, Depth = depth, Shannon = rowMeans(shannon_mat), 
            Shannon_sd = apply(shannon_mat, 1, sd))
    }
    curve_df <- bind_rows(curve_list) %>% left_join(sample_meta, by = "SampleID")
    endpoint_observed_df <- sample_meta %>% transmute(SampleID, Depth = pmin(Total_reads, reference_depth)) %>% 
        anti_join(curve_df %>% select(SampleID, Depth), by = c("SampleID", "Depth"))
    if (nrow(endpoint_observed_df) > 0) {
        endpoint_observed_df <- endpoint_observed_df %>% mutate(Shannon = vegan::diversity(count_mat[SampleID, 
            , drop = FALSE], index = "shannon"), Shannon_sd = NA_real_) %>% left_join(sample_meta, by = "SampleID")
        curve_df <- bind_rows(curve_df, endpoint_observed_df) %>% arrange(Timepoint, SampleID, Depth)
    }
    write.csv(curve_df, repo_output(repo_path(out_dir, "Inputs/rarefaction_curve_points.csv"), "FigureS1A_plot.R"), 
        row.names = FALSE)
} else {
    message("Reading existing Shannon rarefaction curve points...")
    curve_df <- read.csv(repo_input(repo_path(out_dir, "Inputs/rarefaction_curve_points.csv")), check.names = FALSE, 
        stringsAsFactors = FALSE) %>% mutate(Timepoint = factor(Timepoint, levels = c(timepoint_levels, 
        "Other")))
}

summary_df <- curve_df %>% group_by(Timepoint, Depth) %>% summarise(n_samples = n(), Median_shannon = median(Shannon), 
    Q25_shannon = quantile(Shannon, 0.25), Q75_shannon = quantile(Shannon, 0.75), .groups = "drop") %>% 
    filter(n_samples >= 3)

write.csv(summary_df, repo_output(repo_path(out_dir, "rarefaction_curve_group_summary.csv"), "FigureS1A_plot.R"), 
    row.names = FALSE)

endpoint_df <- curve_df %>% filter(Depth == reference_depth) %>% rename(Shannon_at_reference = Shannon)

write.csv(endpoint_df, repo_output(repo_path(out_dir, "rarefaction_reference_depth_shannon.csv"), "FigureS1A_plot.R"), 
    row.names = FALSE)

sample_counts <- sample_meta %>% count(Timepoint, name = "n_samples") %>% mutate(label = paste0(as.character(Timepoint), 
    " (n=", n_samples, ")"))

legend_labels <- setNames(as.character(sample_counts$label), as.character(sample_counts$Timepoint))

curve_depth_breaks <- (scales::breaks_extended(n = 5))(c(0, reference_depth))

depth_breaks <- (scales::breaks_extended(n = 6))(c(0, reference_depth))

p_curve <- ggplot() + geom_line(data = curve_df, aes(Depth, Shannon, group = SampleID, color = Timepoint), 
    linewidth = 0.28, alpha = 0.28) + geom_vline(xintercept = reference_depth, linetype = "22", linewidth = 0.35, 
    color = "grey35") + scale_x_continuous(labels = scales::comma, breaks = curve_depth_breaks, limits = c(0, 
    reference_depth), expand = expansion(mult = c(0.005, 0.02))) + scale_y_continuous(expand = expansion(mult = c(0.01, 
    0.06))) + scale_color_manual(values = timepoint_palette, labels = legend_labels, drop = FALSE) + 
    facet_wrap(~Timepoint, ncol = 3, scales = "free_y") + labs(x = "Subsampled microbial reads (SMRN)", 
    y = "Shannon index", color = NULL, title = "Individual Shannon Rarefaction Curves by Timepoint", 
    subtitle = paste0("No fitted trend is shown; each line is one sample and stops at its observed read depth. ", 
        "Dashed line marks ", scales::comma(reference_depth), " reads.")) + theme_classic(base_size = 9.5) + 
    theme(axis.line = element_line(linewidth = 0.35, color = "grey25"), axis.ticks = element_line(linewidth = 0.3, 
        color = "grey25"), legend.position = "none", legend.justification = "left", legend.key.width = unit(11, 
        "pt"), legend.key.height = unit(10, "pt"), legend.text = element_text(size = 8.2), plot.title = element_text(face = "bold", 
        size = 11.5), plot.subtitle = element_text(color = "grey35", size = 8.5), strip.background = element_rect(fill = "grey94", 
        color = NA), strip.text = element_text(face = "bold", color = "grey20"), panel.grid.major.y = element_line(color = "grey90", 
        linewidth = 0.25), panel.grid.major.x = element_line(color = "grey93", linewidth = 0.22), plot.margin = margin(t = 5.5, 
        r = 18, b = 5.5, l = 5.5))

ggsave(repo_output(repo_path(out_dir, "rarefaction_curve_by_timepoint.pdf"), "FigureS1A_plot.R"), p_curve, 
    width = 6, height = 3)

message("Calculating observed richness rarefaction curves...")

observed_curve_list <- vector("list", nrow(count_mat))

for (i in seq_len(nrow(count_mat))) {
    sample_id <- rownames(count_mat)[i]
    sample_depth <- sum(count_mat[i, ])
    sample_depths <- sort(unique(c(global_depths[global_depths <= min(sample_depth, reference_depth)], 
        min(sample_depth, reference_depth))))
    observed_curve_list[[i]] <- tibble(SampleID = sample_id, Depth = sample_depths, Observed_richness = as.numeric(vegan::rarefy(count_mat[i, 
        , drop = FALSE], sample = sample_depths, se = FALSE)))
}

observed_curve_df <- bind_rows(observed_curve_list) %>% left_join(sample_meta, by = "SampleID")

write.csv(observed_curve_df, repo_output(repo_path(out_dir, "observed_richness_rarefaction_curve_points.csv"), 
    "FigureS1A_plot.R"), row.names = FALSE)

p_observed_curve <- ggplot() + geom_line(data = observed_curve_df, aes(Depth, Observed_richness, group = SampleID, 
    color = Timepoint), linewidth = 0.28, alpha = 0.28) + geom_vline(xintercept = reference_depth, linetype = "22", 
    linewidth = 0.35, color = "grey35") + scale_x_continuous(labels = scales::comma, breaks = curve_depth_breaks, 
    limits = c(0, reference_depth), expand = expansion(mult = c(0.005, 0.02))) + scale_y_continuous(expand = expansion(mult = c(0.01, 
    0.06))) + scale_color_manual(values = timepoint_palette, labels = legend_labels, drop = FALSE) + 
    facet_wrap(~Timepoint, ncol = 3, scales = "free_y") + labs(x = "Subsampled microbial reads (SMRN)", 
    y = "Observed richness", color = NULL, title = "Individual Observed Richness Rarefaction Curves by Timepoint", 
    subtitle = paste0("No fitted trend is shown; each line is one sample and stops at its observed read depth. ", 
        "Dashed line marks ", scales::comma(reference_depth), " reads.")) + theme_classic(base_size = 9.5) + 
    theme(axis.line = element_line(linewidth = 0.35, color = "grey25"), axis.ticks = element_line(linewidth = 0.3, 
        color = "grey25"), legend.position = "none", legend.justification = "left", legend.key.width = unit(11, 
        "pt"), legend.key.height = unit(10, "pt"), legend.text = element_text(size = 8.2), plot.title = element_text(face = "bold", 
        size = 11.5), plot.subtitle = element_text(color = "grey35", size = 8.5), strip.background = element_rect(fill = "grey94", 
        color = NA), strip.text = element_text(face = "bold", color = "grey20"), panel.grid.major.y = element_line(color = "grey90", 
        linewidth = 0.25), panel.grid.major.x = element_line(color = "grey93", linewidth = 0.22), plot.margin = margin(t = 5.5, 
        r = 18, b = 5.5, l = 5.5))

ggsave(repo_output(repo_path(out_dir, "observed_richness_rarefaction_by_timepoint.pdf"), "FigureS1A_plot.R"), 
    p_observed_curve, width = 6, height = 3)

