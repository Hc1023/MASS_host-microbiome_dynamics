# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(ggrepel)

library(patchwork)

if (Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX")) try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

project_dir <- "."

results_dir <- repo_path(project_dir, "Outputs")

out_dir <- repo_path(project_dir, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(day1 = "Inputs/all_GOBP_D1_mortality_divergence.csv", timepoint = "Inputs/all_GOBP_D1_D4_D7_mortality_differences.csv", 
    interaction = "Inputs/all_GOBP_D1_to_D4_D7_extra_divergence.csv", trajectory = "Inputs/all_GOBP_trajectory_divergence_omnibus.csv", 
    trajectory_complete = "Inputs/all_GOBP_trajectory_statistics_complete.csv", group_trajectory = "Inputs/all_GOBP_Survivor_Death_D1_D4_D7_trajectories.csv")

input_paths <- setNames(repo_path(results_dir, unname(input_files)), names(input_files))

if (!all(repo_exists(input_paths))) {
    stop("Missing model result files: ", paste(input_paths[!repo_exists(input_paths)], collapse = ", "))
}

day1_results <- read.csv(repo_input(input_paths[["day1"]]))

timepoint_results <- read.csv(repo_input(input_paths[["timepoint"]]))

interaction_results <- read.csv(repo_input(input_paths[["interaction"]]))

trajectory_results <- read.csv(repo_input(input_paths[["trajectory"]]))

trajectory_complete <- read.csv(repo_input(input_paths[["trajectory_complete"]]))

group_trajectory_results <- read.csv(repo_input(input_paths[["group_trajectory"]]))

clean_label <- function(x, width = 40) {
    x %>% str_remove("^GOBP_") %>% str_replace_all("_", " ") %>% str_to_sentence() %>% str_wrap(width)
}

format_test_value <- function(x) {
    ifelse(x < 0.001, "< 0.001", paste0("= ", format.pval(x, digits = 2)))
}

theme_manuscript <- function(base_size = 10) {
    theme_classic(base_size = base_size) + theme(plot.title = element_text(face = "bold", size = base_size + 
        3), plot.subtitle = element_text(size = base_size - 1, colour = "grey30"), axis.title = element_text(colour = "grey10"), 
        axis.text = element_text(colour = "grey15"), strip.background = element_rect(fill = "grey93", 
            colour = "grey35"), strip.text = element_text(face = "bold"), legend.position = "bottom")
}

save_plot <- function(plot, stem, width, height) {
    ggsave(repo_output(repo_path(out_dir, paste0(stem, ".pdf")), "Figure2A_Figure3A_FigureS5A_plot.R"), plot, 
        width = width, height = height, device = grDevices::pdf, encoding = "WinAnsi.enc", useDingbats = FALSE)
}

up_pattern <- paste(c("CELL_CYCLE", "CELL_DIVISION", "MITOTIC", "CHROMOSOME_SEGREGATION", "CHROMATID_SEGREGATION", 
    "DNA_REPLICATION"), collapse = "|")

down_pattern <- paste(c("T_CELL_ACTIVATION", "T_HELPER", "DENDRITIC_CELL_ANTIGEN_PROCESSING", "ANTIGEN_PROCESSING_AND_PRESENTATION"), 
    collapse = "|")

day1_plot_df <- day1_results %>% mutate(minus_log10_fdr = -log10(pmax(p_adj_BH, .Machine$double.xmin)), 
    Program = case_when(p_adj_BH < 0.05 & estimate > 0 & str_detect(Pathway, up_pattern) ~ "Cell cycle / chromosome / DNA replication", 
        p_adj_BH < 0.05 & estimate < 0 & str_detect(Pathway, down_pattern) ~ "T-cell activation / antigen presentation", 
        p_adj_BH < 0.05 ~ "Other FDR-significant pathway", TRUE ~ "Not significant"), Program = factor(Program, 
        levels = c("Not significant", "Other FDR-significant pathway", "Cell cycle / chromosome / DNA replication", 
            "T-cell activation / antigen presentation")))

day1_label_pathways <- c("GOBP_POSITIVE_REGULATION_OF_CELL_DIVISION", "GOBP_REGULATION_OF_MITOTIC_SISTER_CHROMATID_SEGREGATION", 
    "GOBP_DNA_REPLICATION_INITIATION", "GOBP_T_CELL_ACTIVATION_VIA_T_CELL_RECEPTOR_CONTACT_WITH_ANTIGEN_BOUND_TO_MHC_MOLECULE_ON_ANTIGEN_PRESENTING_CELL", 
    "GOBP_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN")

day1_label_map <- c(GOBP_POSITIVE_REGULATION_OF_CELL_DIVISION = "Positive regulation of\ncell division", 
    GOBP_REGULATION_OF_MITOTIC_SISTER_CHROMATID_SEGREGATION = "Mitotic sister chromatid\nsegregation", 
    GOBP_DNA_REPLICATION_INITIATION = "DNA replication\ninitiation", GOBP_T_CELL_ACTIVATION_VIA_T_CELL_RECEPTOR_CONTACT_WITH_ANTIGEN_BOUND_TO_MHC_MOLECULE_ON_ANTIGEN_PRESENTING_CELL = "T-cell activation via\nTCR\342\200\223MHC interaction", 
    GOBP_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION = "Dendritic-cell antigen\nprocessing/presentation", 
    GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN = "Exogenous antigen\nprocessing/presentation")

day1_label_df <- day1_plot_df %>% filter(Pathway %in% names(day1_label_map)) %>% mutate(label = unname(day1_label_map[Pathway]))

p_day1 <- ggplot(day1_plot_df, aes(estimate, minus_log10_fdr)) + geom_hline(yintercept = -log10(0.05), 
    linetype = 2, linewidth = 0.4, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.35, 
    colour = "grey70") + geom_point(aes(colour = Program), size = 1.25) + scale_colour_manual(values = c(`Not significant` = alpha("grey80", 
    0.6), `Other FDR-significant pathway` = alpha("grey48", 0.6), `Cell cycle / chromosome / DNA replication` = "#D73027", 
    `T-cell activation / antigen presentation` = "#4575B4"), labels = c("Not significant", "Other FDR-significant", 
    "Cell cycle / DNA replication", "T-cell activation / antigen presentation")) + labs(subtitle = "Global Day-1 GO:BP mortality effects", 
    x = "Adjusted mortality effect", y = expression(-log[10](BH - FDR)), colour = NULL) + guides(colour = guide_legend(nrow = 2, 
    byrow = TRUE)) + theme_bw() + theme(legend.position = "bottom", legend.spacing.x = unit(1.5, "mm"), 
    legend.key.width = unit(3, "mm"), legend.key.height = unit(3.5, "mm"), legend.box.spacing = unit(1, 
        "mm"), legend.margin = margin(0, 0, 0, 0))

p_day1

save_plot(p_day1, "Fig_A_global_D1_effect_notext", 4.6, 4.2)

d4_results <- interaction_results %>% filter(Estimand == "D1 to D4 extra divergence")

baseline_dynamic_df <- day1_results %>% dplyr::select(Pathway, D1_effect = estimate, D1_FDR = p_adj_BH) %>% 
    left_join(d4_results %>% dplyr::select(Pathway, D4_interaction = estimate, D4_interaction_FDR = p_adj_BH), 
        by = "Pathway") %>% left_join(trajectory_results %>% dplyr::select(Pathway, trajectory_FDR = p_adj_BH), 
    by = "Pathway") %>% mutate(D1_significant = D1_FDR < 0.05, trajectory_significant = trajectory_FDR < 
    0.05, Evidence = case_when(D1_significant & trajectory_significant ~ "Both", D1_significant ~ "D1 only", 
    trajectory_significant ~ "Trajectory only", TRUE ~ "Neither"), Evidence = factor(Evidence, levels = c("Neither", 
    "D1 only", "Trajectory only", "Both")))

observed_counts <- baseline_dynamic_df %>% count(Evidence, .drop = FALSE)

expected_counts <- c(Neither = 4332L, `D1 only` = 199L, `Trajectory only` = 29L, Both = 1L)

stopifnot(setNames(observed_counts$n, as.character(observed_counts$Evidence)) == expected_counts)

evidence_labels <- setNames(paste0(names(expected_counts), " (n=", expected_counts, ")"), names(expected_counts))

dynamic_label_pathways <- c("GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_DETECTION_OF_MOLECULE_OF_BACTERIAL_ORIGIN", 
    "GOBP_PHAGOSOME_LYSOSOME_FUSION", "GOBP_RIBOSOMAL_LARGE_SUBUNIT_ASSEMBLY")

representative_labels <- setNames(c("MyD88\342\200\223TLR signaling", "Bacterial sensing", "Phagosome\342\200\223lysosome fusion", 
    "Ribosomal assembly"), dynamic_label_pathways)

p_baseline_dynamic <- ggplot(baseline_dynamic_df, aes(D1_effect, D4_interaction)) + geom_hline(yintercept = 0, 
    colour = "grey70", linewidth = 0.35) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) + 
    geom_point(aes(colour = Evidence), size = 1.25) + geom_text_repel(data = filter(baseline_dynamic_df, 
    Pathway %in% dynamic_label_pathways), aes(label = unname(representative_labels[Pathway])), size = 2.5, 
    colour = "grey10", max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0, seed = 260818) + 
    scale_colour_manual(values = c(Neither = alpha("grey80", 0.2), `D1 only` = "#3B78A8", `Trajectory only` = "#D95F02", 
        Both = "#7B3294"), labels = evidence_labels, drop = FALSE) + labs(x = "D1 mortality effect", 
    y = "D1 to D4 interaction effect", colour = NULL) + guides(colour = guide_legend(nrow = 2, byrow = TRUE)) + 
    theme_bw() + theme(legend.position = "bottom", legend.spacing.x = unit(1.5, "mm"), legend.key.width = unit(3, 
    "mm"), legend.key.height = unit(3.5, "mm"), legend.box.spacing = unit(1, "mm"), legend.margin = margin(0, 
    0, 0, 0))

p_baseline_dynamic

save_plot(p_baseline_dynamic, "Fig_B_D1_vs_D4_interaction", 3.5, 3.5)

representative_pathways <- c("GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_DETECTION_OF_MOLECULE_OF_BACTERIAL_ORIGIN", 
    "GOBP_PHAGOSOME_LYSOSOME_FUSION", "GOBP_RIBOSOMAL_LARGE_SUBUNIT_ASSEMBLY")

representative_df <- group_trajectory_results %>% filter(Pathway %in% representative_pathways) %>% mutate(Timepoint = factor(Timepoint, 
    levels = c("D1", "D4", "D7")), Mortality28d = factor(Mortality28d, levels = c("Survivor", "Death")), 
    Pathway_label = factor(representative_labels[Pathway], levels = unname(representative_labels[representative_pathways])))

representative_annotations <- trajectory_complete %>% filter(Pathway %in% representative_pathways) %>% 
    transmute(Pathway, Pathway_label = factor(representative_labels[Pathway], levels = unname(representative_labels[representative_pathways])), 
        label = paste0("Overall trajectory: P ", format_test_value(overall_trajectory_P), "; FDR ", format_test_value(overall_trajectory_FDR)), 
        Timepoint = factor("D1", levels = c("D1", "D4", "D7"))) %>% left_join(representative_df %>% group_by(Pathway) %>% 
    summarise(y = max(emmean + SE) + 0.1, .groups = "drop"), by = "Pathway")

representative_annotations$label = gsub("Overall t", "T", representative_annotations$label)

p_representative <- ggplot(representative_df, aes(Timepoint, emmean, group = Mortality28d, colour = Mortality28d)) + 
    geom_line(linewidth = 0.85) + geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.07, 
    linewidth = 0.48, position = position_dodge(width = 0.04)) + geom_point(size = 2.2) + geom_text(data = representative_annotations, 
    aes(x = Timepoint, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 2.7) + facet_wrap(~Pathway_label, 
    scales = "free_y", ncol = 4) + scale_colour_manual(values = c(Survivor = "#4575B4", Death = "#D73027"), 
    label = c("Survival", "Mortality")) + labs(x = NULL, y = "Adjusted standardized GSVA score", colour = "28-day outcome") + 
    scale_x_discrete(limits = c("D1", "D4", "D7"), expand = expansion(add = c(0.15, 0.15))) + scale_y_continuous(expand = expansion(mult = c(0, 
    0.1))) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.title.y = element_text(size = 9), 
    legend.position = "top")

p_representative

save_plot(p_representative, "Fig_C_representative_group_trajectories", 10.5, 2.7)

trajectory_30 <- trajectory_complete %>% filter(overall_trajectory_FDR < 0.05) %>% arrange(overall_trajectory_FDR)

stopifnot(nrow(trajectory_30) == 30)

effect_order <- c("D1 mortality difference", "D4 mortality difference", "D7 mortality difference", "D1 to D4 extra divergence", 
    "D1 to D7 extra divergence")

effect_labels <- c("D1 effect", "D4 effect", "D7 effect", "D1->D4 interaction", "D1->D7 interaction")

trajectory_heatmap_df <- bind_rows(timepoint_results, interaction_results) %>% filter(Pathway %in% trajectory_30$Pathway, 
    Estimand %in% effect_order) %>% mutate(Column = factor(Estimand, levels = effect_order, labels = effect_labels))

heat_matrix <- trajectory_heatmap_df %>% dplyr::select(Pathway, Column, estimate) %>% pivot_wider(names_from = Column, 
    values_from = estimate) %>% column_to_rownames("Pathway") %>% as.matrix()

row_order <- rownames(heat_matrix)[hclust(dist(heat_matrix))$order]

trajectory_heatmap_df <- trajectory_heatmap_df %>% mutate(Pathway_label = factor(clean_label(Pathway, 
    58), levels = clean_label(rev(row_order), 58)))

trajectory_heatmap_df <- trajectory_heatmap_df %>% mutate(Column = factor(as.character(Column), levels = effect_labels))

heat_limit <- max(abs(trajectory_heatmap_df$estimate), na.rm = TRUE)

p_trajectory_heatmap <- ggplot() + geom_tile(data = trajectory_heatmap_df, aes(Column, Pathway_label, 
    fill = estimate), colour = "white", linewidth = 0.3) + geom_text(data = trajectory_heatmap_df, aes(Column, 
    Pathway_label, label = if_else(p_adj_BH < 0.05, "*", "")), size = 3) + scale_fill_gradient2(low = "#2166AC", 
    mid = "white", high = "#B2182B", midpoint = 0, limits = c(-heat_limit, heat_limit)) + scale_x_discrete(drop = FALSE) + 
    labs(x = NULL, y = NULL, fill = "Adjusted effect\n(standardized beta)") + theme_bw(base_size = 8) + 
    theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 7.5, 
        colour = "grey30"), panel.grid = element_blank(), axis.text.x = element_text(angle = 32, hjust = 1, 
        size = 7.5), axis.text.y = element_text(size = 6.3), legend.position = "right")

save_plot(p_trajectory_heatmap, "Fig_D_30_global_trajectory_pathways", 5.5, 5.3)

combined_figure <- (p_day1 | p_baseline_dynamic)/p_representative + plot_layout(heights = c(1, 1.08))

save_plot(combined_figure, "Figure_GOBP_baseline_and_dynamic_summary", 16, 13)

write.csv(day1_plot_df, repo_output(repo_path(out_dir, "Fig_A_plot_data.csv"), "Figure2A_Figure3A_FigureS5A_plot.R"), 
    row.names = FALSE)

write.csv(baseline_dynamic_df, repo_output(repo_path(out_dir, "Fig_B_plot_data.csv"), "Figure2A_Figure3A_FigureS5A_plot.R"), 
    row.names = FALSE)

write.csv(representative_df, repo_output(repo_path(out_dir, "Fig_C_plot_data.csv"), "Figure2A_Figure3A_FigureS5A_plot.R"), 
    row.names = FALSE)

write.csv(trajectory_heatmap_df, repo_output(repo_path(out_dir, "Fig_D_effect_plot_data.csv"), "Figure2A_Figure3A_FigureS5A_plot.R"), 
    row.names = FALSE)

n_pathways <- nrow(day1_results)

n_d1 <- sum(day1_results$p_adj_BH < 0.05, na.rm = TRUE)

n_d1_up <- sum(day1_results$p_adj_BH < 0.05 & day1_results$estimate > 0, na.rm = TRUE)

n_d1_down <- sum(day1_results$p_adj_BH < 0.05 & day1_results$estimate < 0, na.rm = TRUE)

n_trajectory <- sum(trajectory_results$p_adj_BH < 0.05, na.rm = TRUE)

n_d4_interaction <- sum(d4_results$p_adj_BH < 0.05, na.rm = TRUE)

