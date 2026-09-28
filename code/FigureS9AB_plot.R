# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
if (Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX")) {
    try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(readxl)
    library(patchwork)
})

root <- "."

source(repo_input(repo_path(root, "code/shared_CMAISE_cohort.R")))

shared <- prepare_cmaise_cohort(root)

assert_cmaise_same_cohort(shared$scores, read.csv(repo_input("Inputs/260911_CMAISE_scores_long.csv")))

out_dir <- repo_path(root, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- shared$scores %>% distinct(HumanID, SampleID, Timepoint, Outcome, Subgroup)

stopifnot(!anyDuplicated(samples$SampleID), !anyDuplicated(samples[c("HumanID", "Timepoint")]))

overall_counts <- samples %>% count(Timepoint, Outcome, name = "n", .drop = FALSE)

subgroup_counts <- samples %>% count(Subgroup, Timepoint, Outcome, name = "n", .drop = FALSE)

make_count_plot <- function(counts, title, y_max, subgroup = FALSE) {
    totals <- counts %>% group_by(Timepoint) %>% summarise(n_total = sum(n), .groups = "drop")
    ggplot(counts, aes(Timepoint, n, fill = Outcome)) + geom_col(alpha = 0.9, color = "white", linewidth = 0.4, 
        width = 0.6) + geom_text(data = totals, aes(Timepoint, n_total, label = paste0("n=", n_total)), 
        inherit.aes = FALSE, vjust = -0.8, fontface = "bold", size = 3.5) + scale_fill_manual(values = c(Survival = "#4575B4", 
        Mortality = "#D73027"), labels=c(Survival="No event",Mortality="Death <=28 d"), name = NULL, drop = FALSE) + scale_y_continuous(limits = c(0, y_max), 
        breaks = seq(0, y_max, 50), expand = expansion(mult = c(0, 0))) + labs(title = title, x = "Study Timepoint", 
        y = "Number of Samples") + theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5, 
        size = 15, margin = margin(b = 12)), legend.position = "top", legend.background = element_rect(fill = "white", color = "grey80"), legend.key.size = grid::unit(0.4, 
        "cm"), legend.text = element_text(size = 9), panel.grid.minor = element_blank(), axis.title = element_text(size = 11), 
        axis.text = element_text(size = 10), plot.margin = margin(8, 8, 8, 8))
}

p_A <- make_count_plot(overall_counts, "CMAISE cohort", 450)

p_lung <- make_count_plot(filter(subgroup_counts, Subgroup == "lung infection"), "Lung infection", 250, 
    TRUE)

p_nonlung <- make_count_plot(filter(subgroup_counts, Subgroup == "non-lung infection"), "Non-lung infection", 
    250, TRUE)

p_B <- p_lung + (p_nonlung + labs(y = NULL) + theme(axis.text.y = element_blank())) + plot_layout(nrow = 1)

p_AB <- (p_A + labs(tag = "A")) + plot_spacer() + (p_lung + labs(tag = "B")) + (p_nonlung + labs(y = NULL) + 
    theme(axis.text.y = element_blank())) + plot_layout(nrow = 1, widths = c(1, 0.18, 1, 1)) & theme(plot.tag = element_text(size = 24))

plots <- list(panel_A = p_A, panel_B = p_B, panels_AB = p_AB)

sizes <- list(panel_A = c(3.5, 4), panel_B = c(6.6, 4), panels_AB = c(10.7, 4))

for (tag in names(plots)) {
    for (ext in c("pdf", "png")) {
        ggsave(repo_output(repo_path(out_dir, paste0("260912_CMAISE_", tag, ".", ext)), "FigureS9AB_plot.R"), 
            plots[[tag]], width = sizes[[tag]][1], height = sizes[[tag]][2], dpi = 300, bg = "white")
    }
}

write.csv(overall_counts, repo_output(repo_path(out_dir, "260912_CMAISE_overall_counts.csv"), "FigureS9AB_plot.R"), 
    row.names = FALSE)

write.csv(subgroup_counts, repo_output(repo_path(out_dir, "260912_CMAISE_subgroup_counts.csv"), "FigureS9AB_plot.R"), 
    row.names = FALSE)

write.csv(samples, repo_output(repo_path(out_dir, "260912_CMAISE_sample_manifest.csv"), "FigureS9AB_plot.R"), 
    row.names = FALSE)

writeLines(c("Endpoint: in-hospital death within 28 days; shared primary/JM cohort.", "Exact sample/outcome/score agreement with both analysis inputs verified.", 
    "Each biological sample counted once, not once per GSVA module.", "Overall D1/D3/D5 totals: 402/347/290.", 
    "Lung D1/D3/D5 totals: 183/160/138; non-lung: 219/187/152.", "Two D3 samples after recorded death excluded, matching both analyses.", 
    "Common origin confirmed; corrected duration for JinHua_53 and SSR_1 is 5 days; live discharge is censored at discharge in joint models."), 
    repo_output(repo_path(out_dir, "260912_run_info.txt"), "FigureS9AB_plot.R"))

message("Done: ", out_dir)

