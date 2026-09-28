# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(tidyverse)
})

metadata = read.csv(file = repo_input("Inputs/metadata_table_source_1007.csv"))

out_dir <- repo_path(".")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_df <- metadata %>% transmute(HumanID, StudyCohort417, Mortality28d = factor(Mortality28d, levels = c("0", 
    "1")), MortalityGroup = factor(Mortality28d, levels = c("0", "1"), labels = c("Survival", "Mortality")), 
    LymphocyteCount = suppressWarnings(as.numeric(LymphocyteCount))) %>% filter(!is.na(Mortality28d), 
    !is.na(LymphocyteCount))

plot_lymphocyte_boxplot <- function(df, title) {
    wilcox_res <- wilcox.test(LymphocyteCount ~ Mortality28d, data = df, exact = FALSE)
    p_label <- paste0("Wilcoxon p = ", format.pval(wilcox_res$p.value, digits = 2, eps = 0.001))
    medians <- df %>% group_by(MortalityGroup) %>% summarise(value = median(LymphocyteCount), .groups = "drop")
    y_upper <- quantile(df$LymphocyteCount, 0.99, na.rm = TRUE)
    ggplot(df, aes(x = MortalityGroup, y = LymphocyteCount, fill = Mortality28d, color = Mortality28d)) + 
        stat_boxplot(geom = "errorbar", width = 0.18, linewidth = 0.55) + geom_boxplot(width = 0.6, linewidth = 0.55, 
        outlier.shape = NA) + geom_jitter(position = position_jitter(width = 0.12, height = 0, seed = 260915), 
        size = 1.15, alpha = 0.1) + annotate("text", x = 1.5, y = y_upper * 0.96, label = p_label, fontface = "bold", 
        size = 2.8) + geom_text(data = medians, aes(x = MortalityGroup, y = y_upper * 0.78, label = sprintf("Median\n%.2f", 
        value)), inherit.aes = FALSE, size = 2.8) + scale_fill_manual(values = c(`0` = alpha("#4575B4", 
        0.5), `1` = alpha("#D73027", 0.5)), labels = c(`0` = "Survival", `1` = "Mortality"), guide = "none") + 
        scale_color_manual(values = c(`0` = "#4575B4", `1` = "#D73027"), labels = c(`0` = "Survival", 
            `1` = "Mortality"), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.04, 
        0.12))) + coord_cartesian(ylim = c(0, y_upper)) + labs(x = NULL, y = "Lymphocyte\ncount (10^9/L)") + 
        theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12), panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 10), axis.title = element_text(size = 11))
}

study_df <- plot_df %>% filter(StudyCohort417 == 1)

stopifnot(!anyDuplicated(study_df$HumanID), nrow(study_df) > 0)

p_this_study <- plot_lymphocyte_boxplot(study_df, title = "This study (n=417)")

ggsave(repo_output(repo_path(out_dir, "260915_lymphocyte_count_mortality.pdf"), "FigureS3B_analysis.R"), 
    p_this_study, width = 2.5, height = 1.8)

ggsave(repo_output(repo_path(out_dir, "260915_lymphocyte_count_mortality.png"), "FigureS3B_analysis.R"), 
    p_this_study, width = 2.2, height = 2.8, dpi = 300)

w <- wilcox.test(LymphocyteCount ~ Mortality28d, data = study_df, exact = FALSE)

write_csv(tibble(statistic = unname(w$statistic), p_value = w$p.value, median_survival = median(study_df$LymphocyteCount[study_df$MortalityGroup == 
    "Survival"]), median_mortality = median(study_df$LymphocyteCount[study_df$MortalityGroup == "Mortality"]), 
    n = nrow(study_df)), repo_output(repo_path(out_dir, "260915_wilcoxon_test.csv"), "FigureS3B_analysis.R"))

write_csv(study_df, repo_output(repo_path(out_dir, "260915_lymphocyte_data.csv"), "FigureS3B_analysis.R"))

write_csv(study_df %>% group_by(MortalityGroup) %>% summarise(n = n(), median = median(LymphocyteCount), 
    Q1 = quantile(LymphocyteCount, 0.25), Q3 = quantile(LymphocyteCount, 0.75), .groups = "drop"), repo_output(repo_path(out_dir, 
    "260915_lymphocyte_summary.csv"), "FigureS3B_analysis.R"))

writeLines(c("Source: 2608_ajrccm/260622_lymphocyte.R; manuscript cohort StudyCohort417==1.", "Outcome: 28-day mortality; missing outcome/count excluded.", 
    "Group medians are annotated; log2FC is not used.", "Two-sided Wilcoxon rank-sum, exact=FALSE; all available counts used.", 
    "Display zoomed to the 99th percentile, as in the original; no observations removed from test or box calculations.", 
    "Deterministic horizontal jitter; PDF dimensions 2.2 x 2 inches.", paste("Available patients:", nrow(study_df)), 
    paste("Wilcoxon P:", w$p.value)), repo_output(repo_path(out_dir, "260915_run_info.txt"), "FigureS3B_analysis.R"))

message("Done: ", out_dir)

