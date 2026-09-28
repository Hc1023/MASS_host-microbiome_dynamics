# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggpubr)
})

root <- "."

out_dir <- repo_path(root, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

e <- new.env()

load(repo_input(repo_path(root, "Inputs/1211_metadata.rdata")), envir = e)

load(repo_input(repo_path(root, "Inputs/1616_microbe.rdata")), envir = e)

meta <- e$df_long %>% filter(Timepoint == "D1") %>% mutate(Outcome = factor(Mortality28d, levels = c("0", 
    "1"), labels = c("Survival", "Mortality")))

stopifnot(!anyDuplicated(meta$SampleID), !anyNA(meta$Outcome), all(meta$SampleID %in% colnames(e$data)))

x <- as.matrix(e$data[, meta$SampleID, drop = FALSE])

virus_rows <- e$mapid_vec[rownames(x)] == "Viruses"

stopifnot(!anyNA(virus_rows), all(c("HHV-4", "HCMV") %in% rownames(x)), all(is.finite(x)), all(x >= 0))

dat <- meta %>% transmute(HumanID, SampleID, Outcome, Total = colSums(x), Viruses = colSums(x[virus_rows, 
    , drop = FALSE]), EBV = as.numeric(x["HHV-4", ]), HCMV = as.numeric(x["HCMV", ])) %>% pivot_longer(c(Total, 
    Viruses, EBV, HCMV), names_to = "Feature", values_to = "Mass") %>% mutate(Feature = factor(Feature, 
    levels = c("Total", "Viruses", "EBV", "HCMV")), Mass_log2 = log2(Mass + 1), Positive = Mass > 0)

colors <- c(Survival = "#4575B4", Mortality = "#D73027")

common_theme <- theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank(), legend.position = "top", 
    plot.margin = grid::unit(rep(0.2, 4), "cm"))

format_p <- function(p) ifelse(is.na(p), "Not tested", ifelse(p < 1e-04, "p < 0.0001", paste0("p = ", 
    format.pval(p, digits = 2, eps = 1e-04))))

wilcox_tests <- function(d) d %>% group_by(Feature) %>% group_modify(~{
    ns <- sum(.x$Outcome == "Survival")
    nm <- sum(.x$Outcome == "Mortality")
    p <- if (ns > 0 && nm > 0) 
        wilcox.test(Mass_log2 ~ Outcome, data = .x, exact = FALSE)$p.value
    else NA_real_
    tibble(n_survival = ns, n_mortality = nm, p_value = p)
}) %>% ungroup() %>% mutate(FDR = p.adjust(p_value, "BH"))

all_tests <- wilcox_tests(dat)

pos <- filter(dat, Positive)

pos_tests <- wilcox_tests(pos)

make_abundance <- function(d, tests, positive_only = FALSE) {
    yr <- range(d$Mass_log2)
    span <- max(diff(yr), 1)
    ann <- tests %>% mutate(y = max(d$Mass_log2) + 0.07 * span, label = format_p(p_value))
    p <- ggplot(d, aes(Feature, Mass_log2, fill = Outcome)) + geom_boxplot(aes(color = Outcome), outlier.shape = NA, 
        alpha = 0.75, width = 0.55, position = position_dodge(width = 0.72)) + geom_point(aes(color = Outcome), 
        position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72, seed = 20260914), alpha = 0.65, 
        size = 1.4)
    if (!positive_only) 
        p <- p + stat_compare_means(aes(group = Outcome), method = "wilcox.test", label = "p", hide.ns = TRUE, 
            color = "#c43932", size = 4)
    else p <- p + geom_text(data = ann, aes(Feature, y, label = label), inherit.aes = FALSE, color = "#c43932", 
        size = 4)
    p + scale_color_manual(values = colors, name = NULL) + scale_fill_manual(values = colors, name = NULL) + 
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) + labs(x = NULL, y = "log2(mass + 1)", 
        subtitle = if (positive_only) 
            "Positive samples only (mass > 0)"
        else NULL) + common_theme
}

prev <- dat %>% group_by(Feature, Outcome) %>% summarise(n = n(), n_positive = sum(Positive), prevalence = n_positive/n, 
    .groups = "drop")

prev_tests <- dat %>% group_by(Feature) %>% group_modify(~{
    tab <- table(factor(.x$Outcome, levels = names(colors)), factor(.x$Positive, levels = c(FALSE, TRUE)))
    tibble(p_value = fisher.test(tab)$p.value)
}) %>% ungroup() %>% mutate(FDR = p.adjust(p_value, "BH"), label = format_p(p_value))

prev_ann <- prev %>% group_by(Feature) %>% summarise(y = 1.2, .groups = "drop") %>% left_join(prev_tests, 
    by = "Feature")

p_prev <- ggplot(prev, aes(Feature, prevalence, fill = Outcome)) + geom_col(position = position_dodge(0.72), 
    width = 0.55) + geom_text(aes(x = as.numeric(Feature) + if_else(Outcome == "Survival", -0.255, 0.255), 
    label = sprintf("%d/%d\n(%.1f%%)", n_positive, n, 100 * prevalence)), vjust = -0.3, size = 3) + geom_text(data = prev_ann, 
    aes(Feature, y, label = label), inherit.aes = FALSE, color = "#c43932", size = 4) + scale_fill_manual(values = colors, 
    name = NULL) + scale_y_continuous(labels = c("0", "25", "50", "75", "100"), breaks = seq(0, 1, 0.25), 
    limits = c(0, 1.3), expand = expansion(mult = c(0, 0))) + labs(x = NULL, y = "Prevalence (mass > 0)") + 
    common_theme

plots <- list(abundance_all = make_abundance(dat, all_tests), prevalence = p_prev, abundance_positive_only = make_abundance(pos, 
    pos_tests, TRUE))

for (tag in names(plots)) ggsave(repo_output(repo_path(out_dir, paste0("260914_D1_", tag, ".pdf")), "Figure2B_analysis.R"), 
    plots[[tag]], width = 5, height = 3.3)

tag = names(plots)[2]

ggsave(repo_output(repo_path(out_dir, paste0("260914_D1_", tag, ".pdf")), "Figure2B_analysis.R"), 
    plots[[tag]], width = 4.9, height = 3)

write_csv(dat, repo_output(repo_path(out_dir, "260914_D1_sample_data.csv"), "Figure2B_analysis.R"))

write_csv(prev, repo_output(repo_path(out_dir, "260914_D1_prevalence_counts.csv"), "Figure2B_analysis.R"))

write_csv(all_tests, repo_output(repo_path(out_dir, "260914_D1_abundance_all_tests.csv"), "Figure2B_analysis.R"))

write_csv(pos_tests, repo_output(repo_path(out_dir, "260914_D1_abundance_positive_tests.csv"), "Figure2B_analysis.R"))

write_csv(prev_tests, repo_output(repo_path(out_dir, "260914_D1_prevalence_tests.csv"), "Figure2B_analysis.R"))

writeLines(c("Source: 2608_ajrccm/260621_biomass_comparison.R; Day 1, 28-day outcome.", "Total = sum of all taxa; Viruses = sum of taxa mapped to Viruses; EBV = HHV-4.", 
    "EBV/HCMV are retained in Total and Viruses, matching the original biomass analysis.", "Positive means mass > 0; no additional detection threshold is applied.", 
    "Prevalence denominators: all Day-1 samples in each outcome group.", "Positive-only abundance: independently restrict each feature to mass > 0; log2(mass + 1).", 
    "Abundance: two-sided Wilcoxon rank-sum; prevalence: two-sided Fisher exact test.", "Plots show raw P values. CSVs also include BH FDR across four features separately for each analysis.", 
    "Each PDF is 5 x 3.3 inches; boxes show median/IQR with default 1.5-IQR whiskers; points are samples."), 
    repo_output(repo_path(out_dir, "260914_run_info.txt"), "Figure2B_analysis.R"))

message("Done: ", out_dir)

