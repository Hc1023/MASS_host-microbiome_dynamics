# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(tidyverse)
    library(lmerTest)
    library(patchwork)
})

root <- "."

input_dir <- repo_path(root, "Outputs")

out_dir <- repo_path(root, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "NK cells", "Monocytes", "Neutrophils")

model_data <- read.csv(repo_input(repo_path(input_dir, "Inputs/composition_adjusted_model_data.csv")), 
    check.names = FALSE)

cache <- readRDS(repo_input(repo_path(input_dir, "Inputs/MASS_D1_D4_D7_xCell_FPKM_scores.rds")))

stopifnot(!anyDuplicated(model_data$SampleID), all(cells %in% rownames(cache)), all(model_data$SampleID %in% 
    colnames(cache)))

for (cell in cells) stopifnot(isTRUE(all.equal(as.numeric(model_data[[cell]]), as.numeric(cache[cell, 
    model_data$SampleID]), tolerance = 1e-12)))

dat <- model_data %>% select(HumanID, SampleID, Timepoint, Mortality28d, all_of(cells)) %>% pivot_longer(all_of(cells), 
    names_to = "Cell_type", values_to = "Score") %>% mutate(HumanID = factor(HumanID), Cell_type = factor(Cell_type, 
    levels = cells), Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7")), Outcome = factor(Mortality28d, 
    levels = c("0", "1"), labels = c("Survival", "Mortality")))

stopifnot(!anyNA(dat), all(is.finite(dat$Score)), !anyDuplicated(dat[c("SampleID", "Cell_type")]))

summary_df <- dat %>% group_by(Cell_type, Timepoint, Outcome) %>% summarise(n = n(), mean = mean(Score), 
    SD = sd(Score), SE = SD/sqrt(n), .groups = "drop") %>% mutate(Lower = pmax(0, mean - SE), Upper = mean + 
    SE)

timepoint_tests <- dat %>% group_by(Cell_type, Timepoint) %>% group_modify(~{
    w <- wilcox.test(Score ~ Outcome, data = .x, exact = FALSE)
    tibble(n_survival = sum(.x$Outcome == "Survival"), n_mortality = sum(.x$Outcome == "Mortality"), 
        statistic = unname(w$statistic), p_value = w$p.value)
}) %>% ungroup() %>% mutate(p_adj_BH = p.adjust(p_value, "BH"))

timepoint_tests <- timepoint_tests %>% mutate(label = case_when(p_adj_BH < 0.001 ~ "***", p_adj_BH < 
    0.01 ~ "**", p_adj_BH < 0.05 ~ "*", TRUE ~ "ns"), hjust = 0.5)

fits <- list()

fit_interaction <- function(cell) {
    d <- filter(dat, Cell_type == cell)
    fit <- lmerTest::lmer(Score ~ Outcome * Timepoint + (1 | HumanID), data = d, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa", 
        optCtrl = list(maxfun = 1e+05)))
    fits[[cell]] <<- fit
    cn <- names(lme4::fixef(fit))
    idx <- grep("^OutcomeMortality:Timepoint(D4|D7)$", cn)
    stopifnot(length(idx) == 2L, nobs(fit) == nrow(d))
    L <- matrix(0, 2, length(cn), dimnames = list(NULL, cn))
    L[cbind(1:2, idx)] <- 1
    test <- lmerTest::contestMD(fit, L, ddf = "Satterthwaite")
    tibble(Cell_type = cell, NumDF = test[["NumDF"]], DenDF = test[["DenDF"]], F_value = test[["F value"]], 
        interaction_p = test[["Pr(>F)"]], n_samples = nrow(d), n_patients = n_distinct(d$HumanID), singular_fit = lme4::isSingular(fit), 
        convergence_message = paste(fit@optinfo$conv$lme4$messages, collapse = "; "))
}

trajectory_tests <- map_dfr(cells, fit_interaction) %>% mutate(interaction_fdr = p.adjust(interaction_p, 
    "BH"), label = paste0("Trajectory interaction FDR ", if_else(interaction_fdr < 0.001, "<0.001", formatC(interaction_fdr, 
    format = "f", digits = 3))))

outcome_colors <- c(Survival = "#4575B4", Mortality = "#D73027")

make_plot <- function(module_id) {
    s <- filter(summary_df, Cell_type == module_id)
    ann <- filter(trajectory_tests, Cell_type == module_id)
    tp_ann <- filter(timepoint_tests, Cell_type == module_id)
    observed_range <- diff(range(c(s$Lower, s$Upper), na.rm = TRUE))
    if (!is.finite(observed_range) || observed_range == 0) 
        observed_range <- 1
    timepoint_y <- max(s$Upper, na.rm = TRUE) + 0.08 * observed_range
    trajectory_y <- timepoint_y + 0.25 * observed_range
    ggplot(s, aes(Timepoint, mean, color = Outcome, group = Outcome)) + annotate("rect", xmin = 1, xmax = 2, 
        ymin = -Inf, ymax = Inf, fill = "grey92", color = NA) + geom_line(linewidth = 1.05) + geom_point(size = 2.7) + 
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.12, linewidth = 0.55) + geom_text(data = tp_ann, 
        aes(x = Timepoint, y = timepoint_y, label = label, hjust = hjust), inherit.aes = FALSE, vjust = 0, 
        size = 2.75, lineheight = 0.92) + geom_label(data = ann, aes(x = 1, y = trajectory_y, label = label), 
        inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3.1, fontface = "bold", linewidth = 0, fill = scales::alpha("white", 
            0.82)) + scale_color_manual(values = outcome_colors, name = NULL, drop = FALSE) + scale_x_discrete(expand = expansion(add = c(0.15, 
        0.15))) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) + labs(title = unname(module_id), 
        x = NULL, y = "xCell score (mean <U+00B1> SE)") + theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", 
        size = 12), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_text(face = "bold"))
}

plots <- setNames(map(cells, make_plot), cells)

pdf(repo_output(repo_path(out_dir, "260915_xcell_trajectory_by_cell.pdf"), "FigureS3F_analysis.R"), width = 3.1, 
    height = 2.9, onefile = TRUE)

walk(plots, print)

dev.off()

combined <- wrap_plots(plots, ncol = 3, guides = "collect", axis_titles = "collect_y") & theme(legend.position = "right")

ggsave(repo_output(repo_path(out_dir, "260915_xcell_trajectories.pdf"), "FigureS3F_analysis.R"), combined, 
    width = 10, height = 4.5)

ggsave(repo_output(repo_path(out_dir, "260915_xcell_trajectories.png"), "FigureS3F_analysis.R"), combined, 
    width = 10.8, height = 5.8, dpi = 300)

write_csv(summary_df, repo_output(repo_path(out_dir, "260915_xcell_mean_SE.csv"), "FigureS3F_analysis.R"))

write_csv(timepoint_tests, repo_output(repo_path(out_dir, "260915_xcell_timepoint_tests.csv"), "FigureS3F_analysis.R"))

write_csv(trajectory_tests, repo_output(repo_path(out_dir, "260915_xcell_interaction_tests.csv"), "FigureS3F_analysis.R"))

write_csv(dat, repo_output(repo_path(out_dir, "260915_xcell_scores_long.csv"), "FigureS3F_analysis.R"))

saveRDS(fits, repo_output(repo_path(out_dir, "260915_xcell_mixed_models.rds"), "FigureS3F_analysis.R"))

writeLines(c("Original plot: 260822_3_modgene_cell_plot.R (p_xcell).", "Existing composition_adjusted_model_data.csv reused, with scores verified against cached xCell RDS.", 
    "Scores are enrichment scores, not cell fractions; no xCell rerun or score transformation.", "Timepoint: two-sided Wilcoxon rank-sum, exact=FALSE; BH across all 18 tests.", 
    "Timepoint labels use BH FDR: *** <0.001; ** <0.01; * <0.05; ns >=0.05.", "Trajectory: Score ~ Outcome * categorical Timepoint + (1 | HumanID), ML.", 
    "Joint 2-df Satterthwaite F test for D4/D7 interaction; BH across six cell types.", "Outcome: 28-day mortality. Original sample set retained; all observations used.", 
    "Error bars: mean +/- SE, lower bound truncated at zero as in the original figure.", "Gray band: D1-D4 interval. Colors: Survival #4575B4; Mortality #D73027."), 
    repo_output(repo_path(out_dir, "260915_run_info.txt"), "FigureS3F_analysis.R"))

capture.output(sessionInfo(), file = repo_path(out_dir, "260915_sessionInfo.txt"))

message("Done: ", out_dir)

