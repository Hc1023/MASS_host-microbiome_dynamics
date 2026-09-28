# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(tidyverse)
    library(lmerTest)
})

root <- "."

input_file <- "Outputs/Figure3C_prepare_MASS_GSVA_scores_long_with_metadata.csv"

out_dir <- repo_path(root, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(repo_exists(input_file))

old_to_id <- c(`PRR/TLR/TNF/NF-kB signaling` = "M1", `Phagolysosome/autophagy` = "M2", `Type I interferon response` = "M3", 
    `Ribosome biogenesis` = "M4")

module_labels <- c(M1 = "PRR/TLR/NF-kB signaling", M2 = "Phagolysosomal/autophagy", M3 = "Type I IFN response", 
    M4 = "Ribosome biogenesis")

dat <- read.csv(repo_input(input_file), check.names = FALSE, stringsAsFactors = FALSE) %>% transmute(HumanID = factor(HumanID), 
    SampleID, Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7")), Outcome = factor(as.character(Mortality28d), 
        levels = c("0", "1"), labels = c("Survival", "Mortality")), Module_ID = unname(old_to_id[Module]), 
    Score = as.numeric(Score)) %>% filter(!is.na(Timepoint), !is.na(Outcome), !is.na(Module_ID), !is.na(Score)) %>% 
    mutate(Module_ID = factor(Module_ID, levels = names(module_labels)))

stopifnot(!anyDuplicated(dat[c("SampleID", "Module_ID")]), nlevels(droplevels(dat$Module_ID)) == 4L)

fit_interaction <- function(module_id) {
    d <- filter(dat, Module_ID == module_id)
    fit <- lmerTest::lmer(Score ~ Outcome * Timepoint + (1 | HumanID), data = d, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa", 
        optCtrl = list(maxfun = 1e+05)))
    cn <- names(lme4::fixef(fit))
    interaction_columns <- grep("^OutcomeMortality:Timepoint(D4|D7)$", cn)
    if (length(interaction_columns) != 2L) {
        stop("Could not identify both interaction coefficients for ", module_id)
    }
    L <- matrix(0, nrow = 2, ncol = length(cn), dimnames = list(NULL, cn))
    L[cbind(seq_len(2), interaction_columns)] <- 1
    joint <- lmerTest::contestMD(fit, L, ddf = "Satterthwaite")
    tibble(Module_ID = module_id, NumDF = joint[["NumDF"]], DenDF = joint[["DenDF"]], F_value = joint[["F value"]], 
        interaction_p = joint[["Pr(>F)"]], n_samples = nrow(d), n_patients = n_distinct(d$HumanID), singular_fit = lme4::isSingular(fit))
}

trajectory_tests <- map_dfr(names(module_labels), fit_interaction) %>% mutate(interaction_fdr = p.adjust(interaction_p, 
    method = "BH"), Module = unname(module_labels[Module_ID]))

format_p <- function(x) {
    if_else(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
}

trajectory_tests <- trajectory_tests %>% mutate(significance = case_when(interaction_fdr < 0.001 ~ "***", 
    interaction_fdr < 0.01 ~ "**", interaction_fdr < 0.05 ~ "*", TRUE ~ "ns"), label = paste0("Trajectory interaction FDR ", 
    format_p(interaction_fdr)))

summary_df <- dat %>% group_by(Module_ID, Timepoint, Outcome) %>% summarise(n = n(), mean = mean(Score), 
    SD = sd(Score), SE = SD/sqrt(n), .groups = "drop")

timepoint_tests <- dat %>% group_by(Module_ID, Timepoint) %>% summarise(n_survival = sum(Outcome == "Survival"), 
    n_mortality = sum(Outcome == "Mortality"), statistic = unname(wilcox.test(Score ~ Outcome, exact = FALSE)$statistic), 
    p_value = wilcox.test(Score ~ Outcome, exact = FALSE)$p.value, .groups = "drop") %>% mutate(p_adj_BH = p.adjust(p_value, 
    method = "BH"), significance = case_when(p_adj_BH < 0.001 ~ "***", p_adj_BH < 0.01 ~ "**", p_adj_BH < 
    0.05 ~ "*", TRUE ~ "ns"), label = significance)

outcome_colors <- c(Survival = "#4575B4", Mortality = "#D73027")

make_plot <- function(module_id) {
    s <- filter(summary_df, Module_ID == module_id)
    ann <- filter(trajectory_tests, Module_ID == module_id)
    tp_ann <- filter(timepoint_tests, Module_ID == module_id)
    observed_range <- diff(range(c(s$mean - s$SE, s$mean + s$SE), na.rm = TRUE))
    if (!is.finite(observed_range) || observed_range == 0) 
        observed_range <- 1
    timepoint_y <- max(s$mean + s$SE, na.rm = TRUE) + 0.08 * observed_range
    trajectory_y <- timepoint_y + 0.14 * observed_range
    ggplot(s, aes(Timepoint, mean, color = Outcome, group = Outcome)) + annotate("rect", xmin = 1, xmax = 2, 
        ymin = -Inf, ymax = Inf, fill = "grey92", color = NA) + geom_line(linewidth = 1.05) + geom_point(size = 2.7) + 
        geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.12, linewidth = 0.55) + geom_text(data = tp_ann, 
        aes(x = Timepoint, y = timepoint_y, label = label), inherit.aes = FALSE, hjust = 0.5, vjust = 0, 
        size = 2.75, lineheight = 0.92) + geom_label(data = ann, aes(x = 1, y = trajectory_y, label = label), 
        inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3.1, fontface = "bold", linewidth = 0, fill = scales::alpha("white", 
            0.82)) + scale_color_manual(values = outcome_colors, name = NULL, drop = FALSE) + scale_x_discrete(expand = expansion(add = c(0.15, 
        0.15))) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) + labs(title = unname(module_labels[module_id]), 
        x = NULL, y = "GSVA score (mean +/- SE)") + theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", 
        size = 12), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_text(face = "bold"))
}

plots <- setNames(map(names(module_labels), make_plot), names(module_labels))

pdf(repo_output(repo_path(out_dir, "260907_four_modules_trajectory.pdf"), "Figure3C_plot.R"), width = 3.1, 
    height = 2.9, onefile = TRUE)

walk(plots, print)

dev.off()

write.csv(trajectory_tests %>% select(Module_ID, Module, NumDF, DenDF, F_value, interaction_p, interaction_fdr, 
    n_samples, n_patients, singular_fit), repo_output(repo_path(out_dir, "260907_trajectory_interaction_tests.csv"), 
    "Figure3C_plot.R"), row.names = FALSE)

write.csv(summary_df %>% mutate(Module = unname(module_labels[as.character(Module_ID)])) %>% select(Module_ID, 
    Module, Timepoint, Outcome, n, mean, SD, SE), repo_output(repo_path(out_dir, "260907_GSVA_trajectory_mean_SE.csv"), 
    "Figure3C_plot.R"), row.names = FALSE)

write.csv(timepoint_tests %>% mutate(Module = unname(module_labels[as.character(Module_ID)])) %>% select(Module_ID, 
    Module, Timepoint, n_survival, n_mortality, statistic, p_value, p_adj_BH, significance), repo_output(repo_path(out_dir, 
    "260907_GSVA_wilcoxon_D1_D4_D7.csv"), "Figure3C_plot.R"), row.names = FALSE)

writeLines(c(paste("GSVA input:", input_file), "Plot: observed mean +/- SE; one module per PDF page", 
    "Model: Score ~ Outcome * categorical Timepoint + (1 | HumanID)", "Interaction P: joint 2-df Satterthwaite F test of D4 and D7 interactions", 
    "Interaction FDR: Benjamini-Hochberg adjustment across the four modules", "Timepoint tests: two-sided Wilcoxon rank-sum, Survival versus Mortality", 
    "Timepoint FDR: Benjamini-Hochberg adjustment across all 12 comparisons"), repo_output(repo_path(out_dir, 
    "260907_run_info.txt"), "Figure3C_plot.R"))

message("Done. Results written to: ", out_dir)

