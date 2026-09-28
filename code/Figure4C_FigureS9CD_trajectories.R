# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
if (Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX")) {
    try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE)
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(lmerTest)
    library(readxl)
})

root <- "."

source(repo_input(repo_path(root, "code/shared_CMAISE_cohort.R")))

shared <- prepare_cmaise_cohort(root)

scores <- shared$scores

module_labels <- shared$module_labels

out_dir <- repo_path(root, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(shared$audit, repo_output(repo_path(out_dir, "260911_CMAISE_input_audit.csv"), "Figure4C_FigureS9CD_trajectories.R"), 
    row.names = FALSE)

wilcox_one <- function(d) {
    wt <- wilcox.test(Score ~ Outcome, data = d, exact = FALSE)
    tibble(n_survival = sum(d$Outcome == "Survival"), n_mortality = sum(d$Outcome == "Mortality"), statistic = unname(wt$statistic), 
        p_value = wt$p.value)
}

overall_tests <- scores %>% group_by(Module_ID, Timepoint) %>% group_modify(~wilcox_one(.x)) %>% ungroup() %>% 
    mutate(p_adj_BH = p.adjust(p_value, "BH"))

subgroup_tests <- scores %>% group_by(Subgroup, Module_ID, Timepoint) %>% group_modify(~wilcox_one(.x)) %>% 
    ungroup() %>% mutate(p_adj_BH = p.adjust(p_value, "BH"))

format_p <- function(x) ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))

format_tp <- function(x) ifelse(x < 0.001, "<0.001", paste0("=", formatC(x, format = "f", digits = 3)))

label_tests <- function(d) d %>% mutate(label = if_else(!is.na(p_value) & p_value < 0.05, paste0("P", 
    format_tp(p_value), "\nFDR", format_tp(p_adj_BH)), "ns"), hjust = case_when(Timepoint == "D1" ~ 0, 
    Timepoint == "D5" ~ 1, TRUE ~ 0.5))

outcome_colors <- c(Survival = "#4575B4", Mortality = "#D73027")

fit_interaction <- function(module_id) {
    d <- filter(dat, Module_ID == module_id)
    fit <- lmerTest::lmer(Score ~ Outcome * Timepoint + (1 | HumanID), data = d, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa", 
        optCtrl = list(maxfun = 1e+05)))
    cn <- names(lme4::fixef(fit))
    interaction_columns <- grep("^OutcomeMortality:Timepoint(D3|D5)$", cn)
    if (length(interaction_columns) != 2L) {
        stop("Could not identify both interaction coefficients for ", module_id)
    }
    L <- matrix(0, nrow = 2, ncol = length(cn), dimnames = list(NULL, cn))
    L[cbind(seq_len(2), interaction_columns)] <- 1
    joint <- lmerTest::contestMD(fit, L, ddf = "Satterthwaite")
    tibble(Module_ID = module_id, NumDF = joint[["NumDF"]], DenDF = joint[["DenDF"]], F_value = joint[["F value"]], 
        interaction_p = joint[["Pr(>F)"]], n_samples = nrow(d), n_patients = n_distinct(d$HumanID), singular_fit = lme4::isSingular(fit))
}

make_plot <- function(module_id) {
    s <- filter(summary_df, Module_ID == module_id)
    ann <- filter(trajectory_tests, Module_ID == module_id)
    tp_ann <- filter(timepoint_tests, Module_ID == module_id)
    observed_range <- diff(range(c(s$mean - s$SE, s$mean + s$SE), na.rm = TRUE))
    if (!is.finite(observed_range) || observed_range == 0) 
        observed_range <- 1
    timepoint_y <- max(s$mean + s$SE, na.rm = TRUE) + 0.08 * observed_range
    trajectory_y <- timepoint_y + 0.25 * observed_range
    ggplot(s, aes(Timepoint, mean, color = Outcome, group = Outcome)) + annotate("rect", xmin = 1, xmax = 2, 
        ymin = -Inf, ymax = Inf, fill = "grey92", color = NA) + geom_line(linewidth = 1.05) + geom_point(size = 2.7) + 
        geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.12, linewidth = 0.55) + geom_text(data = tp_ann, 
        aes(x = Timepoint, y = timepoint_y, label = label, hjust = hjust), inherit.aes = FALSE, vjust = 0, 
        size = 2.75, lineheight = 0.92) + geom_label(data = ann, aes(x = 1, y = trajectory_y, label = label), 
        inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3.1, fontface = "bold", linewidth = 0, fill = scales::alpha("white", 
            0.82)) + scale_color_manual(values = outcome_colors, labels=c(Survival="No event",Mortality="Death <=28 d"), name = NULL, drop = FALSE) + scale_x_discrete(expand = expansion(add = c(0.15, 
        0.15))) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) + labs(title = unname(module_labels[module_id]), 
        x = NULL, y = "GSVA score (mean +/- SE)") + theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", 
        size = 12), panel.grid.minor = element_blank(), legend.position = "top", legend.title=element_text(size=8),legend.text=element_text(size=8), axis.text.x = element_text(face = "bold"))
}

cohorts <- c(overall = "Overall CMAISE", lung = "lung infection", nonlung = "non-lung infection")

for (tag in names(cohorts)) {
    dat <- if (tag == "overall") 
        scores
    else filter(scores, Subgroup == cohorts[[tag]])
    dat <- mutate(dat, HumanID = factor(HumanID))
    timepoint_tests <- label_tests(if (tag == "overall") 
        overall_tests
    else filter(subgroup_tests, Subgroup == cohorts[[tag]]))
    trajectory_tests <- map_dfr(names(module_labels), fit_interaction) %>% mutate(interaction_fdr = p.adjust(interaction_p, 
        "BH"), Module = unname(module_labels[Module_ID]), label = paste0("Trajectory interaction FDR ", 
        format_p(interaction_fdr)))
    summary_df <- dat %>% group_by(Module_ID, Timepoint, Outcome) %>% summarise(n = n(), mean = mean(Score), 
        SD = sd(Score), SE = SD/sqrt(n), .groups = "drop")
    plots <- setNames(map(names(module_labels), make_plot), names(module_labels))
    prefix <- repo_path(out_dir, paste0("260911_CMAISE_", tag))
    pdf(repo_output(paste0(prefix, "_trajectory.pdf"), "Figure4C_FigureS9CD_trajectories.R"), width = 3.1, height = 2.9, 
        onefile = TRUE)
    walk(plots, print)
    dev.off()
    write.csv(timepoint_tests, repo_output(paste0(prefix, "_timepoint_tests.csv"), "Figure4C_FigureS9CD_trajectories.R"), 
        row.names = FALSE)
    write.csv(trajectory_tests, repo_output(paste0(prefix, "_interaction_tests.csv"), "Figure4C_FigureS9CD_trajectories.R"), 
        row.names = FALSE)
    write.csv(summary_df, repo_output(paste0(prefix, "_mean_SE.csv"), "Figure4C_FigureS9CD_trajectories.R"), row.names = FALSE)
}

write.csv(scores, repo_output(repo_path(out_dir, "260911_CMAISE_scores_long.csv"), "Figure4C_FigureS9CD_trajectories.R"), 
    row.names = FALSE)

writeLines(c("Shared cohort: shared_CMAISE_cohort.R; endpoint 28-day mortality", paste("GSVA input:", 
    shared$score_file), paste("Clinical input:", shared$clinical_file), "GSVA scores reused; mean/SE, Wilcoxon tests and interaction models recomputed.", 
    "Timepoint P/FDR values recalculated on shared cohort; original label format preserved.", "Overall timepoint BH-FDR: 12 tests; subgroup timepoint BH-FDR: 24 tests together.", 
    "Timepoint labels: raw P < 0.05 shows P and FDR; otherwise ns (original rule).", "Interaction: Score ~ Outcome * categorical Timepoint + (1 | HumanID), ML.", 
    "Joint 2-df Satterthwaite F test of D3/D5 interactions; BH across 4 modules per cohort.", "PDF: one module per page, 3.1 x 2.9 inches, matching Figure3C_plot.R.", 
    "Extra vertical space and edge alignment retain two-line P/FDR labels."), repo_output(repo_path(out_dir, 
    "260911_run_info.txt"), "Figure4C_FigureS9CD_trajectories.R"))

capture.output(sessionInfo(), file = repo_path(out_dir, "260911_sessionInfo.txt"))

message("Done. Results written to: ", out_dir)

