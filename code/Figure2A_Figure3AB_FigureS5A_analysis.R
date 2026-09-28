# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(lme4)

library(lmerTest)

library(emmeans)

library(parallel)

project_dir <- "."

invisible(NULL)

out_dir <- repo_path(project_dir, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

load(repo_input("Inputs/1211_metadata.rdata"))

vars <- c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

metadata_analysis <- df_long %>% left_join(meta %>% dplyr::select(HumanID, all_of(vars)), by = "HumanID")

gsva_file <- repo_path(project_dir, "Outputs", "Inputs/MSigDB_2025.1.Hs_C5_GO_BP_GSVA_all_samples.rds")

gsva_scores <- readRDS(repo_input(gsva_file))

pathways <- rownames(gsva_scores)

analysis_meta <- metadata_analysis %>% filter(Timepoint %in% c("D1", "D4", "D7")) %>% mutate(Timepoint = factor(Timepoint, 
    levels = c("D1", "D4", "D7")), Mortality28d = factor(Mortality28d, c(0, 1), c("Survivor", "Death")), 
    HumanID = factor(HumanID), across(c(Gender, CenterGroup, PneumoniaTypeGroup, Immunosuppression, MV), 
        factor))

model_vars <- c("Mortality28d", "Timepoint", vars, "HumanID")

stopifnot(!anyDuplicated(analysis_meta$SampleID), all(analysis_meta$SampleID %in% colnames(gsva_scores)), 
    !anyNA(analysis_meta[, model_vars]))

gsva_analysis <- gsva_scores[, analysis_meta$SampleID, drop = FALSE]

model_formula <- Score_z ~ Mortality28d * Timepoint + Gender + Age + CenterGroup + PneumoniaTypeGroup + 
    CCI + SOFA_24h + Immunosuppression + MV + (1 | HumanID)

contrast_names <- c("D1 mortality difference", "D4 mortality difference", "D7 mortality difference", 
    "D1 to D4 extra divergence", "D1 to D7 extra divergence")

fit_one_pathway <- function(pathway_name) {
    dat <- analysis_meta
    dat$Score_z <- as.numeric(scale(gsva_analysis[pathway_name, ]))
    fit <- tryCatch(lmerTest::lmer(model_formula, dat, REML = FALSE, control = lmerControl(optimizer = "bobyqa", 
        optCtrl = list(maxfun = 1e+05))), error = identity)
    if (inherits(fit, "error")) 
        return(list(contrasts = tibble(Pathway = pathway_name, Estimand = contrast_names, estimate = NA_real_, 
            SE = NA_real_, df = NA_real_, lower.CL = NA_real_, upper.CL = NA_real_, p.value = NA_real_), 
            omnibus = tibble(Pathway = pathway_name, NumDF = NA_real_, DenDF = NA_real_, F.value = NA_real_, 
                p.value = NA_real_), diagnostics = tibble(Pathway = pathway_name, n_samples = nrow(dat), 
                n_patients = n_distinct(dat$HumanID), singular = NA, convergence_message = conditionMessage(fit)), 
            adjusted_trajectory = tidyr::expand_grid(Pathway = pathway_name, Mortality28d = factor(c("Survivor", 
                "Death"), levels = c("Survivor", "Death")), Timepoint = factor(c("D1", "D4", "D7"), levels = c("D1", 
                "D4", "D7"))) %>% mutate(emmean = NA_real_, SE = NA_real_, df = NA_real_, lower.CL = NA_real_, 
                upper.CL = NA_real_), observed_trajectory = tidyr::expand_grid(Pathway = pathway_name, 
                Mortality28d = factor(c("Survivor", "Death"), levels = c("Survivor", "Death")), Timepoint = factor(c("D1", 
                  "D4", "D7"), levels = c("D1", "D4", "D7"))) %>% mutate(n = NA_integer_, observed_mean = NA_real_, 
                observed_SD = NA_real_)))
    cn <- names(fixef(fit))
    i1 <- match("Mortality28dDeath", cn)
    i4 <- match("Mortality28dDeath:TimepointD4", cn)
    i7 <- match("Mortality28dDeath:TimepointD7", cn)
    if (anyNA(c(i1, i4, i7))) 
        stop("Missing expected coefficient: ", pathway_name)
    L <- matrix(0, 5, length(cn), dimnames = list(contrast_names, cn))
    L[1, i1] <- 1
    L[2, c(i1, i4)] <- 1
    L[3, c(i1, i7)] <- 1
    L[4, i4] <- 1
    L[5, i7] <- 1
    contrast_table <- map_dfr(seq_len(5), function(i) {
        z <- contest1D(fit, L[i, ], ddf = "Satterthwaite", confint = TRUE)
        tibble(Pathway = pathway_name, Estimand = rownames(L)[i], estimate = z[["Estimate"]], SE = z[["Std. Error"]], 
            df = z[["df"]], lower.CL = z[["lower"]], upper.CL = z[["upper"]], p.value = z[["Pr(>|t|)"]])
    })
    joint <- contestMD(fit, L[c(4, 5), , drop = FALSE], ddf = "Satterthwaite")
    omnibus <- tibble(Pathway = pathway_name, NumDF = joint[["NumDF"]], DenDF = joint[["DenDF"]], F.value = joint[["F value"]], 
        p.value = joint[["Pr(>F)"]])
    diagnostics <- tibble(Pathway = pathway_name, n_samples = nobs(fit), n_patients = n_distinct(dat$HumanID), 
        singular = isSingular(fit), convergence_message = paste(fit@optinfo$conv$lme4$messages %||% "", 
            collapse = "; "))
    adjusted_trajectory <- emmeans(fit, ~Mortality28d * Timepoint, lmer.df = "satterthwaite") %>% confint(level = 0.95) %>% 
        as.data.frame() %>% as_tibble() %>% mutate(Pathway = pathway_name, .before = 1)
    observed_trajectory <- dat %>% group_by(Mortality28d, Timepoint) %>% summarise(n = n(), observed_mean = mean(Score_z), 
        observed_SD = sd(Score_z), .groups = "drop") %>% mutate(Pathway = pathway_name, .before = 1)
    list(contrasts = contrast_table, omnibus = omnibus, diagnostics = diagnostics, adjusted_trajectory = adjusted_trajectory, 
        observed_trajectory = observed_trajectory)
}

detected_cores <- detectCores()

if (is.na(detected_cores)) detected_cores <- 4L

n_cores <- min(4L, max(1L, detected_cores - 1L))

message("Fitting ", length(pathways), " GO:BP models using ", n_cores, " cores")

model_results <- mclapply(pathways, fit_one_pathway, mc.cores = n_cores, mc.preschedule = TRUE)

contrast_results <- map_dfr(model_results, "contrasts") %>% group_by(Estimand) %>% mutate(p_adj_BH = p.adjust(p.value, 
    "BH")) %>% ungroup()

trajectory_results <- map_dfr(model_results, "omnibus") %>% mutate(p_adj_BH = p.adjust(p.value, "BH"))

model_diagnostics <- map_dfr(model_results, "diagnostics")

adjusted_group_trajectories <- map_dfr(model_results, "adjusted_trajectory")

observed_group_trajectories <- map_dfr(model_results, "observed_trajectory")

day1_results <- contrast_results %>% filter(Estimand == "D1 mortality difference") %>% arrange(p_adj_BH)

interaction_results <- contrast_results %>% filter(str_detect(Estimand, "extra divergence")) %>% arrange(Estimand, 
    p_adj_BH)

timepoint_results <- contrast_results %>% filter(str_detect(Estimand, "mortality difference")) %>% arrange(Estimand, 
    p_adj_BH)

trajectory_statistics <- trajectory_results %>% transmute(Pathway, overall_trajectory_NumDF = NumDF, 
    overall_trajectory_DenDF = DenDF, overall_trajectory_F = F.value, overall_trajectory_P = p.value, 
    overall_trajectory_FDR = p_adj_BH) %>% left_join(interaction_results %>% filter(Estimand == "D1 to D4 extra divergence") %>% 
    transmute(Pathway, D1_to_D4_beta = estimate, D1_to_D4_SE = SE, D1_to_D4_df = df, D1_to_D4_lower95 = lower.CL, 
        D1_to_D4_upper95 = upper.CL, D1_to_D4_P = p.value, D1_to_D4_FDR = p_adj_BH), by = "Pathway") %>% 
    left_join(interaction_results %>% filter(Estimand == "D1 to D7 extra divergence") %>% transmute(Pathway, 
        D1_to_D7_beta = estimate, D1_to_D7_SE = SE, D1_to_D7_df = df, D1_to_D7_lower95 = lower.CL, D1_to_D7_upper95 = upper.CL, 
        D1_to_D7_P = p.value, D1_to_D7_FDR = p_adj_BH), by = "Pathway") %>% arrange(overall_trajectory_FDR)

all_group_trajectories <- adjusted_group_trajectories %>% left_join(observed_group_trajectories, by = c("Pathway", 
    "Mortality28d", "Timepoint")) %>% left_join(trajectory_statistics %>% dplyr::select(Pathway, overall_trajectory_P, 
    overall_trajectory_FDR, D1_to_D4_P, D1_to_D4_FDR, D1_to_D7_P, D1_to_D7_FDR), by = "Pathway") %>% 
    arrange(Pathway, Mortality28d, Timepoint)

write.csv(arrange(trajectory_results, p_adj_BH), repo_output(repo_path(out_dir, "all_GOBP_trajectory_divergence_omnibus.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

write.csv(day1_results, repo_output(repo_path(out_dir, "all_GOBP_D1_mortality_divergence.csv"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    row.names = FALSE)

write.csv(interaction_results, repo_output(repo_path(out_dir, "all_GOBP_D1_to_D4_D7_extra_divergence.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

write.csv(timepoint_results, repo_output(repo_path(out_dir, "all_GOBP_D1_D4_D7_mortality_differences.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

write.csv(model_diagnostics, repo_output(repo_path(out_dir, "all_GOBP_model_diagnostics.csv"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    row.names = FALSE)

write.csv(trajectory_statistics, repo_output(repo_path(out_dir, "all_GOBP_trajectory_statistics_complete.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

write.csv(all_group_trajectories, repo_output(repo_path(out_dir, "all_GOBP_Survivor_Death_D1_D4_D7_trajectories.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

write.csv(adjusted_group_trajectories, repo_output(repo_path(out_dir, "all_GOBP_adjusted_group_means_95CI.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

write.csv(observed_group_trajectories, repo_output(repo_path(out_dir, "all_GOBP_observed_group_mean_SD.csv"), 
    "Figure2A_Figure3AB_FigureS5A_analysis.R"), row.names = FALSE)

saveRDS(list(trajectory = trajectory_results, contrasts = contrast_results, adjusted_group_trajectories = adjusted_group_trajectories, 
    observed_group_trajectories = observed_group_trajectories, trajectory_statistics = trajectory_statistics, 
    diagnostics = model_diagnostics), repo_output(repo_path(out_dir, "all_GOBP_lmer_results.rds"), "Figure2A_Figure3AB_FigureS5A_analysis.R"))

clean_label <- function(x, width = 45) x %>% str_remove("^GOBP_") %>% str_replace_all("_", " ") %>% str_to_sentence() %>% 
    str_wrap(width)

d4_results <- interaction_results %>% filter(Estimand == "D1 to D4 extra divergence")

overview_df <- bind_rows(mutate(day1_results, Panel = "D1 mortality difference"), mutate(d4_results, 
    Panel = "D1 to D4 extra divergence")) %>% mutate(significant = p_adj_BH < 0.05, minus_log10_fdr = -log10(pmax(p_adj_BH, 
    .Machine$double.xmin)))

p_overview <- ggplot(overview_df, aes(estimate, minus_log10_fdr, colour = significant)) + geom_hline(yintercept = -log10(0.05), 
    linetype = 2, colour = "grey55") + geom_vline(xintercept = 0, colour = "grey75") + geom_point(alpha = 0.65, 
    size = 0.8) + facet_wrap(~Panel, scales = "free_x") + scale_colour_manual(values = c(`FALSE` = "grey65", 
    `TRUE` = "#C51B7D")) + labs(x = "Adjusted difference in standardized GSVA score", y = expression(-log[10](BH - 
    FDR)), colour = "FDR < 0.05") + theme_bw(base_size = 10) + theme(legend.position = "top")

ggsave(repo_output(repo_path(out_dir, "all_GOBP_D1_and_D4_divergence_overview.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    p_overview, width = 10, height = 4.8)

top_n <- 30L

forest_df <- bind_rows(day1_results %>% slice_head(n = top_n) %>% mutate(Analysis = "D1 mortality difference"), 
    d4_results %>% arrange(p_adj_BH) %>% slice_head(n = top_n) %>% mutate(Analysis = "D1 to D4 extra divergence")) %>% 
    group_by(Analysis) %>% arrange(estimate, .by_group = TRUE) %>% mutate(Pathway_label = factor(clean_label(Pathway), 
    levels = clean_label(Pathway))) %>% ungroup()

p_forest <- ggplot(forest_df, aes(estimate, Pathway_label, xmin = lower.CL, xmax = upper.CL)) + geom_vline(xintercept = 0, 
    colour = "grey65", linetype = 2) + geom_errorbar(orientation = "y", width = 0, linewidth = 0.45) + 
    geom_point(aes(colour = p_adj_BH < 0.05), size = 1.5) + facet_wrap(~Analysis, scales = "free", ncol = 2) + 
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#C51B7D")) + labs(x = "Adjusted standardized difference (95% CI)", 
    y = NULL, colour = "FDR < 0.05") + theme_bw(base_size = 8) + theme(legend.position = "top")

ggsave(repo_output(repo_path(out_dir, "top_GOBP_D1_and_D4_divergence_forest.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    p_forest, width = 12, height = 12)

top_trajectory_pathways <- trajectory_results %>% arrange(p_adj_BH) %>% slice_head(n = 24) %>% pull(Pathway)

fit_means <- function(pw) {
    dat <- analysis_meta
    dat$Score_z <- as.numeric(scale(gsva_analysis[pw, ]))
    fit <- lmerTest::lmer(model_formula, dat, REML = FALSE)
    emmeans(fit, ~Mortality28d * Timepoint, lmer.df = "satterthwaite") %>% confint() %>% as.data.frame() %>% 
        mutate(Pathway = pw)
}

adjusted_means <- map_dfr(top_trajectory_pathways, fit_means) %>% mutate(Pathway_label = clean_label(Pathway, 
    32))

write.csv(adjusted_means, repo_output(repo_path(out_dir, "top_trajectory_adjusted_means.csv"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    row.names = FALSE)

p_trajectory <- ggplot(adjusted_means, aes(Timepoint, emmean, group = Mortality28d, colour = Mortality28d, 
    fill = Mortality28d)) + geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_ribbon(aes(ymin = lower.CL, 
    ymax = upper.CL), alpha = 0.12, colour = NA) + geom_line(linewidth = 0.7) + geom_point(size = 1.4) + 
    facet_wrap(~Pathway_label, scales = "free_y", ncol = 4) + scale_colour_manual(values = c(Survivor = "#2C7BB6", 
    Death = "#D7191C")) + scale_fill_manual(values = c(Survivor = "#2C7BB6", Death = "#D7191C")) + labs(x = NULL, 
    y = "Adjusted standardized GSVA score", colour = "28-day outcome", fill = "28-day outcome") + theme_bw(base_size = 8) + 
    theme(strip.text = element_text(size = 6.5), legend.position = "top")

ggsave(repo_output(repo_path(out_dir, "top_GOBP_adjusted_longitudinal_trajectories.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    p_trajectory, width = 12, height = 14)

library(ggrepel)

library(patchwork)

panel_a_df <- day1_results %>% dplyr::select(Pathway, D1_effect = estimate, D1_FDR = p_adj_BH) %>% left_join(d4_results %>% 
    dplyr::select(Pathway, D4_interaction = estimate, D4_interaction_FDR = p_adj_BH), by = "Pathway") %>% 
    left_join(trajectory_results %>% dplyr::select(Pathway, trajectory_FDR = p_adj_BH), by = "Pathway") %>% 
    mutate(D1_significant = D1_FDR < 0.05, trajectory_significant = trajectory_FDR < 0.05, Evidence = case_when(D1_significant & 
        trajectory_significant ~ "Both", D1_significant ~ "D1 only", trajectory_significant ~ "Trajectory only", 
        TRUE ~ "Neither"), Evidence = factor(Evidence, levels = c("Neither", "D1 only", "Trajectory only", 
        "Both")))

evidence_counts <- panel_a_df %>% count(Evidence) %>% mutate(Legend = paste0(Evidence, " (n=", n, ")"))

evidence_labels <- setNames(evidence_counts$Legend, as.character(evidence_counts$Evidence))

panel_a_labels <- c("GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_DETECTION_OF_MOLECULE_OF_BACTERIAL_ORIGIN", 
    "GOBP_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE", "GOBP_RESPONSE_TO_MINERALOCORTICOID")

panel_a <- ggplot(panel_a_df, aes(D1_effect, D4_interaction)) + geom_hline(yintercept = 0, colour = "grey75", 
    linewidth = 0.35) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) + geom_point(aes(colour = Evidence), 
    alpha = 0.72, size = 1.15) + ggrepel::geom_text_repel(data = filter(panel_a_df, Pathway %in% panel_a_labels), 
    aes(label = clean_label(Pathway, 25)), size = 2.4, colour = "grey15", max.overlaps = Inf, box.padding = 0.35, 
    min.segment.length = 0) + scale_colour_manual(values = c(Neither = "grey78", `D1 only` = "#3B78A8", 
    `Trajectory only` = "#D95F02", Both = "#7B3294"), labels = evidence_labels, drop = FALSE) + labs(title = "A", 
    x = "D1 mortality effect, Death - Survivor (SD)", y = "D1 to D4 interaction effect (SD)", colour = NULL) + 
    theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 13), legend.position = "bottom", 
    legend.text = element_text(size = 7))

effect_levels <- c("D1 mortality difference", "D4 mortality difference", "D7 mortality difference", "D1 to D4 extra divergence", 
    "D1 to D7 extra divergence")

effect_labels <- c("D1 effect", "D4 effect", "D7 effect", "D1->D4 interaction", "D1->D7 interaction")

trajectory_30 <- trajectory_results %>% filter(p_adj_BH < 0.05) %>% arrange(p_adj_BH)

panel_b_df <- contrast_results %>% filter(Pathway %in% trajectory_30$Pathway, Estimand %in% effect_levels) %>% 
    mutate(Effect = factor(Estimand, levels = effect_levels, labels = effect_labels), significant = p_adj_BH < 
        0.05)

effect_matrix <- panel_b_df %>% dplyr::select(Pathway, Effect, estimate) %>% pivot_wider(names_from = Effect, 
    values_from = estimate) %>% column_to_rownames("Pathway") %>% as.matrix()

row_order <- rownames(effect_matrix)[hclust(dist(effect_matrix))$order]

panel_b_df <- panel_b_df %>% mutate(Pathway_label = factor(clean_label(Pathway, 50), levels = clean_label(rev(row_order), 
    50)))

heat_limit <- max(abs(panel_b_df$estimate), na.rm = TRUE)

panel_b <- ggplot(panel_b_df, aes(Effect, Pathway_label, fill = estimate)) + geom_tile(colour = "white", 
    linewidth = 0.25) + geom_text(aes(label = if_else(significant, "*", "")), size = 3) + scale_fill_gradient2(low = "#2166AC", 
    mid = "white", high = "#B2182B", midpoint = 0, limits = c(-heat_limit, heat_limit)) + labs(title = "B", 
    x = NULL, y = NULL, fill = "Adjusted effect\n(SD)") + theme_bw(base_size = 8) + theme(plot.title = element_text(face = "bold", 
    size = 13), panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.5), 
    legend.position = "right")

representative_pathways <- c("GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_DETECTION_OF_MOLECULE_OF_BACTERIAL_ORIGIN", 
    "GOBP_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE", "GOBP_NEGATIVE_REGULATION_OF_COMPLEMENT_ACTIVATION", 
    "GOBP_RESPONSE_TO_MINERALOCORTICOID", "GOBP_REGULATION_OF_RUFFLE_ASSEMBLY")

panel_c_means <- map_dfr(representative_pathways, fit_means) %>% mutate(Pathway_label = clean_label(Pathway, 
    30))

panel_c_annotations <- d4_results %>% filter(Pathway %in% representative_pathways) %>% transmute(Pathway, 
    Pathway_label = clean_label(Pathway, 30), label = paste0("D1->D4 interaction FDR = ", format.pval(p_adj_BH, 
        2)), Timepoint = factor("D1", levels = c("D1", "D4", "D7"))) %>% left_join(panel_c_means %>% 
    group_by(Pathway) %>% summarise(y = max(upper.CL) + 0.12, .groups = "drop"), by = "Pathway")

panel_c <- ggplot(panel_c_means, aes(Timepoint, emmean, group = Mortality28d, colour = Mortality28d, 
    fill = Mortality28d)) + geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_ribbon(aes(ymin = lower.CL, 
    ymax = upper.CL), alpha = 0.12, colour = NA) + geom_line(linewidth = 0.7) + geom_point(size = 1.5) + 
    geom_text(data = panel_c_annotations, aes(x = Timepoint, y = y, label = label), inherit.aes = FALSE, 
        hjust = 0, size = 2.3) + facet_wrap(~Pathway_label, scales = "free_y", ncol = 2) + scale_colour_manual(values = c(Survivor = "#2C7BB6", 
    Death = "#D7191C")) + scale_fill_manual(values = c(Survivor = "#2C7BB6", Death = "#D7191C")) + labs(title = "C", 
    x = NULL, y = "Adjusted standardized GSVA score", colour = "28-day outcome", fill = "28-day outcome") + 
    theme_bw(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 13), strip.text = element_text(size = 6.5), 
    legend.position = "bottom")

ggsave(repo_output(repo_path(out_dir, "Fig3A_baseline_vs_D4_interaction.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    panel_a, width = 7, height = 6)

ggsave(repo_output(repo_path(out_dir, "Fig3B_trajectory_pathway_effect_heatmap.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    panel_b, width = 8.2, height = 9.5)

ggsave(repo_output(repo_path(out_dir, "Fig3C_representative_adjusted_trajectories.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    panel_c, width = 8.5, height = 8)

figure_3 <- (panel_a | panel_c)/panel_b + plot_layout(heights = c(0.9, 1.35), widths = c(0.9, 1.1), guides = "collect") & 
    theme(legend.position = "bottom")

ggsave(repo_output(repo_path(out_dir, "Figure3_GOBP_baseline_and_trajectory_divergence.pdf"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    figure_3, width = 16, height = 17)

ggsave(repo_output(repo_path(out_dir, "Figure3_GOBP_baseline_and_trajectory_divergence.png"), "Figure2A_Figure3AB_FigureS5A_analysis.R"), 
    figure_3, width = 16, height = 17, dpi = 300)

run_summary <- c(paste("GSVA input:", normalizePath(repo_input(gsva_file))), paste("Available D1/D4/D7 samples:", 
    nrow(analysis_meta)), paste("Patients:", n_distinct(analysis_meta$HumanID)), paste("GO:BP pathways tested:", 
    length(pathways)), paste("Parallel cores:", n_cores), "Standardization: pathway-specific z score across all D1/D4/D7 samples", 
    paste("Model:", paste(deparse(model_formula), collapse = " ")), "Reference levels: Mortality28d=Survivor; Timepoint=D1", 
    "Trajectory test: joint 2-df test of Mortality x D4 and Mortality x D7", paste("Failed models:", 
        sum(is.na(model_diagnostics$singular))), paste("Singular models:", sum(model_diagnostics$singular, 
        na.rm = TRUE)), capture.output(sessionInfo()))

writeLines(run_summary, repo_output(repo_path(out_dir, "all_GOBP_model_run_info.txt"), "Figure2A_Figure3AB_FigureS5A_analysis.R"))

