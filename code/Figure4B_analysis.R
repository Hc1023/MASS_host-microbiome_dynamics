# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(tidyverse)
    library(nlme)
    library(survival)
    library(JM)
    library(emmeans)
})

root <- "."

score_file <- "Outputs/Figure3C_prepare_MASS_GSVA_scores_long_with_metadata.csv"

metadata_file <- repo_path(dirname(root), "Inputs/1211_metadata.rdata")

original_file <- repo_path(root, "Outputs", "Inputs/260907_trajectory_interaction_tests.csv")

out_dir <- repo_path(root, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(repo_exists(score_file), repo_exists(metadata_file), repo_exists(original_file))

old_to_id <- c(`PRR/TLR/TNF/NF-kB signaling` = "M1", `Phagolysosome/autophagy` = "M2", `Type I interferon response` = "M3", 
    `Ribosome biogenesis` = "M4")

module_labels <- c(M1 = "PRR/TLR/NF-kB signaling", M2 = "Phagolysosomal/autophagy", M3 = "Type I IFN response", 
    M4 = "Ribosome biogenesis")

e <- new.env(parent = emptyenv())

load(repo_input(metadata_file), envir = e)

if (!exists("meta", envir = e, inherits = FALSE)) stop("Missing metadata object: meta")

clinical <- e$meta %>% transmute(HumanID, Mortality28d = as.numeric(as.character(Mortality28d)), event_time = pmin(as.numeric(SurvivalTimeWithin28Days), 
    28)) %>% filter(Mortality28d %in% 0:1, !is.na(event_time), event_time > 0) %>% distinct(HumanID, 
    .keep_all = TRUE)

if (any(clinical$Mortality28d == 0 & clinical$event_time != 28)) stop("Survivors are not consistently censored at Day 28")

dat <- read.csv(repo_input(score_file), check.names = FALSE, stringsAsFactors = FALSE) %>% transmute(HumanID, 
    SampleID, Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7")), TimeDays = recode(as.character(Timepoint), 
        D1 = 1, D4 = 4, D7 = 7) %>% as.numeric(), TimeSinceD1 = TimeDays - 1, Module_ID = unname(old_to_id[Module]), 
    Score) %>% inner_join(clinical, by = "HumanID") %>% filter(!is.na(Score), !is.na(TimeDays), !is.na(Module_ID)) %>% 
    mutate(Module_ID = factor(Module_ID, levels = names(module_labels))) %>% arrange(Module_ID, HumanID, 
    TimeDays)

stopifnot(!anyDuplicated(dat[c("SampleID", "Module_ID")]))

fit_one <- function(id) {
    message("Fitting joint model: ", id)
    d <- filter(dat, Module_ID == id) %>% mutate(Mortality28d = factor(Mortality28d, levels = c(0, 1))) %>% 
        arrange(HumanID, TimeDays)
    lf <- try(nlme::lme(Score ~ TimeSinceD1 * Mortality28d, random = ~TimeSinceD1 | HumanID, data = d, 
        method = "REML", na.action = na.omit, control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)), 
        silent = TRUE)
    random_structure <- "random intercept + TimeSinceD1 slope"
    if (inherits(lf, "try-error")) {
        lf <- nlme::lme(Score ~ TimeSinceD1 * Mortality28d, random = ~1 | HumanID, data = d, method = "REML", 
            na.action = na.omit, control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
        random_structure <- "random intercept (slope model did not converge)"
    }
    sd <- d %>% transmute(HumanID, time = event_time, event = as.numeric(as.character(Mortality28d))) %>% 
        distinct(HumanID, .keep_all = TRUE)
    cf <- coxph(Surv(time, event) ~ 1, data = sd, x = TRUE)
    jf <- JM::jointModel(lmeObject = lf, survObject = cf, timeVar = "TimeSinceD1", method = "weibull-PH-aGH")
    long_tab <- summary(jf)$`CoefTable-Long`
    ir <- grep("TimeSinceD1:Mortality28d", rownames(long_tab))
    if (length(ir) != 1) 
        stop("Cannot identify JM interaction for ", id)
    event_tab <- summary(jf)$`CoefTable-Event`
    ar <- grep("Assoct|association", rownames(event_tab), ignore.case = TRUE)
    if (length(ar) != 1) 
        ar <- nrow(event_tab)
    rf <- if (grepl("slope", random_structure)) 
        ~TimeSinceD1 | HumanID
    else ~1 | HumanID
    full <- nlme::lme(Score ~ Timepoint * Mortality28d, random = rf, data = d, method = "ML", na.action = na.omit, 
        control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
    reduced <- update(full, fixed = Score ~ Timepoint + Mortality28d)
    omnibus_p <- anova(reduced, full)$`p-value`[2]
    cbg <- emmeans::contrast(emmeans(full, ~Timepoint | Mortality28d), method = "trt.vs.ctrl", ref = 1)
    contrasts <- as.data.frame(emmeans::contrast(cbg, method = "revpairwise", by = "contrast", adjust = "none")) %>% 
        transmute(Module_ID = id, interval = contrast, estimate, SE, df, t_ratio = t.ratio, p_value = p.value)
    pg <- expand_grid(TimeDays = c(1, 4, 7), Mortality28d = factor(c(0, 1), levels = c(0, 1))) %>% mutate(TimeSinceD1 = TimeDays - 
        1, HumanID = d$HumanID[1], Score = 0, Timepoint = factor(paste0("D", TimeDays), levels = c("D1", 
        "D4", "D7")))
    pred <- as.data.frame(predict(jf, newdata = pg, type = "Marginal", interval = "confidence", returnData = TRUE)) %>% 
        transmute(Module_ID = id, TimeDays, Timepoint, Mortality28d, Predicted = pred, Lower = low, Upper = upp)
    list(jm = jf, pred = pred, contrasts = contrasts, random_structure = random_structure, slope = unname(long_tab[ir, 
        "Value"]), slope_se = unname(long_tab[ir, "Std.Err"]), slope_p = unname(long_tab[ir, "p-value"]), 
        alpha = unname(event_tab[ar, "Value"]), alpha_se = unname(event_tab[ar, "Std.Err"]), alpha_p = unname(event_tab[ar, 
            "p-value"]), omnibus_p = omnibus_p)
}

fits <- setNames(lapply(names(module_labels), fit_one), names(module_labels))

pred <- bind_rows(lapply(fits, `[[`, "pred")) %>% mutate(Module = factor(unname(module_labels[Module_ID]), 
    levels = unname(module_labels)), Outcome = factor(Mortality28d, 0:1, c("Survival", "Mortality")))

contrasts <- bind_rows(lapply(fits, `[[`, "contrasts")) %>% mutate(fdr = p.adjust(p_value, "BH"))

tests <- tibble(Module_ID = names(module_labels), Module = unname(module_labels), JM_slope_difference = vapply(fits, 
    `[[`, numeric(1), "slope"), JM_slope_SE = vapply(fits, `[[`, numeric(1), "slope_se"), JM_slope_P = vapply(fits, 
    `[[`, numeric(1), "slope_p"), alpha = vapply(fits, `[[`, numeric(1), "alpha"), alpha_SE = vapply(fits, 
    `[[`, numeric(1), "alpha_se"), alpha_P = vapply(fits, `[[`, numeric(1), "alpha_p"), categorical_trajectory_P = vapply(fits, 
    `[[`, numeric(1), "omnibus_p")) %>% mutate(JM_slope_FDR = p.adjust(JM_slope_P, "BH"), categorical_trajectory_FDR = p.adjust(categorical_trajectory_P, 
    "BH"), alpha_lower95 = alpha - 1.96 * alpha_SE, alpha_upper95 = alpha + 1.96 * alpha_SE)

write.csv(tests, repo_output(repo_path(out_dir, "260907_JM_estimates.csv"), "Figure4B_analysis.R"), row.names = FALSE)

write.csv(contrasts, repo_output(repo_path(out_dir, "260907_JM_trajectory_change_contrasts.csv"), "Figure4B_analysis.R"), 
    row.names = FALSE)

write.csv(pred, repo_output(repo_path(out_dir, "260907_JM_marginal_predictions.csv"), "Figure4B_analysis.R"), 
    row.names = FALSE)

ann <- tests %>% mutate(Module = factor(Module, levels = unname(module_labels)), label = paste0("JM trajectory interaction\nFDR ", 
    if_else(JM_slope_FDR < 0.001, "<0.001", paste0("= ", formatC(JM_slope_FDR, digits = 3, format = "f")))))

p <- ggplot(pred, aes(TimeDays, Predicted, color = Outcome, fill = Outcome, group = Outcome)) + geom_ribbon(aes(ymin = Lower, 
    ymax = Upper), alpha = 0.16, linewidth = 0, color = NA) + geom_line(linewidth = 0.65) + geom_point(size = 1.6) + 
    geom_label(data = ann, aes(x = 1, y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1.1, 
        size = 2.2, lineheight = 0.95, fontface = "bold", linewidth = 0, fill = scales::alpha("white", 
            0.8)) + scale_color_manual(values = c(Survival = "#4575B4", Mortality = "#D73027"), name = NULL) + 
    scale_fill_manual(values = c(Survival = "#4575B4", Mortality = "#D73027"), name = NULL) + scale_x_continuous(breaks = c(1, 
    4, 7), labels = c("D1", "D4", "D7")) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) + 
    facet_wrap(~Module, scales = "free_y", ncol = 2) + labs(x = NULL, y = "Joint-model estimated GSVA score") + 
    theme_bw(base_size = 8) + theme(legend.position = "top", axis.text.x = element_text(face = "bold"), 
    panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 7), legend.key.size = grid::unit(3, 
        "mm"), legend.margin = margin(0, 0, 0, 0), panel.spacing = grid::unit(2, "mm"), plot.margin = margin(3, 
        3, 3, 3))

p

ggsave(repo_output(repo_path(out_dir, "260907_four_modules_trajectory_JM.pdf"), "Figure4B_analysis.R"), 
    p, width = 4, height = 3)

original <- read.csv(repo_input(original_file), check.names = FALSE)

cc <- contrasts %>% mutate(day = if_else(str_detect(interval, "D4"), "D4", "D7")) %>% dplyr::select(Module_ID, 
    day, JM_effect = estimate) %>% pivot_wider(names_from = day, values_from = JM_effect, names_prefix = "JM_")

comparison <- original %>% dplyr::select(Module_ID, Original_overall_P = interaction_p, Original_overall_FDR = interaction_fdr) %>% 
    left_join(cc, by = "Module_ID") %>% mutate(Module = unname(module_labels[Module_ID]), .after = Module_ID)

write.csv(comparison, repo_output(repo_path(out_dir, "260907_mixed_vs_JM_comparison.csv"), "Figure4B_analysis.R"), 
    row.names = FALSE)

diagnostics <- c("Joint longitudinal-survival model specification", "Longitudinal: Score ~ TimeSinceD1 * Mortality28d", 
    "Survival: Surv(true 28-day event time, event) ~ 1", "JM: JM::jointModel, weibull-PH-aGH; default current-value association", 
    "", unlist(map(names(module_labels), ~c(paste0("[", .x, "] ", module_labels[.x]), paste("random structure:", 
        fits[[.x]]$random_structure), paste("JM convergence code:", fits[[.x]]$jm$convergence), ""))))

writeLines(diagnostics, repo_output(repo_path(out_dir, "260907_JM_model_diagnostics.txt"), "Figure4B_analysis.R"))

summary_lines <- tests %>% transmute(x = paste0(Module_ID, " (", Module, "): JM slope interaction ", 
    signif(JM_slope_difference, 3), " (P=", signif(JM_slope_P, 3), "); alpha ", signif(alpha, 3), " [", 
    signif(alpha_lower95, 3), ", ", signif(alpha_upper95, 3), "].")) %>% pull(x)

writeLines(c(summary_lines, "", "This analysis uses the established JM framework without modification."), 
    repo_output(repo_path(out_dir, "260907_JM_key_results.txt"), "Figure4B_analysis.R"))

capture.output(sessionInfo(), file = repo_path(out_dir, "260907_JM_sessionInfo.txt"))

message("Done. Joint-model results: ", out_dir)

