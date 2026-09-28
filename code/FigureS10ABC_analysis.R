# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(emmeans)
    library(patchwork)
})

root <- "."

host_root <- repo_path(dirname(root), ".")

input_file <- repo_path(host_root, "Inputs/1616_meta_model.rdata")

reference_script <- repo_path(host_root, "code/shared_CTS_transition_groups.R")

out <- repo_path(root, "Outputs")

repo_dir(out, recursive = TRUE, showWarnings = FALSE)

path <- function(x) repo_path(out, paste0("260915_", x))

write_table <- function(x, name) write.csv(x, repo_output(path(paste0(name, ".csv")), "FigureS10ABC_analysis.R"), 
    row.names = FALSE)

save_plot <- function(p, name, width, height) {
    ggsave(repo_output(path(paste0(name, ".pdf")), "FigureS10ABC_analysis.R"), p, width = width, height = height, 
        device = grDevices::pdf, family = "Helvetica", useDingbats = FALSE)
    ggsave(repo_output(path(paste0(name, ".png")), "FigureS10ABC_analysis.R"), p, width = width, height = height, 
        dpi = 300, bg = "white")
}

theme_pub <- function() theme_bw(base_size = 8, base_family = "Helvetica") + theme(panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill = "grey94", colour = "grey50", linewidth = 0.3), strip.text = element_text(face = "bold", 
        size = 7), legend.position = "top", legend.key.size = grid::unit(3, "mm"), legend.margin = margin(0, 
        0, 0, 0), axis.text = element_text(size = 6.5), plot.caption = element_text(size = 6), plot.title = element_text(face = "bold", 
        hjust = 0.5, size = 8), plot.margin = margin(3, 3, 3, 3), panel.spacing = grid::unit(2, "mm"))

cts_cols <- c(CTS3 = "#777777", CTS2 = "#8DA0CB", CTS1 = "#4575B4")

pool_cols <- c(`Stable CTS3` = "#777777", `CTS3 to CTS1/2` = "#8DA0CB")

features <- c(Total = "Total", Viral = "Viruses", BacteriaFungi = "Bacf")

feature_labels <- c(Total = "Total", Viral = "Viruses", BacteriaFungi = "Bacteria/Fungi")

covariates <- c("Age", "Gender", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

input_env <- new.env()

load("Inputs/1616_meta_model.rdata", input_env)

meta <- input_env$meta_model

stopifnot(!anyDuplicated(meta$HumanID))

required <- c("HumanID", "Mortality28d", "D1_CTS", "D4_CTS", "D7_CTS")

meta$complete_CTS <- complete.cases(meta[, required])

meta$trajectory_group <- NA_character_

meta$trajectory_group[meta$complete_CTS] <- with(meta[meta$complete_CTS, ], case_when(D1_CTS == 3 & D7_CTS != 
    3 ~ "Deteriorating", D1_CTS != 3 & D7_CTS == 3 ~ "Improving", D1_CTS == 3 & D4_CTS == 3 & D7_CTS == 
    3 ~ "Stable low-risk", D1_CTS != 3 & D4_CTS != 3 & D7_CTS != 3 ~ "Stable high-risk", TRUE ~ "Mixed"))

reference_lines <- readLines(repo_input(reference_script), warn = FALSE)

start <- grep("^df1 = meta_model", reference_lines)

end <- grep("^df_sum <-", reference_lines)

end <- end[end > start][1]

stopifnot(length(start) == 1, length(end) == 1, end > start)

reference_env <- new.env()

reference_env$meta_model <- input_env$meta_model

reference_env$`%<>%` <- magrittr::`%<>%`

eval(parse(text = reference_lines[start:(end - 1)]), reference_env)

reference_groups <- reference_env$df1

new_groups <- meta[meta$complete_CTS, ]

stopifnot(setequal(new_groups$HumanID, reference_groups$HumanID), identical(new_groups$trajectory_group, 
    reference_groups$grp[match(new_groups$HumanID, reference_groups$HumanID)]))

meta$included <- meta$complete_CTS & !is.na(meta$trajectory_group) & meta$trajectory_group != "Mixed"

meta$exclusion_reason <- ifelse(!meta$complete_CTS, "Missing CTS or mortality required by main figure", 
    ifelse(meta$trajectory_group == "Mixed", "Mixed trajectory excluded by main figure", "Included"))

write_table(meta %>% select(all_of(required), complete_CTS, trajectory_group, included, exclusion_reason), 
    "cohort_audit")

main_counts <- meta %>% filter(complete_CTS) %>% group_by(trajectory_group) %>% summarise(N = n(), deaths = sum(Mortality28d == 
    1), survivors = sum(Mortality28d == 0), .groups = "drop")

write_table(main_counts, "main_figure_group_counts")

wide <- meta %>% filter(included)

for (o in names(features)) for (day in c("D1", "D7")) wide[[paste0(o, "_", day)]] <- wide[[paste0(day, 
    "_mi2_", features[o])]]

measurements <- unlist(lapply(names(features), function(o) paste0(o, "_", c("D1", "D7"))))

stopifnot(all(complete.cases(wide[, c(measurements, covariates)])), all(vapply(wide[measurements], function(x) all(is.finite(x)), 
    logical(1))))

wide$CTS_D1 <- factor(paste0("CTS", wide$D1_CTS), c("CTS1", "CTS2", "CTS3"))

wide$CTS_D7 <- factor(paste0("CTS", wide$D7_CTS), c("CTS1", "CTS2", "CTS3"))

write_table(wide %>% select(HumanID, all_of(required[-1]), trajectory_group, all_of(covariates), all_of(measurements)), 
    "analysis_dataset")

paired <- bind_rows(lapply(names(features), function(o) wide %>% transmute(HumanID, trajectory_group, 
    CTS_D1, CTS_D7, outcome = o, baseline = .data[[paste0(o, "_D1")]], D7 = .data[[paste0(o, "_D7")]], 
    change = D7 - baseline)))

transition_stats <- paired %>% group_by(outcome, CTS_D1, CTS_D7) %>% summarise(N = n(), median_delta = median(change), 
    mean_delta = mean(change), P = if (all(change == 0)) 1 else wilcox.test(D7, baseline, paired = TRUE, 
        exact = FALSE)$p.value, .groups = "drop") %>% mutate(FDR = p.adjust(P, "BH"))

write_table(transition_stats, "01_transition_statistics")

heat <- lapply(names(features), function(o) {
    right_title <- c(Total = "Total", Viral = "Viruses", BacteriaFungi = "Bacteria/\nFungi")
    ggplot(filter(transition_stats, outcome == o), aes(CTS_D7, CTS_D1, fill = median_delta)) + geom_tile(colour = "white", 
        linewidth = 0.8) + geom_text(aes(label = sprintf("%.2f\nP = %.3f", median_delta, P)), size = 2) + 
        scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) + scale_x_discrete(drop = FALSE) + 
        scale_y_discrete(drop = FALSE) + labs(x = if (o == "BacteriaFungi") 
        "D7 CTS"
    else NULL, y = "D1 CTS", fill = paste0(right_title[o], "\n\nMedian")) + theme_pub() + theme(legend.position = "right", 
        panel.grid = element_blank(), legend.title = element_text(size = 8, lineheight = 0.9), axis.text.x = if (o == 
            "BacteriaFungi") 
            element_text(size = 6.5)
        else element_blank()) + guides(fill = guide_colourbar(title.position = "top", barwidth = grid::unit(2, 
        "mm"), barheight = grid::unit(13, "mm")))
})

p_heat <- wrap_plots(heat, ncol = 1, guides = "keep")

save_plot(p_heat, "01_CTS_transition_microbes", 3.2, 4.5)

pdat <- wide %>% filter(trajectory_group %in% c("Stable low-risk", "Deteriorating")) %>% mutate(CTS_D7 = factor(CTS_D7, 
    c("CTS3", "CTS2", "CTS1")), pooled_group = factor(ifelse(trajectory_group == "Stable low-risk", "Stable CTS3", 
    "CTS3 to CTS1/2"), levels = names(pool_cols)))

stopifnot(!anyDuplicated(pdat$HumanID), all(table(pdat$CTS_D7)>0))

write_table(pdat %>% count(pooled_group, CTS_D7, name = "N"), "CTS3_model_group_counts")

models <- list()

emm_tables <- list()

contrasts <- list()

diagnostics <- list()

coefficients <- list()

for (group in c("CTS_D7", "pooled_group")) for (o in names(features)) {
    key <- paste(group, o, sep = "_")
    formula <- reformulate(c(paste0(o, "_D1"), group, covariates), response = paste0(o, "_D7"))
    fit <- lm(formula, data = pdat, na.action = na.fail)
    stopifnot(nobs(fit) == nrow(pdat), fit$rank == ncol(model.matrix(fit)), !anyNA(coef(fit)))
    models[[key]] <- fit
    em <- emmeans(fit, specs = group, weights = "equal")
    emm_tables[[key]] <- as.data.frame(confint(em)) %>% mutate(outcome = o, analysis = group, N = nobs(fit))
    ct <- if (group == "CTS_D7") 
        contrast(em, list(CTS2_vs_stable_CTS3 = c(-1, 1, 0), CTS1_vs_stable_CTS3 = c(-1, 0, 1)), adjust = "none")
    else contrast(em, list(Deteriorating_vs_stable_CTS3 = c(-1, 1)), adjust = "none")
    contrasts[[key]] <- as.data.frame(summary(ct, infer = c(TRUE, TRUE))) %>% mutate(outcome = o, analysis = group, 
        N = nobs(fit))
    diagnostics[[key]] <- data.frame(analysis = group, outcome = o, N = nobs(fit), rank = fit$rank, parameters = ncol(model.matrix(fit)), 
        residual_df = df.residual(fit), max_cooks_distance = max(cooks.distance(fit)))
    coefficients[[key]] <- broom::tidy(fit, conf.int = TRUE) %>% mutate(outcome = o, analysis = group)
}

emms <- bind_rows(emm_tables)

cts <- bind_rows(contrasts) %>% group_by(analysis) %>% mutate(FDR = p.adjust(p.value, "BH")) %>% ungroup()

write_table(emms, "adjusted_emmeans")

write_table(cts, "adjusted_contrasts")

write_table(bind_rows(coefficients), "model_coefficients")

write_table(bind_rows(diagnostics), "model_diagnostics")

saveRDS(models, repo_output(path("models.rds"), "FigureS10ABC_analysis.R"))

n_dest <- table(pdat$CTS_D7)

dest_ann <- cts %>% filter(analysis == "CTS_D7") %>% mutate(CTS_D7 = ifelse(contrast == "CTS2_vs_stable_CTS3", 
    "CTS2", "CTS1"), outcome = factor(outcome, names(features)), label = sprintf("beta = %.2f\nP = %.3f", 
    estimate, p.value))

pool_ann <- cts %>% filter(analysis == "pooled_group") %>% mutate(outcome = factor(outcome, names(features)), 
    label = sprintf("beta = %.2f; P = %.3f", estimate, p.value))

dest_labels <- c(CTS3 = paste0("Stable CTS3\nN=", n_dest["CTS3"]), CTS2 = paste0("CTS3 -> CTS2\nN=", 
    n_dest["CTS2"]), CTS1 = paste0("CTS3 -> CTS1\nN=", n_dest["CTS1"]))

p_dest <- ggplot(filter(emms, analysis == "CTS_D7") %>% mutate(outcome = factor(outcome, names(features))), 
    aes(CTS_D7, emmean, colour = CTS_D7)) + geom_line(aes(group = 1), colour = "grey70", linewidth = 0.45) + 
    geom_point(size = 1.6) + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.08, linewidth = 0.6) + 
    geom_text(data = dest_ann, aes(x = CTS_D7, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.1, 
        size = 2.8) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.27))) + facet_wrap(~outcome, 
    scales = "free_y", labeller = as_labeller(feature_labels)) + scale_x_discrete(limits = c("CTS3", 
    "CTS2", "CTS1"), labels = dest_labels) + scale_colour_manual(values = cts_cols) + labs(x = NULL, 
    y = "Adjusted D7 burden", colour = "D7 CTS", title = NULL, caption = "Beta vs stable CTS3; unadjusted P.") + 
    theme_pub() + theme(axis.text.x = element_text(size = 6))

save_plot(p_dest, "02_CTS3_destination_adjusted_D7", 6.7, 2.6)

pool_n <- table(pdat$pooled_group)

pool_labels <- setNames(paste0(c("Stable CTS3", "CTS3 -> CTS1/2"), "\nN=", as.integer(pool_n)), names(pool_n))

p_pool <- ggplot(filter(emms, analysis == "pooled_group") %>% mutate(outcome = factor(outcome, names(features))), 
    aes(pooled_group, emmean, colour = pooled_group)) + geom_point(size = 1.6) + geom_errorbar(aes(ymin = lower.CL, 
    ymax = upper.CL), width = 0.08) + geom_text(data = pool_ann, aes(x = 1.5, y = Inf, label = label), 
    inherit.aes = FALSE, vjust = 1.2, size = 2.8) + scale_y_continuous(expand = expansion(mult = c(0.05, 
    0.18))) + facet_wrap(~outcome, scales = "free_y", labeller = as_labeller(feature_labels)) + scale_x_discrete(labels = pool_labels) + 
    scale_colour_manual(values = pool_cols) + labs(x = NULL, y = "Adjusted D7 burden", colour = NULL, 
    title = NULL, caption = "Beta: deteriorating minus stable CTS3; unadjusted P.") + theme_pub()

save_plot(p_pool, "03_CTS3_pooled_adjusted_D7", 5.5, 2.4)

write_table(paired, "paired_microbe_data")

writeLines(c("Sources: 260904.R (heatmap), 260904_2.R (destination ANCOVA), 260904_3.R (pooled ANCOVA).", 
    paste("Grouping reference:", reference_script), paste("Input:", input_file), "Cohort and every group assignment verified against the executable grouping block of the main-figure script.", 
    paste("Complete D1/D4/D7 CTS and mortality:", sum(meta$complete_CTS)), paste("Mixed excluded:", sum(meta$trajectory_group == 
        "Mixed", na.rm = TRUE)), paste("Non-Mixed cohort for heatmap:", nrow(wide)), paste("Rebuilt CTS3-origin model counts:",paste(capture.output(table(pdat$CTS_D7)),collapse=" ")), 
    "Historical extra deterioration patients ZZRM_091 and TZEZ_023 have missing D4 CTS.", "Use existing D1/D7_mi2_Total, mi2_Viruses and mi2_Bacf without another transformation.", 
    "Heatmap: median paired D7-D1 change; per-cell N and paired Wilcoxon P/FDR in statistics CSV.", paste("ANCOVA: D7 burden ~ D1 burden + transition group +", 
        paste(covariates, collapse = " + ")), "Points/intervals: emmeans with equal factor weighting and 95% confidence intervals.", 
    "BH correction separately across 6 destination contrasts and 3 pooled contrasts; 27 heatmap tests.", 
    "All six ANCOVAs full rank; no patients dropped during model fitting."), repo_output(path("run_info.txt"), 
    "FigureS10ABC_analysis.R"))

writeLines(capture.output(sessionInfo()), repo_output(path("sessionInfo.txt"), "FigureS10ABC_analysis.R"))

message("Completed three figures. Main-figure-aligned cohort N=", nrow(wide), "; CTS3-origin models N=", 
    nrow(pdat))

