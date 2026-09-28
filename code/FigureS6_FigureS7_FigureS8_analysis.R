# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(tidyverse)
    library(lmerTest)
})

args <- commandArgs(trailingOnly = FALSE)

script_arg <- grep("^--file=", args, value = TRUE)

root <- "."

stopifnot(repo_exists(repo_path(root, "code/Figure3C_plot.R")))

input_root <- repo_path(dirname(root), ".")

out_root <- repo_path(root, "Outputs")

repo_dir(out_root, recursive = TRUE, showWarnings = FALSE)

meta_file <- repo_path(input_root, "Inputs/1211_metadata.rdata")

biomass_file <- repo_path(input_root, "Inputs/1616_microbe.rdata")

gsva_file <- "Outputs/Module_scores.csv"

pathogen_file <- repo_path(root, "Inputs/meta_pathogen_binary_types.csv")

me <- new.env()

load(repo_input(meta_file), envir = me)

mi <- new.env()

load(repo_input(biomass_file), envir = mi)

module_labels <- c(M1 = "PRR/TLR/NF-kB signaling", M2 = "Phagolysosomal/autophagy", M3 = "Type I IFN response", 
    M4 = "Ribosome biogenesis")

module_ids <- c(up1 = "M1", up2 = "M2", up3 = "M3", dw1 = "M4")

outcome_colors <- c(Survival = "#4575B4", Mortality = "#D73027")

stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

format_p <- function(x) if_else(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))

meta <- me$meta %>% select(HumanID, Immunosuppression, PneumoniaTypeGroup) %>% distinct()

pathogen <- read.csv(repo_input(pathogen_file), check.names = FALSE) %>% select(HumanID, Bacteria_Infection, 
    Fungus_Infection, Virus_Infection)

stopifnot(!anyDuplicated(meta$HumanID), !anyDuplicated(pathogen$HumanID))

sample_meta <- me$df_long %>% select(HumanID, SampleID, Timepoint, Mortality28d) %>% distinct() %>% left_join(meta, 
    by = "HumanID") %>% left_join(pathogen, by = "HumanID") %>% mutate(Outcome = factor(as.character(Mortality28d), 
    levels = c("0", "1"), labels = c("Survival", "Mortality")), Timepoint = factor(Timepoint, levels = c("D1", 
    "D4", "D7", "D14", "D21")))

stopifnot(!anyDuplicated(sample_meta$SampleID))

gsva <- as.matrix(read.csv(repo_input(gsva_file), row.names = 1, check.names = FALSE))

stopifnot(all(names(module_ids) %in% rownames(gsva)))

host <- as.data.frame(t(gsva[names(module_ids), , drop = FALSE])) %>% rownames_to_column("SampleID") %>% 
    inner_join(sample_meta, by = "SampleID") %>% pivot_longer(all_of(names(module_ids)), names_to = "Original_module", 
    values_to = "Score") %>% mutate(Module_ID = factor(unname(module_ids[Original_module]), levels = names(module_labels)), 
    Module = unname(module_labels[as.character(Module_ID)])) %>% filter(!is.na(Outcome), is.finite(Score))

stopifnot(!anyDuplicated(host[c("SampleID", "Module_ID")]))

d1 <- sample_meta %>% filter(Timepoint == "D1", !is.na(Outcome))

stopifnot(all(d1$SampleID %in% colnames(mi$data)))

bm <- mi$data[, d1$SampleID, drop = FALSE]

virus_rows <- mi$mapid_vec[rownames(bm)] == "Viruses"

stopifnot(!anyNA(virus_rows), all(c("HHV-4", "HCMV") %in% rownames(bm)))

biomass <- d1 %>% mutate(Total = colSums(bm), Viruses = colSums(bm[virus_rows, , drop = FALSE]), EBV = as.numeric(bm["HHV-4", 
    ]), HCMV = as.numeric(bm["HCMV", ])) %>% pivot_longer(c(Total, Viruses, EBV, HCMV), names_to = "Feature", 
    values_to = "Mass") %>% mutate(Feature = factor(Feature, levels = c("Total", "Viruses", "EBV", "HCMV")), 
    Mass_log2 = log2(Mass + 1))

stopifnot(all(is.finite(biomass$Mass_log2)))

select_subgroup <- function(d, key) {
    switch(key, No_immunosuppression = filter(d, as.character(Immunosuppression) == "0"), Immunosuppression = filter(d, 
        as.character(Immunosuppression) == "1"), Bacteria_fungus_infection = filter(d, Bacteria_Infection == 
        1 | Fungus_Infection == 1), Viral_infection = filter(d, Virus_Infection == 1), CAP = filter(d, 
        PneumoniaTypeGroup == "CAP"), NP = filter(d, PneumoniaTypeGroup == "non-CAP"))
}

wilcox_summary <- function(d, value) {
    x <- d[[value]]
    stopifnot(all(c("Survival", "Mortality") %in% as.character(d$Outcome)))
    w <- wilcox.test(x ~ d$Outcome, exact = FALSE)
    tibble(n_survival = sum(d$Outcome == "Survival"), n_mortality = sum(d$Outcome == "Mortality"), statistic = unname(w$statistic), 
        p_value = w$p.value)
}

fit_interaction <- function(module_id, dat) {
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

make_plot <- function(module_id, summary_df, trajectory_tests, timepoint_tests) {
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

make_d1_host <- function(d, tests, label) {
    ann <- d %>% group_by(Module_ID) %>% summarise(y = max(Score) + 0.06 * diff(range(Score)), .groups = "drop") %>% 
        left_join(tests, by = "Module_ID")
    counts <- d %>% distinct(HumanID, Outcome)
    subtitle <- sprintf("%s (n = %d; Survival %d vs Mortality %d)", label, nrow(counts), sum(counts$Outcome == 
        "Survival"), sum(counts$Outcome == "Mortality"))
    ggplot(d, aes(Module_ID, Score, fill = Outcome, color = Outcome)) + stat_boxplot(geom = "errorbar", 
        width = 0.28, linewidth = 0.55, position = position_dodge(0.72)) + geom_boxplot(width = 0.65, 
        linewidth = 0.55, outlier.shape = 21, outlier.size = 1.8, position = position_dodge(0.72), alpha = 0.9) + 
        geom_text(data = ann, aes(Module_ID, y, label = stars(p_value)), inherit.aes = FALSE, fontface = "bold", 
            size = 3.2) + scale_fill_manual(values = outcome_colors, name = NULL) + scale_color_manual(values = outcome_colors, 
        name = NULL) + scale_x_discrete(labels = c(M1 = "PRR/TLR/NF-kB\nsignaling", M2 = "Phagolysosomal/\nautophagy", 
        M3 = "Type I IFN\nresponse", M4 = "Ribosome\nbiogenesis")) + scale_y_continuous(expand = expansion(mult = c(0.08, 
        0.12))) + labs(x = NULL, y = "GSVA score", title = "D1 host transcriptomic GSVA scores", subtitle = subtitle) + 
        theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), 
        legend.position = "top")
}

make_d1_biomass <- function(d, tests, label) {
    ann <- tests %>% mutate(y = max(d$Mass_log2) * 1.07, label = paste0("p = ", format.pval(p_value, 
        digits = 3, eps = 1e-05)), label_color = if_else(p_value < 0.05, "p_lt_0.05", "p_ge_0.05"))
    ggplot(d, aes(Feature, Mass_log2, fill = Outcome)) + geom_boxplot(aes(color = Outcome), outlier.shape = NA, 
        alpha = 0.75, width = 0.55, position = position_dodge(0.72)) + geom_point(aes(color = Outcome), 
        size = 1.4, alpha = 0.65, position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72, 
            seed = 260916)) + geom_text(data = ann, aes(Feature, y, label = label, color = label_color), 
        inherit.aes = FALSE, show.legend = FALSE, size = 3.5) + scale_color_manual(values = c(outcome_colors, 
        p_lt_0.05 = "#c43932", p_ge_0.05 = "black"), breaks = names(outcome_colors), name = NULL) + scale_fill_manual(values = outcome_colors, 
        name = NULL) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) + labs(x = NULL, 
        y = "log2(mass + 1)", title = "D1 microbe biomass comparison", subtitle = label) + theme_bw(base_size = 11) + 
        theme(panel.grid.minor = element_blank(), legend.position = "top", plot.title = element_text(face = "bold"))
}

analyses <- list(S8_Immunosuppression = c(No_immunosuppression = "No immunosuppression", Immunosuppression = "Immunosuppression"), 
    S9_Infection_type = c(Bacteria_fungus_infection = "Bacteria/fungus infection", Viral_infection = "Viral infection"), 
    S10_Pneumonia_type = c(CAP = "CAP", NP = "NP"))

for (analysis in names(analyses)) {
    out_dir <- repo_path(out_root, analysis)
    repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)
    counts <- list()
    for (key in names(analyses[[analysis]])) {
        label <- analyses[[analysis]][[key]]
        dat <- select_subgroup(host, key) %>% filter(Timepoint %in% c("D1", "D4", "D7")) %>% mutate(HumanID = factor(HumanID), 
            Timepoint = droplevels(Timepoint))
        a <- filter(dat, Timepoint == "D1")
        b <- select_subgroup(biomass, key)
        stopifnot(nrow(a) > 0, nrow(b) > 0)
        trajectory_tests <- map_dfr(names(module_labels), fit_interaction, dat = dat) %>% mutate(interaction_fdr = p.adjust(interaction_p, 
            method = "BH"), Module = unname(module_labels[Module_ID]), label = paste0("Trajectory interaction FDR ", 
            format_p(interaction_fdr)))
        timepoint_tests <- dat %>% group_by(Module_ID, Timepoint) %>% group_modify(~wilcox_summary(.x, 
            "Score")) %>% ungroup() %>% mutate(p_adj_BH = p.adjust(p_value, "BH"), significance = stars(p_adj_BH), 
            label = significance)
        summary_df <- dat %>% group_by(Module_ID, Timepoint, Outcome) %>% summarise(n = n(), mean = mean(Score), 
            SD = sd(Score), SE = SD/sqrt(n), .groups = "drop")
        a_tests <- a %>% group_by(Module_ID) %>% group_modify(~wilcox_summary(.x, "Score")) %>% ungroup() %>% 
            mutate(p_adj_BH = p.adjust(p_value, "BH"))
        b_tests <- b %>% group_by(Feature) %>% group_modify(~wilcox_summary(.x, "Mass_log2")) %>% ungroup() %>% 
            mutate(p_adj_BH = p.adjust(p_value, "BH"))
        save_pdf <- function(plot, suffix, width, height) {
            ggsave(repo_output(repo_path(out_dir, paste0(key, "_", suffix, ".pdf")), "FigureS6_FigureS7_FigureS8_analysis.R"), 
                plot = plot, device = grDevices::pdf, width = width, height = height, useDingbats = FALSE)
        }
        save_pdf(make_d1_host(a, a_tests, label), "A_D1_GSVA", 5.7, 3.7)
        save_pdf(make_d1_biomass(b, b_tests, label), "B_D1_biomass", 5.2, 3.5)
        for (id in names(module_labels)) {
            save_pdf(make_plot(id, summary_df, trajectory_tests, timepoint_tests), paste0("C_trajectory_", 
                id), 3.1, 2.9)
        }
        tables <- list(trajectory_interaction_tests = trajectory_tests, trajectory_mean_SE = summary_df, 
            trajectory_timepoint_tests = timepoint_tests, D1_GSVA_tests = a_tests, D1_biomass_tests = b_tests, 
            host_plot_data = dat, biomass_plot_data = b)
        for (nm in names(tables)) write_csv(tables[[nm]], repo_output(repo_path(out_dir, paste0(key, 
            "_", nm, ".csv")), "FigureS6_FigureS7_FigureS8_analysis.R"))
        counts[[key]] <- bind_rows(a %>% distinct(HumanID, Outcome) %>% count(Outcome) %>% mutate(panel = "A_D1_GSVA"), 
            b %>% distinct(HumanID, Outcome) %>% count(Outcome) %>% mutate(panel = "B_D1_biomass"), dat %>% 
                distinct(HumanID, Outcome, Timepoint) %>% count(Timepoint, Outcome) %>% mutate(panel = "C_trajectory")) %>% 
            mutate(subgroup = key, .before = 1)
        message(analysis, " / ", key, ": 6 single-page PDFs written")
    }
    write_csv(bind_rows(counts), repo_output(repo_path(out_dir, "patient_counts.csv"), "FigureS6_FigureS7_FigureS8_analysis.R"))
}

writeLines(c("Source scripts: see the header of FigureS6_FigureS7_FigureS8_analysis.R.", paste("Original subgroup GSVA input:", 
    gsva_file), paste("Metadata:", meta_file), paste("Biomass:", biomass_file), paste("Infection strata:", 
    pathogen_file), "Subgroup analyses use the shared MASS module-score matrix; scores are subset without re-normalization.", 
    "Module display names, trajectory model, style and FDR annotation follow Figure3C_plot.R.", "C: observed mean +/- SE; Score ~ Outcome * categorical Timepoint + (1 | HumanID), ML fit.", 
    "C: joint 2-df Satterthwaite interaction F test; BH across 4 modules within each subgroup.", "C: timepoint two-sided Wilcoxon rank-sum; BH across 4 modules x 3 visits within each subgroup.", 
    "C: *** FDR<0.001; ** FDR<0.01; * FDR<0.05; ns otherwise.", "A/B: original unadjusted Wilcoxon P annotations retained; BH values also exported in CSV.", 
    "Infection subgroups are nonexclusive; bacterial/fungal-virus coinfections enter both.", "NP follows original code: PneumoniaTypeGroup == non-CAP.", 
    "Counts are recomputed per panel from available observations, never copied from screenshot labels.", 
    "Each folder contains 12 separate, single-page PDF panels; no PNG files generated.", "Trajectory subgroup identity is in the filename; module title/layout matches the main figure."), 
    repo_output(repo_path(out_root, "run_info.txt"), "FigureS6_FigureS7_FigureS8_analysis.R"))

write_csv(tibble(path = c(meta_file, biomass_file, gsva_file, pathogen_file), md5 = unname(tools::md5sum(repo_input(c(meta_file, 
    biomass_file, gsva_file, pathogen_file))))), repo_output(repo_path(out_root, "input_manifest.csv"), 
    "FigureS6_FigureS7_FigureS8_analysis.R"))

capture.output(sessionInfo(), file = repo_path(out_root, "sessionInfo.txt"))

message("Done: ", out_root)

