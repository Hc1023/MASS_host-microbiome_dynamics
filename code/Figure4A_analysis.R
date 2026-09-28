# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(tidyverse)
    library(lmerTest)
})

root <- "."

input <- repo_path(root, "Outputs")

out <- repo_path(root, "Outputs")

repo_dir(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(repo_input(repo_path(input, "Inputs/composition_adjusted_model_data.csv")), check.names = FALSE)

keys <- c("Module1_PRR_TLR_TNF_NFkB", "Module2_Phagolysosome_Autophagy", "Module3_IFN_I_Response", "Module4_Ribosome_Biogenesis")

labels <- setNames(c("PRR/TLR/NF-kB signaling", "Phagolysosomal/autophagy", "Type I IFN response", "Ribosome biogenesis"), 
    keys)

clinical <- c("Gender", "Age_z", "CenterGroup", "PneumoniaTypeGroup", "CCI_z", "SOFA_24h_z", "Immunosuppression", 
    "MV")

pcs <- paste0("xCell_PC", 1:3)

required <- c("HumanID", "SampleID", "Mortality28d", "Timepoint", clinical, pcs, paste0(keys, "_z"))

stopifnot(all(required %in% names(d)), !anyNA(d[, required]), !anyDuplicated(d$SampleID))

d <- d[, required]

d$Mortality28d <- factor(d$Mortality28d, levels = c(0, 1))

d$Timepoint <- factor(d$Timepoint, levels = c("D1", "D4", "D7"))

for (nm in c("HumanID", "Gender", "CenterGroup", "PneumoniaTypeGroup", "Immunosuppression", "MV")) d[[nm]] <- factor(d[[nm]])

control <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

fits <- list()

coefficients <- list()

omnibus <- list()

diagnostics <- list()

for (key in keys) for (model in c("Base", "xCell_PC_adjusted")) {
    extra <- c(clinical, if (model == "xCell_PC_adjusted") pcs)
    ff <- as.formula(paste(paste0(key, "_z"), "~ Mortality28d * Timepoint +", paste(extra, collapse = " + "), 
        "+ (1 | HumanID)"))
    rf <- as.formula(paste(paste0(key, "_z"), "~ Mortality28d + Timepoint +", paste(extra, collapse = " + "), 
        "+ (1 | HumanID)"))
    message("Fitting ", key, ": ", model)
    full <- lmerTest::lmer(ff, data = d, REML = FALSE, control = control)
    reduced <- lmerTest::lmer(rf, data = d, REML = FALSE, control = control)
    id <- paste(key, model, sep = ":")
    fits[[id]] <- list(full = full, reduced = reduced)
    tab <- as.data.frame(summary(full)$coefficients) %>% rownames_to_column("Term")
    coefficients[[id]] <- tab %>% filter(grepl("Mortality28d1:TimepointD[47]", Term)) %>% transmute(Module_key = key, 
        Module = unname(labels[key]), Model = model, Interval = if_else(grepl("D4", Term), "D4-D1", "D7-D1"), 
        Beta_SD = Estimate, SE = `Std. Error`, Lower95 = Estimate - 1.96 * SE, Upper95 = Estimate + 1.96 * 
            SE, p_value = `Pr(>|t|)`)
    lr <- anova(reduced, full)
    omnibus[[id]] <- tibble(Module_key = key, Module = unname(labels[key]), Model = model, LRT_chisq = lr$Chisq[2], 
        df = lr$Df[2], p_value = lr$`Pr(>Chisq)`[2])
    diagnostics[[id]] <- tibble(Module_key = key, Model = model, n_samples = nobs(full), n_patients = n_distinct(d$HumanID), 
        singular = isSingular(full), convergence_message = paste(full@optinfo$conv$lme4$messages, collapse = "; "))
}

coefficients <- bind_rows(coefficients) %>% group_by(Model, Interval) %>% mutate(FDR = p.adjust(p_value, 
    "BH")) %>% ungroup()

omnibus <- bind_rows(omnibus) %>% group_by(Model) %>% mutate(FDR = p.adjust(p_value, "BH")) %>% ungroup()

write.csv(coefficients, repo_output(repo_path(out, "mixed_model_interaction_coefficients.csv"), "Figure4A_analysis.R"), 
    row.names = FALSE)

write.csv(omnibus, repo_output(repo_path(out, "mixed_model_omnibus_interaction_tests.csv"), "Figure4A_analysis.R"), 
    row.names = FALSE)

write.csv(bind_rows(diagnostics), repo_output(repo_path(out, "model_diagnostics.csv"), "Figure4A_analysis.R"), 
    row.names = FALSE)

saveRDS(fits, repo_output(repo_path(out, "xCell_adjusted_mixed_models.rds"), "Figure4A_analysis.R"))

plot_data <- coefficients %>% mutate(Module = factor(Module, levels = unname(labels)), Interval = factor(Interval, 
    levels = c("D4-D1", "D7-D1")), Model = factor(Model, levels = c("Base", "xCell_PC_adjusted")))

p <- ggplot(plot_data, aes(Beta_SD, Model, color = Model)) + geom_vline(xintercept = 0, color = "grey55", 
    linetype = 2) + geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = 0.18, 
    linewidth = 0.6) + geom_point(size = 1.9) + facet_grid(Module ~ Interval, scales = "free_x", labeller = labeller(Module = c(`PRR/TLR/NF-kB signaling` = "PRR/TLR/NF-kB\nsignaling", 
    `Phagolysosomal/autophagy` = "Phagolysosomal/\nautophagy", `Type I IFN response` = "Type I IFN\nresponse", 
    `Ribosome biogenesis` = "Ribosome\nbiogenesis"))) + scale_color_manual(values = c(Base = "grey35", 
    xCell_PC_adjusted = "#0072B2"), labels = c(Base = "Base", xCell_PC_adjusted = "xCell PC-adjusted")) + 
    scale_y_discrete(labels = c(Base = "Base", xCell_PC_adjusted = "xCell PC-adjusted")) + labs(x = "Interaction coefficient (module-score SD)", 
    y = NULL, subtitle = paste0("Mortality-by-time difference-in-differences; 95% CI; n = ", nrow(d))) + 
    theme_bw(base_size = 8) + theme(legend.position = "none", strip.text = element_text(face = "bold"), 
    strip.text.y = element_text(angle = 270, size = 6), plot.subtitle = element_text(size = 7), panel.spacing = grid::unit(2, 
        "mm"), plot.margin = margin(3, 3, 3, 3), panel.grid.minor = element_blank())

p

ggsave(repo_output(repo_path(out, "03_mixed_model_interaction_forest.pdf"), "Figure4A_analysis.R"), 
    p, width = 5, height = 3.7)

writeLines(c("Source: saved composition_adjusted_model_data.csv from 260822_2_mod_genes_cell.", "Existing xCell PC1-PC3 and standardized module scores reused; no deconvolution or GSVA rerun.", 
    "Base: mortality * categorical time + eight clinical covariates + patient random intercept.", "xCell-adjusted: Base + existing standardized xCell PC1-PC3.", 
    "ML fitting; coefficient P: Satterthwaite; intervals: estimate +/- 1.96 SE (original method).", "Omnibus: full versus no-interaction likelihood ratio test, 2 df.", 
    "BH across four modules separately for each model/interval; omnibus BH separately per model."), repo_output(repo_path(out, 
    "methods.txt"), "Figure4A_analysis.R"))

capture.output(sessionInfo(), file = repo_path(out, "sessionInfo.txt"))

