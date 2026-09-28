# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(lmerTest)
    library(emmeans)
    library(RColorBrewer)
    library(scales)
})

root <- "."

host_root <- repo_path(dirname(root), ".")

out <- repo_path(root, "Outputs")

repo_dir(out, recursive = TRUE, showWarnings = FALSE)

input_env <- new.env()

load(repo_input(repo_path(root, "Inputs/260823_meta_model.rdata")), input_env)

microbe_env <- new.env()

load(repo_input(repo_path(host_root, "Inputs/1616_microbe.rdata")), microbe_env)

source("code/shared_sources.R")
wide <- module_overlay(input_env$meta_model)

microbe_mass <- as.matrix(microbe_env$data)

modules <- c("up1", "up2", "up3", "dw1")

module_labels <- setNames(c("PRR/TLR/NF-kB signaling", "Phagolysosomal/autophagy", "Type I IFN response", 
    "Ribosome biogenesis"), modules)

clinical <- c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

other_microbes <- c("EBV", "TTV", "Streptococcus", "Staphylococcus", "Prevotella", "Veillonella", "HPgV", 
    "Aspergillus", "Mycobacterium")

model_microbes <- c("HCMV", other_microbes)

raw_names <- setNames(ifelse(model_microbes == "EBV", "HHV-4", model_microbes), model_microbes)

timepoints <- c("D1", "D4", "D7")

long <- bind_rows(lapply(timepoints, function(tp) {
    z <- wide[, c("HumanID", clinical), drop = FALSE]
    z$SampleID <- wide[[tp]]
    z$Timepoint <- tp
    for (m in modules) z[[m]] <- wide[[paste0(tp, "_mod_", m)]]
    for (m in model_microbes) {
        current_name <- paste0(tp, "_mi_", ifelse(m == "EBV", "HHV.4", m))
        stopifnot(current_name %in% names(wide), raw_names[m] %in% rownames(microbe_mass))
        z[[m]] <- wide[[current_name]]
        reference <- log2(microbe_mass[raw_names[m], match(z$SampleID, colnames(microbe_mass))] + 1)
        ok <- complete.cases(z[[m]], reference)
        stopifnot(sum(ok) > 0, max(abs(z[[m]][ok] - reference[ok])) < 1e-10)
    }
    z[!is.na(z$SampleID) & nzchar(z$SampleID), ]
}))

long$Timepoint <- factor(long$Timepoint, levels = timepoints)

long$HumanID <- factor(long$HumanID)

stopifnot(!anyDuplicated(long[c("HumanID", "Timepoint")]))

fdr_stars <- function(x) ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", "")))

common_theme <- theme_bw(base_size = 10, base_family = "Helvetica") + theme(panel.grid = element_blank(), 
    axis.title = element_blank(), plot.title = element_text(face = "bold", size = 11), legend.title = element_text(size = 9), 
    legend.text = element_text(size = 8))

save_pdf <- function(plot, filename, width, height) {
    grDevices::pdf(repo_output(repo_path(out, filename), "Figure3DE_FigureS5C_analysis.R"), width = width, height = height, 
        family = "Helvetica", useDingbats = FALSE)
    print(plot)
    grDevices::dev.off()
}

prevalent <- rownames(microbe_mass)[rowSums(microbe_mass > 0) > 0.05 * ncol(microbe_mass)]

prevalent <- prevalent[order(rowSums(microbe_mass[prevalent, , drop = FALSE] > 0), decreasing = TRUE)]

all_microbes <- log2(microbe_mass[prevalent, , drop = FALSE] + 1)

cor_data <- long[, c("SampleID", modules)]

for (m in prevalent) cor_data[[m]] <- all_microbes[m, match(cor_data$SampleID, colnames(all_microbes))]

cor_results <- expand_grid(host = modules, microbe_raw = prevalent) %>% rowwise() %>% mutate(test = list(cor.test(cor_data[[host]], 
    cor_data[[microbe_raw]], method = "spearman", exact = FALSE)), rho = unname(test$estimate), P_raw = test$p.value) %>% 
    ungroup() %>% select(-test) %>% mutate(FDR_BH = p.adjust(P_raw, "BH"), logFDR = -log10(pmax(FDR_BH, 
    .Machine$double.xmin)), star = fdr_stars(FDR_BH), host = factor(host, levels = modules), microbe = ifelse(microbe_raw == 
    "HHV-4", "EBV", microbe_raw), microbe = factor(microbe, levels = rev(ifelse(prevalent == "HHV-4", 
    "EBV", prevalent))))

cor_limit <- max(0.3, max(abs(cor_results$rho), na.rm = TRUE))

p_cor <- ggplot(cor_results, aes(microbe, host)) + geom_point(aes(fill = rho, size = logFDR), shape = 21, 
    color = "black", stroke = 0.35, alpha = 0.9) + geom_text(aes(label = star), size = 2.7) + scale_x_discrete(limits = rev(levels(cor_results$microbe))) + 
    scale_y_discrete(limits = rev(modules), labels = module_labels) + scale_fill_gradientn(colors = brewer.pal(10, 
    "PuOr"), limits = c(-cor_limit, cor_limit), name = "Spearman rho") + scale_size_area(max_size = 9, 
    name = "-log10(FDR)") + labs(title = "Microbe-host module correlations") + common_theme + theme(axis.text.x = element_text(angle = 40, 
    hjust = 1), legend.position = "right", legend.box = "horizontal", plot.margin = margin(8, 18, 8, 
    8))

save_pdf(p_cor, "03_microbe_mod_cor.pdf", 7, 2.5)

emm_options(lmer.df = "kenward-roger")

lmm_results <- list()
models_current <- list()
rows_current <- list()

for (m in modules) {
    dat <- droplevels(long[complete.cases(long[, c(m, clinical, model_microbes)]), ])
    formula_lmm <- as.formula(paste(m, "~ HCMV * Timepoint +", paste(c(other_microbes, clinical), collapse = " + "), 
        "+ (1 | HumanID)"))
    fit <- lmer(formula_lmm, data = dat, REML = FALSE)
    models_current[[paste0("LMM_",m)]] <- fit
    rows_current[[paste0("LMM_",m)]] <- dat
    stopifnot(!lme4::isSingular(fit), is.null(fit@optinfo$conv$lme4$messages))
    em <- as.data.frame(summary(emtrends(fit, specs = "Timepoint", var = "HCMV", lmer.df = "kenward-roger"), 
        infer = c(TRUE, TRUE), adjust = "none"))
    lmm_results[[m]] <- data.frame(host = m, Timepoint = as.character(em$Timepoint), beta = em$HCMV.trend, 
        SE = em$SE, CI_low = em$lower.CL, CI_high = em$upper.CL, P_raw = em$p.value, N_samples=nrow(dat),N_patients=dplyr::n_distinct(dat$HumanID),N_D1=sum(dat$Timepoint=="D1"),N_D4=sum(dat$Timepoint=="D4"),N_D7=sum(dat$Timepoint=="D7"))
}

lmm_results <- bind_rows(lmm_results) %>% mutate(FDR_BH = p.adjust(P_raw, "BH"), logFDR = -log10(pmax(FDR_BH, 
    .Machine$double.xmin)), star = fdr_stars(FDR_BH), host = factor(host, levels = rev(modules)), Timepoint = factor(Timepoint, 
    levels = timepoints))

lmm_limit <- max(0.08, max(abs(lmm_results$beta)))

p_lmm <- ggplot(lmm_results, aes(Timepoint, host)) + geom_point(aes(fill = beta, size = logFDR), shape = 21, 
    color = "black", stroke = 0.35, alpha = 0.9) + geom_text(aes(label = star), size = 3) + scale_y_discrete(labels = module_labels) + 
    scale_fill_gradientn(colors = rev(brewer.pal(10, "PuOr")), limits = c(-lmm_limit, lmm_limit), name = "HCMV effect (beta)", 
        guide = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(1.8, "cm"))) + scale_size_area(max_size = 10, 
    name = "-log10(FDR)") + labs(title = "Adjusted contemporaneous HCMV-host coupling") + common_theme + 
    theme(plot.margin = margin(8, 18, 8, 8))

save_pdf(p_lmm, "03_microbe_mod_cor_CMV.pdf", 4.5, 2)

gradient_legend <- cowplot::get_legend(p_lmm)

pdf(repo_output("Outputs/260905_HCMV_host/HCMV_gradient_legend.pdf", "Figure3DE_FigureS5C_analysis.R"), width = 2, 
    height = 4)

grid::grid.draw(gradient_legend)

dev.off()

d1 <- long[long$Timepoint == "D1", c("HumanID", clinical, model_microbes, modules)]

d4 <- long[long$Timepoint == "D4", c("HumanID", modules)]

d4_models <- list()

for (m in modules) {
    baseline <- d1[, c("HumanID", clinical, model_microbes, m)]
    names(baseline)[names(baseline) == m] <- "baseline_module"
    future <- d4[, c("HumanID", m)]
    names(future)[2] <- "D4_module"
    dat <- droplevels(inner_join(baseline, future, by = "HumanID"))
    dat <- dat[complete.cases(dat), ]
    fit <- lm(as.formula(paste("D4_module ~ HCMV + baseline_module +", paste(c(other_microbes, clinical), 
        collapse = " + "))), data = dat)
    models_current[[paste0("AR_",m)]] <- fit
    rows_current[[paste0("AR_",m)]] <- dat
    co <- summary(fit)$coefficients["HCMV", ]
    ci <- confint(fit, "HCMV", level = 0.95)
    d4_models[[m]] <- data.frame(host = m, N = nobs(fit), beta = co["Estimate"], SE = co["Std. Error"], 
        CI_low = ci[1], CI_high = ci[2], P_raw = co["Pr(>|t|)"])
}

d4_results <- bind_rows(d4_models) %>% mutate(FDR_BH = p.adjust(P_raw, "BH"), star = fdr_stars(FDR_BH), 
    host = factor(host, levels = rev(modules)))

span <- diff(range(c(d4_results$CI_low, d4_results$CI_high, 0)))

p_d4 <- ggplot(d4_results, aes(beta, host)) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", 
    linewidth = 0.4) + geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.14, 
    linewidth = 0.6) + geom_point(size = 2.5) + geom_text(aes(x = CI_high + 0.04 * span, label = star), 
    hjust = 0, size = 3.5) + scale_y_discrete(labels = module_labels) + scale_x_continuous(expand = expansion(mult = c(0.1, 
    0.2))) + labs(title = "D1 HCMV and D4 host-response state", subtitle = "Adjusted for D1 module state, clinical covariates, and other D1 microbes", 
    x = "Adjusted beta for D1 HCMV burden", y = NULL, caption = "Points: adjusted beta; bars: 95% CI; stars: BH FDR across four module tests.") + 
    theme_bw() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.minor = element_blank(), 
    plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9), plot.caption = element_text(size = 8, 
        hjust = 0), plot.margin = margin(8, 18, 8, 8))

save_pdf(p_d4, "D4_module_D1_HCMV_autoregressive.pdf", 4, 2.8)

cat("\nMicrobe-module Spearman tests (one BH family):", nrow(cor_results), "\n")

cat("Contemporaneous HCMV LMM tests (one BH family):", nrow(lmm_results), "\n")

cat("\nD1 HCMV -> D4 autoregressive results (BH across four modules):\n")

print(d4_results %>% transmute(Module = module_labels[as.character(host)], N, beta, SE, CI_low, CI_high, 
    P_raw, FDR_BH), row.names = FALSE)

stopifnot(nrow(lmm_results)==12, nrow(d4_results)==4)
write.csv(cor_results,"Outputs/FigureS5C_correlations.csv",row.names=FALSE)
write.csv(lmm_results,"Outputs/Figure3D_HCMV_effects.csv",row.names=FALSE)
write.csv(d4_results,"Outputs/Figure3E_HCMV_effects.csv",row.names=FALSE)

write.csv(lmm_results,"Outputs/Figure3D_HCMV_effects_complete.csv",row.names=FALSE)
saveRDS(list(models=models_current,rows=rows_current),"Outputs/HCMV_models.rds")
