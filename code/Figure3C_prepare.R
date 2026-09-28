# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(edgeR)

library(limma)

library(GSVA)

library(ggpubr)

library(rstatix)

library(scales)

project_dir <- "."

main_dir <- "."

metadata_file <- repo_path(project_dir, "Inputs/1211_metadata.rdata")

transcriptome_file <- repo_path(project_dir, "Inputs/1211_transcriptome.rdata")

module_file <- repo_path(main_dir, "Inputs/selected_genes_by_module_wide.csv")

out_dir <- repo_path(main_dir, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(metadata_file, transcriptome_file, module_file)

if (!all(repo_exists(required_files))) {
    stop("Missing input file(s): ", paste(required_files[!repo_exists(required_files)], collapse = ", "))
}

module_wide <- read.csv(repo_input(module_file), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", 
    "NA"))

expected_modules <- c("Module1_PRR_TLR_TNF_NFkB", "Module2_Phagolysosome_Autophagy", "Module3_IFN_I_Response", 
    "Module4_Ribosome_Biogenesis")

if (!identical(names(module_wide), expected_modules)) {
    stop("Expected these four module columns in order: ", paste(expected_modules, collapse = ", "))
}

fixed_gene_sets <- lapply(module_wide, function(x) {
    unique(trimws(x[!is.na(x) & trimws(x) != ""]))
})

if (any(lengths(fixed_gene_sets) == 0L)) stop("At least one module is empty")

gs_list <- setNames(fixed_gene_sets, c("up1", "up2", "up3", "dw1"))

MASS_env <- new.env(parent = globalenv())

load(repo_input(metadata_file), envir = MASS_env)

load(repo_input(transcriptome_file), envir = MASS_env)

required_objects <- c("df_long", "counts", "gene_attr")

if (!all(required_objects %in% ls(MASS_env))) {
    stop("MASS cohort inputs must contain: ", paste(required_objects, collapse = ", "))
}

metadata_analysis <- MASS_env$df_long %>% filter(Timepoint %in% c("D1", "D4", "D7"), !is.na(Mortality28d)) %>% 
    mutate(Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7")), Mortality28d = factor(as.character(Mortality28d), 
        levels = c("0", "1"))) %>% arrange(Timepoint, HumanID)

if (anyDuplicated(metadata_analysis$SampleID)) stop("Duplicated SampleID")

if (!all(metadata_analysis$SampleID %in% colnames(MASS_env$counts))) {
    stop("Some MASS cohort samples are absent from the count matrix")
}

counts_analysis <- as.matrix(MASS_env$counts[, metadata_analysis$SampleID, drop = FALSE])

stopifnot(identical(metadata_analysis$SampleID, colnames(counts_analysis)))

dge <- DGEList(counts = counts_analysis)

dge <- calcNormFactors(dge, method = "TMM")

expr <- cpm(dge, log = TRUE, prior.count = 1)

gene_attr <- MASS_env$gene_attr

sym <- gene_attr[rownames(expr), "SYMBOL"]

is_valid_symbol <- !is.na(sym) & sym != "" & !grepl("^ENSG", sym)

expr_sym <- limma::avereps(expr[is_valid_symbol, , drop = FALSE], ID = sym[is_valid_symbol])

stopifnot(!anyDuplicated(rownames(expr_sym)), identical(colnames(expr_sym), metadata_analysis$SampleID))

module_coverage <- imap_dfr(gs_list, function(genes, module) {
    tibble(Module = module, gene = genes, detected_in_MASS = genes %in% rownames(expr_sym))
})

coverage_summary <- module_coverage %>% group_by(Module) %>% summarise(n_defined = n(), n_detected = sum(detected_in_MASS), 
    fraction_detected = mean(detected_in_MASS), .groups = "drop")

if (any(coverage_summary$n_detected < 10L)) {
    stop("Fewer than 10 detected genes in module(s): ", paste(coverage_summary$Module[coverage_summary$n_detected < 
        10L], collapse = ", "))
}

gsva_param <- gsvaParam(expr = as.matrix(expr_sym), geneSets = gs_list, kcdf = "Gaussian", minSize = 10, 
    maxSize = 500)

gsva_mat <- gsva(gsva_param, verbose = FALSE)

stopifnot(identical(colnames(gsva_mat), metadata_analysis$SampleID))

meta_gsva <- metadata_analysis %>% dplyr::select(HumanID, Timepoint, SampleID, Mortality28d)

plot_df <- bind_cols(meta_gsva, as_tibble(t(gsva_mat), .name_repair = "minimal"))

plot_long <- plot_df %>% pivot_longer(cols = -c(HumanID, Timepoint, SampleID, Mortality28d), names_to = "Module", 
    values_to = "Score") %>% mutate(Module = factor(Module, levels = c("up1", "up2", "up3", "dw1"), labels = c(Module1_PRR_TLR_TNF_NFkB = "PRR/TLR/TNF/NF-kB signaling", 
    Module2_Phagolysosome_Autophagy = "Phagolysosome/autophagy", Module3_IFN_I_Response = "Type I interferon response", 
    Module4_Ribosome_Biogenesis = "Ribosome biogenesis")))

wilcox_df <- plot_long %>% group_by(Module, Timepoint) %>% wilcox_test(Score ~ Mortality28d) %>% ungroup() %>% 
    mutate(p_value = p, p_adj_BH = p.adjust(p_value, method = "BH"))

summary_df <- plot_long %>% group_by(Module, Timepoint, Mortality28d) %>% summarise(n = sum(!is.na(Score)), 
    mean = mean(Score, na.rm = TRUE), SD = sd(Score, na.rm = TRUE), SE = SD/sqrt(n), .groups = "drop") %>% 
    mutate(Outcome = factor(Mortality28d, levels = c("0", "1"), labels = c("Survival", "Mortality")))

format_test_value <- function(x) {
    ifelse(x < 0.001, "<0.001", paste0("=", formatC(x, digits = 3, format = "f")))
}

annotation_y <- summary_df %>% group_by(Module, Timepoint) %>% summarise(y_base = max(mean + SE, na.rm = TRUE), 
    .groups = "drop") %>% left_join(plot_long %>% group_by(Module) %>% summarise(score_range = diff(range(Score, 
    na.rm = TRUE)), .groups = "drop"), by = "Module") %>% mutate(y.position = y_base + pmax(0.08 * score_range, 
    0.012)) %>% dplyr::select(Module, Timepoint, y.position)

wilcox_df <- wilcox_df %>% left_join(annotation_y, by = c("Module", "Timepoint")) %>% mutate(annotation = if_else(!is.na(p_value) & 
    p_value < 0.05, paste0("P", format_test_value(p_value), "\nFDR", format_test_value(p_adj_BH)), ""))

wilcox_output <- wilcox_df %>% transmute(Module, Timepoint, n_survival = n1, n_mortality = n2, statistic, 
    p_value, p_adj_BH)

p <- ggplot(summary_df, aes(x = Timepoint, y = mean, color = Outcome, group = Outcome)) + annotate("rect", 
    xmin = 1, xmax = 2, ymin = -Inf, ymax = Inf, fill = "grey92", color = NA) + geom_line(linewidth = 1.05) + 
    geom_point(size = 2.7) + geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.12, linewidth = 0.55) + 
    stat_pvalue_manual(filter(wilcox_df, annotation != "", Timepoint != "D7"), label = "annotation", 
        x = "Timepoint", y.position = "y.position", tip.length = 0, size = 2.7, inherit.aes = FALSE) + 
    stat_pvalue_manual(filter(wilcox_df, annotation != "", Timepoint == "D7"), label = "annotation", 
        x = "Timepoint", y.position = "y.position", tip.length = 0, size = 2.7, hjust = 1, inherit.aes = FALSE) + 
    scale_color_manual(values = c(Survival = "#4575B4", Mortality = "#D73027"), name = NULL, drop = FALSE) + 
    facet_wrap(~Module, scales = "free_y", ncol = 2) + scale_x_discrete(expand = expansion(add = c(0.12, 
    0.12))) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) + labs(x = NULL, y = "GSVA score (mean +/- SE)", 
    caption = "Labels show two-sided Wilcoxon P values and BH-FDR across 12 module-by-timepoint tests; P>=0.05 is not shown.") + 
    theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(), 
    legend.position = "top", axis.text.x = element_text(face = "bold"), plot.caption = element_text(size = 7.5, 
        hjust = 0))

write.csv(module_coverage, repo_output(repo_path(out_dir, "MASS_module_gene_coverage.csv"), "Figure3C_prepare.R"), 
    row.names = FALSE)

write.csv(coverage_summary, repo_output(repo_path(out_dir, "MASS_module_gene_coverage_summary.csv"), 
    "Figure3C_prepare.R"), row.names = FALSE)

write.csv(tibble(Module = rownames(gsva_mat)) %>% bind_cols(as_tibble(gsva_mat, .name_repair = "minimal")), 
    repo_output(repo_path(out_dir, "MASS_GSVA_scores_wide.csv"), "Figure3C_prepare.R"), row.names = FALSE)

write.csv(plot_long, repo_output(repo_path(out_dir, "MASS_GSVA_scores_long_with_metadata.csv"), "Figure3C_prepare.R"), 
    row.names = FALSE)

write.csv(summary_df, repo_output(repo_path(out_dir, "MASS_GSVA_trajectory_mean_SE.csv"), "Figure3C_prepare.R"), 
    row.names = FALSE)

write.csv(wilcox_output, repo_output(repo_path(out_dir, "MASS_GSVA_wilcoxon_D1_D4_D7.csv"), "Figure3C_prepare.R"), 
    row.names = FALSE)

saveRDS(gsva_mat, repo_output(repo_path(out_dir, "MASS_GSVA_scores.rds"), "Figure3C_prepare.R"))

pdf(repo_output(repo_path(out_dir, "MASS_GSVA_Survival_Mortality_trajectories.pdf"), "Figure3C_prepare.R"), 
    width = 6.8, height = 5.4)

print(p)

dev.off()

ggsave(repo_output(repo_path(out_dir, "MASS_GSVA_Survival_Mortality_trajectories.png"), "Figure3C_prepare.R"), 
    p, width = 6.8, height = 5.4, dpi = 300)

writeLines(c(paste("Metadata input:", metadata_file), paste("Transcriptome input:", transcriptome_file), 
    paste("Fixed module input:", module_file), paste("Samples:", nrow(metadata_analysis)), paste("Patients:", 
        n_distinct(metadata_analysis$HumanID)), "Expression: TMM-normalized log2-CPM (prior.count=1)", 
    "Duplicate symbols: averaged with limma::avereps", "GSVA: Gaussian kernel; minSize=10; maxSize=500", 
    "Tests: two-sided Wilcoxon rank-sum, Survival versus Mortality at each time point", "BH-FDR: adjusted across all 12 module-by-timepoint tests"), 
    repo_output(repo_path(out_dir, "MASS_GSVA_run_info.txt"), "Figure3C_prepare.R"))

message("Done. MASS cohort GSVA results are in: ", out_dir)

