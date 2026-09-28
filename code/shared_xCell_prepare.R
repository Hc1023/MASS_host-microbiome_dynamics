# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
project_root <- "."

suppressPackageStartupMessages({
    library(tidyverse)
})

host_dir <- repo_path(project_root, ".")

mortality_dir <- repo_path(project_root, ".")

module_file <- repo_path(mortality_dir, "Inputs/selected_genes_by_module_wide.csv")

metadata_file <- repo_path(host_dir, "Inputs/1211_metadata.rdata")


frozen_score_file <- "Outputs/Figure3C_prepare_MASS_GSVA_scores_long_with_metadata.csv"

score_module_file <- repo_path(mortality_dir, "Inputs/selected_genes_by_module_wide.csv")

out_dir <- repo_path(mortality_dir, "Outputs")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(repo_exists(module_file), repo_exists(metadata_file), repo_exists(frozen_score_file), 
    repo_exists(score_module_file))

module_labels <- c(Module1_PRR_TLR_TNF_NFkB = "M1 PRR/TLR/TNF/NF-kB", Module2_Phagolysosome_Autophagy = "M2 Phagolysosome/autophagy", 
    Module3_IFN_I_Response = "M3 Type I interferon", Module4_Ribosome_Biogenesis = "M4 Ribosome biogenesis")

module_wide <- read.csv(repo_input(module_file), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", 
    "NA"))

if (!identical(names(module_wide), names(module_labels))) {
    stop("Unexpected fixed-module columns")
}

module_sets <- lapply(module_wide, function(x) {
    unique(trimws(x[!is.na(x) & trimws(x) != ""]))
})

score_module_wide <- read.csv(repo_input(score_module_file), check.names = FALSE, stringsAsFactors = FALSE, 
    na.strings = c("", "NA"))

score_module_sets <- lapply(score_module_wide, function(x) {
    unique(trimws(x[!is.na(x) & trimws(x) != ""]))
})

if (!identical(names(score_module_sets), names(module_sets)) || !all(mapply(setequal, score_module_sets, 
    module_sets))) {
    stop("Frozen MASS GSVA scores were generated from different module genes")
}

metadata_env <- new.env(parent = emptyenv())

load(repo_input(metadata_file), envir = metadata_env)



stopifnot(all(c("df_long", "meta") %in% ls(metadata_env)))

clinical_covariates <- c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

sample_meta <- metadata_env$df_long %>% filter(Timepoint %in% c("D1", "D4", "D7"), as.character(Mortality28d) %in% 
    c("0", "1")) %>% transmute(HumanID, SampleID, Timepoint = factor(as.character(Timepoint), levels = c("D1", 
    "D4", "D7")), Mortality28d = factor(as.character(Mortality28d), levels = c("0", "1"))) %>% distinct(SampleID, 
    .keep_all = TRUE) %>% left_join(metadata_env$meta %>% dplyr::select(HumanID, all_of(clinical_covariates)) %>% 
    distinct(HumanID, .keep_all = TRUE), by = "HumanID") %>% arrange(Timepoint, Mortality28d, HumanID)

if (anyDuplicated(sample_meta$SampleID)) stop("Duplicated SampleID")

if (anyNA(sample_meta[, clinical_covariates])) {
    stop("Missing clinical covariates in the analysis cohort")
}

xcell_cache <- "Inputs/MASS_D1_D4_D7_xCell_FPKM_scores.rds"

xcell_scores <- readRDS(repo_input(xcell_cache))

stopifnot(setequal(colnames(xcell_scores), sample_meta$SampleID))

xcell_scores <- xcell_scores[, sample_meta$SampleID, drop = FALSE]

xcell_cells <- c("Neutrophils", "Monocytes", "CD4+ T-cells", "CD8+ T-cells", "B-cells", "NK cells")

if (length(setdiff(xcell_cells, rownames(xcell_scores)))) {
    stop("Missing required xCell types: ", paste(setdiff(xcell_cells, rownames(xcell_scores)), collapse = ", "))
}

xcell_df <- as.data.frame(t(xcell_scores[xcell_cells, , drop = FALSE]), check.names = FALSE) %>% rownames_to_column("SampleID")

xcell_pca <- prcomp(xcell_df[xcell_cells], center = TRUE, scale. = TRUE)

n_xcell_pcs <- min(3L, ncol(xcell_pca$x))

xcell_pc_names <- paste0("xCell_PC", seq_len(n_xcell_pcs))

xcell_pc_scores <- as_tibble(xcell_pca$x[, seq_len(n_xcell_pcs), drop = FALSE]) %>% setNames(xcell_pc_names) %>% 
    mutate(SampleID = xcell_df$SampleID, .before = 1)

xcell_pc_variance <- 100 * xcell_pca$sdev^2/sum(xcell_pca$sdev^2)

xcell_pc_variance_table <- tibble(PC = paste0("PC", seq_along(xcell_pc_variance)), Variance_percent = xcell_pc_variance, 
    Cumulative_percent = cumsum(xcell_pc_variance)) %>% slice_head(n = n_xcell_pcs)

xcell_pc_loading_table <- as.data.frame(xcell_pca$rotation[, seq_len(n_xcell_pcs), drop = FALSE], check.names = FALSE) %>% 
    rownames_to_column("Cell_type")

frozen_module_names <- c(Module1_PRR_TLR_TNF_NFkB = "PRR/TLR/TNF/NF-kB signaling", Module2_Phagolysosome_Autophagy = "Phagolysosome/autophagy", 
    Module3_IFN_I_Response = "Type I interferon response", Module4_Ribosome_Biogenesis = "Ribosome biogenesis")

frozen_scores_long <- read.csv(repo_input(frozen_score_file), check.names = FALSE, stringsAsFactors = FALSE)

if (!setequal(unique(frozen_scores_long$Module), unname(frozen_module_names))) {
    stop("Frozen GSVA module names do not match the expected four modules")
}

module_score_df <- frozen_scores_long %>% dplyr::select(SampleID, Module, Score) %>% distinct() %>% pivot_wider(names_from = Module, 
    values_from = Score)

names(module_score_df)[match(unname(frozen_module_names), names(module_score_df))] <- names(frozen_module_names)

if (!setequal(module_score_df$SampleID, sample_meta$SampleID)) {
    stop("Frozen GSVA samples do not match the analysis cohort")
}

analysis_data <- sample_meta %>% left_join(xcell_df, by = "SampleID") %>% left_join(xcell_pc_scores, 
    by = "SampleID") %>% left_join(module_score_df, by = "SampleID") %>% mutate(Gender = factor(Gender), 
    CenterGroup = factor(CenterGroup), PneumoniaTypeGroup = factor(PneumoniaTypeGroup), Immunosuppression = factor(Immunosuppression), 
    MV = factor(MV), Age_z = as.numeric(scale(Age)), CCI_z = as.numeric(scale(CCI)), SOFA_24h_z = as.numeric(scale(SOFA_24h)))

for (pc in xcell_pc_names) {
    analysis_data[[pc]] <- as.numeric(scale(analysis_data[[pc]]))
}

module_z_names <- paste0(names(module_sets), "_z")

for (i in seq_along(module_sets)) {
    analysis_data[[module_z_names[i]]] <- as.numeric(scale(analysis_data[[names(module_sets)[i]]]))
}

write.csv(analysis_data, repo_output(repo_path(out_dir, "composition_adjusted_model_data.csv"), "shared_xCell_prepare.R"), 
    row.names = FALSE)

