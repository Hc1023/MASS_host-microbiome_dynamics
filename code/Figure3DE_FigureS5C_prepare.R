# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(edgeR)
    library(limma)
    library(ConsensusTranscriptomicSubtype)
})

project_dir <- "."

main_dir <- "."

meta_model_file <- repo_path(project_dir, "Inputs/1616_meta_model.rdata")

sofa_file <- paste0(".", "Inputs/7.Sofa_apa_info.csv")

metadata_file <- repo_path(project_dir, "Inputs/1211_metadata.rdata")

transcriptome_file <- repo_path(project_dir, "Inputs/1211_transcriptome.rdata")

module_file <- repo_path(main_dir, "Inputs/selected_genes_by_module_wide.csv")

output_file <- repo_path(main_dir, "Inputs/260823_meta_model.rdata")

required_files <- c(meta_model_file, sofa_file, metadata_file, transcriptome_file, module_file)

if (!all(repo_exists(required_files))) {
    stop("Missing input file(s): ", paste(required_files[!repo_exists(required_files)], collapse = ", "))
}

model_env <- new.env(parent = emptyenv())

load(repo_input(meta_model_file), envir = model_env)

if (!"meta_model" %in% ls(model_env)) {
    stop("Object 'meta_model' is absent from ", meta_model_file)
}

meta_model <- model_env$meta_model

if (!"HumanID" %in% names(meta_model) || anyDuplicated(meta_model$HumanID)) {
    stop("meta_model must contain one unique row per HumanID")
}

sofa_96h <- readr::read_csv(repo_input(sofa_file), show_col_types = FALSE, name_repair = "minimal", locale = readr::locale(encoding = "UTF-8")) %>% 
    dplyr::select(HumanID, SOFA_96h) %>% distinct()

if (anyDuplicated(sofa_96h$HumanID)) {
    stop("SOFA input contains more than one SOFA_96h value for a HumanID")
}

meta_model <- meta_model %>% dplyr::select(-any_of("SOFA_96h")) %>% left_join(sofa_96h, by = "HumanID")

module_wide <- read.csv(repo_input(module_file), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", 
    "NA"))

module_genes_all <- module_wide %>% pivot_longer(everything(), names_to = "Module", values_to = "Gene") %>% 
    transmute(Gene = trimws(Gene)) %>% filter(!is.na(Gene), Gene != "") %>% distinct(Gene) %>% pull(Gene)

data("exp_core_g", package = "ConsensusTranscriptomicSubtype", envir = environment())

cts_core_ensembl <- unique(rownames(exp_core_g))

input_env <- new.env(parent = emptyenv())

load(repo_input(metadata_file), envir = input_env)

load(repo_input(transcriptome_file), envir = input_env)

required_objects <- c("df_long", "counts", "gene_attr")

if (!all(required_objects %in% ls(input_env))) {
    stop("Metadata/transcriptome inputs must contain: ", paste(required_objects, collapse = ", "))
}

sample_meta <- input_env$df_long %>% filter(Timepoint %in% c("D1", "D4", "D7"), !is.na(SampleID)) %>% 
    transmute(HumanID = as.character(HumanID), SampleID = as.character(SampleID), Timepoint = factor(as.character(Timepoint), 
        levels = c("D1", "D4", "D7"))) %>% distinct(HumanID, Timepoint, .keep_all = TRUE) %>% arrange(HumanID, 
    Timepoint)

if (anyDuplicated(sample_meta$SampleID)) stop("Duplicated SampleID in metadata")

if (!all(sample_meta$SampleID %in% colnames(input_env$counts))) {
    stop("Some D1/D4/D7 samples are absent from the count matrix")
}

counts_analysis <- as.matrix(input_env$counts[, sample_meta$SampleID, drop = FALSE])

dge <- DGEList(counts = counts_analysis)

dge <- calcNormFactors(dge, method = "TMM")

log_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

gene_attr <- input_env$gene_attr

if (!"SYMBOL" %in% colnames(gene_attr) || !all(rownames(log_cpm) %in% rownames(gene_attr))) {
    stop("gene_attr does not provide SYMBOL for every count-matrix row")
}

gene_symbol <- as.character(gene_attr[rownames(log_cpm), "SYMBOL"])

valid_symbol <- !is.na(gene_symbol) & gene_symbol != "" & !str_detect(gene_symbol, "^ENSG")

expr_symbol <- limma::avereps(log_cpm[valid_symbol, , drop = FALSE], ID = gene_symbol[valid_symbol])

cts_core_symbols <- gene_attr[intersect(cts_core_ensembl, rownames(gene_attr)), "SYMBOL"] %>% as.character() %>% 
    .[!is.na(.) & . != "" & !str_detect(., "^ENSG")] %>% unique()

cts_genes_detected <- intersect(cts_core_symbols, rownames(expr_symbol))

missing_cts_genes <- setdiff(cts_core_symbols, rownames(expr_symbol))

if (length(cts_genes_detected) == 0L) {
    stop("None of the ConsensusTranscriptomicSubtype core genes was detected")
}

module_genes <- setdiff(module_genes_all, cts_core_symbols)

missing_module_genes <- setdiff(module_genes, rownames(expr_symbol))

module_genes_detected <- intersect(module_genes, rownames(expr_symbol))

if (length(module_genes_detected) == 0L) {
    stop("None of the selected module genes was detected")
}

selected_genes <- c(cts_genes_detected, module_genes_detected)

gene_group <- c(setNames(rep("CTSg", length(cts_genes_detected)), cts_genes_detected), setNames(rep("DYNg", 
    length(module_genes_detected)), module_genes_detected))

gene_expr_wide <- as.data.frame(t(expr_symbol[selected_genes, sample_meta$SampleID, drop = FALSE]), check.names = FALSE) %>% 
    rownames_to_column("SampleID") %>% left_join(sample_meta, by = "SampleID") %>% pivot_longer(cols = all_of(selected_genes), 
    names_to = "Gene", values_to = "Expression") %>% mutate(Feature = paste0(Timepoint, "_", gene_group[Gene], 
    "_", Gene)) %>% dplyr::select(HumanID, Feature, Expression) %>% pivot_wider(names_from = Feature, 
    values_from = Expression)

if (anyDuplicated(gene_expr_wide$HumanID)) {
    stop("Gene-expression table contains duplicated HumanID rows")
}

meta_model <- meta_model %>% dplyr::select(-matches("^D[147]_(CTSg|DYNg)_")) %>% left_join(gene_expr_wide, 
    by = "HumanID")

model_prepare_info <- list(expression_scale = "TMM-normalized log2-CPM (prior.count = 1)", timepoints = c("D1", 
    "D4", "D7"), cts_core_ensembl_defined = cts_core_ensembl, cts_genes_defined = cts_core_symbols, cts_genes_detected = cts_genes_detected, 
    cts_genes_missing = missing_cts_genes, module_genes_defined = module_genes_all, module_genes_detected = module_genes_detected, 
    module_genes_missing = missing_module_genes, module_genes_excluded_as_CTS = intersect(module_genes_all, 
        cts_core_symbols), source_meta_model = meta_model_file, source_module_file = module_file)

save(meta_model, model_prepare_info, file = repo_output(output_file, "Figure3DE_FigureS5C_prepare.R"))

