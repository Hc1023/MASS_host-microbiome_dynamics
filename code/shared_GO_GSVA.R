# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(magrittr)

library(clusterProfiler)

library(limma)

library(edgeR)

library(AnnotationDbi)

library(org.Hs.eg.db)

library(ggplot2)

library(ggpubr)

library(scales)

library(stringr)

library(ComplexHeatmap)

library(circlize)

project_dir <- "."

invisible(NULL)

out_dir <- repo_path(".")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

load(repo_input("Inputs/1211_metadata.rdata"))

load(repo_input("Inputs/1211_transcriptome.rdata"))

metadata_analysis = df_long

vars = c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))

metadata_analysis %<>% left_join(tmp, by = "HumanID")

counts_analysis = counts %>% dplyr::select(all_of(metadata_analysis$SampleID))

stopifnot(identical(metadata_analysis$SampleID, colnames(counts_analysis)))

library(GSVA)

library(BiocParallel)

msigdb_expected_version <- "2025.1.Hs"

geneset_file <- normalizePath(repo_input(repo_path(project_dir, "Inputs/c5.go.bp.v2025.1.Hs.symbols.gmt")), 
    mustWork = TRUE)

if (!grepl("c5\\.go\\.bp\\.v2025\\.1\\.Hs\\.symbols\\.gmt$", geneset_file)) {
    stop("The selected GMT is not the pinned MSigDB C5 GO:BP 2025.1.Hs file")
}

bp_tbl <- clusterProfiler::read.gmt(repo_input(geneset_file))

bp_tbl <- bp_tbl %>% distinct(term, gene)

bp_sets <- split(bp_tbl$gene, bp_tbl$term)

dge_all <- DGEList(counts = as.matrix(counts))

dge_all <- calcNormFactors(dge_all, method = "TMM")

expr_logcpm <- cpm(dge_all, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

gene_annot <- gene_attr %>% as.data.frame() %>% rownames_to_column("row_id") %>% transmute(ENSEMBL = if_else(is.na(gene_id) | 
    gene_id == "", row_id, gene_id), SYMBOL = SYMBOL) %>% distinct(ENSEMBL, .keep_all = TRUE)

symbol <- gene_annot$SYMBOL[match(rownames(expr_logcpm), gene_annot$ENSEMBL)]

keep_gene <- !is.na(symbol) & symbol != "" & !grepl("^ENSG", symbol)

expr_symbol <- expr_logcpm[keep_gene, , drop = FALSE]

rownames(expr_symbol) <- symbol[keep_gene]

expr_symbol <- limma::avereps(expr_symbol, ID = rownames(expr_symbol))

stopifnot(identical(colnames(expr_symbol), colnames(counts)), !anyDuplicated(rownames(expr_symbol)), 
    all(is.finite(expr_symbol)))

gsva_param <- gsvaParam(exprData = expr_symbol, geneSets = bp_sets, kcdf = "Gaussian", minSize = 10, 
    maxSize = 500)

gsva_scores <- gsva(gsva_param, BPPARAM = SerialParam(progressbar = TRUE), verbose = TRUE)

stopifnot(identical(colnames(gsva_scores), colnames(counts)))

gsva_by_sample <- as.data.frame(t(gsva_scores), check.names = FALSE) %>% rownames_to_column("SampleID") %>% 
    left_join(metadata_analysis, by = "SampleID")

saveRDS(gsva_scores, repo_output(repo_path(out_dir, "MSigDB_2025.1.Hs_C5_GO_BP_GSVA_all_samples.rds"), 
    "shared_GO_GSVA.R"))

write.csv(gsva_scores, repo_output(repo_path(out_dir, "MSigDB_2025.1.Hs_C5_GO_BP_GSVA_all_samples.csv"), 
    "shared_GO_GSVA.R"), quote = FALSE)

write.csv(gsva_by_sample, repo_output(repo_path(out_dir, "MSigDB_2025.1.Hs_C5_GO_BP_GSVA_all_samples_with_metadata.csv"), 
    "shared_GO_GSVA.R"), row.names = FALSE, quote = FALSE)

write.csv(bp_tbl, repo_output(repo_path(out_dir, "MSigDB_2025.1.Hs_C5_GO_BP_gene_sets.csv"), "shared_GO_GSVA.R"), 
    row.names = FALSE, quote = FALSE)

file.copy(geneset_file, repo_path(out_dir, basename(geneset_file)), overwrite = TRUE)

run_info <- c(paste("MSigDB version:", msigdb_expected_version), paste("GMT source:", geneset_file), 
    paste("GMT MD5:", unname(tools::md5sum(repo_input(geneset_file)))), paste("GSVA version:", packageVersion("GSVA")), 
    paste("Number of samples:", ncol(gsva_scores)), paste("Number of scored gene sets:", nrow(gsva_scores)), 
    paste("Number of expression genes after SYMBOL mapping:", nrow(expr_symbol)), capture.output(sessionInfo()))

writeLines(run_info, repo_output(repo_path(out_dir, "GSVA_run_info.txt"), "shared_GO_GSVA.R"))

