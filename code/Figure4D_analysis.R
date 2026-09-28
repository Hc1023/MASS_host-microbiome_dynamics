# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(tidyverse)
    library(edgeR)
    library(limma)
    library(fgsea)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

root <- "."

earli <- "."

out <- repo_path(root, "Outputs")

repo_dir(out, recursive = TRUE, showWarnings = FALSE)

gmt <- "Inputs/c5.go.bp.v2025.1.Hs.symbols.gmt"

micro <- read.csv(repo_input(repo_path(earli, "Inputs/microbe.csv")), check.names = FALSE)

counts <- read.csv(repo_input(repo_path(earli, "Inputs/earli_counts_kallisto_mortality.csv")), row.names = 1, 
    check.names = FALSE)

meta <- read.csv(repo_input(repo_path(earli, "Inputs/Cleaned_Metadata_070324.csv")))

meta$Barcode <- paste0("EARLI_", meta$Barcode)

ids <- intersect(intersect(micro$Barcode, colnames(counts)), meta$Barcode)

stopifnot(!anyDuplicated(micro$Barcode), !anyDuplicated(meta$Barcode))

md <- meta[match(ids, meta$Barcode), ]

rownames(md) <- ids

abundance <- t(as.matrix(micro[match(ids, micro$Barcode), -(1:3)]))

colnames(abundance) <- ids

prevalence <- tibble(Pathogen = rownames(abundance), n_detected = rowSums(abundance > 0), n_samples = length(ids)) %>% 
    mutate(prevalence = n_detected/n_samples, retained = prevalence > 0.05)

model_taxa <- prevalence %>% filter(retained) %>% arrange(desc(n_detected)) %>% pull(Pathogen)

taxa <- c("HHV5", "HHV4")

stopifnot(all(taxa %in% model_taxa))

abundance <- abundance[model_taxa, , drop = FALSE]

vars <- c("Gender", "Age", "Group", "APACHEIII", "Immunocomp_manuscript", "Intubated")

md[c("Gender", "Group", "Immunocomp_manuscript", "Intubated")] <- lapply(md[c("Gender", "Group", "Immunocomp_manuscript", 
    "Intubated")], factor)

stopifnot(!anyNA(md[, vars]), all(is.finite(abundance)))

md <- cbind(md, t(log2(abundance + 1)))

design <- model.matrix(reformulate(c(model_taxa, vars)), md)

stopifnot(nrow(design) == length(ids), qr(design)$rank == ncol(design))

cts <- as.matrix(counts[, ids, drop = FALSE])

cts <- matrix(as.integer(cts), nrow = nrow(cts), dimnames = dimnames(cts))

stopifnot(identical(colnames(cts), rownames(md)), !anyNA(cts))

y <- calcNormFactors(DGEList(counts = cts))

v <- voom(y, design, plot = FALSE)

fit <- eBayes(lmFit(v, design), robust = TRUE)

de <- map_dfr(taxa, function(tx) {
    topTable(fit, coef = tx, number = Inf, sort.by = "none") %>% rownames_to_column("gene_name") %>% 
        mutate(Pathogen = tx, .before = 1)
})

write.csv(de, repo_output(repo_path(out, "gene_associations_all_taxa.csv"), "Figure4D_analysis.R"), row.names = FALSE)

mapping <- bitr(rownames(cts), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

sets <- read.gmt(repo_input(gmt))

sets <- split(sets$gene, sets$term)

gsea <- map_dfr(taxa, function(tx) {
    message("GSEA: ", tx)
    ranked <- de %>% filter(Pathogen == tx) %>% left_join(mapping, by = c(gene_name = "ENSEMBL")) %>% 
        na.omit() %>% arrange(desc(abs(t))) %>% distinct(SYMBOL, .keep_all = TRUE)
    ranks <- sort(setNames(ranked$t, ranked$SYMBOL), decreasing = TRUE)
    set.seed(1)
    fgsea(pathways = sets, stats = ranks, minSize = 15, maxSize = 500, nproc = 1) %>% as_tibble() %>% 
        mutate(Pathogen = tx, n_ranked_genes = length(ranks), .before = 1)
}) %>% mutate(FDR_across_all_taxa = p.adjust(pval, "BH"), leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

write.csv(gsea, repo_output(repo_path(out, "GSEA_all_taxa.csv"), "Figure4D_analysis.R"), row.names = FALSE)

write.csv(filter(gsea, padj < 0.05), repo_output(repo_path(out, "GSEA_significant_per_taxon_FDR.csv"), 
    "Figure4D_analysis.R"), row.names = FALSE)

summary <- de %>% group_by(Pathogen) %>% summarise(n_genes = n(), n_gene_FDR05 = sum(adj.P.Val < 0.05), 
    .groups = "drop") %>% left_join(gsea %>% group_by(Pathogen) %>% summarise(n_pathways = n(), n_positive_FDR05 = sum(padj < 
    0.05 & NES > 0), n_negative_FDR05 = sum(padj < 0.05 & NES < 0), n_pathways_global_FDR05 = sum(FDR_across_all_taxa < 
    0.05), .groups = "drop"), by = "Pathogen") %>% left_join(prevalence, by = "Pathogen")

write.csv(summary, repo_output(repo_path(out, "microbial_association_summary.csv"), "Figure4D_analysis.R"), 
    row.names = FALSE)

write.csv(prevalence, repo_output(repo_path(out, "microbial_prevalence.csv"), "Figure4D_analysis.R"), row.names = FALSE)

cohort <- md %>% count(Group, Hospital_Death, name = "n")

write.csv(cohort, repo_output(repo_path(out, "cohort_summary.csv"), "Figure4D_analysis.R"), row.names = FALSE)

old <- read.csv(repo_input(repo_path(earli, "Inputs/gsea.csv")))

check <- filter(gsea, Pathogen == "HHV5") %>% dplyr::select(pathway, new_NES = NES, new_FDR = padj) %>% 
    full_join(old %>% dplyr::select(pathway, old_NES = NES, old_FDR = padj), by = "pathway") %>% mutate(NES_difference = new_NES - 
    old_NES, significance_agrees = (new_FDR < 0.05) == (old_FDR < 0.05), direction_agrees = sign(new_NES) == 
    sign(old_NES))

write.csv(check, repo_output(repo_path(out, "HCMV_previous_GSEA_comparison.csv"), "Figure4D_analysis.R"), 
    row.names = FALSE)

plots <- map(taxa, function(tx) {
    d <- gsea %>% filter(Pathogen == tx, padj < 0.05, !grepl("^GOBP_(SENSORY_|DETECTION_)", pathway)) %>% 
        mutate(Direction = if_else(NES > 0, "Positive", "Negative")) %>% group_by(Direction) %>% arrange(padj, 
        desc(abs(NES)), .by_group = TRUE) %>% slice_head(n = 15) %>% ungroup() %>% mutate(signed_FDR = sign(NES) * 
        -log10(pmax(padj, .Machine$double.xmin)), label = str_to_lower(str_replace_all(str_remove(pathway, 
        "^GOBP_"), "_", " ")), label = str_replace_all(label, c(trna = "tRNA", rrna = "rRNA", rna = "RNA", 
        dna = "DNA", `\\batp\\b` = "ATP")), label = case_when(pathway == "GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS" ~ 
        "Somatic recombination-based adaptive immunity", pathway == "GOBP_ADENYLATE_CYCLASE_MODULATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY" ~ 
        "GPCR-cAMP signaling", TRUE ~ label))
    ggplot(d, aes(signed_FDR, reorder(label, signed_FDR), fill = Direction)) + geom_col(width = 0.8) + 
        geom_vline(xintercept = 0) + geom_vline(xintercept = c(-1, 1) * -log10(0.05), linetype = 2, color = "grey50") + 
        scale_fill_manual(values = c(Positive = "#D73027", Negative = "#4575B4")) + labs(title = paste("EARLI", 
        if (tx == "HHV5") 
            "HCMV (HHV5)"
        else tx), subtitle = if (nrow(d)) 
        "Top 15 terms per direction; within-taxon BH-FDR <0.05"
    else "No pathways with within-taxon BH-FDR <0.05", x = "-log10(FDR) \303\227 direction", y = NULL) + 
        theme_bw(base_size = 10) + theme(legend.position = "none", panel.grid.minor = element_blank())
})

pdf(repo_output(repo_path(out, "EARLI_all_taxa_GSEA.pdf"), "Figure4D_analysis.R"), width = 4, height = 4.3)

walk(plots, print)

dev.off()

writeLines(c("Reference: EARLI_mortality-main/1.R", paste("Samples:", length(ids), "Unique barcodes:", 
    n_distinct(md$Barcode)), paste("Hospital deaths:", sum(md$Hospital_Death == 1)), "Primary FDR: BH within each taxon across pathways; additional global BH across all taxon-pathway tests.", 
    "GSEA uses the original minSize=15, maxSize=500, seed=1 and package defaults.", "Figure excludes SENSORY_/DETECTION_ terms only for display; all tests remain in FDR correction.", 
    "Previous gsea.csv provenance is uncertain because original export references undefined gsea_.", 
    capture.output(print(summary, width = Inf))), repo_output(repo_path(out, "analysis_summary.txt"), 
    "Figure4D_analysis.R"))

capture.output(sessionInfo(), file = repo_path(out, "sessionInfo.txt"))

message("Done: ", out)

