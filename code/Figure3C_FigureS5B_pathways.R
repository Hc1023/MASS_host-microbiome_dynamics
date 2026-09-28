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

load(repo_input("Inputs/1211_metadata.rdata"))

load(repo_input("Inputs/1211_transcriptome.rdata"))

metadata_analysis = df_long %>% filter(Timepoint %in% c("D1", "D4", "D7")) %>% droplevels()

vars = c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))

metadata_analysis %<>% left_join(tmp, by = "HumanID")

counts_analysis = counts %>% dplyr::select(all_of(metadata_analysis$SampleID))

stopifnot(identical(metadata_analysis$SampleID, colnames(counts_analysis)))

design <- model.matrix(~Mortality28d * Timepoint + Gender + Age + CenterGroup + PneumoniaTypeGroup + 
    CCI + SOFA_24h + Immunosuppression + MV, data = metadata_analysis)

qr(design)$rank

ncol(design)

colnames(design) <- gsub("Mortality28d", "M", colnames(design))

colnames(design) <- gsub("Timepoint", "", colnames(design))

colnames(design) <- make.names(colnames(design))

dge <- DGEList(counts = as.matrix(counts_analysis))

dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = FALSE)

fit <- lmFit(v, design)

fit <- eBayes(fit)

colnames(design)

contrast_mat <- makeContrasts(D1 = M1, D4 = M1 + M1.D4, D7 = M1 + M1.D7, Dyn_D4vsD1 = M1.D4, Dyn_D7vsD1 = M1.D7, 
    Overall = (M1 + (M1 + M1.D4) + (M1 + M1.D7))/3, levels = design)

fit2 <- contrasts.fit(fit, contrast_mat)

fit2 <- eBayes(fit2)

res_D1 <- topTable(fit2, coef = "D1", number = Inf, adjust.method = "BH")

res_D4 <- topTable(fit2, coef = "D4", number = Inf, adjust.method = "BH")

res_D7 <- topTable(fit2, coef = "D7", number = Inf, adjust.method = "BH")

res_Dyn41 <- topTable(fit2, coef = "Dyn_D4vsD1", number = Inf, adjust.method = "BH")

res_Dyn71 <- topTable(fit2, coef = "Dyn_D7vsD1", number = Inf, adjust.method = "BH")

res_Overall <- topTable(fit2, coef = "Overall", number = Inf, adjust.method = "BH")

head(res_Dyn41)

head(res_Dyn71)

library(msigdbr)

library(fgsea)

library(clusterProfiler)

geneset <- read.gmt(repo_input("Inputs/c5.go.bp.v2025.1.Hs.symbols.gmt"))

bp_sets <- split(geneset$gene, geneset$term)

make_ranks_symbol <- function(res_tt, gene_annot, id_col = "ENSEMBL", sym_col = "SYMBOL", stat = c("t", 
    "logFC")) {
    stat <- match.arg(stat)
    df <- res_tt %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% left_join(gene_annot, by = "ENSEMBL") %>% 
        filter(!is.na(SYMBOL))
    df2 <- df %>% arrange(desc(abs(.data[[stat]]))) %>% distinct(SYMBOL, .keep_all = TRUE)
    ranks <- df2[[stat]]
    names(ranks) <- df2$SYMBOL
    ranks <- sort(ranks, decreasing = TRUE)
    return(ranks)
}

run_fgsea_bp <- function(ranks, bp_sets, minSize = 15, maxSize = 500) {
    fgseaRes <- fgsea(pathways = bp_sets, stats = ranks, minSize = minSize, maxSize = maxSize) %>% as_tibble() %>% 
        arrange(padj, desc(abs(NES)))
    return(fgseaRes)
}

gene_annot = gene_attr

colnames(gene_annot) = c("ENSEMBL", "SYMBOL")

gene_annot = gene_annot[!grepl("ENSG", gene_annot$SYMBOL), ]

ranks_Dyn41 <- make_ranks_symbol(res_Dyn41, gene_annot, stat = "t")

set.seed(1)

gsea_Dyn41 <- run_fgsea_bp(ranks_Dyn41, bp_sets)

ranks_Dyn71 <- make_ranks_symbol(res_Dyn71, gene_annot, stat = "t")

set.seed(1)

gsea_Dyn71 <- run_fgsea_bp(ranks_Dyn71, bp_sets)

plot_fgsea_bidir <- function(gseaRes, top_n = 15, title = "", padj_cut = NULL, wrap_width = 60, thr = NULL) {
    df <- gseaRes %>% filter(!is.na(padj), !is.na(NES))
    df$pathway = gsub("GOBP_", "", df$pathway)
    drop_pat <- "^(SENSORY_|DETECTION_)"
    df %<>% filter(!grepl(drop_pat, pathway))
    if (!is.null(padj_cut)) {
        df <- df %>% filter(padj <= padj_cut)
    }
    up <- df %>% filter(NES > 0) %>% arrange(padj, desc(NES)) %>% slice_head(n = top_n) %>% mutate(direction = "Up")
    down <- df %>% filter(NES < 0) %>% arrange(padj, NES) %>% slice_head(n = top_n) %>% mutate(direction = "Down")
    top_combined <- bind_rows(up, down) %>% mutate(pathway2 = str_wrap(pathway, width = wrap_width)) %>% 
        arrange(NES) %>% mutate(pathway2 = factor(pathway2, levels = pathway2))
    p <- ggplot(top_combined, aes(x = NES, y = pathway2, fill = direction)) + geom_col(width = 0.8) + 
        geom_vline(xintercept = 0) + scale_fill_manual(values = c(Up = "#D73027", Down = "#4575B4")) + 
        labs(title = title, x = "NES", y = NULL) + theme_bw() + theme(panel.grid.minor = element_blank())
    if (!is.null(thr)) {
        p <- p + geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) + geom_vline(xintercept = -thr, 
            linetype = "dashed", alpha = 0.5)
    }
    return(p)
}

p1 <- plot_fgsea_bidir(gsea_Dyn41, top_n = 15, title = "Dyn D4 vs D1 (interaction)", padj_cut = 0.05)

p2 <- plot_fgsea_bidir(gsea_Dyn71, top_n = 15, title = "Dyn D7 vs D1 (interaction)", padj_cut = 0.1)

p1

p2

pdf(repo_output(paste0("Outputs/1222_DynGSEA.pdf"), "Figure3C_FigureS5B_pathways.R"), width = 7.3, height = 4.7)

print(p1)

dev.off()

pdf(repo_output(paste0("Outputs/1222_DynGSEA_2.pdf"), "Figure3C_FigureS5B_pathways.R"), width = 8.1, height = 4.7)

print(p2)

dev.off()

get_pos <- function(gsea_res, tp, padj_cut = 0.1) {
    out <- gsea_res %>% filter(!is.na(NES), !is.na(padj), NES > 0, padj < padj_cut) %>% arrange(padj) %>% 
        mutate(Timepoint = tp)
    out$pathway = gsub("GOBP_", "", out$pathway)
    drop_pat <- "^(SENSORY_|DETECTION_)"
    if (!is.null(drop_pat)) {
        out <- out %>% filter(!grepl(drop_pat, pathway))
    }
    out
}

gsea_pos_D7 <- get_pos(gsea_Dyn71, "D7", padj_cut = 1)

gsea_pos_D4 <- get_pos(gsea_Dyn41, "D4", padj_cut = 1)

top10_union <- union(gsea_pos_D4 %>% slice_head(n = 10) %>% pull(pathway), gsea_pos_D7 %>% slice_head(n = 10) %>% 
    pull(pathway))

df_long <- bind_rows(gsea_pos_D4 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D4"), 
    gsea_pos_D7 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D7")) %>% filter(pathway %in% 
    top10_union) %>% distinct(pathway, Timepoint, .keep_all = TRUE)

mat_nes <- df_long %>% dplyr::select(pathway, Timepoint, NES) %>% pivot_wider(names_from = Timepoint, 
    values_from = NES) %>% tibble::column_to_rownames("pathway") %>% as.matrix()

mat_padj <- df_long %>% dplyr::select(pathway, Timepoint, padj) %>% pivot_wider(names_from = Timepoint, 
    values_from = padj) %>% tibble::column_to_rownames("pathway") %>% as.matrix()

ord <- order(rowSums(mat_nes, na.rm = TRUE), decreasing = TRUE)

mat_nes <- mat_nes[ord, , drop = FALSE]

mat_padj <- mat_padj[ord, , drop = FALSE]

rng <- max(abs(mat_nes), na.rm = TRUE)

col_fun <- colorRamp2(c(0, 3), c("white", "#D73027"))

sig_symbol <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 0.05, "*", ifelse(p < 0.1, "\302\267", "")))
}

go_label_map_pos <- c(MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY = "MyD88-dependent TLR signaling", 
    REGULATION_OF_LYSOSOMAL_LUMEN_PH = "regulation of lysosomal lumen PH", RESPONSE_TO_TYPE_I_INTERFERON = "response to IFN-1", 
    PHAGOLYSOSOME_ASSEMBLY = "phagolysosome assembly", VACUOLAR_ACIDIFICATION = "vacuolar acidification", 
    PHAGOSOME_MATURATION = "phagosome maturation", TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY = "TNF-mediated signaling", 
    REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS = "regulation of response to cytokine stimulus", REGULATORY_NCRNA_MEDIATED_GENE_SILENCING = "regulatory ncRNA-mediated gene silencing", 
    ENDOLYSOSOMAL_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY = "endolysosomal TLR signaling", REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS = "regulation of epithelial cell apoptotic process", 
    CYTOKINE_MEDIATED_SIGNALING_PATHWAY = "cytokine-mediated signaling pathway", CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY = "cytoplasmic PRR signaling", 
    CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION = "canonical NF-kB signal transduction", NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS = "negative regulation of epithelial apoptosis", 
    REGULATION_OF_HUMORAL_IMMUNE_RESPONSE = "regulation of humoral immune response", POSITIVE_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS = "positive regulation of response to biotic stimulus", 
    REGULATION_OF_INNATE_IMMUNE_RESPONSE = "regulation of innate immune response", POSITIVE_REGULATION_OF_DEFENSE_RESPONSE = "positive regulation of defense response", 
    VACUOLE_ORGANIZATION = "vacuole organization", REGULATION_OF_INFLAMMATORY_RESPONSE = "regulation of inflammatory response", 
    MACROAUTOPHAGY = "macroautophagy", VESICLE_ORGANIZATION = "vesicle organization")

rownames(mat_nes) <- go_label_map_pos[rownames(mat_nes)]

rownames(mat_padj) <- rownames(mat_nes)

rn_w <- max_text_width(rownames(mat_nes), gp = gpar(fontsize = 8))

ht1 = Heatmap(mat_nes, name = "NES", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "grey90", 
    column_names_rot = 0, row_names_side = "left", row_names_gp = gpar(fontsize = 8), row_names_max_width = rn_w + 
        unit(6, "mm"), column_title = "Trajectory up", cell_fun = function(j, i, x, y, width, height, 
        fill) {
        grid.text(sig_symbol(mat_padj[i, j]), x, y, gp = gpar(fontsize = 14))
    })

ht1

pdf(repo_output(paste0("Outputs/1222_DynTop_up.pdf"), "Figure3C_FigureS5B_pathways.R"), width = 3.6, height = 4)

print(ht1)

dev.off()

get_neg <- function(gsea_res, tp, padj_cut = 0.1) {
    out <- gsea_res %>% filter(!is.na(NES), !is.na(padj), NES < 0, padj < padj_cut) %>% arrange(padj) %>% 
        mutate(Timepoint = tp)
    out$pathway = gsub("GOBP_", "", out$pathway)
    drop_pat <- "^(SENSORY_|DETECTION_)"
    if (!is.null(drop_pat)) {
        out <- out %>% filter(!grepl(drop_pat, pathway))
    }
    out
}

gsea_neg_D7 <- get_neg(gsea_Dyn71, "D7", padj_cut = 1)

gsea_neg_D4 <- get_neg(gsea_Dyn41, "D4", padj_cut = 1)

top10_union <- union(gsea_neg_D4 %>% slice_head(n = 10) %>% pull(pathway), gsea_neg_D7 %>% slice_head(n = 10) %>% 
    pull(pathway))

df_long <- bind_rows(gsea_neg_D4 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D4"), 
    gsea_neg_D7 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D7")) %>% filter(pathway %in% 
    top10_union) %>% distinct(pathway, Timepoint, .keep_all = TRUE)

mat_nes <- df_long %>% dplyr::select(pathway, Timepoint, NES) %>% pivot_wider(names_from = Timepoint, 
    values_from = NES) %>% tibble::column_to_rownames("pathway") %>% as.matrix() %>% na.omit()

mat_padj <- df_long %>% dplyr::select(pathway, Timepoint, padj) %>% pivot_wider(names_from = Timepoint, 
    values_from = padj) %>% tibble::column_to_rownames("pathway") %>% as.matrix() %>% na.omit()

ord <- order(rowSums(mat_nes, na.rm = TRUE))

mat_nes <- mat_nes[ord, , drop = FALSE]

mat_padj <- mat_padj[ord, , drop = FALSE]

rng <- max(abs(mat_nes), na.rm = TRUE)

col_fun <- colorRamp2(c(-rng, 0), c("#4575B4", "white"))

sig_symbol <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 0.05, "*", ifelse(p < 0.1, "\302\267", "")))
}

go_label_map <- c(RIBOSOME_BIOGENESIS = "ribosome biogenesis", RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS = "ribosomal small subunit biogenesis", 
    RRNA_PROCESSING = "rRNA processing", RRNA_METABOLIC_PROCESS = "rRNA metabolic process", INNATE_IMMUNE_RESPONSE_IN_MUCOSA = "innate immune response in mucosa", 
    NUCLEOSOME_ORGANIZATION = "nucleosome organization", ANTIBACTERIAL_HUMORAL_RESPONSE = "antibacterial humoral response", 
    ANTIMICROBIAL_HUMORAL_RESPONSE = "antimicrobial humoral response", ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE = "AMP-mediated humoral immune response", 
    RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS = "ribonucleoprotein complex biogenesis", TRNA_METABOLIC_PROCESS = "tRNA metabolic process", 
    MITOCHONDRIAL_GENE_EXPRESSION = "mitochondrial gene expression", MITOCHONDRIAL_TRANSLATION = "mitochondrial translation", 
    RNA_MODIFICATION = "RNA modification", TRNA_PROCESSING = "tRNA processing", DISRUPTION_OF_ANATOMICAL_STRUCTURE_IN_ANOTHER_ORGANISM = "disruption of host tissue structure", 
    PROTEIN_DNA_COMPLEX_ORGANIZATION = "protein-DNA complex organization", HUMORAL_IMMUNE_RESPONSE = "humoral immune response", 
    MITOTIC_NUCLEAR_DIVISION = "mitotic nuclear division", ORGANELLE_FISSION = "organelle fission")

rownames(mat_nes) <- go_label_map[rownames(mat_nes)]

rownames(mat_padj) <- rownames(mat_nes)

rn_w <- max_text_width(rownames(mat_nes), gp = gpar(fontsize = 8))

ht2 = Heatmap(mat_nes, name = "NES", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "grey90", 
    column_names_rot = 0, row_names_side = "left", row_names_gp = gpar(fontsize = 8), row_names_max_width = rn_w + 
        unit(6, "mm"), column_title = "Trajectory down", cell_fun = function(j, i, x, y, width, height, 
        fill) {
        grid.text(sig_symbol(mat_padj[i, j]), x, y, gp = gpar(fontsize = 14))
    })

pdf(repo_output(paste0("Outputs/1222_DynTop_down.pdf"), "Figure3C_FigureS5B_pathways.R"), width = 3.5, height = 4)

print(ht2)

dev.off()

library(GSVA)

library(GSEABase)

meta_gsva = metadata_analysis %>% dplyr::select(HumanID, Timepoint, SampleID, Mortality28d)

expr <- v$E

sym <- gene_attr[rownames(expr), "SYMBOL"]

is_valid_symbol <- !is.na(sym) & !grepl("^ENSG", sym)

expr_sym <- expr[is_valid_symbol, , drop = FALSE]

rownames(expr_sym) <- sym[is_valid_symbol]

sym_tab <- table(rownames(expr_sym))

unique_sym <- names(sym_tab[sym_tab == 1])

expr_sym <- expr_sym[rownames(expr_sym) %in% unique_sym, , drop = FALSE]

up1_pathways <- c("GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION", 
    "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY", 
    "GOBP_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY")

up2_pathways <- c("GOBP_PHAGOLYSOSOME_ASSEMBLY", "GOBP_PHAGOSOME_MATURATION", "GOBP_VACUOLE_ORGANIZATION", 
    "GOBP_MACROAUTOPHAGY")

up3_pathways <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON")

dw1_pathways <- c("GOBP_RIBOSOME_BIOGENESIS", "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS", "GOBP_RRNA_PROCESSING", 
    "GOBP_RRNA_METABOLIC_PROCESS")

dw2_pathways <- c("GOBP_NUCLEOSOME_ORGANIZATION", "GOBP_PROTEIN_DNA_COMPLEX_ORGANIZATION", "GOBP_MITOTIC_NUCLEAR_DIVISION", 
    "GOBP_ORGANELLE_FISSION")

extract_leading_edge <- function(gsea_res, pathways) {
    gsea_res %>% filter(pathway %in% pathways) %>% dplyr::select(pathway, leadingEdge) %>% tidyr::unnest(leadingEdge) %>% 
        distinct(pathway, leadingEdge)
}

get_leading_edge_intersect <- function(gsea_D4, gsea_D7, pathways) {
    le_D4 <- extract_leading_edge(gsea_D4, pathways)$leadingEdge
    le_D7 <- extract_leading_edge(gsea_D7, pathways)$leadingEdge
    intersect(le_D4, le_D7)
}

modules <- list(up1 = up1_pathways, up2 = up2_pathways, up3 = up3_pathways, dw1 = dw1_pathways)

gs_list <- lapply(modules, function(pw) get_leading_edge_intersect(gsea_Dyn41, gsea_Dyn71, pw))

gsva_param <- gsvaParam(expr = as.matrix(expr_sym), geneSets = gs_list, kcdf = "Gaussian", minSize = 10, 
    maxSize = 500)

gsva_mat <- gsva(gsva_param)

identical(colnames(gsva_mat), meta_gsva$SampleID)

plot_df = bind_cols(meta_gsva[, 1:4], t(gsva_mat))

plot_long <- plot_df %>% pivot_longer(cols = -c(HumanID, Timepoint, SampleID, Mortality28d), names_to = "Module", 
    values_to = "Score")

table(plot_long$Module)

plot_long <- plot_long %>% mutate(Module = factor(Module, levels = c("up1", "up2", "up3", "dw1"), labels = c("Inflammatory signaling", 
    "Phagolysosome function", "IFN signaling", "Ribosome biogenesis")))

library(rstatix)

wilcox_df <- plot_long %>% group_by(Module, Timepoint) %>% wilcox_test(Score ~ Mortality28d) %>% ungroup() %>% 
    mutate(p.label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "\302\267", 
        TRUE ~ "ns"))

y_pos_df <- plot_long %>% group_by(Module, Timepoint, Mortality28d) %>% summarise(m = mean(Score, na.rm = TRUE), 
    se = sd(Score, na.rm = TRUE)/sqrt(sum(!is.na(Score))), .groups = "drop") %>% group_by(Module, Timepoint) %>% 
    summarise(y.position = max(m + se, na.rm = TRUE), .groups = "drop") %>% left_join(plot_long %>% group_by(Module, 
    Timepoint) %>% summarise(rng = diff(range(Score, na.rm = TRUE)), .groups = "drop"), by = c("Module", 
    "Timepoint")) %>% mutate(y.position = y.position + 0.01 * rng) %>% select(Module, Timepoint, y.position)

wilcox_df <- wilcox_df %>% left_join(y_pos_df, by = c("Module", "Timepoint"))

p = ggplot(plot_long, aes(x = Timepoint, y = Score, color = factor(Mortality28d), group = interaction(HumanID, 
    Mortality28d))) + stat_summary(aes(group = Mortality28d), fun = mean, geom = "line", linewidth = 1.3) + 
    stat_summary(aes(group = Mortality28d), fun = mean, geom = "point", size = 3) + stat_summary(aes(group = Mortality28d), 
    fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.6) + stat_pvalue_manual(wilcox_df, 
    label = "p.label", x = "Timepoint", y.position = "y.position", tip.length = 0, size = 3) + theme_bw() + 
    scale_color_manual(values = alpha(c(`0` = "#4575B4", `1` = "#D73027"), 0.9), name = "", labels = c("Survival", 
        "Mortality")) + facet_wrap(~Module, scales = "free_y", ncol = 4) + theme(strip.text = element_text(face = "bold"), 
    panel.grid.minor = element_blank()) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

pdf(repo_output(paste0("Outputs/1222_traj.pdf"), "Figure3C_FigureS5B_pathways.R"), width = 9, height = 2.2)

print(p)

dev.off()

