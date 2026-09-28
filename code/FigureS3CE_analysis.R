# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(magrittr)

library(clusterProfiler)

library(AnnotationDbi)

library(org.Hs.eg.db)

library(ggplot2)

library(scales)

library(limma)

library(edgeR)

library(tidyverse)

library(ggpubr)

project_dir <- "."

invisible(NULL)

out_dir <- repo_path(".")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

load(repo_input("Inputs/1211_metadata.rdata"))

load(repo_input("Inputs/1211_transcriptome.rdata"))

get_res = function(tp) {
    metadata_analysis = df_long %>% filter(Timepoint == tp) %>% droplevels()
    vars = c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
        "MV")
    tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))
    metadata_analysis %<>% left_join(tmp, by = "HumanID")
    counts_analysis = counts %>% dplyr::select(all_of(metadata_analysis$SampleID))
    stopifnot(identical(metadata_analysis$SampleID, colnames(counts_analysis)))
    design <- model.matrix(~Mortality28d + Gender + Age + CenterGroup + PneumoniaTypeGroup + CCI + SOFA_24h + 
        Immunosuppression + MV, data = metadata_analysis)
    qr(design)$rank
    ncol(design)
    colnames(design) <- gsub("Mortality28d", "M", colnames(design))
    colnames(design) <- make.names(colnames(design))
    dge <- DGEList(counts = as.matrix(counts_analysis))
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    colnames(design)
    contrast_mat <- makeContrasts(M1 = M1, levels = design)
    fit2 <- contrasts.fit(fit, contrast_mat)
    fit2 <- eBayes(fit2)
    res <- topTable(fit2, coef = "M1", number = Inf, adjust.method = "BH")
    res$gene_symbol <- as.character(gene_attr[rownames(res), "SYMBOL"])
    return(res)
}

res_D1 = get_res("D1")

res_D4 = get_res("D4")

res_D7 = get_res("D7")

write.csv(res_D1, file = repo_output(repo_path(out_dir, "deg_D1.csv"), "FigureS3CE_analysis.R"))

write.csv(res_D4, file = repo_output(repo_path(out_dir, "deg_D4.csv"), "FigureS3CE_analysis.R"))

write.csv(res_D7, file = repo_output(repo_path(out_dir, "deg_D7.csv"), "FigureS3CE_analysis.R"))

plot_fun = function(res) {
    res <- res[order(res$adj.P.Val, decreasing = F), ]
    res$adj.P.Val[is.na(res$adj.P.Val)] <- 1
    res_table <- data.frame(res)
    res_table$sig <- res_table$adj.P.Val < 0.1
    res_table$sig[res_table$adj.P.Val < 0.1 & res_table$logFC > 0] <- "1"
    res_table$sig[res_table$adj.P.Val < 0.1 & res_table$logFC < 0] <- "0"
    {
        table(res_table$sig)
        deg_up = res_table %>% filter(sig == "1")
        gene_up <- rownames(deg_up)
        gene_up_df <- bitr(gene_up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        ego_up <- enrichGO(gene = gene_up_df$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", 
            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
        ego_up_simpl_0.7 <- simplify(ego_up, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
        deg_down = res_table %>% filter(sig == "0")
        gene_down <- rownames(deg_down)
        gene_down_df <- bitr(gene_down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        ego_down <- enrichGO(gene = gene_down_df$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
            ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
        ego_down_simpl_0.7 <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
        top_combined <- bind_rows(as.data.frame(ego_up_simpl_0.7) %>% mutate(direction = "Up"), as.data.frame(ego_down_simpl_0.7) %>% 
            mutate(direction = "Down")) %>% mutate(p_signed = ifelse(direction == "Down", log10(p.adjust), 
            -log10(p.adjust))) %>% group_by(direction) %>% slice_min(abs(p.adjust), n = 15) %>% ungroup()
        if (F) {
            top_combined$Description[2] = "TCR-mediated T cell activation"
            top_combined$Description[7] = "Antigen receptor\342\200\223mediated adaptive immunity"
        }
        thr <- -log10(0.05)
        p2 = ggplot(top_combined, aes(x = p_signed, y = reorder(Description, p_signed), fill = direction)) + 
            geom_col() + scale_fill_manual(values = c(Up = "#D73027", Down = "#4575B4")) + labs(x = "-log10(adj.P)*Direction", 
            y = "") + theme_bw() + theme(panel.grid.minor = element_blank()) + geom_vline(xintercept = 0) + 
            geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) + geom_vline(xintercept = -thr, 
            linetype = "dashed", alpha = 0.5)
    }
    p2
    return(p2)
}

go_D1 = plot_fun(res_D1)

go_D4 = plot_fun(res_D4)

go_D7 = plot_fun(res_D7)

pdf(repo_output(repo_path(out_dir, "go_D1.pdf"), "FigureS3CE_analysis.R"), width = 6.5, height = 5)

print(go_D1)

dev.off()

pdf(repo_output(repo_path(out_dir, "go_D4.pdf"), "FigureS3CE_analysis.R"), width = 6.5, height = 5)

print(go_D4)

dev.off()

pdf(repo_output(repo_path(out_dir, "go_D7.pdf"), "FigureS3CE_analysis.R"), width = 6.5, height = 3)

print(go_D7)

dev.off()

