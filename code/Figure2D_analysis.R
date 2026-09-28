# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(magrittr)

library(data.table)

library(ggpubr)

library(scales)

library(edgeR)
library(limma)

load(repo_input("Inputs/1211_metadata.rdata"))

load(repo_input("Inputs/1211_microbe.rdata"))

load(repo_input("Inputs/1211_transcriptome.rdata"))

geneID = setNames(gene_attr$SYMBOL, gene_attr$gene_id)

data_filtered = data[rowSums(data > 0) > 0.05 * ncol(data), ]

data_filtered = data_filtered[order(rowSums(data_filtered > 0), decreasing = T), ]

data_filtered_log2 = log2(data_filtered + 1)

rownames(data_filtered_log2) %<>% make.names()

pathogen_vars = rownames(data_filtered_log2)

pathogen_str = paste(pathogen_vars, collapse = " + ")

df_tp <- df_long

counts_tp <- counts[, df_tp$SampleID]

data_filtered_log2_tp = data_filtered_log2[, df_tp$SampleID]

df_tp = bind_cols(df_tp, t(data_filtered_log2_tp))

{
    df_tp %<>% left_join(meta[, -2], by = "HumanID")
    vars = c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
        "MV")
    vars_str = paste(vars, collapse = " + ")
}

design <- model.matrix(as.formula(paste("~ ", pathogen_str, " + ", vars_str, "+ Timepoint")), data = df_tp)

stopifnot(identical(colnames(counts_tp), df_tp$SampleID))

y <- DGEList(counts = counts_tp)

y <- calcNormFactors(y)

v <- voom(y, design, plot = F)

fit <- lmFit(v, design)

fit <- eBayes(fit, robust = TRUE)

res_list <- lapply(pathogen_vars, function(p) {
    res = topTable(fit, coef = p, number = Inf, sort.by = "none")
    res <- res[order(res$adj.P.Val, decreasing = F), ]
    res$gene_name <- geneID[rownames(res)]
    res$Pathogen = p
    return(res)
})

res_pathogen <- do.call(rbind, res_list)

sig_list_pos <- res_pathogen %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% group_by(Pathogen) %>% summarise(Genes = list(gene_name))

sig_list_neg <- res_pathogen %>% filter(adj.P.Val < 0.05 & logFC < 0) %>% group_by(Pathogen) %>% summarise(Genes = list(gene_name))

library(clusterProfiler)

library(AnnotationDbi)

library(org.Hs.eg.db)

ego_fun = function(gene_vec) {
    entrez_ids <- bitr(gene_vec, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", 
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    return(ego)
}

ego_pos_list = list()

ego_neg_list = list()

for (i in 1:2) {
    print(i)
    ego_pos_list[[i]] = ego_fun(sig_list_pos$Genes[[i]])
    ego_neg_list[[i]] = ego_fun(sig_list_neg$Genes[[i]])
}

p3_list = list()

for (i in 1:2) {
    ego_up = data.frame(ego_pos_list[[i]])
    ego_down = data.frame(ego_neg_list[[i]])
    top_combined <- bind_rows(as.data.frame(ego_up) %>% mutate(direction = "Up"), as.data.frame(ego_down) %>% 
        mutate(direction = "Down")) %>% mutate(p_signed = ifelse(direction == "Down", log10(p.adjust), 
        -log10(p.adjust))) %>% group_by(direction) %>% slice_min(abs(p.adjust), n = 15) %>% ungroup()
    if (i == 2) {
        idx = top_combined$Description == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"
        top_combined$Description[idx] = "Somatic-recombination\342\200\223based adaptive immune response (Ig-SF)"
    }
    thr <- -log10(0.05)
    p3 = ggplot(top_combined, aes(x = p_signed, y = reorder(Description, p_signed), fill = direction)) + 
        geom_col() + scale_fill_manual(values = c(Up = "#D73027", Down = "#4575B4")) + labs(x = "-log10(FDR)*Direction", 
        y = "") + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
        geom_vline(xintercept = 0) + geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) + 
        geom_vline(xintercept = -thr, linetype = "dashed", alpha = 0.5)
    p3_list[[i]] = p3
}

p3_list[[1]]

p3_list[[2]]

invisible(NULL)

invisible(NULL)

pdf(repo_output(paste0("Outputs/1423_microbe_DEGO.pdf"), "Figure2D_analysis.R"), width = 5.5, height = 4.3)

print(p3_list[[1]])

print(p3_list[[2]])

dev.off()


# Export only HCMV/EBV results; all tested pathways are retained for their FDR families.
write.csv(res_pathogen %>% filter(Pathogen %in% c("HCMV", "HHV.4")), "Outputs/Figure2D_gene_associations.csv", row.names=FALSE)
enrichment_results <- bind_rows(lapply(1:2, function(i) bind_rows(
  transform(ego_pos_list[[i]]@result, Virus=c("HCMV","EBV")[i], Direction="Up"),
  transform(ego_neg_list[[i]]@result, Virus=c("HCMV","EBV")[i], Direction="Down"))))
write.csv(enrichment_results, "Outputs/Figure2D_enrichment.csv", row.names=FALSE)
write.csv(bind_rows(lapply(1:2, function(i) transform(p3_list[[i]]$data, Virus=c("HCMV","EBV")[i]))), "Outputs/Figure2D_plot.csv", row.names=FALSE)
