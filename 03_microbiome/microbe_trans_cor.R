rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(scales)
library(edgeR) 

load('Inputs/1211_metadata.rdata')
load('Inputs/1211_microbe.rdata')
load('Inputs/1211_transcriptome.rdata')
geneID = setNames(gene_attr$SYMBOL, gene_attr$gene_id)

data_filtered = data[rowSums(data>0) > 0.05*ncol(data),]
data_filtered = data_filtered[order(rowSums(data_filtered>0), decreasing = T),]
data_filtered_log2 = log2(data_filtered+1)
rownames(data_filtered_log2) %<>% make.names()
pathogen_vars = rownames(data_filtered_log2)
pathogen_str = paste(pathogen_vars, collapse = " + ")

df_tp <- df_long
counts_tp <- counts[, df_tp$SampleID]
data_filtered_log2_tp = data_filtered_log2[,df_tp$SampleID]
df_tp = bind_cols(df_tp, t(data_filtered_log2_tp))

{
  df_tp %<>% left_join(meta[,-2], by = 'HumanID')
  vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
           'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
  vars_str = paste(vars, collapse = " + ")
}  

###Create design model for differential expression
design <- model.matrix(as.formula(paste("~ ", pathogen_str,
                                        " + ", vars_str, 
                                        "+ Timepoint")),
                       data = df_tp)


# 确保顺序一致
stopifnot(identical(colnames(counts_tp), df_tp$SampleID))

# 构建 DGE 对象
y <- DGEList(counts = counts_tp)
y <- calcNormFactors(y)

# voom + limma
v <- voom(y, design, plot = F)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust=TRUE)

res_list <- lapply(pathogen_vars, function(p) {
  res = topTable(fit, coef = p, number = Inf, sort.by = "none")
  res <- res[order(res$adj.P.Val, decreasing = F), ]
  res$gene_name <- geneID[rownames(res)]
  res$Pathogen = p
  return(res)
})

res_pathogen <- do.call(rbind, res_list)

if(F){
  res_out = res_pathogen %>%
    filter(adj.P.Val < 0.05 & logFC>0)
  res_out %<>% filter(Pathogen %in% c('HHV.4', 'HCMV', 'HPgV', 'TTV'))
  res_out$Pathogen[res_out$Pathogen == 'HHV.4'] = 'EBV'
  write.csv(res_out, 'Outputs/Supplementary_data_microbe_trans_gene.csv')
  
}
sig_list_pos <- res_pathogen %>%
  filter(adj.P.Val < 0.05 & logFC>0) %>%
  group_by(Pathogen) %>%
  summarise(Genes = list(gene_name))
# sig_list_pos
sig_list_neg <- res_pathogen %>%
  filter(adj.P.Val < 0.05 & logFC<0) %>%
  group_by(Pathogen) %>%
  summarise(Genes = list(gene_name))
# sig_list_neg


library(ggVennDiagram)
sig_list = sig_list_pos
gene_sets <- setNames(sig_list$Genes, sig_list$Pathogen)
gene_sets <- gene_sets[lengths(gene_sets) >= 50]
names(gene_sets)[2] = 'EBV'

p1 = ggVennDiagram(
  gene_sets,
  label_alpha = 0
) +
  scale_fill_gradient(low = "white", high = "steelblue", name = 'Genes') +
  theme_void()
p1
sig_list = sig_list_neg
gene_sets <- setNames(sig_list$Genes, sig_list$Pathogen)
gene_sets <- gene_sets[lengths(gene_sets) >= 50]
names(gene_sets)[2] = 'EBV'

# names(gene_sets) = NA
p2 <- ggVennDiagram(gene_sets, label_alpha = 0,
                    set_label_geom = "none" ) +
  scale_fill_gradient(low = "white", high = "steelblue", name = 'Genes') +
  theme_void()
if(F){
  pdf('Outputs/03_venn_all.pdf', height = 3, width = 5)
  print(p1)
  print(p2)
  dev.off()
}


library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

ego_fun = function(gene_vec){
  # 基因ID映射（如果是基因符号 -> ENTREZID）
  entrez_ids <- bitr(gene_vec, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  # 富集分析
  ego <- enrichGO(
    gene          = entrez_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",         # 生物过程
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE          # 转回基因符号
  )
  return(ego)
}

ego_pos_list = list()
ego_neg_list = list()
for (i in 1:2) {
  print(i)
  ego_pos_list[[i]] = ego_fun(sig_list_pos$Genes[[i]])
  ego_neg_list[[i]] = ego_fun(sig_list_neg$Genes[[i]])
  
}

if(F){
  microbe_trans_ego = list()
  tmp = c('HCMV','EBV')
  for (i in 1:2) {
    ego_up = data.frame(ego_pos_list[[i]])
    ego_down = data.frame(ego_neg_list[[i]])
    ego_up$Direction = 'Mortality'
    ego_down$Direction = 'Survival'
    microbe_trans_ego[[tmp[i]]] = bind_rows(ego_up, ego_down)
  }
  library(openxlsx)
  # write to excel
  write.xlsx(
    microbe_trans_ego,
    file = "Outputs/Supplementary_data_microbe_trans_GO.xlsx",
    asTable = TRUE
  )
}

p3_list = list()

for (i in 1:2) {
  ego_up = data.frame(ego_pos_list[[i]])
  ego_down = data.frame(ego_neg_list[[i]])
  
  top_combined <- bind_rows(
    as.data.frame(ego_up) %>% mutate(direction = "Up"),
    as.data.frame(ego_down) %>% mutate(direction = "Down")
  ) %>%
    mutate(p_signed = ifelse(direction == "Down", log10(p.adjust), -log10(p.adjust))) %>%
    group_by(direction) %>%
    slice_min(abs(p.adjust), n = 15) %>%
    ungroup()
  if(i == 2){
    idx = top_combined$Description == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"
    top_combined$Description[idx] = 'Somatic-recombination–based adaptive immune response (Ig-SF)'
  }
  # 绘图
  thr <- -log10(0.05)
  p3 = ggplot(top_combined,
              aes(x = p_signed,
                  y = reorder(Description, p_signed),
                  fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("Up"="#D73027", "Down"="#4575B4")) +
    labs(x="-log10(FDR)*Direction", y="") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept =  thr, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = -thr, linetype = "dashed", alpha = 0.5)
  
  p3_list[[i]] = p3
}


pdf(paste0("Outputs/03_microbe_trans_GO.pdf"), width = 5.5, height = 4.3)
print(p3_list[[1]])
print(p3_list[[2]])
dev.off()




