rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(scales)
library(edgeR) 

# Read data
data <-read.csv("Inputs/EARLI_microbe.csv")
bc1 = data$Barcode

###Gene count data
genecountspc<-read.csv("Inputs/earli_counts_kallisto_mortality.csv", header=TRUE)
genecountspc<-as.matrix(genecountspc)
rownames(genecountspc)<-genecountspc[,1]
genecountspc<-genecountspc[,-1]
bc2 = intersect(bc1, colnames(genecountspc))

###Metadata
metadata<-read.csv("Inputs/Cleaned_Metadata_070324.csv")
metadata$Barcode = paste0("EARLI_", metadata$Barcode)
bc3 = intersect(bc2, metadata$Barcode)


#### prepare microbe data
rownames(data) = data$Barcode
data = data[,-c(1:3)]
data = t(data)
rownames(data)

data = data.frame(data) %>% dplyr::select(all_of(bc3))

# > rowSums(data>0)[order(rowSums(data>0), decreasing = T)][1:5]
# HHV4   Escherichia    Klebsiella Streptococcus 
# 20            17            12            12 
# HHV5 
# 12 

data_filtered = data[rowSums(data>0) > 0.05*ncol(data),]
data_filtered = data_filtered[order(rowSums(data_filtered>0), decreasing = T),]
data_filtered_log2 = log2(data_filtered+1)

pathogen_vars = rownames(data_filtered_log2)
pathogen_str = paste(pathogen_vars, collapse = " + ")

counts <- data.frame(genecountspc) %>% dplyr::select(all_of(bc3))
counts <- apply(counts, c(1,2), function(x) { (as.integer(x))})

rownames(metadata) = metadata$Barcode
df = metadata[bc3,]

table(df$Hospital_Death)
# 0  1 
# 81 42

vars = c('Gender', 'Age', 'Group', 'APACHEIII',
         'Immunocomp_manuscript', 'Intubated')
vars_str = paste(vars, collapse = " + ")

df = df %>% dplyr::select(Barcode, Hospital_Death, 
                          daysofsurvival, all_of(vars))
df$Gender %<>% factor()
df$Group %<>% factor()
df$Immunocomp_manuscript %<>% factor()
df$Intubated %<>% factor()
str(df)

df = bind_cols(df, t(data_filtered_log2))


design <- model.matrix(as.formula(paste("~ ", pathogen_str,
                                        " + ", vars_str)),
                       data = df)

stopifnot(identical(colnames(counts), rownames(df)))

# 构建 DGE 对象
y <- DGEList(counts = counts)
y <- calcNormFactors(y)

# voom + limma
v <- voom(y, design, plot = F)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust=TRUE)

res_list <- lapply(pathogen_vars, function(p) {
  res = topTable(fit, coef = p, number = Inf, sort.by = "none")
  res <- res[order(res$adj.P.Val, decreasing = F), ]
  res$gene_name <- rownames(res)
  res$Pathogen = p
  return(res)
})


res_pathogen <- do.call(rbind, res_list)

res = res_pathogen %>% filter(Pathogen == 'HHV5')

library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

# 准备 GO BP 基因集（MSigDB C5:BP）去GSEA官网下载
geneset <- read.gmt("../../AI_MASS/Utility/c5.go.bp.v2025.1.Hs.symbols.gmt")
# gene sets: list(term -> gene symbols)
bp_sets <- split(geneset$gene, geneset$term)

gene_df <- bitr(
  res$gene_name,
  fromType = "ENSEMBL",
  toType = c("SYMBOL"),
  OrgDb = org.Hs.eg.db
)


make_ranks_symbol <- function(res, gene_df,
                              id_col = "ENSEMBL", sym_col = "SYMBOL") {
  stat <- "t"
  
  df <- res %>%
    base::as.data.frame() %>%
    dplyr::left_join(gene_df, by = c("gene_name"= "ENSEMBL")) %>%
    na.omit()
  
  # 一个 SYMBOL 对应多个 ENSEMBL 时，保留 |stat| 最大的那条
  df2 <- df %>%
    arrange(desc(abs(.data[[stat]]))) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  
  ranks <- df2[[stat]]
  names(ranks) <- df2$SYMBOL
  ranks <- sort(ranks, decreasing = TRUE)
  return(ranks)
}

## 跑 fgsea


run_fgsea_bp <- function(ranks, bp_sets,
                         minSize = 15, maxSize = 500) {
  fgseaRes <- fgsea(
    pathways = bp_sets,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize
  ) %>% as_tibble() %>%
    arrange(padj, desc(abs(NES)))
  return(fgseaRes)
}


ranks <- make_ranks_symbol(res, gene_df)
set.seed(1)
gsea <- run_fgsea_bp(ranks, bp_sets)

if(F){
  write.csv(gsea[,-ncol(gsea)], 'Outputs/Supplementary_data_EARLI_gsea.csv', row.names = F)
}
plot_fgsea_bidir <- function(gseaRes,
                             top_n = 15,
                             title = "",
                             padj_cut = NULL,      # 可选：比如 0.05
                             wrap_width = 60,      # 过长通路名换行
                             thr = NULL            # 可选：画 NES 阈值虚线，比如 1.5
) {
  
  df <- gseaRes %>%
    filter(!is.na(padj), !is.na(NES))
  df$pathway = gsub('GOBP_', '', df$pathway)
  ## 可选：过滤掉“感知刺激/味觉嗅觉”等常见伪信号 GO
  drop_pat <- "^(SENSORY_|DETECTION_)"
  df %<>% filter(!grepl(drop_pat, pathway))
  
  if (!is.null(padj_cut)) {
    df <- df %>% filter(padj <= padj_cut)
  }
  
  up <- df %>%
    filter(NES > 0) %>%
    arrange(padj, desc(NES)) %>%
    slice_head(n = top_n) %>%
    mutate(direction = "Up")
  
  down <- df %>%
    filter(NES < 0) %>%
    arrange(padj, NES) %>%
    slice_head(n = top_n) %>%
    mutate(direction = "Down")
  
  top_combined <- bind_rows(up, down)
  
  pathway2 <- top_combined$pathway
  
  pathway2 <- pathway2 |>
    tolower() |>
    gsub("_", " ", x = _)
  
  # 恢复常见缩写
  pathway2 <- gsub("trna", "tRNA", pathway2)
  pathway2 <- gsub("rrna", "rRNA", pathway2)
  pathway2 <- gsub("rna", "RNA", pathway2)
  pathway2 <- gsub("dna", "DNA", pathway2)
  pathway2 <- gsub("gpcr", "GPCR", pathway2)
  
  
  pathway2[3] = "GPCR-cAMP signaling"
  top_combined$pathway2 <- pathway2
  
  top_combined %<>% 
    mutate(p_signed = ifelse(direction == "Down", 
                             log10(padj), -log10(padj))) %>%
    arrange(p_signed)
    
  
  p <- ggplot(top_combined,
              aes(x = p_signed, y = reorder(pathway2, p_signed), 
                  fill = direction)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
    labs(title = title, x = "", y = NULL) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = c(-4:1)*5, limits = c(-20,5))
  # p
  if (!is.null(thr)) {
    p <- p +
      geom_vline(xintercept =  thr, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = -thr, linetype = "dashed", alpha = 0.5)
  }

  return(p)
}

p = plot_fgsea_bidir(gsea, top_n = 15,
                 title = "",
                 padj_cut = 0.05,
                 thr = -log10(0.05))

pdf("Outputs/03_EARLI_gsea.pdf", height = 4.3, width = 5.3)
print(p)
dev.off()


