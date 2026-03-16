rm(list = ls())
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

load('Inputs/1211_metadata.rdata')
load('Inputs/1211_transcriptome.rdata')


metadata_analysis = df_long %>% 
  filter(Timepoint %in% c('D1','D4','D7')) %>%
  droplevels()

vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
         'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))
metadata_analysis %<>% left_join(tmp, by = 'HumanID')
counts_analysis = counts %>% dplyr::select(all_of(metadata_analysis$SampleID))

stopifnot(identical(metadata_analysis$SampleID, colnames(counts_analysis)))


## 设计矩阵：主效应+交互
design <- model.matrix(
  ~ Mortality28d * Timepoint +
    Gender + Age + CenterGroup + PneumoniaTypeGroup +
    CCI + SOFA_24h + Immunosuppression + MV,
  data = metadata_analysis
)

# 立刻检查：一定要满秩
qr(design)$rank
ncol(design)

# 让列名更好读：M0:D1 / M1:D1 ...
colnames(design) <- gsub("Mortality28d", "M", colnames(design))
colnames(design) <- gsub("Timepoint", "", colnames(design))

## 把列名变成 R 语法合法名字（处理 - 和 : 等）
colnames(design) <- make.names(colnames(design))

## calcNormFactors + voom
dge <- DGEList(counts = as.matrix(counts_analysis))
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = FALSE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

## 构建 contrasts（用合法后的列名）
colnames(design)
# 变成：M1.D4、M1.D7、PneumoniaTypeGroupnon.CAP 等

contrast_mat <- makeContrasts(
  ## 每个时间点：死亡 vs 存活
  D1 = M1,
  D4 = M1 + M1.D4,
  D7 = M1 + M1.D7,
  
  ## 动态差异（交互项本身）
  Dyn_D4vsD1 = M1.D4,
  Dyn_D7vsD1 = M1.D7,
  
  ## 早期总体差异（D1/D4/D7 等权平均）
  Overall = ( M1 + (M1 + M1.D4) + (M1 + M1.D7) ) / 3,
  
  levels = design
)


fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

res_D1 <- topTable(fit2, coef="D1", number=Inf, adjust.method="BH")
res_D4 <- topTable(fit2, coef="D4", number=Inf, adjust.method="BH")
res_D7 <- topTable(fit2, coef="D7", number=Inf, adjust.method="BH")
res_Dyn41 <- topTable(fit2, coef="Dyn_D4vsD1", number=Inf, adjust.method="BH")
res_Dyn71 <- topTable(fit2, coef="Dyn_D7vsD1", number=Inf, adjust.method="BH")
res_Overall <- topTable(fit2, coef="Overall", number=Inf, adjust.method="BH")

head(res_Dyn41)
head(res_Dyn71)

res = res_Overall
res = res_D7
plot_fun = function(res){
  res <- res[order(res$adj.P.Val, decreasing = F), ]
  res$gene_name <- as.character(gene_attr[rownames(res),'SYMBOL'])
  
  x = res[res$gene_name %in% c('BPI','DEFA','MMP8','CD3D','LDHA'),c(1,4,5,7)]
  x
  # logFC      P.Value    adj.P.Val gene_name
  # ENSG00000134333  0.3506174 3.657284e-06 0.0001043044      LDHA
  # ENSG00000118113  0.7079262 2.315865e-04 0.0023852160      MMP8
  # ENSG00000167286 -0.3515506 7.810301e-04 0.0059316362      CD3D
  # ENSG00000101425  0.2919730 3.238781e-02 0.0966444843       BPI
  
  # replace NA values with 1s and keep only significant genes
  res$adj.P.Val[is.na(res$adj.P.Val)] <- 1
  
  res_table <- data.frame(res)
  res_table$sig <- res_table$adj.P.Val < 0.1
  res_table$sig[res_table$adj.P.Val < 0.1 & res_table$logFC>0] <- "1"
  res_table$sig[res_table$adj.P.Val < 0.1 & res_table$logFC<0] <- "0"
  
  if (sum(res_table$sig!=FALSE)==0){
    label_data <- res_table[res_table$sig, ]
  }else{
    label_data <- res_table[1:min(c(25, nrow(res_table[res_table$sig, ]))), ]
    
    label_data_logfc <- res_table[order(abs(res_table$logFC), decreasing=T), ]
    label_data_logfc <- label_data_logfc[1:25, ][label_data_logfc[1:25, "sig"], ]
    
    label_data <- rbind(label_data, label_data_logfc)
    label_data <- unique(label_data)
  }
  
  # mortality_col_fun = c("0" = "#4575B4", "1" = "#D73027")
  
  p = ggplot(res_table, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = sig)) +
    scale_color_manual(breaks=c("1", "0", "FALSE"), values = alpha(c("#D73027", "#4575B4", "black"), 0.5)) + 
    ggrepel::geom_text_repel(data = label_data,
                             aes(label = gene_name),
                             max.overlaps=10) +
    labs(x="Log2 Fold Change", y="-Log10 adj.P") +
    theme_bw() +
    theme(legend.position = "none")
  
  {
    table(res_table$sig)
    # 上调
    deg_up = res_table %>% filter(sig == '1')
    gene_up <- rownames(deg_up)
    gene_up_df <- bitr(gene_up, 
                       fromType="ENSEMBL", 
                       toType="ENTREZID", 
                       OrgDb=org.Hs.eg.db)
    
    ego_up <- enrichGO(gene=gene_up_df$ENTREZID,
                       OrgDb=org.Hs.eg.db,
                       keyType="ENTREZID",
                       ont="BP",
                       pAdjustMethod="BH",
                       pvalueCutoff=0.05,
                       qvalueCutoff=0.05,
                       readable=TRUE)
    
    ego_up_simpl_0.7 <- simplify(
      ego_up,
      cutoff     = 0.7,          # 0.5~0.9 可调；0.7 常用
      by         = "p.adjust",
      select_fun = min,
      measure    = "Wang"        # BP 推荐 Wang
    )
    
    # 下调
    deg_down = res_table %>% filter(sig == '0')
    gene_down <- rownames(deg_down)
    gene_down_df <- bitr(gene_down, 
                         fromType="ENSEMBL", 
                         toType="ENTREZID", 
                         OrgDb=org.Hs.eg.db)
    ego_down <- enrichGO(gene=gene_down_df$ENTREZID,
                         OrgDb=org.Hs.eg.db,
                         keyType="ENTREZID",
                         ont="BP",
                         pAdjustMethod="BH",
                         pvalueCutoff=0.05,
                         qvalueCutoff=0.05,
                         readable=TRUE)
    
    ego_down_simpl_0.7 <- simplify(
      ego_down,
      cutoff     = 0.7,          # 0.5~0.9 可调；0.7 常用
      by         = "p.adjust",
      select_fun = min,
      measure    = "Wang"        # BP 推荐 Wang
    )

    
    top_combined <- bind_rows(
      as.data.frame(ego_up_simpl_0.7) %>% mutate(direction = "Up"),
      as.data.frame(ego_down_simpl_0.7) %>% mutate(direction = "Down")
    ) %>%
      mutate(p_signed = ifelse(direction == "Down", log10(p.adjust), -log10(p.adjust))) %>%
      group_by(direction) %>%
      slice_min(abs(p.adjust), n = 15) %>%
      ungroup()
    
    if(F){
      # D1
      top_combined$Description[2] = "TCR-mediated T cell activation"
      top_combined$Description[7] = "Antigen receptor–mediated adaptive immunity"
      
    }

    thr <- -log10(0.05)
    p2 = ggplot(top_combined,
                aes(x = p_signed,
                    y = reorder(Description, p_signed),
                    fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c("Up"="#D73027", "Down"="#4575B4")) +
      labs(x="-log10(adj.P)*Direction", y="") +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      geom_vline(xintercept = 0) +
      geom_vline(xintercept =  thr, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = -thr, linetype = "dashed", alpha = 0.5)
    
    }
  if(F){
    plist_res_overall = list()
    plist_res_overall[[1]] = p; plist_res_overall[[2]] = p2;
    plist1 = list(p,p2)
    plist2 = list(p,p2)
    plist3 = list(p,p2)
    
  }
  return(list(p, p2))
  
}  


pdf(paste0("Outputs/1215_DE.pdf"), width = 4.5, height = 5)
print(plist1[[1]])
print(plist2[[1]])
print(plist3[[1]])
print(plist_res_overall[[1]])
dev.off()

pdf(paste0("Outputs/1215_DEGO.pdf"), width = 8, height = 5)
print(plist1[[2]])
print(plist2[[2]])
print(plist3[[2]])
print(plist_res_overall[[2]])
dev.off()


#### 展现selected genes ####

res = res_Overall
res <- res[order(res$adj.P.Val, decreasing = F), ]
res$gene_name <- as.character(gene_attr[rownames(res),'SYMBOL'])

x = res[res$gene_name %in% c('BPI','DEFA','MMP8','CD3D','LDHA','CD4'),c(1,4,5,7)]
x
# logFC      P.Value    adj.P.Val gene_name
# ENSG00000134333  0.3506174 3.657284e-06 0.0001043044      LDHA
# ENSG00000118113  0.7079262 2.315865e-04 0.0023852160      MMP8
# ENSG00000167286 -0.3515506 7.810301e-04 0.0059316362      CD3D
# ENSG00000101425  0.2919730 3.238781e-02 0.0966444843       BPI

tmp = v$E
df_genes = tmp[rownames(x)[1:4],]
identical(colnames(df_genes), metadata_analysis$SampleID)
tmp = metadata_analysis %>% dplyr::select(SampleID, Mortality28d)
df_genes  = bind_cols(tmp, t(df_genes))
colnames(df_genes)[3:6] = x$gene_name[1:4]


df_long <- df_genes %>%
  pivot_longer(
    cols = -c(SampleID, Mortality28d),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    Mortality28d = factor(
      Mortality28d,
      levels = c(0, 1),
      labels = c("Survival", "Mortality")
    )
  )

df_long %<>% filter(Expression>0)
p <- ggplot(df_long, aes(x = Mortality28d, y = Expression, 
                         color = Mortality28d, fill = Mortality28d)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = 'black') +
  geom_jitter(width = 0.15, alpha = 0.1, size = 1) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("#4575B4", "#D73027")) +
  scale_fill_manual(values = alpha(c("#4575B4", "#D73027"), 0.5)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_blank()
  ) +
  ylab("Expression (logCPM)")

library(ggpubr)

pval_df <- x[1:4, ] %>%
  as.data.frame() %>%
  dplyr::select(gene_name, adj.P.Val) %>%
  dplyr::rename(Gene = gene_name) %>%
  mutate(
    group1 = "Survival",
    group2 = "Mortality",
    label = case_when(
      adj.P.Val < 0.001 ~ "FDR < 0.001",
      TRUE ~ paste0("FDR = ", sprintf("%.3f", adj.P.Val))
    )
  )


y_pos <- df_long %>%
  group_by(Gene) %>%
  summarise(
    y.position = max(Expression, na.rm = TRUE) * 1.05
  )

pval_df <- left_join(pval_df, y_pos, by = "Gene")

p +
  stat_pvalue_manual(
    pval_df,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    vjust = -0.3,
    size = 3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.12))
  ) -> p2


pdf(paste0("Outputs/1221_DE_selected_genes.pdf"), width = 6, height = 2)
print(p2)
dev.off()

