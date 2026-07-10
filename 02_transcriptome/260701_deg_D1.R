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

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)

out_dir <- file.path(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260701_deg_D1"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


load('Inputs/1211_metadata.rdata')
load('Inputs/1211_transcriptome.rdata')

# metadata_analysis = df_long %>% 
#   filter(Timepoint %in% c('D1','D4','D7')) %>%
#   droplevels()



get_res = function(tp){
  metadata_analysis = df_long %>% 
    filter(Timepoint == tp) %>%
    droplevels()
  vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
           'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
  tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))
  metadata_analysis %<>% left_join(tmp, by = 'HumanID')
  counts_analysis = counts %>% dplyr::select(all_of(metadata_analysis$SampleID))
  
  stopifnot(identical(metadata_analysis$SampleID, colnames(counts_analysis)))
  
  
  ## 设计矩阵：主效应+交互
  design <- model.matrix(
    ~ Mortality28d +
      Gender + Age + CenterGroup + PneumoniaTypeGroup +
      CCI + SOFA_24h + Immunosuppression + MV,
    data = metadata_analysis
  )
  
  # 立刻检查：一定要满秩
  qr(design)$rank
  ncol(design)
  
  # 让列名更好读：M0:D1 / M1:D1 ...
  colnames(design) <- gsub("Mortality28d", "M", colnames(design))
  
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
    M1 = M1,
    levels = design
  )
  
  
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)
  
  res <- topTable(fit2, coef="M1", number=Inf, adjust.method="BH")
  res$gene_symbol <- as.character(gene_attr[rownames(res),"SYMBOL"])
  return(res)
}

res_D1 = get_res("D1")
res_D4 = get_res("D4")
res_D7 = get_res("D7")

write.csv(res_D1, file = file.path(out_dir, "deg_D1.csv"))
write.csv(res_D4, file = file.path(out_dir, "deg_D4.csv"))
write.csv(res_D7, file = file.path(out_dir, "deg_D7.csv"))

plot_fun = function(res){
  res <- res[order(res$adj.P.Val, decreasing = F), ]
  # genes_interest = c('BPI','DEFA','MMP8','CD3D','CD4','LDHA',
  #                    'S100A8', 'S100A9', 'IL1B', 'CXCL8',
  #                    'NFKBIA', 'TNFAIP3', 'IFI27')
  # x = res[res$gene_symbol %in% genes_interest,c(1,4,5,7)]
  # x
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
  p2

  return(p2)
  
}  


go_D1 = plot_fun(res_D1)
go_D4 = plot_fun(res_D4)
go_D7 = plot_fun(res_D7)

pdf(file.path(out_dir, "go.pdf"), width = 6.5, height = 5)
print(go_D1)
print(go_D4)
print(go_D7)
dev.off()


