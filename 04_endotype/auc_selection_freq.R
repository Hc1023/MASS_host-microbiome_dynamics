rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(e1071)    # SVM 模型
library(caret)    # 数据切分 + 评价
library(pROC)     # ROC AUC
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(ConsensusTranscriptomicSubtype)
library("SepstratifieR")
library(ComplexHeatmap)
library(circlize)
load('Inputs/1616_meta_model.rdata')

#### selection prob heatmap ####

## heatmap



load('Inputs/1211_metadata.rdata')
load('Inputs/1211_transcriptome.rdata')

# new_expr_data must be a data.frame or matrix
# First column = gene identifiers (optional), rownames = Ensembl IDs
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE, log =T)

## CTS
result <- run_subtype_classifier(new_expr_data = cpm_norm)
new_mat <- result$expression_corrected$new   # genes x samples
rownames(new_mat) = gene_attr[rownames(new_mat),'SYMBOL']
ann_col <- result$predictions
ann_col$Timepoint = sub("^([^_]+_[^_]+)_(.*)$", "\\2", ann_col$Sample)
ann_col$Timepoint %<>% factor(., levels = paste0("D",c(1,4,7,14,21)))
ann_col %<>% dplyr::select(-Sample)

## SRS
srs_preds <- stratifyPatients(t(cpm_norm), gene_set = "extended")
tmp = srs_preds@SRS_probs
tmp$SRS <- max.col(tmp[,c(1,2)], ties.method = "first")
identical(rownames(tmp), rownames(ann_col))
ann_col$SRS = factor(tmp$SRS)
ann_col = ann_col[,c(1,3,2)]
col_fun = setNames(RColorBrewer::brewer.pal(6,'Set3')[-2], levels(ann_col$Timepoint))
scales::show_col(col_fun)

ann_colors <- list(
  CTS = c(
    `1` = "royalblue",
    `2` = "#B2DF8A",
    `3` = "orange"
  ),
  SRS = c(`1` = "#FDC086", `2` = "#BEAED4"),
  Timepoint = col_fun
)


ann_col %<>% arrange(CTS, SRS, Timepoint)
new_mat_ord <- new_mat[, rownames(ann_col)]

rn = c(
  "ACER3",
  "SERPINB1",
  "HK3",
  "TDRD9",
  "NLRC4",
  "PGD",
  "UBE2H",
  "METTL9",
  "STOM",
  "SNX3",
  "GADD45A",
  "BTN3A3",
  "BPGM",
  "CA1",
  "SLC4A1",
  "EPB42",
  "FECH",
  "GLRX5"
)

new_mat_ord = new_mat_ord[rn,]



pheatmap(
  new_mat_ord,
  color = colorRampPalette(rev(brewer.pal(11, "PuOr")))(50),
  scale = "row",
  cluster_rows = F,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  border_color = NA,
  treeheight_row = 0
)

## sel_CTSg for mortality28d D1/4/7
# run pre_fun function

## varialble selection 
### CTSg
selfun = function(d = "D1"){
  df_G1 = meta_model %>% dplyr::select(Mortality28d, SurvivalTimeWithin28Days, APACHEII_24h, SOFA_24h,
                                       starts_with(paste0(d,"_")))
  df_G1 %<>% na.omit()
  table(df_G1$Mortality28d)
  X = df_G1 %>% dplyr::select(starts_with(paste0(d,"_CTSg")))
  y = as.integer(as.character(df_G1$Mortality28d))
  B <- 100
  sel_mat <- matrix(0, nrow = ncol(X), ncol = B)
  set.seed(123)
  for (b in 1:B) {
    idx <- sample(1:nrow(X), replace = TRUE)
    fit <- cv.glmnet(as.matrix(X[idx, ]), y[idx], family = "binomial", alpha = 1)
    sel_mat[, b] <- as.numeric(coef(fit, s = "lambda.1se")[-1] != 0) 
  }
  rownames(sel_mat) = sub(paste0(d, '_CTSg_'), "", colnames(X))
  tmp = data.frame(SYMBOL = rownames(sel_mat), Freq = rowMeans(sel_mat)) 
  return(tmp)
}
seldf_D1 = selfun(d = "D1")
seldf_D4 = selfun(d = "D4")
seldf_D7 = selfun(d = "D7")
seldf_all = bind_cols(seldf_D1,seldf_D4$Freq, seldf_D7$Freq)
colnames(seldf_all)[2:4] = paste0("D",c(1,4,7))
# 对齐到 heatmap 的行（基因）顺序
freq_row <- seldf_all[,-1] %>%
  .[rownames(new_mat_ord), , drop = FALSE]




# 1) freq 矩阵：确保是 numeric matrix，行顺序对齐 new_mat_ord
freq_mat <- as.matrix(freq_row)   # 3 cols: D1 D4 D7

# 2) 颜色函数（统一配色：一个 legend）
col_freq <- colorRamp2(c(0, 1), c("white", "darkred"))

# 3) 左边表达热图颜色（你原来 PuOr + scale=row 的感觉）
# 4) 如果你需要“row scale”，先自己对 new_mat_ord 做 z-score
mat_z <- t(scale(t(new_mat_ord)))
mat_z[is.na(mat_z)] <- 0

pal_expr <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(11, "PuOr"))
)(50)
rng <- range(mat_z, na.rm = TRUE)
col_expr <- colorRamp2(
  seq(-5, 5, length.out = length(pal_expr)),
  pal_expr
)

top_ha <- HeatmapAnnotation(
  df  = ann_col,
  col = ann_colors,
  annotation_name_side = "left",
  # annotation_name_gp = gpar(fontsize = 9),
  gap = unit(1.5, "mm")
)


ht_expr <- Heatmap(
  mat_z,
  name = "Z-score",
  col = col_expr,
  top_annotation = top_ha, 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  border = FALSE
)

ht_freq <- Heatmap(
  freq_mat,
  name = "Freq",          # <- 只有一个 legend，覆盖 D1/D4/D7 三列
  col = col_freq,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T,
  column_title = NULL,
  border = FALSE,
  width = unit(12, "mm")  # 右侧热图窄一点（可调）
)

p = draw(ht_expr + ht_freq, heatmap_legend_side = "right", annotation_legend_side = "right")


pdf("Outputs/04_CTSg_heatmap.pdf", width = 6.6, height = 4)
print(p)
dev.off()

if(F){
  write.csv(seldf_all, 'Outputs/Supplementary_data_selection_ctsg.csv',
            row.names = F)
}

#### microbe selection prob heatmap ####
selfun = function(d = "D1"){
  df_G1 = meta_model %>% dplyr::select(Mortality28d, SurvivalTimeWithin28Days, APACHEII_24h, SOFA_24h,
                                       starts_with(paste0(d,"_")))
  df_G1 %<>% na.omit()
  table(df_G1$Mortality28d)
  X = df_G1 %>% dplyr::select(starts_with(paste0(d,"_mi")))
  y = as.integer(as.character(df_G1$Mortality28d))
  B <- 100
  sel_mat <- matrix(0, nrow = ncol(X), ncol = B)
  set.seed(123)
  for (b in 1:B) {
    idx <- sample(1:nrow(X), replace = TRUE)
    fit <- cv.glmnet(as.matrix(X[idx, ]), y[idx], family = "binomial", alpha = 1)
    sel_mat[, b] <- as.numeric(coef(fit, s = "lambda.1se")[-1] != 0) 
  }
  rownames(sel_mat) = sub(paste0(d, '_mi_'), "", colnames(X))
  tmp = data.frame(SYMBOL = rownames(sel_mat), Freq = rowMeans(sel_mat)) 
  return(tmp)
}
seldf_D1 = selfun(d = "D1")
seldf_D4 = selfun(d = "D4")
seldf_D7 = selfun(d = "D7")
seldf_all = bind_cols(seldf_D1,seldf_D4$Freq, seldf_D7$Freq)
rownames(seldf_all) = gsub('D1_mi2_','',rownames(seldf_all))
colnames(seldf_all)[2:4] = paste0("D",c(1,4,7))
seldf_all = seldf_all[,-1]
freq_mat = as.matrix(seldf_all)
rownames(freq_mat) = gsub('\\.','-',rownames(freq_mat))
rownames(freq_mat)[rownames(freq_mat) == 'Influenza-A'] = 'Influenza A'
rownames(freq_mat)[rownames(freq_mat) == 'HHV-4'] = 'EBV'
rownames(freq_mat)[rownames(freq_mat) == 'Bacf'] = 'Bacteria/Fungi'
col_freq <- colorRamp2(c(0, 1), c("white", "darkred"))

ht_freq <- Heatmap(
  freq_mat,
  name = "Freq",          # <- 只有一个 legend，覆盖 D1/D4/D7 三列
  col = col_freq,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = T,
  column_title = NULL,
  border = FALSE,
  width = unit(12, "mm")  # 右侧热图窄一点（可调）
)

pdf("Outputs/04_mi_heatmap.pdf", width = 3, height = 4)
print(ht_freq)
dev.off()

if(F){
  write.csv(seldf_all, 'Outputs/Supplementary_data_selection_microbe.csv')
}




