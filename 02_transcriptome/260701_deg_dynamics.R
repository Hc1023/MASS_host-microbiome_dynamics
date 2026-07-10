rm(list = ls())
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


project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)

out_dir <- file.path(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260701_deg_dyn"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

## using the outputs in deg analysis
sheets = sisiUtils::read_excel_allsheets('Outputs/Supplementary_data_1_4_DE.xlsx')
res_Dyn41 = sheets[['Dyn_D4vsD1']]
rownames(res_Dyn41) = res_Dyn41$rowname
res_Dyn41 = res_Dyn41[,-1]
res_Dyn71 = sheets[['Dyn_D7vsD1']]
rownames(res_Dyn71) = res_Dyn71$rowname
res_Dyn71 = res_Dyn71[,-1]

#### GSEA ####
library(msigdbr)
library(fgsea)
library(clusterProfiler)

# 准备 GO BP 基因集（MSigDB C5:BP）去GSEA官网下载
geneset <- read.gmt("../../AI_MASS/Utility/c5.go.bp.v2025.1.Hs.symbols.gmt")
# gene sets: list(term -> gene symbols)
bp_sets <- split(geneset$gene, geneset$term)

# 转换基因名字为SYMBOL

make_ranks_symbol <- function(res_tt, gene_annot,
                              id_col = "ENSEMBL", sym_col = "SYMBOL",
                              stat = c("t", "logFC")) {
  stat <- match.arg(stat)
  
  df <- res_tt %>%
    as.data.frame() %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(gene_annot, by = "ENSEMBL") %>%
    filter(!is.na(SYMBOL))
  
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


gene_annot = gene_attr
colnames(gene_annot) =  c("ENSEMBL", "SYMBOL")
gene_annot = gene_annot[!grepl('ENSG', gene_annot$SYMBOL),]

ranks_Dyn41 <- make_ranks_symbol(res_Dyn41, gene_annot, stat = "t")
set.seed(1)
gsea_Dyn41  <- run_fgsea_bp(ranks_Dyn41, bp_sets)
ranks_Dyn71 <- make_ranks_symbol(res_Dyn71, gene_annot, stat = "t")
set.seed(1)
gsea_Dyn71  <- run_fgsea_bp(ranks_Dyn71, bp_sets)

if(F){
  library(openxlsx)
  
  # put results into a list
  res_list <- list(
    gsea_Dyn41 = gsea_Dyn41[,-ncol(gsea_Dyn41)],
    gsea_Dyn71 = gsea_Dyn71[,-ncol(gsea_Dyn71)]
  )
  
  # write to excel
  write.xlsx(
    res_list,
    file = file.path(out_dir, "Supplementary_data_GSEA_dynamics.xlsx"),
    asTable = TRUE
  )
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
    arrange(padj, NES) %>%          # NES 越负越靠前
    slice_head(n = top_n) %>%
    mutate(direction = "Down")
  
  top_combined <- bind_rows(up, down) %>%
    mutate(pathway2 = str_wrap(pathway, width = wrap_width)) %>%
    # 让排序按 NES 从小到大（画出来更直观：Down 在下面，Up 在上面）
    arrange(NES) %>%
    mutate(pathway2 = factor(pathway2, levels = pathway2))
  
  p <- ggplot(top_combined,
              aes(x = NES, y = pathway2, fill = direction)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
    labs(title = title, x = "NES", y = NULL) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if (!is.null(thr)) {
    p <- p +
      geom_vline(xintercept =  thr, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = -thr, linetype = "dashed", alpha = 0.5)
  }
  
  return(p)
}

p1 <- plot_fgsea_bidir(gsea_Dyn41, top_n = 15,
                       title = "Dyn D4 vs D1 (interaction)",
                       padj_cut = 0.05)

p2 <- plot_fgsea_bidir(gsea_Dyn71, top_n = 15,
                       title = "Dyn D7 vs D1 (interaction)",
                       padj_cut = 0.10)

p1
p2


get_pos <- function(gsea_res, tp, padj_cut = 0.10) {
  out <- gsea_res %>%
    filter(!is.na(NES), !is.na(padj), NES > 0, padj < padj_cut) %>%
    arrange(padj) %>%
    mutate(Timepoint = tp)
  out$pathway = gsub('GOBP_', '', out$pathway)
  drop_pat <- "^(SENSORY_|DETECTION_)"
  if (!is.null(drop_pat)) {
    out <- out %>% filter(!grepl(drop_pat, pathway))
  }
  out
}

gsea_pos_D7 <- get_pos(gsea_Dyn71, "D7", padj_cut = 1)
gsea_pos_D4 <- get_pos(gsea_Dyn41, "D4", padj_cut = 1)

# union of top10
top10_union <- union(
  gsea_pos_D4 %>% slice_head(n = 10) %>% pull(pathway),
  gsea_pos_D7 %>% slice_head(n = 10) %>% pull(pathway)
)

df_long <- bind_rows(
  gsea_pos_D4 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D4"),
  gsea_pos_D7 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D7")
) %>%
  filter(pathway %in% top10_union) %>%
  distinct(pathway, Timepoint, .keep_all = TRUE)

# NES matrix
mat_nes <- df_long %>%
  dplyr::select(pathway, Timepoint, NES) %>%
  pivot_wider(names_from = Timepoint, values_from = NES) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# padj matrix
mat_padj <- df_long %>%
  dplyr::select(pathway, Timepoint, padj) %>%
  pivot_wider(names_from = Timepoint, values_from = padj) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# 行顺序（按 NES 总体强度）
ord <- order(rowSums(mat_nes, na.rm = TRUE), decreasing = TRUE)
mat_nes  <- mat_nes[ord, , drop = FALSE]
mat_padj <- mat_padj[ord, , drop = FALSE]


# NES 发散色
rng <- max(abs(mat_nes), na.rm = TRUE)
col_fun <- colorRamp2(
  c(0, 3),
  c("white", "#D73027")
)

# padj → symbol
sig_symbol <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.05, "*",
                ifelse(p < 0.10, "·", "")))
}


go_label_map_pos <- c(
  MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY =
    "MyD88-dependent TLR signaling",
  REGULATION_OF_LYSOSOMAL_LUMEN_PH =
    "regulation of lysosomal lumen PH",
  RESPONSE_TO_TYPE_I_INTERFERON = 
    "response to IFN-1",
  PHAGOLYSOSOME_ASSEMBLY =
    "phagolysosome assembly",
  VACUOLAR_ACIDIFICATION =
    "vacuolar acidification",
  PHAGOSOME_MATURATION =
    "phagosome maturation",
  TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY =
    "TNF-mediated signaling",
  REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS =
    "regulation of response to cytokine stimulus",
  REGULATORY_NCRNA_MEDIATED_GENE_SILENCING =
    "regulatory ncRNA-mediated gene silencing",
  ENDOLYSOSOMAL_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY =
    "endolysosomal TLR signaling",
  REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS =
    "regulation of epithelial cell apoptotic process",
  CYTOKINE_MEDIATED_SIGNALING_PATHWAY =
    "cytokine-mediated signaling pathway",
  CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY =
    "cytoplasmic PRR signaling",
  CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION =
    "canonical NF-kB signal transduction",
  NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS =
    "negative regulation of epithelial apoptosis",
  REGULATION_OF_HUMORAL_IMMUNE_RESPONSE = 
    "regulation of humoral immune response",
  POSITIVE_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS =
    "positive regulation of response to biotic stimulus",
  REGULATION_OF_INNATE_IMMUNE_RESPONSE =
    "regulation of innate immune response",
  POSITIVE_REGULATION_OF_DEFENSE_RESPONSE =
    "positive regulation of defense response",
  VACUOLE_ORGANIZATION =
    "vacuole organization",
  REGULATION_OF_INFLAMMATORY_RESPONSE =
    "regulation of inflammatory response",
  MACROAUTOPHAGY =
    "macroautophagy",
  VESICLE_ORGANIZATION =
    "vesicle organization"
)

rownames(mat_nes)  <- go_label_map_pos[rownames(mat_nes)]
rownames(mat_padj) <- rownames(mat_nes)

rn_w <- max_text_width(rownames(mat_nes), gp = gpar(fontsize = 8))

ht1 = Heatmap(
  mat_nes,
  name = "NES",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  na_col = "grey90",
  column_names_rot = 0,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),              # 字体小一点
  row_names_max_width = rn_w + unit(6, "mm"),     # 给 rowname 留足空间（关键）
  column_title = "Trajectory up",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sig_symbol(mat_padj[i, j]),
      x, y,
      gp = gpar(fontsize = 14)
    )
  }
)

ht1

pdf(paste0("Outputs/02_DynTop_up.pdf"), width = 3.6, height = 4)
print(ht1)
dev.off()

get_neg <- function(gsea_res, tp, padj_cut = 0.10) {
  out <- gsea_res %>%
    filter(!is.na(NES), !is.na(padj), NES < 0, padj < padj_cut) %>%
    arrange(padj) %>%
    mutate(Timepoint = tp)
  out$pathway = gsub('GOBP_', '', out$pathway)
  drop_pat <- "^(SENSORY_|DETECTION_)"
  if (!is.null(drop_pat)) {
    out <- out %>% filter(!grepl(drop_pat, pathway))
  }
  out
}

gsea_neg_D7 <- get_neg(gsea_Dyn71, "D7", padj_cut = 1)
gsea_neg_D4 <- get_neg(gsea_Dyn41, "D4", padj_cut = 1)

# union of top10
top10_union <- union(
  gsea_neg_D4 %>% slice_head(n = 10) %>% pull(pathway),
  gsea_neg_D7 %>% slice_head(n = 10) %>% pull(pathway)
)

df_long <- bind_rows(
  gsea_neg_D4 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D4"),
  gsea_neg_D7 %>% dplyr::select(pathway, NES, padj) %>% mutate(Timepoint = "D7")
) %>%
  filter(pathway %in% top10_union) %>%
  distinct(pathway, Timepoint, .keep_all = TRUE)

# NES matrix
mat_nes <- df_long %>%
  dplyr::select(pathway, Timepoint, NES) %>%
  pivot_wider(names_from = Timepoint, values_from = NES) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix() %>% na.omit()

# padj matrix
mat_padj <- df_long %>%
  dplyr::select(pathway, Timepoint, padj) %>%
  pivot_wider(names_from = Timepoint, values_from = padj) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix() %>% na.omit()

# 行顺序（按 NES 总体强度）
ord <- order(rowSums(mat_nes, na.rm = TRUE))
mat_nes  <- mat_nes[ord, , drop = FALSE]
mat_padj <- mat_padj[ord, , drop = FALSE]


# NES 发散色
rng <- max(abs(mat_nes), na.rm = TRUE)
col_fun <- colorRamp2(
  c(-rng, 0),
  c("#4575B4", "white")
)

# padj → symbol
sig_symbol <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.05, "*",
                ifelse(p < 0.10, "·", "")))
}


go_label_map <- c(
  RIBOSOME_BIOGENESIS = "ribosome biogenesis",
  RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS = "ribosomal small subunit biogenesis",
  RRNA_PROCESSING = "rRNA processing",
  RRNA_METABOLIC_PROCESS = "rRNA metabolic process",
  INNATE_IMMUNE_RESPONSE_IN_MUCOSA = "innate immune response in mucosa",
  NUCLEOSOME_ORGANIZATION = "nucleosome organization",
  ANTIBACTERIAL_HUMORAL_RESPONSE = "antibacterial humoral response",
  ANTIMICROBIAL_HUMORAL_RESPONSE = "antimicrobial humoral response",
  ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE = "AMP-mediated humoral immune response",
  RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS = "ribonucleoprotein complex biogenesis",
  TRNA_METABOLIC_PROCESS = "tRNA metabolic process",
  MITOCHONDRIAL_GENE_EXPRESSION = 'mitochondrial gene expression',
  MITOCHONDRIAL_TRANSLATION = "mitochondrial translation",
  RNA_MODIFICATION = "RNA modification",
  TRNA_PROCESSING = "tRNA processing",
  DISRUPTION_OF_ANATOMICAL_STRUCTURE_IN_ANOTHER_ORGANISM =
    "disruption of host tissue structure",
  PROTEIN_DNA_COMPLEX_ORGANIZATION = "protein-DNA complex organization",
  HUMORAL_IMMUNE_RESPONSE = "humoral immune response",
  MITOTIC_NUCLEAR_DIVISION = "mitotic nuclear division",
  ORGANELLE_FISSION = "organelle fission"
)

rownames(mat_nes)  <- go_label_map[rownames(mat_nes)]
rownames(mat_padj) <- rownames(mat_nes)

rn_w <- max_text_width(rownames(mat_nes), gp = gpar(fontsize = 8))

ht2 = Heatmap(
  mat_nes,
  name = "NES",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  na_col = "grey90",
  column_names_rot = 0,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),              # 字体小一点
  row_names_max_width = rn_w + unit(6, "mm"),     # 给 rowname 留足空间（关键）
  column_title = "Trajectory down",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sig_symbol(mat_padj[i, j]),
      x, y,
      gp = gpar(fontsize = 14)
    )
  }
)

pdf(paste0("Outputs/02_DynTop_down.pdf"), width = 3.5, height = 4)
print(ht2)
dev.off()


#### trajectory ####
library(GSVA)
library(GSEABase)

meta_gsva = metadata_analysis %>% dplyr::select(HumanID,Timepoint,SampleID,Mortality28d)

## v$E: genes × samples（行名是 ENSEMBL 或 SYMBOL？这里用 SYMBOL）
expr <- v$E
# 行：gene (Ensembl ID)
# 列：SampleID
sym <- gene_attr[rownames(expr), "SYMBOL"]
is_valid_symbol <- !is.na(sym) & !grepl("^ENSG", sym)
expr_sym <- expr[is_valid_symbol, , drop = FALSE]
rownames(expr_sym) <- sym[is_valid_symbol]
sym_tab <- table(rownames(expr_sym))
unique_sym <- names(sym_tab[sym_tab == 1])
expr_sym <- expr_sym[rownames(expr_sym) %in% unique_sym, , drop = FALSE]
## 如果 rownames 还是 ENSEMBL，需要先映射到 SYMBOL
## 假设你之前已经在 GSEA 时做过 ENSG -> SYMBOL，这里直接用 SYMBOL
## 确保 rownames(expr) 是 SYMBOL

## 选择你要展示的 GO term
### up1: Inflammatory_signaling
up1_pathways <- c(
  "GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION",
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY"
)

### up2: Phagolysosome_function
up2_pathways <- c(
  "GOBP_PHAGOLYSOSOME_ASSEMBLY",
  "GOBP_PHAGOSOME_MATURATION",
  "GOBP_VACUOLE_ORGANIZATION",
  "GOBP_MACROAUTOPHAGY"
)
up3_pathways <- c(
  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON"
)

### down1: Ribosome_biogenesis
dw1_pathways <- c(
  "GOBP_RIBOSOME_BIOGENESIS",
  "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
  "GOBP_RRNA_PROCESSING",
  "GOBP_RRNA_METABOLIC_PROCESS"
)
### down2: Ribosome_biogenesis
dw2_pathways <- c(
  "GOBP_NUCLEOSOME_ORGANIZATION",
  "GOBP_PROTEIN_DNA_COMPLEX_ORGANIZATION",
  "GOBP_MITOTIC_NUCLEAR_DIVISION",
  "GOBP_ORGANELLE_FISSION"
)


## 从 gsea_Dyn41 / gsea_Dyn71 里取 leading-edge
extract_leading_edge <- function(gsea_res, pathways) {
  gsea_res %>%
    filter(pathway %in% pathways) %>%
    dplyr::select(pathway, leadingEdge) %>%
    tidyr::unnest(leadingEdge) %>%
    distinct(pathway, leadingEdge)
}
get_leading_edge_intersect <- function(gsea_D4, gsea_D7, pathways) {
  le_D4 <- extract_leading_edge(gsea_D4, pathways)$leadingEdge
  le_D7 <- extract_leading_edge(gsea_D7, pathways)$leadingEdge
  intersect(le_D4, le_D7)
}

modules <- list(
  up1 = up1_pathways,
  up2 = up2_pathways,
  up3 = up3_pathways,
  dw1 = dw1_pathways
)

gs_list <- lapply(
  modules,
  function(pw) get_leading_edge_intersect(gsea_Dyn41, gsea_Dyn71, pw)
)


gsva_param <- gsvaParam(
  expr = as.matrix(expr_sym),
  geneSets = gs_list,
  kcdf = "Gaussian",
  minSize = 10,
  maxSize = 500
)
gsva_mat <- gsva(gsva_param)

if(F){
  write.csv(gsva_mat, 'Outputs/Supplementary_data_module_gsva.csv')
}
identical(colnames(gsva_mat), meta_gsva$SampleID)

plot_df = bind_cols(meta_gsva[,1:4], t(gsva_mat))

plot_long <- plot_df %>%
  pivot_longer(
    cols = -c(HumanID, Timepoint, SampleID, Mortality28d),
    names_to = "Module",
    values_to = "Score"
  )

table(plot_long$Module)
plot_long <- plot_long %>%
  mutate(
    Module = factor(
      Module,
      levels = c("up1", "up2", "up3", "dw1"),
      labels = c(
        "Inflammatory signaling",
        "Phagolysosome function",
        "IFN signaling",
        "Ribosome biogenesis"
      )
    )
  )
library(rstatix)
wilcox_df <- plot_long %>%
  group_by(Module, Timepoint) %>%
  wilcox_test(
    Score ~ Mortality28d
  ) %>%
  ungroup() %>%
  mutate(
    p.label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ "·",
      TRUE      ~ "ns"
    )
  )
y_pos_df <- plot_long %>%
  group_by(Module, Timepoint, Mortality28d) %>%
  summarise(
    m  = mean(Score, na.rm = TRUE),
    se = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
    .groups = "drop"
  ) %>%
  group_by(Module, Timepoint) %>%
  summarise(
    y.position = max(m + se, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 避免贴着 error bar
  left_join(
    plot_long %>%
      group_by(Module, Timepoint) %>%
      summarise(rng = diff(range(Score, na.rm = TRUE)), .groups = "drop"),
    by = c("Module", "Timepoint")
  ) %>%
  mutate(y.position = y.position + 0.01 * rng) %>%
  select(Module, Timepoint, y.position)

wilcox_df <- wilcox_df %>%
  left_join(y_pos_df, by = c("Module", "Timepoint"))




p = ggplot(plot_long, aes(
  x = Timepoint,
  y = Score,
  color = factor(Mortality28d),
  group = interaction(HumanID, Mortality28d)
)) +
  # geom_line(alpha = 0.25) +
  stat_summary(
    aes(group = Mortality28d),
    fun = mean,
    geom = "line",
    linewidth = 1.3
  ) +
  stat_summary(
    aes(group = Mortality28d),
    fun = mean,
    geom = "point",
    size = 3
  ) +
  stat_summary(
    aes(group = Mortality28d),
    fun.data = mean_se,      
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.6
  ) +
  stat_pvalue_manual(
    wilcox_df,
    label = "p.label",
    x = "Timepoint",
    y.position = "y.position",
    tip.length = 0,
    size = 3
  ) +
  theme_bw() +
  scale_color_manual(values = alpha(c("0" = "#4575B4", "1" = "#D73027"), 0.9),
                     name   = "",
                     labels = c("Survival", "Mortality")) +
  facet_wrap(~ Module, scales = "free_y", ncol = 4) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) 


pdf(paste0("Outputs/02_Dyntraj.pdf"), width = 9, height = 2.2)
print(p)
dev.off()


#### table output ####

library(tidyverse)

## input: gs_list
## gs_list$up1, gs_list$up2, gs_list$up3, gs_list$dw1

module_names <- c(
  up1 = "Module1_PRR_TLR_TNF_NFkB",
  up2 = "Module2_Phagolysosome_Autophagy",
  up3 = "Module3_IFN_I_Response",
  dw1 = "Module4_Ribosome_Biogenesis"
)

## clean gene lists, preserving original order
gs_clean <- gs_list[names(module_names)] %>%
  lapply(function(x) unique(na.omit(trimws(as.character(x)))))

names(gs_clean) <- module_names[names(gs_clean)]

## convert to wide-format table
max_len <- max(lengths(gs_clean))

module_gene_wide <- gs_clean %>%
  lapply(function(x) {
    length(x) <- max_len
    x
  }) %>%
  as.data.frame(check.names = FALSE)

## save
write.csv(
  module_gene_wide,
  file = file.path(out_dir, "selected_genes_by_module_wide.csv"),
  row.names = FALSE,
  na = ""
)
library(ggVennDiagram)
p_venn <- ggVennDiagram(gs_clean, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Overlap of selected genes across modules",
    fill = "Gene count"
  ) 
p_venn

pdf(file.path(out_dir, "venn.pdf"), width = 5, height = 5)
print(p_venn)
dev.off()
