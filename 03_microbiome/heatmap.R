rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)
library(ComplexHeatmap)
library(circlize)

load('Inputs/1616_microbe.rdata')
load('Inputs/1211_metadata.rdata')

data_log2 = log2(data + 1)
# name_ = 'Survival'; virus = 1
myplot = function(){
  df_long1 = df_long
  genus_sum2 = data_log2[, df_long1$SampleID]

  {
    # 1. 给每行定义分类
    group_vec <- mapid_vec[rownames(genus_sum2)]
    # 2. 计算 prevalence（出现次数）
    prev <- rowSums(genus_sum2 > 0)
    # 3. 拆分
    bac_df <- genus_sum2[group_vec == "Bacteria", , drop = FALSE]
    fun_df <- genus_sum2[group_vec == "Fungi", , drop = FALSE]
    virus_df <- genus_sum2[group_vec == "Viruses", , drop = FALSE]
    # 4. 对两个组分别按 prevalence 排序
    virus_ordered <- virus_df[order(prev[rownames(virus_df)], decreasing = TRUE), ]
    bac_ordered <- bac_df[order(prev[rownames(bac_df)], decreasing = TRUE), ]
    # 5. 各取前 10（如果不足 10 自动取所有）
    virus_top10 <- head(virus_ordered, 10)
    bac_top10 <- head(bac_ordered, 10)
    # 6. 合并（保持顺序）
    genus_sum3 <- rbind(bac_top10, fun_df, virus_top10)
  }
  
  df_long2 = cbind(df_long1, t(genus_sum3))
  df_long2$Timepoint = factor(df_long2$Timepoint, 
                              levels = c('D1','D4','D7','D14','D21'))
  df_long2 <- df_long2 %>%
    arrange(Timepoint, Mortality28d, desc(across(5:ncol(df_long2))))
  genus_sum3 = df_long2[,5:ncol(df_long2)]
  # rownames(genus_sum3) = genus_sum3[,1]
  genus_sum3 = t(genus_sum3)
  
  rownames(genus_sum3)[rownames(genus_sum3) == 'HHV-4'] = 'EBV'
  names(mapid_vec)[names(mapid_vec) == 'HHV-4'] = 'EBV'
  expr_col_fun <- colorRamp2(
    c(0, median(genus_sum3[genus_sum3>0])/2, max(genus_sum3)),   # centered at 0
    c("white", "#FF8080", "red")
  )
  # RColorBrewer::display.brewer.pal(6,'Set3')
  # RColorBrewer::brewer.pal(4,'Set2')
  col_fun = setNames(RColorBrewer::brewer.pal(6,'Set3')[-2], levels(df_long2$Timepoint))

  mortality_col_fun = c("0" = "#4575B4", "1" = "#D73027")
  
  top_ha <- HeatmapAnnotation(
    Timepoint = df_long2$Timepoint,
    Mortality28d  = df_long2$Mortality28d,
    col = list(Timepoint = col_fun,
               Mortality28d = mortality_col_fun),
    show_legend = FALSE
  )
  row_group <- mapid_vec[rownames(genus_sum3)] %>% tools::toTitleCase()
  ht = Heatmap(
    as.matrix(genus_sum3),
    col = expr_col_fun,
    name = paste0("\nlog2(mass+1)"),
    cluster_rows = F,       # cluster samples
    cluster_columns = FALSE,
    show_column_names = FALSE,
    top_annotation = top_ha,
    column_split = df_long2$Timepoint,
    row_split = row_group,
    border = T,
    row_names_gp = gpar(fontsize = 10))
  # ht
  return(ht)
}

p1 = myplot()

lgd_mort <- Legend(
  title = "Mortality28d",
  labels = c("Survival", "Mortality"),
  legend_gp = gpar(fill = c("0" = "#4575B4", "1" = "#D73027"))
)

pdf("Outputs/03_microbe_heatmap.pdf", width = 6.6, height = 4)
print(p1)
grid::grid.newpage()
grid::grid.draw(lgd_mort)
dev.off()


