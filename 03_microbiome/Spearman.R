rm(list = ls())
library(tidyverse)
library(data.table)
library(Hmisc)

load('Inputs/1616_microbe.rdata')

{
  genus_sum2 = data
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


rownames(genus_sum3)[rownames(genus_sum3) == 'HHV-4'] = 'EBV'
# 计算菌属之间的 Spearman 相关
overallres = rcorr(t(genus_sum3), type = "spearman")


diag(overallres$P) <- 0


pdf("Outputs/03_microbe_cor.pdf", width = 6, height = 4.5)
corrplot::corrplot(overallres$r,  
                   p.mat = overallres$P, 
                   sig.level = 0.05, insig = "blank", 
                   tl.cex=0.8, type="lower", 
                   method = "color", outline=TRUE, 
                   col= rev(RColorBrewer::brewer.pal(n=10, name="PuOr")), 
                   tl.srt=45,  tl.offset = 0.3, tl.col="black", 
                   mar=c(0.1,0.1,0.1,0.1))

{
  # 2. 添加水平线
  # 找到 Aspergillus 对应的 y 坐标
  # corrplot 用的是 1:n 从上到下排列行名
  species <- rownames(overallres$r)
  n <- length(species)
  row_num <- which(species == "Aspergillus")
  
  # 横坐标范围
  x_left <- 0.5
  x_right <- row_num + 0.5
  y_pos_upper <- n - row_num + 1 + 0.5  # 上线
  y_pos_lower <- n - row_num + 1 - 0.5  # 下线
  
  # 画线
  rect(x_left, y_pos_lower, x_right, y_pos_upper, 
       border = "black", lwd=0.3)
  
}

dev.off()


