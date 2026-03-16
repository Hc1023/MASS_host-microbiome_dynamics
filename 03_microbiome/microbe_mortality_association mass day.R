rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(gridExtra)

load('Inputs/1211_microbe.rdata')
load('Inputs/1211_metadata.rdata')

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
  genus_sum = genus_sum3
  pathogens = rownames(genus_sum3)
  # save(data, df_long, metadata, mapid_vec, pathogens, file = 'Inputs/1124_microbe_data.rdata')
  
}

d = 'D1'
myplot = function(d){
  # 1. select samples and pathogens
  df_long1 = df_long[df_long$Timepoint == d,]
  genus_sum3 = genus_sum[, df_long1$SampleID]
  
  genus_long <- genus_sum3 %>%
    tibble::rownames_to_column(var = "Pathogen") %>%   # move rownames into a column
    pivot_longer(
      cols = -Pathogen,            # all other columns
      names_to = "SampleID",
      values_to = "Abundance"
    )
  genus_long %<>% filter(Abundance>0)
  genus_long$Pathogen %<>% factor(.,levels = pathogens)
  genus_long %<>% left_join(df_long %>% select(Mortality28d, SampleID), by = 'SampleID')
  genus_long$Abundance = log2(genus_long$Abundance+1)
  genus_long %<>%
    mutate(Abundance = ifelse(Mortality28d == "1", -Abundance, Abundance))
  
  results <- genus_long %>%
    group_by(Pathogen) %>%
    summarise(
      n0 = sum(Mortality28d == 0),
      n1 = sum(Mortality28d == 1),
      P_value = ifelse(
        length(unique(Mortality28d)) < 2 | min(n0, n1) < 3,  # 样本太少则跳过
        NA,
        wilcox.test(abs(Abundance) ~ Mortality28d, exact = FALSE)$p.value
      )
    ) %>% as.data.frame()

  # 添加标记列
  results <- results %>%
    mutate(sig_label = case_when(
      P_value < 0.05 ~ "*",
      P_value < 0.1 ~ ".",
      TRUE ~ ""
    ))
  
  # 找每个 Pathogen 的最大条形高度，用于放置标记
  max_y <- genus_long %>%
    group_by(Pathogen) %>%
    summarise(max_y = max(Abundance))
  
  # 合并显著性标记
  annot_df <- max_y %>%
    left_join(results %>% select(Pathogen, sig_label), by = "Pathogen")
  
  annot_df$max_y = ifelse(annot_df$sig_label == '*',
                          annot_df$max_y * 0.9,
                          annot_df$max_y * 1.1)
  
  # 2. 绘图
  df_longplot = genus_long %>% filter(Pathogen %in% pathogens[1:11])
  p1 = ggplot(df_longplot, aes(x = Pathogen, y = Abundance, fill = Mortality28d)) +
    geom_boxplot(position = position_identity(),
                 width = 0.7, color = "black", linewidth = 0.3,
                 outlier.alpha = 0.2) +
    # geom_jitter(width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = c("0" = "#4575B4", "1" = "#D73027")) +
    labs(x = "", y = "", fill = "") +
    theme_bw() + labs(title = d) +
    scale_y_continuous(labels = function(x) abs(x), # y轴显示正数
                       limits = c(-6,6)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = 'None',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background  = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.margin = margin(l = 0.5, r = 0.1, t = 0.1, unit = "cm")) +
    geom_text(data = annot_df[1:11,], aes(x = Pathogen, y = max_y + 0.01, label = sig_label),
              inherit.aes = FALSE, vjust = 0, size = 5)
  # p1
  df_longplot = genus_long %>% filter(Pathogen %in% pathogens[12:21])
  levels(df_longplot$Pathogen)[levels(df_longplot$Pathogen) == 'HHV-4'] = 'EBV'
  levels(annot_df$Pathogen)[levels(annot_df$Pathogen) == 'HHV-4'] = 'EBV'
  p2 = ggplot(df_longplot, aes(x = Pathogen, y = Abundance, fill = Mortality28d)) +
    geom_boxplot(position = position_identity(),
                 width = 0.7, color = "black", linewidth = 0.3,
                 outlier.alpha = 0.2) +
    # geom_jitter(width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = c("0" = "#4575B4", "1" = "#D73027")) +
    labs(x = "", y = "", fill = "") +
    theme_bw() + labs(title = d) +
    scale_y_continuous(labels = function(x) abs(x), # y轴显示正数
                       limits = c(-15,15)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = 'None',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background  = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          plot.margin = margin(l = 0.5, r = 0.1, t = 0.1, unit = "cm")) +
    geom_text(data = annot_df[12:21,], aes(x = Pathogen, y = max_y + 0.01, label = sig_label),
              inherit.aes = FALSE, vjust = 0, size = 5)
  # p2
  return(list(p1,p2))
}


plist1 = myplot(d = 'D1')
plist2 = myplot(d = 'D4')
plist3 = myplot(d = 'D7')
p1 = plist1[[1]]; p2 = plist1[[2]];
p3 = plist2[[1]]; p4 = plist2[[2]];
p5 = plist3[[1]]; p6 = plist3[[2]];
pdf(paste0("Outputs/03_microbe_mortality_association mass.pdf"), width = 2.6, height = 2.4)
print(p1);print(p2);print(p3);print(p4);print(p5);print(p6)
dev.off()
