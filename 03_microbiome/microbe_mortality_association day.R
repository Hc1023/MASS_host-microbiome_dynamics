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
}

d = 'D1';virus=0
myplot = function(d){
  # 1. select samples and pathogens
  df_long1 = df_long[df_long$Timepoint == d,]
  genus_sum3 = genus_sum[, df_long1$SampleID]
  
  results <- data.frame(
    Pathogen = pathogens,
    Percent_mortality = NA,
    Percent_survival = NA,
    P_value = NA
  )
  
  # 遍历每个属
  for (i in seq_along(pathogens)) {
    bug <- pathogens[i]
    
    # --- 死亡组 ---
    dead_samples <- df_long1 %>% filter(Mortality28d == 1) %>% pull(SampleID)
    dead_data <- genus_sum2[bug, dead_samples, drop = FALSE] %>% unlist()
    n_dead <- length(dead_samples)
    n_dead_with <- sum(dead_data > 0)
    
    # --- 存活组 ---
    surv_samples <- df_long1 %>% filter(Mortality28d == 0) %>% pull(SampleID)
    surv_data <- genus_sum2[bug, surv_samples, drop = FALSE] %>% unlist()
    n_surv <- length(surv_samples)
    n_surv_with <- sum(surv_data > 0)
    
    # --- 计算百分比 ---
    results$Percent_mortality[i] <- n_dead_with / n_dead
    results$Percent_survival[i] <- n_surv_with / n_surv
    
    # --- 构建 2×2 表并做 Fisher’s exact test ---
    contingency <- matrix(
      c(n_dead_with, n_dead - n_dead_with,
        n_surv_with, n_surv - n_surv_with),
      nrow = 2,
      byrow = TRUE
    )
    
    ft <- fisher.test(contingency)
    results$P_value[i] <- ft$p.value
    
  }

  # 1. 将数据转换为长格式，方便 ggplot
  df_longplot <- results %>%
    select(Pathogen, Percent_mortality, Percent_survival) %>%
    pivot_longer(
      cols = c(Percent_mortality, Percent_survival),
      names_to = "Outcome",
      values_to = "Percent"
    ) %>%
    mutate(
      Percent = 100*ifelse(Outcome == "Percent_mortality", -Percent, Percent),
      Outcome = ifelse(Outcome == "Percent_mortality", "Mortality", "Survival")
    )
  df_longplot$Pathogen <- factor(df_longplot$Pathogen, levels = pathogens)
  
  
  # 添加标记列
  results <- results %>%
    mutate(sig_label = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      P_value < 0.1 ~ ".",
      TRUE ~ ""
    ))
  
  # 找每个 Pathogen 的最大条形高度，用于放置标记
  max_percent <- df_longplot %>%
    group_by(Pathogen) %>%
    summarise(max_y = max(Percent))
  
  # 合并显著性标记
  annot_df <- max_percent %>%
    left_join(results %>% select(Pathogen, sig_label), by = "Pathogen")
  annot_df$max_y = ifelse(annot_df$sig_label %in% c("*","**","***"), 
                          annot_df$max_y * 0.9,
                          annot_df$max_y * 1.1)
  
  levels(df_longplot$Pathogen)[levels(df_longplot$Pathogen) == 'HHV-4'] = 'EBV'
  annot_df$Pathogen[annot_df$Pathogen == 'HHV-4'] = 'EBV'
  # 2. 绘图
  # mortality_col_fun = c("0" = "#4575B4", "1" = "#D73027")
  p1 = ggplot(df_longplot[1:22,], aes(x = Pathogen, y = Percent, fill = Outcome)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Survival" = "#4575B4", "Mortality" = "#D73027")) +
    labs(x = "", y = "", fill = "Outcome") +
    theme_bw() + labs(title = d) +
    scale_y_continuous(labels = function(x) abs(x), # y轴显示正数
                       expand = expansion(mult = c(0.05, 0.12))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = 'None',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background  = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA)) +
    geom_text(data = annot_df[1:11,], aes(x = Pathogen, y = max_y + 0.01, label = sig_label),
              inherit.aes = FALSE, vjust = 0, size = 5)
  p2 = ggplot(df_longplot[23:42,], aes(x = Pathogen, y = Percent, fill = Outcome)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Survival" = "#4575B4", "Mortality" = "#D73027")) +
    labs(x = "", y = "", fill = "Outcome") +
    theme_bw() + labs(title = d) +
    scale_y_continuous(labels = function(x) abs(x), # y轴显示正数
                       expand = expansion(mult = c(0.05, 0.12))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = 'None',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background  = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA)) +
    geom_text(data = annot_df[12:21,], aes(x = Pathogen, y = max_y + 0.01, label = sig_label),
              inherit.aes = FALSE, vjust = 0, size = 5)
  return(list(p1,p2))
}


plist1 = myplot(d = 'D1')
plist2 = myplot(d = 'D4')
plist3 = myplot(d = 'D7')
p1 = plist1[[1]]; p2 = plist1[[2]];
p3 = plist2[[1]]; p4 = plist2[[2]];
p5 = plist3[[1]]; p6 = plist3[[2]];
pdf(paste0("Outputs/03_microbe_mortality_association.pdf"), width = 2.6, height = 2.5)
print(p1);print(p2);print(p3);print(p4);print(p5);print(p6)
dev.off()
