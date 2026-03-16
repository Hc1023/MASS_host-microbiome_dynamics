rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(aplot)

load('Inputs/1211_microbe.rdata')
load('Inputs/1211_metadata.rdata')

df_long = df_long %>% 
  left_join(meta %>% select(HumanID, SurvivalTimeWithin28Days), 
            by = 'HumanID')

tmp = data.frame(x = rev(c(4,7,10,14,29)), 
                 y = rev(c('1-3','4-6','7-9','10-13','14-28')),
                 y_num = 1:5)
df_long$SurvTime = 'Survival'
df_long$SurvTime_num = 0
for (i in 1:nrow(tmp)) {
  idx = (df_long$SurvivalTimeWithin28Days < tmp[i,1] & df_long$Mortality28d == '1')
  df_long$SurvTime[idx] =  tmp[i,2]
  df_long$SurvTime_num[idx] = tmp[i,3]
}

## 确定HumanID order
{
  df_summary <- df_long %>%
    group_by(HumanID) %>%
    summarise(
      Mortality28d = first(Mortality28d),   # 每个患者的死亡状态相同
      SurvivalTime = first(SurvTime_num),
      n_samples = n()
    ) %>%
    arrange(Mortality28d, SurvivalTime, desc(n_samples)) 
  
  levels_ = df_summary$HumanID
  
  # 1个D14样本是D13
  df_long[df_long$Timepoint == 'D14' & df_long$SurvTime =='10-13',]
  # 5个D4样本是D3
  df_long[df_long$Timepoint == 'D4' & df_long$SurvTime =='1-3',]
  
  which(levels_ %in% df_long[df_long$Timepoint == 'D4' & df_long$SurvTime =='1-3','HumanID'])
  # [1] 375 376 377 415
  tmp = levels_[378]
  levels_[378] = levels_[415]
  levels_[415] = tmp
  
  df_long <- arrange(df_long, match(HumanID, levels_))
  df_long$HumanID = factor(df_long$HumanID, levels = levels_)
  
}



df_long$SurvTime = factor(df_long$SurvTime, levels = rev(c('1-3','4-6','7-9','10-13','14-28','Survival')))
color_SurvTime = c("#4575B4", "#a5c4dc", "#D9EF8B","#FEE090", "#FDAE61", "#c64032")
names(color_SurvTime) = rev(c('1-3','4-6','7-9','10-13','14-28','Survival'))

df_long$Timepoint %<>% factor(., levels = rev(levels(.)))
p = ggplot(df_long, aes(y = Timepoint, x = HumanID, fill = SurvTime)) +
  geom_tile(color = NA, linewidth = 0) +  # color 可以画每格边框
  scale_fill_manual(values = color_SurvTime) +
  labs(y = "Timepoint", x = "HumanID", fill = "28-day Mortality") +
  theme_bw() +
  geom_hline(yintercept = c(1.5,2.5,3.5), 
             color = "grey80", linewidth = 0.2) +
  geom_vline(xintercept = table(df_summary$Mortality28d)[1] + 0.5, 
             color = "black", linewidth = 0.2) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()  # 去掉背景网格
  ) +
  scale_y_discrete(expand = c(0,0)) 

p


df_top <- df_long %>%
  group_by(HumanID) %>%
  summarise(Mortality28d = first(Mortality28d)) %>%
  mutate(Timepoint = "Mortality")
p_anno <- ggplot(df_top, aes(y = 1, x = HumanID, fill = factor(Mortality28d))) +
  geom_tile(color = NA) +
  scale_fill_manual(values = c("0" = "#4575B4", "1" = "#D73027")) + # 可自定义颜色
  theme_void() +
  theme(
    legend.position = "none"
  ) 



p %>% insert_top(p_anno, height = 0.1)

df_right = df_long %>%
  count(Timepoint, Mortality28d) 
p_right = ggplot(df_right, aes(y = Timepoint, x = n, fill = Mortality28d)) +
  geom_col(position = "dodge", alpha = 0.9, 
           color = 'black', linewidth = 0.4) +
  geom_text(aes(label = n), 
            position = position_dodge(width = 1),  # 与柱子对齐
            hjust = -0.1,                            # 数值在柱子右侧
            size = 3) +
  scale_fill_manual(values = c("0" = "#4575B4", "1" = "#D73027")) +
  labs(x = 'Samples', y = NULL, fill = "Mortality28d") +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text.x.top = element_text(),        # secondary axis 上方文字显示
    axis.ticks.x.top = element_line(),          # 隐藏主轴刻度
    axis.text.x = element_blank(),           # 隐藏主轴文字
    axis.ticks.x = element_blank(),          # 隐藏主轴刻度
    axis.title.y = element_blank(),      # 去掉 y 轴标题
    axis.text.y = element_blank(),       # 去掉 y 轴文字
    axis.ticks.y = element_blank(),      # 去掉 y 轴刻度
    # panel.grid.major.y = element_blank(),# 去掉 y 轴网格线（建议）
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(
    expand = expansion(add = c(0, 80)),
    sec.axis = dup_axis(name = NULL)  # 复制主轴作为 secondary axis
  ) 

p_right
p_spacer <- ggplot() + theme_void()
p_final = p %>% 
  insert_top(p_anno, height = 0.1) %>% 
  insert_right(p_spacer, width = 0.02) %>%
  insert_right(p_right, width = 0.35)

pdf(file = 'Outputs/01_overview.pdf', width = 6, height = 2.4)
print(p_final)
dev.off()
