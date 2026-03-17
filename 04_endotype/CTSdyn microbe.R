rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(purrr)
library(ggplot2)
library(stringr)
library(forcats)

load('Inputs/1315_meta_model.rdata')

df1 = meta_model %>% dplyr::select(HumanID, Mortality28d, D1_CTS, D4_CTS, D7_CTS)
df1 %<>% na.omit()

# 3 as good group; 1/2 as bad group
# deterioration: good -> bad -> bad; good -> good -> bad
# improving: bad -> good -> good; bad -> bad -> good
# stable good: good -> good -> good
# stable bad: bad -> bad -> bad
df1$grp = "Mixed"
idx = df1$D1_CTS == 3 & df1$D7_CTS != 3
df1$grp[idx] = 'Deteriorating'
idx = df1$D1_CTS != 3 & df1$D7_CTS == 3
df1$grp[idx] = 'Improving'
idx = df1$D1_CTS == 3 & df1$D4_CTS == 3 & df1$D7_CTS == 3
df1$grp[idx] = 'Stable low-risk'
idx = df1$D1_CTS != 3 & df1$D4_CTS != 3 & df1$D7_CTS != 3
df1$grp[idx] = 'Stable high-risk'



df_sum <- df1 %>%
  filter(!is.na(grp), !is.na(Mortality28d)) %>%
  group_by(grp) %>%
  summarise(
    n = n(),
    death_rate = mean(Mortality28d == 1),
    .groups = "drop"
  )

df_sum <- df_sum %>%
  mutate(
    death_pct = death_rate * 100,
    n_label = paste0("n = ", n)
  )
unique(df_sum$grp)
df_sum$grp = factor(df_sum$grp, 
                    levels = c("Stable high-risk", "Stable low-risk",
                               "Deteriorating", "Improving", "Mixed"))
p = ggplot(df_sum, aes(x = grp, y = death_pct)) +
  geom_col(
    fill = "#D55E00",
    color = 'black',
    width = 0.6,
    alpha = 0.85
  ) +
  
  geom_text(
    aes(label = n_label),
    # vjust = -0.5,
    hjust = -0.1,
    size = 3,
    color = "black"
  ) +
  
  scale_y_continuous(
    name = "Mortality (%)",
    limits = c(0, max(df_sum$death_pct) * 1.3),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 11),
    # plot.margin = margin(5.5, 5.5, 10, 5.5),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL) +
  coord_flip()


pdf(file = 'Outputs/04_dyngrp_samplen.pdf', width = 3, height = 1.6)
print(p)
dev.off()

#### microbe changes ####

tmp = meta_model %>% dplyr::select(HumanID, matches('D[147]_mi'))
df1_mi = df1 %>% left_join(tmp, by = 'HumanID')

names(df1_mi)
delta_long <- df1_mi %>%
  pivot_longer(
    cols = matches("^D(1|4|7)_(mi2?|mi)_"),
    names_to = c("Day", "feature"),
    names_pattern = "^D(1|4|7)_(.+)$",
    values_to = "value"
  ) %>%
  mutate(Day = paste0("D", Day)) %>%
  pivot_wider(names_from = Day, values_from = value) %>%
  mutate(delta = D7 - D1)

delta_long$feature <- gsub("^(mi2|mi)_", "", delta_long$feature)
delta_long %<>% filter(grp != 'Mixed')

head(data.frame(delta_long))


library(rstatix)

res_delta <- delta_long %>%
  filter(!is.na(D1), !is.na(D7)) %>%
  group_by(grp, feature) %>%
  summarise(
    n = n(),
    median_delta = median(delta, na.rm = TRUE),
    mean_delta   = mean(delta, na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(D7, D1, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

mean_delta_mat <- res_delta %>%
  dplyr::select(feature, grp, mean_delta) %>%
  tidyr::pivot_wider(
    names_from  = grp,
    values_from = mean_delta
  )

pval_mat <- res_delta %>%
  dplyr::select(feature, grp, p_value) %>%
  tidyr::pivot_wider(
    names_from  = grp,
    values_from = p_value
  )

pval_mat_2 = pval_mat %>% filter(!is.na(Deteriorating)) %>% filter(Deteriorating<1)
pval_mat_2 %<>% arrange(Deteriorating)


bubble_df <- pval_mat_2 %>%
  pivot_longer(
    cols = -feature,
    names_to = "grp",
    values_to = "p"
  ) %>%
  left_join(
    mean_delta_mat %>%
      pivot_longer(
        cols = -feature,
        names_to = "grp",
        values_to = "delta"
      ),
    by = c("feature", "grp")
  ) %>%
  # filter(!is.na(p)) %>%
  mutate(
    neg_log10_p = -log10(p)
  )

bubble_df$feature = factor(bubble_df$feature, levels = rev(pval_mat_2$feature))

bubble_df <- bubble_df %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1  ~ "\u00B7",
      TRUE ~ ""
    )
  )

levels(bubble_df$feature)[levels(bubble_df$feature) == "HHV.4"] <- "EBV"
levels(bubble_df$feature)[levels(bubble_df$feature) == "Influenza.A"] <- "Influenza A"
levels(bubble_df$feature)[levels(bubble_df$feature) == "Bacf"] <- "Bacteria/Fungi"

p = ggplot(bubble_df, aes(x = grp, y = feature)) +
  geom_tile(
    aes(fill = delta),
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", delta), sig)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name =  expression(Delta~"(D7-D1)")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c('CTS3 -> 1/2',
                              'CTS1/2 -> 3',
                              'Stable CTS1/2',
                              'Stable CTS3'))

pdf(file = 'Outputs/04_dyngrp_mi.pdf', width = 5, height = 3)
print(p)
dev.off()


#### 折线图  ####

feature_x = 'Klebsiella'
feature_x = 'Total'

delta_long$feature[delta_long$feature == 'HHV.4'] = 'EBV'
delta_long$feature = gsub('\\.',' ',delta_long$feature)
delta_long$feature[delta_long$feature == 'SARS CoV 2'] = 'SARS-CoV-2'
delta_long$feature[delta_long$feature == 'HAdV B'] = 'HAdV-B'
delta_long$feature[delta_long$feature == 'HHV 6B'] = 'HHV-6B'
delta_long$feature[delta_long$feature == 'Bacf'] = 'Bacteria/Fungi'
unique(delta_long$feature)
plist = list()
for (feature_x in unique(delta_long$feature)) {
  dfx = delta_long %>% filter(feature == feature_x)
  head(data.frame(dfx))
  
  df_long <- dfx %>%
    pivot_longer(
      cols = c(D1, D4, D7),
      names_to = "Day",
      values_to = "value"
    ) %>%
    mutate(
      Day = factor(Day, levels = c("D1", "D4", "D7"))
    )
  
  stat_day <- df_long %>%
    group_by(Day) %>%
    kruskal_test(value ~ grp) %>%
    ungroup() %>%
    mutate(
      label = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        p < 0.1   ~ "\u00B7",
        TRUE ~ ""
      )
    )
  y_pos <- df_long %>%
    group_by(Day, grp) %>%
    summarise(
      m  = mean(value, na.rm = TRUE),
      se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    ) %>%
    mutate(upper = m + se) %>%
    group_by(Day) %>%
    summarise(
      y.position = max(upper, na.rm = TRUE) * 1.05,
      .groups = "drop"
    )
  
  stat_day <- stat_day %>%
    left_join(y_pos, by = "Day")
  
  df_long <- df_long %>%
    mutate(
      line_type = ifelse(grp %in% c("Deteriorating", "Improving"),
                         "focus", "other")
    )
  
  p = ggplot(df_long,
             aes(x = Day, y = value,
                 color = grp, group = grp)) +
    
    ## mean ± SE
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.15,
      linewidth = 0.6
    ) +
    
    ## 组均值线（加 linetype）
    stat_summary(
      aes(linetype = line_type),
      fun = mean,
      geom = "line",
      linewidth = 0.6
    ) +
    
    ## 组均值点
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2
    ) +
    
    scale_linetype_manual(
      values = c(
        focus = "solid",
        other = "dashed"
      ),
      guide = "none"   # 不在图例里显示 linetype（推荐）
    ) +
    
    theme_bw() +
    labs(
      x = NULL,
      y = "log2(mass+1)",
      color = "Trajectory group"
    ) +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'none') +
    ggtitle(feature_x) +
    geom_text(
      data = stat_day,
      aes(x = Day, y = y.position, label = label),
      inherit.aes = FALSE,
      size = 5
    ) +
    scale_color_manual(
      values = c(
        "Deteriorating"      = "#D73027",  # red
        "Stable high-risk"   = "#FC8D59",  # orange
        "Improving"          = "#1A9850",  # green
        "Stable low-risk"    = "#4575B4",  # blue
        "Mixed"              = "grey60"    # optional
      )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  plist[[feature_x]] = p
}


p

library(gridExtra)

pdf(paste0("Outputs/04_CTSdyn microbe.pdf"), width = 8, height = 10)
grid.arrange(grobs = plist, ncol = 4)
dev.off()

p1 = plist[['Total']]
p2 = plist[['Klebsiella']]

p_leg = p1 + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(legend.position = 'top') +scale_color_manual(
    values = c(
      "Deteriorating"    = "#D73027",
      "Stable high-risk" = "#FC8D59",
      "Improving"        = "#1A9850",
      "Stable low-risk"  = "#4575B4"
    ),
    labels = c(
      "Deteriorating"    = 'CTS3->1/2',
      "Stable high-risk" = 'Stable CTS1/2',
      "Improving"        = 'CTS1/2->3',
      "Stable low-risk"  = 'Stable CTS3'
    )
  )

legend_only <- get_legend(p_leg)
grid::grid.newpage()
grid::grid.draw(legend_only)


p1 = p1 + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_x_discrete(expand = expansion(mult = c(0.1,0.1)))
p2 = p2 + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_x_discrete(expand = expansion(mult = c(0.1,0.1)))

pdf(paste0("Outputs/04_CTSdyn microbe_2.pdf"), width = 2, height = 2)
print(p1)
print(p2)
dev.off()

pdf(paste0("Outputs/04_CTSdyn microbe_2leg.pdf"), width = 7, height = 3)
grid::grid.draw(legend_only)
dev.off()

