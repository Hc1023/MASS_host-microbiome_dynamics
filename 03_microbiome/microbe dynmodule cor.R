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

load('Inputs/1211_metadata.rdata')
load('Inputs/1211_transcriptome.rdata')


metadata_analysis = df_long %>% 
  filter(Timepoint %in% c('D1','D4','D7')) %>%
  droplevels()

vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
         'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))
metadata_analysis %<>% left_join(tmp, by = 'HumanID')
meta_gsva = metadata_analysis %>% dplyr::select(HumanID,Timepoint,SampleID,Mortality28d)

#### microbe module cor ####
# load('Inputs/1315_module.rdata')
gsva_mat = read.csv('Outputs/Supplementary_data_module_gsva.csv',
                     row.names = 1)
load('Inputs/1211_microbe.rdata')
data_filtered = data[rowSums(data>0) > 0.05*ncol(data),]
data_filtered = data_filtered[order(rowSums(data_filtered>0), decreasing = T),]
data_filtered_log2 = log2(data_filtered+1)
data_filtered_log2 %<>% dplyr::select(all_of(colnames(gsva_mat)))


identical(meta_gsva$SampleID, colnames(gsva_mat))
dfmm = cbind(meta_gsva, t(gsva_mat), t(data_filtered_log2))

rownames(data_filtered_log2) %<>% make.names()
microbe_vars = rownames(data_filtered_log2)
host_vars = c("up1", "up2", "up3", "dw1")
microbe_vars[microbe_vars == 'HHV.4'] = 'EBV'
colnames(dfmm)[colnames(dfmm) == 'HHV-4'] = 'EBV'
cor_df <- expand.grid(
  host = host_vars,
  microbe = microbe_vars,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    test = list(
      cor.test(
        dfmm[[host]],
        dfmm[[microbe]],
        method = "spearman",
        exact = FALSE
      )
    ),
    rho = test$estimate,
    p.value = test$p.value
  ) %>%
  ungroup() %>%
  dplyr::select(-test) %>%
  mutate(
    adj.p = p.adjust(p.value, method = "BH")
  )

cor_plot <- cor_df %>%
  filter(
    adj.p < 0.1,
  )

cor_df2 <- cor_df %>%
  mutate(
    logFDR = -log10(adj.p),
    sig = case_when(
      adj.p < 0.001 ~ "***",
      adj.p < 0.01  ~ "**",
      adj.p < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

host_label_map = c(up1 = 'Inflammatory signaling', 
                   up2 = 'Phagolysosome function', 
                   up3 = 'IFN signaling', 
                   dw1 = 'Ribosome biogenesis')

cor_df2 %<>%
  mutate(
    host = factor(host, levels = host_vars),
    microbe = factor(microbe, levels = microbe_vars)
  )

library(scales)
library(RColorBrewer)


p_bubble_all <- cor_df2 %>%
  mutate(
    host = factor(host, levels = host_vars),
    microbe = factor(microbe, levels = microbe_vars)
  ) %>%
  ggplot(aes(x = host, y = microbe)) +
  geom_point(
    aes(fill = rho, size = logFDR),
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sig),
    color = "black",
    size = 3
  ) +
  scale_fill_gradientn(
    colours = brewer.pal(10, "PuOr"),
    values = rescale(c(0.3, 0, -0.3)),
    limits = c(-0.3, 0.3),
    name = expression(Spearman~rho)
  ) +
  scale_size(
    range = c(1.2, 10),
    name = expression(-log[10]*"(FDR)")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = host_label_map)  +
  theme(
    plot.margin = margin(t = 30, r = 30, b = 5, l = 5)
  ) 

pdf(paste0("Outputs/03_microbe_mod_cor.pdf"), width = 4, height = 4)
print(p_bubble_all)
dev.off()

#### microbe HCMV module D1/4/7 glm ####

# 评估 HCMV 在不同 Timepoint（D1 / D4 / D7）上的回归效应，
# 同时 控制一组临床协变量。

vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
         'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
tmp = meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days, all_of(vars))
dfmm %<>% left_join(tmp, by = 'HumanID')

# 临床协变量
vars_str <- paste(vars, collapse = " + ")

# 已经 log 的微生物协变量（排除 HCMV）
microbe_vars <- c(
  "EBV", "HCMV", "TTV",
  "Streptococcus", "Lactobacillus", "Staphylococcus",
  "Prevotella", "Veillonella", "HPgV", "Aspergillus",
  "Mycobacterium"
)
microbe_covars <- setdiff(microbe_vars, "HCMV")
microbe_str <- paste(microbe_covars, collapse = " + ")

# Timepoint 设定参考水平（确保 D1 是 reference）
dfmm$Timepoint <- factor(dfmm$Timepoint, levels = c("D1","D4","D7"))

# design：HCMV × Timepoint + 其他微生物 + 临床协变量
library(lme4)
library(lmerTest)
library(broom.mixed)
library(purrr)
outcomes <- c("up1", "up2", "up3", "dw1")

fit_list <- map(
  outcomes,
  ~ lmer(
    as.formula(
      paste(
        .x, "~ HCMV * Timepoint +",
        microbe_str, "+",
        vars_str,
        "+ (1 | HumanID)"
      )
    ),
    data = dfmm,
    REML = FALSE
  )
)

names(fit_list) <- outcomes

library(emmeans)

emm_res <- map2_dfr(
  fit_list,
  names(fit_list),
  ~ {
    em <- emtrends(.x, specs = "Timepoint", var = "HCMV")
    test(em) %>%
      as.data.frame() %>%
      mutate(outcome = .y)
  }
)

emm_res <- emm_res %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  )


pathway_map <- c(
  up1 = "Inflammatory signaling",
  up2 = "Phagolysosome function",
  up3 = "IFN signaling",
  dw1 = "Ribosome biogenesis"
)

plot_df <- emm_res %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7")),
    Pathway   = factor(outcome,
                       levels = c("up1","up2","up3","dw1"),
                       labels = pathway_map[c("up1","up2","up3","dw1")]),
    logFDR = -log10(p_adj)
  )


library(ggplot2)
library(RColorBrewer)
library(scales)

plot_df <- plot_df %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

p_bubble <- ggplot(
  plot_df,
  aes(x = Timepoint, y = Pathway)
) +
  geom_point(
    aes(fill = HCMV.trend, size = logFDR),
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_gradientn(
    colours = rev(brewer.pal(10, "PuOr")),
    values  = rescale(c(-0.08, 0, 0.08)),  # 覆盖你现在的系数范围
    limits  = c(-0.08, 0.08),
    name    = expression("HCMV effect (" * beta * ")")
  ) +
  scale_size(
    range = c(2, 10),
    name  = expression(-log[10] * "(FDR)")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.margin  = margin(5, 30, 5, 5)
  )

p_bubble <- p_bubble +
  geom_text(
    data = subset(plot_df, sig != ""),
    aes(label = sig),
    color = "black",
    size = 3
  )

p_bubble

pdf(paste0("Outputs/03_microbe_mod_cor_CMV.pdf"), width = 4.5, height = 3.2)
print(p_bubble +
        theme(
          plot.margin = margin(t = 30, r = 30, b = 5, l = 5)
        )
)
dev.off()




