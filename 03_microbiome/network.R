rm(list = ls())
library(tidyverse)
library(magrittr)
library(Hmisc)
library(scales)

df_tab = read.csv("Outputs/Source_data_table1.csv")
df_tab %<>% mutate(SampleID = paste0(HumanID, "_D1"))

load("Inputs/1616_microbe.rdata")
load("Inputs/1211_metadata.rdata")

sel_D1 = df_long %>% filter(Timepoint == "D1")

microbe_sel = data[pathogens, ] %>% dplyr::select(sel_D1$SampleID)
df_tab_sel = df_tab %>% filter(SampleID %in% sel_D1$SampleID)

microbe_sel_log2 = log2(microbe_sel + 1)
microbe_sel_log2 %<>% t()
microbe_sel_log2 %<>% data.frame() %>% rownames_to_column(var = "SampleID")

df_combined = df_tab_sel %>% left_join(microbe_sel_log2, by = "SampleID")

vars = colnames(df_combined)
vars = vars[-1]
vars = vars[vars != "SampleID"]
names(vars)[1:24] = rep("Clinical_var", 24)
names(vars)[34:43] = "Viruses"
names(vars)[24:33] = "Bacteria/Fungi"

df_combined_sel = df_combined %>% dplyr::select(all_of(unname(vars)))

## -----------------------------
## 1. 准备 node 信息
## -----------------------------
node_df = tibble(
  id = names(df_combined_sel),
  label = names(df_combined_sel),
  group = names(vars)
)

## 可选：如果你想让 Cytoscape 显示更漂亮，可加一个 type 列
node_df = node_df %>%
  mutate(
    type = case_when(
      group == "Clinical_var" ~ "Clinical",
      group == "Viruses" ~ "Virus",
      group == "Bacteria/Fungi" ~ "Bacteria/Fungi"
    )
  )

## -----------------------------
## 2. 准备相关性矩阵
## -----------------------------
dat_cor = df_combined_sel

## 把字符型变量转成 factor -> numeric
## 注意：这对 Gender / PneumoniaType 能跑，但解释需谨慎
dat_cor = dat_cor %>%
  mutate(across(where(is.character), ~ as.numeric(factor(.x))))

## 如果你不想让 PneumoniaType 进网络，可以打开这一行
dat_cor = dat_cor %>% dplyr::select(-PneumoniaType)

## 确保都是 numeric
dat_cor = dat_cor %>% mutate(across(everything(), as.numeric))

## 计算 Spearman correlation
cor_res = Hmisc::rcorr(as.matrix(dat_cor), type = "spearman")

cor_mat = cor_res$r
p_mat   = cor_res$P

## -----------------------------
## 3. 整理成 edge table
## -----------------------------
var_names = colnames(dat_cor)

edge_df = expand.grid(
  source = var_names,
  target = var_names,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(source != target) %>%
  rowwise() %>%
  mutate(
    rho = cor_mat[source, target],
    pvalue = p_mat[source, target]
  ) %>%
  ungroup()

## 只保留上三角，避免重复边
edge_df = edge_df %>%
  mutate(
    pair = map2_chr(source, target, ~ paste(sort(c(.x, .y)), collapse = "__"))
  ) %>%
  distinct(pair, .keep_all = TRUE) %>%
  dplyr::select(-pair)

## FDR 校正
edge_df = edge_df %>%
  mutate(
    fdr = p.adjust(pvalue, method = "BH"),
    abs_rho = abs(rho),
    sign = ifelse(rho > 0, "positive", "negative")
  )

## 加入 node group 信息
edge_df = edge_df %>%
  left_join(node_df %>% dplyr::select(id, source_group = group), by = c("source" = "id")) %>%
  left_join(node_df %>% dplyr::select(id, target_group = group), by = c("target" = "id")) %>%
  mutate(
    edge_type = paste(source_group, target_group, sep = "_")
  )

## -----------------------------
## 4. 筛选边
## -----------------------------
## pvalue < 0.05
edge_sig = edge_df %>%
  dplyr::filter(!is.na(pvalue), pvalue < 0.05)

## 不保留Clinical_var_Clinical_var组内边，打开这一行
edge_sig %<>% dplyr::filter(edge_type != "Clinical_var_Clinical_var")

## 去掉negative边
edge_sig %<>% filter(sign == 'positive')

## 给 Cytoscape 一个可直接映射的宽度
## 这里把 abs(rho) 映射到 1~8
if (nrow(edge_sig) > 0) {
  edge_sig = edge_sig %>%
    mutate(
      edge_width = scales::rescale(abs_rho, to = c(1, 8))
    )
} else {
  edge_sig = edge_sig %>%
    mutate(edge_width = numeric(0))
}

## -----------------------------
## 5. 只保留网络中实际出现的节点
## -----------------------------
node_use = node_df %>%
  dplyr::filter(id %in% unique(c(edge_sig$source, edge_sig$target)))

## 可加 node size，按 degree
if (nrow(edge_sig) > 0) {
  degree_df = tibble(id = c(edge_sig$source, edge_sig$target)) %>%
    count(id, name = "degree")
  
  node_use = node_use %>%
    left_join(degree_df, by = "id") %>%
    mutate(
      degree = replace_na(degree, 0),
      node_size = scales::rescale(degree, to = c(20, 80))
    )
} else {
  node_use = node_use %>%
    mutate(
      degree = 0,
      node_size = 20
    )
}

## -----------------------------
## 6. 导出 Cytoscape 用 csv
## -----------------------------
write.csv(node_use, "03_microbiome/network/network_nodes.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(edge_sig, "03_microbiome/network/network_edges.csv", 
          row.names = FALSE, quote = FALSE)

## 看一下结果
cat("Number of nodes:", nrow(node_use), "\n")
cat("Number of edges:", nrow(edge_sig), "\n")

## -----------------------------
## 7. Cytoscape legend
## -----------------------------
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(grid)

rho_vals = c("0.1", "0.2", "0.3", "0.4")

df_dummy = data.frame(
  x = 1:4,
  y = 1,
  rho = factor(rho_vals, levels = rho_vals)
)

p_dummy = ggplot(df_dummy, aes(x = x, y = y, color = rho, linewidth = rho)) +
  geom_segment(
    aes(x = x - 0.35, xend = x + 0.35, y = y, yend = y)
  ) +
  scale_color_manual(
    values = brewer.pal(9, "YlOrRd")[c(2, 4, 6, 8)],
    name = "Spearman correlation"
  ) +
  scale_linewidth_manual(
    values = c(1.2, 2.2, 3.2, 4.2),
    guide = "none"
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      # byrow = TRUE,
      override.aes = list(
        linewidth = c(1.2, 2.2, 3.2, 4.2)
      )
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.background = element_rect(fill = "white", colour = NA)
  )

legend_only = cowplot::get_legend(p_dummy)

pdf("03_microbiome/network/spearman_legend_only.pdf", width = 10, height = 1.8)
grid.newpage()
grid.draw(legend_only)
dev.off()

