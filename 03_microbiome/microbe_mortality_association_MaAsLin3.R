rm(list = ls())
library(tidyverse)
library(data.table)
library(maaslin3)
library(ggplot2)
library(magrittr)

load('Inputs/1616_microbe.rdata')
load('Inputs/1211_metadata.rdata')

# 物种丰度表
taxa_table = data.frame(t(data[pathogens,]))
## 进行log2(m+1)处理
taxa_table = log2(taxa_table+1)
colnames(taxa_table) = pathogens

# 元数据表
# 选择协变量Center
d = 1

df_long1 = df_long %>% filter(Timepoint == paste0('D',d)) %>% 
  left_join(meta[,-2], by = 'HumanID')
df_long1$PneumoniaType = ifelse(df_long1$PneumoniaType == '1.CAP',
                                'CAP','Others') %>% as.factor()
df_long1 %<>%
  mutate(
    Center = ifelse(Center %in% names(which(table(Center) < 10)), "Others", Center)
  )

# cov = colnames(df_long1)[c(6,10:21)]
# # Center作为随机/分层变量
# # 其他13个变量作为固定变量
vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
         'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
var_str = paste(vars, collapse = '+')
df_long1$Gender %<>% as.factor()
rownames(df_long1) = df_long1$SampleID

taxa1 = taxa_table[df_long1$SampleID,]
fit_out <- maaslin3(
  input_data = taxa1,
  input_metadata = df_long1,
  output = 'Outputs/03_mas_output',
  formula = paste('~ Mortality28d +', var_str),
  normalization = 'NONE',       # 总丰度标准化
  transform = 'NONE',           # 对数转换
  augment = TRUE,              # 启用数据增强（解决逻辑回归线性可分问题）
  standardize = F,             # 连续变量Z-score标准化
  max_significance = 0.1,      # FDR阈值设为0.1
  median_comparison_abundance = F,  # 丰度系数与中位数比较（应对组成性效应）
  median_comparison_prevalence = FALSE,# 存在率系数与0比较
  cores = 1,                    # 单核运行
  warn_prevalence = F
)


all_mas = read.delim2('Outputs/03_mas_output/all_results.tsv')
all_mas %<>% filter(metadata == 'Mortality28d')
all_mas$coef %<>% as.numeric()
all_mas$null_hypothesis %<>% as.numeric()
all_mas$xmin = all_mas$coef - as.numeric(all_mas$stderr)
all_mas$xmax = all_mas$coef + as.numeric(all_mas$stderr)
all_mas$pval_individual %<>% as.numeric()
head(all_mas)

median_df = data.frame(
  model = c('Prevalence','Abundance'),
  null = c(all_mas$null_hypothesis[all_mas$model == 'prevalence'][1],
           all_mas$null_hypothesis[all_mas$model == 'abundance'][1])
)

all_mas$feature %<>% factor(., levels = rev(pathogens))

mytheme = theme_bw() +
  theme(axis.text = element_text(color = 'black', size = 13),
        axis.ticks.length = unit(0.2, "lines"),
        axis.ticks = element_line(size = 0.5, color = 'black'),
        title = element_text(color = 'black', size = 13),
        panel.border = element_rect(linewidth = 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = 'right',
        plot.margin = margin(l = 0.1, r = 0.1, t = 1, b = 0.1, unit = "cm"),
        legend.spacing.y = unit(2, "pt")
  )

levels(all_mas$feature)[levels(all_mas$feature) == 'HHV-4'] = 'EBV'
pp = ggplot(all_mas, aes(x = coef, y = feature)) +
  ggplot2::guides(linetype = ggplot2::guide_legend(title = "Null hypothesis", 
                                                   order = 1), ) + 
  ggplot2::geom_vline(data = median_df,ggplot2::aes(xintercept = null, linetype = model), color = "darkgray", size = 0.3) + 
  ggplot2::scale_linetype_manual(values = c(Prevalence = "dashed", Abundance = "solid")) +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0.2, linewidth = 0.4) +
  geom_point(aes(shape =model), size = 3.2) +
  ggplot2::scale_shape_manual(name = "Association", 
                              values = c(21, 24)) + ggplot2::guides(shape = ggplot2::guide_legend(order = 2), 
                              ) + ggplot2::labs(x = expression(paste(beta, " coefficient")), 
                                                y = "") +
  
  geom_point(all_mas[all_mas$model == 'abundance', ], 
             mapping = aes(shape = model, fill = pval_individual), size = 3.2,
             shape = 21, stroke = 0.3) +
  ggplot2::scale_fill_gradient(low = "#8B008B", high = "white",
                               limits = c(1e-5, 1), breaks = c(1e-5, 0.05, 1), 
                               labels =  c(1e-5, 0.05, 1), transform = scales::pseudo_log_trans(sigma = 0.001), 
                               name = "Abundance P") +
  ggnewscale::new_scale_fill() +
  geom_point(all_mas[all_mas$model == 'prevalence', ], 
             mapping = aes(shape = model, fill = pval_individual), size = 3.2,
             shape = 24, stroke = 0.3) +
  ggplot2::scale_fill_gradient(low = "#008B8B", high = "white",
                               limits = c(1e-5, 1), breaks = c(1e-5, 0.05, 1), 
                               labels =  c(1e-5, 0.05, 1), transform = scales::pseudo_log_trans(sigma = 0.001), 
                               name = "Prevalence P") +
  # facet_grid(Type ~ ., space = 'free', scales = 'free') +
  mytheme + theme(strip.text = element_text(size = 13), panel.grid = element_line(linewidth = 0.3)) 
pp

pdf(file = 'Outputs/03_maaslin.pdf', width = 6, height = 5.5)
print(pp)
dev.off()

tmp = all_mas %>% filter(pval_individual < 0.05)

