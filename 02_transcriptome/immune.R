rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(magrittr)

### 查看IL-6
mysheets <- sisiUtils::read_excel_allsheets("../../AI_MASS/MASS_metadata_database_20250607/免疫学数据汇总 2025.6.14.xlsx")
df1 = mysheets[[1]]

colnames(df1) = df1[2,]
df1 = df1[-c(1,2),]

hospital = read.csv('../../AI_MASS/MASS_metadata_database_20250607/1.Hospital_info.csv')
colnames(df1)[1:2] = c('id','Name')
df1 = left_join(df1, hospital, by = 'Name')
df1$id_padded <- sprintf("%03d", as.numeric(df1$id))
df1$HumanID = paste0(df1$Code,'_',df1$id_padded)


load('Inputs/1211_metadata.rdata')
table(meta$HumanID %in% df1$HumanID)

df2 = meta[,c(1,2)] %>% 
  left_join(df1[,-c(1,2)], by = 'HumanID')

df2 = df2[1:20]
df2[3:20] <- lapply(df2[3:20], function(x) as.numeric(as.character(x)))

colnames(df2)[3:20] = c('C4','C3','IgM','IgG','IgA','IL2','IL4','IL6',
                        'IL8','IL1b','IL10','TNFa','IFNg','CD3T','CD3CD4T','CD3CD8T',
                        'CD19B','CD3nCD16CD56NK')

df2$CD4CD8ratio = df2$CD3CD4T/df2$CD3CD8T


experiment_df = read.csv('../../AI_MASS/MASS_metadata_database_20250607/6.Experiment_info.csv')
df2 %<>% left_join(experiment_df %>% dplyr::select(HumanID, hsCRP, LymphocyteCount), by = 'HumanID')

df2$CD4abs = df2$LymphocyteCount * df2$CD3CD4T/100
df2$CD3Tabs = df2$LymphocyteCount * df2$CD3T/100
colnames(df2)

if(F){
  write.csv(df2, 'Outputs/Source_data_immune_assays.csv',
            row.names = F)
}

vars_to_test <- c('C4','C3','IL6','IL8','TNFa','CD3Tabs','CD4abs','hsCRP')


res_baseline <- map_dfr(vars_to_test, function(v){
  
  x_surv <- df2 %>% filter(Mortality28d == 0) %>% pull(!!sym(v))
  x_dead <- df2 %>% filter(Mortality28d == 1) %>% pull(!!sym(v))
  
  # 去掉 NA
  x_surv <- x_surv[!is.na(x_surv)]
  x_dead <- x_dead[!is.na(x_dead)]
  
  # 如果某一组样本太少，直接跳过
  if (length(x_surv) < 3 | length(x_dead) < 3) {
    return(tibble(
      variable = v,
      median_survivor = NA,
      median_death = NA,
      fold_change = NA,
      log2FC = NA,
      p_value = NA
    ))
  }
  
  med_surv <- median(x_surv)
  med_dead <- median(x_dead)
  
  fc <- med_dead / med_surv
  log2fc <- log2(fc)
  
  pval <- wilcox.test(x_dead, x_surv)$p.value
  
  tibble(
    variable = v,
    median_survivor = med_surv,
    median_death = med_dead,
    fold_change = fc,
    log2FC = log2fc,
    p_value = pval
  )
})

res_baseline <- res_baseline %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value)


library(ggplot2)

res_baseline$variable[3] = 'CD3CD4T'
res_baseline$variable[6] = 'CD3T'
dfh <- res_baseline %>%
  mutate(
    variable = factor(variable, levels = variable[order(log2FC)]),
    sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )


breaks_fc <- c(-2, -1,  0,
               1, 2)

fc_colors <- c(
  "#313695", "#4575B4",  
  "white",
 "#F46D43", "#D73027"
)


p = ggplot(dfh, aes(x = "Baseline", y = variable, fill = log2FC)) +
  geom_tile(color = "black", linewidth = 0.8) +
  geom_text(aes(label = sig)) +
  scale_fill_gradientn(
    colors = fc_colors,
    breaks = breaks_fc,
    limits = range(breaks_fc),
    name = "log2FC"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold")
  )

pdf("Outputs/02_immune_markers.pdf", width = 2.6, height = 4)
print(p)
dev.off()

if(F){
  write.csv(dfh, 'Outputs/Supplementary_data_immune_markers.csv',
            row.names = F)
}

