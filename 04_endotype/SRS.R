rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(edgeR)
library("SepstratifieR")

load('Inputs/1211_metadata.rdata')
load('Inputs/1211_transcriptome.rdata')
load('Inputs/1616_microbe.rdata')

dge <- DGEList(counts)
dge <- calcNormFactors(dge)
cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE, log = T)
srs_preds <- stratifyPatients(t(cpm_norm), gene_set = "extended")
identical(df_long$SampleID, colnames(cpm_norm))
tmp = srs_preds@SRS_probs
max_col_index <- max.col(tmp[,c(1,2)], ties.method = "first")
df_long$SRS = max_col_index
table(df_long$Mortality28d, df_long$SRS)

# D1/4/7 have samples or death
meta_analysis = meta %>% filter(!is.na(D1) & !is.na(D4.death) & !is.na(D7.death))
table(meta_analysis$Mortality28d)

meta_analysis %<>%
  dplyr::select(
    HumanID,
    Mortality28d,
    D1,
    D4 = D4.death,
    D7 = D7.death
  ) %>%
  pivot_longer(
    cols = c(D1, D4, D7),
    names_to = "Timepoint",
    values_to = "SampleID"
  )
tmp = df_long %>% dplyr::select(SampleID, SRS)

meta_analysis %<>% left_join(tmp, by = 'SampleID')

meta_analysis$SRS = factor(meta_analysis$SRS, 
                           levels = c("Death","Alive","1", "2"))
meta_analysis$SRS[meta_analysis$SampleID == 'Death'] = 'Death'

d28_rows <- meta_analysis %>%
  distinct(HumanID, Mortality28d) %>%   # 保留每个人的结局信息
  mutate(
    Timepoint = "D28",
    SampleID  = NA_character_,
    SRS       = ifelse(Mortality28d == 1, 'Death', 'Alive')
  )

meta_analysis %<>%
  bind_rows(d28_rows) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D1","D4","D7","D28"))
  ) %>%
  arrange(HumanID, Timepoint)


df_plot = meta_analysis %>% dplyr::select(HumanID, Timepoint, State = SRS)
df_plot$State[df_plot$State == 'Death'] = 'Mortality'
df_plot$State[df_plot$State == 'Alive'] = 'Survival'
df_plot$State %<>% factor(., levels = c("Mortality","Survival","1", "2"))


library(ggalluvial)

p = ggplot(df_plot,
       aes(x = Timepoint,
           stratum = State,
           alluvium = HumanID,
           fill = State,
           label = State)) +
  
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            alpha = .8, width = 1/12) +
  
  geom_stratum(width = 1/12, color = "black") +
  geom_text(stat = "stratum") +
  
  scale_x_discrete(expand = c(.1, .1)) +
  
  labs(title = "SRS Trajectory Alluvial Plot",
       x = "Timepoint",
       y = "Number of Patients",
       fill = "State") +
  
  theme_classic() +
  scale_fill_manual(values = c(
    "Mortality" = "#d73027",   
    "1"     = "#FDC086",
    "2"     = "#BEAED4" ,   
    "Survival" = "#4575b4"    
  ),
  label = c('1' = 'SRS1','2' = 'SRS2')) 

if(F){
  pdf(paste0("Outputs/04_SRS.pdf"), width = 5, height = 3)
  print(p)
  dev.off()
}


#### State mortality at each timepoint ####

df_mort_prop <- meta_analysis %>%
  filter(Timepoint %in% c("D1", "D4", "D7"),
         SRS %in% c("1", "2")) %>%
  mutate(
    SRS = factor(SRS, levels = c("1", "2")),
    Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7"))
  )

mort_summary <- df_mort_prop %>%
  group_by(Timepoint, SRS) %>%
  summarise(
    n_total = n(),
    n_dead  = sum(Mortality28d == 1),
    prop_dead = n_dead / n_total,
    .groups = "drop"
  )

mort_summary

# 转成 matrix：行 = Timepoint，列 = CTS
mat_prob <- mort_summary %>%
  dplyr::select(Timepoint, SRS, prop_dead) %>%
  pivot_wider(
    names_from  = SRS,
    values_from = prop_dead
  ) %>%
  column_to_rownames("Timepoint") %>%
  as.matrix()

library(ComplexHeatmap)
library(circlize)
library(grid)
library(scales)
library(RColorBrewer)

col_fun <- colorRamp2(
  c(0, 0.3, 0.6),
  c("#F7FBFF", "#FDAE61", "#B2182B")
)

p = Heatmap(
  mat_prob,
  name = "P(death)",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title = "Timepoint",
  column_title = "SRS",
  
  row_names_side = "left",
  column_names_side = "top",
  
  rect_gp = gpar(col = "white", lwd = 1),
  
  # 在格子里写百分比
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(
      percent(mat_prob[i, j], accuracy = 0.1),
      x = x,
      y = y,
      gp = gpar(
        fontsize = 10,
        col = "black"
      )
    )
  },
  
  heatmap_legend_param = list(
    at = c(0, 0.3, 0.6),
    labels = percent(c(0, 0.3, 0.6)),
    title = "Mortality\nprobability"
  )
)

pdf(paste0("Outputs/04_SRS_mort.pdf"), width = 2.6, height = 2.2)
print(p)
dev.off()
