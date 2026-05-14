rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)
library(vegan)
library(aplot)

load('Inputs/1211_metadata.rdata')
load('Inputs/1616_microbe.rdata')

#### beta diversity ####

plot_beta = function(data2, d, title_){
  df_long1 = df_long[df_long$Timepoint == d,]
  data2 = data2[,df_long1$SampleID]
  df_long2 = df_long1[colSums(data2)>0,]
  data2 = data2[,df_long2$SampleID]
  
  identical(colnames(data2), df_long2$SampleID)
  table(df_long2$Mortality28d)
  
  method = "bray"; #method = "jaccard"
  dist <- vegdist(t(data2), method = method)
  pcoa_res <- cmdscale(dist, eig = TRUE, k = 2)
  var_exp <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)
  var_exp[1:2]  # PC1 和 PC2
  pcoa_df <- data.frame(
    PC1 = pcoa_res$points[,1],
    PC2 = pcoa_res$points[,2],
    Mortality28d = factor(df_long2$Mortality28d)
  )
  set.seed(123) 
  adonis_res <- adonis2(dist ~ Mortality28d, data = df_long2, permutations = 1000)

  p_val <- adonis_res$`Pr(>F)`[1]  # 取模型行的 p 值
  p_text <- paste0("PERMANOVA p = ", signif(p_val, 2))

  
  p = ggplot(pcoa_df, aes(PC1, PC2, color = factor(Mortality28d))) +
    geom_point(size = 2, alpha = 0.9) +
    stat_ellipse(level = 0.95, linetype = 2) +
    scale_color_manual(values = c("0" = "#68a0db", "1" = "#c43932"),
                       labels = c("0" = "Survival", "1" = "Mortality")) +
    xlab(paste0("PC1 (", var_exp[1], "%)")) +
    ylab(paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw() +
    labs(title = title_) +
    theme(legend.title = element_blank(),
          legend.position = 'none',
          plot.margin = unit(c(0.1, 0.15, 0.1, 0.1), "cm"),
          panel.grid.minor = element_blank()) +
    annotate("text", x = min(pcoa_df$PC1)*1.5, y = max(pcoa_df$PC2)*1.2,
             label = p_text, hjust = 0, vjust = 1, size = 4)
  return(p)
}

data_log2 = log2(data+1)


plist = list()
for (d in c('D1','D4','D7')) {
  print(d)
  p1 = plot_beta(data2 = data_log2, 
                 d = d, title_ = 'Total')
  p2 = plot_beta(data2 = data_log2[mapid_vec[rownames(data_log2)] != 'Viruses',], 
                 d = d, title_ = 'Bacteria/Fungi')
  p3 = plot_beta(data2 = data_log2[mapid_vec[rownames(data_log2)] == 'Viruses',], 
                 d = d, title_ = 'Viruses')
  plist[[d]] = ggarrange(p1,p2,p3, ncol = 3, nrow = 1) 
}

pdf(paste0("Outputs/03_biomass_beta.pdf"), width = 7.2, height = 2.2)
print(plist[[1]])
print(plist[[2]])
print(plist[[3]])
dev.off()

