rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)

load('Inputs/1211_metadata.rdata')
load('Inputs/1211_microbe.rdata')

# !choose
d = 'D1'
df_long1 = df_long[df_long$Timepoint == d,]


#### ERCC comparison ####
myplot = function(df_long1){
  data1 = data[,df_long1$SampleID]
  
  df_long1$Microbial_mass = colSums(data1)
  df_long1$Microbial_mass_nonvirus = colSums(data1[mapid_vec[rownames(data1)] != 'Viruses',])
  df_long1$Microbial_mass_virus = colSums(data1[mapid_vec[rownames(data1)] == 'Viruses',])
  
  # table(df_long1$Microbial_mass == 0)
  table(df_long1$Microbial_mass_nonvirus == 0)
  
  
  p1 = ggplot(df_long1, aes(x = Mortality28d, y = Microbial_mass, 
                            group = Mortality28d, color = Mortality28d)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Mortality28d), alpha = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    scale_color_manual(values = c("0" = "#68a0db", "1" = "#c43932")) +
    scale_fill_manual(values = c("0" = "#68a0db", "1" = "#c43932")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format",
                       comparisons = list(c("0", "1")),
                       hide.ns = F) +
    coord_cartesian(ylim = 10^3*c(0.5e-5, 35)) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.2)),
                  labels = trans_format("log10", math_format(10^.x)),
                  breaks = 10^3*c(1e1,1e-2,1e-5)) +
    scale_x_discrete(labels = c("0" = "Survival", "1" = "Mortality")) +
    xlab('') + ylab('Microbial mass (fg)') + labs(title = 'Total')
  
  p2 = ggplot(df_long1, aes(x = Mortality28d, y = Microbial_mass_nonvirus, 
                            group = Mortality28d, color = Mortality28d)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Mortality28d), alpha = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    scale_color_manual(values = c("0" = "#68a0db", "1" = "#c43932")) +
    scale_fill_manual(values = c("0" = "#68a0db", "1" = "#c43932")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format",
                       comparisons = list(c("0", "1")),
                       hide.ns = F) +
    coord_cartesian(ylim = 10^3*c(0.5e-5, 1)) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.2)),
                  labels = trans_format("log10", math_format(10^.x)),
                  breaks = 10^3*c(1e0,1e-2,1e-4)) +
    scale_x_discrete(labels = c("0" = "Survival", "1" = "Mortality")) +
    xlab('') + ylab('') + labs(title = 'Bacteria/Fungi')

  p3 = ggplot(df_long1, aes(x = Mortality28d, y = Microbial_mass_virus, 
                            group = Mortality28d, color = Mortality28d)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Mortality28d), alpha = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    scale_color_manual(values = c("0" = "#68a0db", "1" = "#c43932")) +
    scale_fill_manual(values = c("0" = "#68a0db", "1" = "#c43932")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.format",
                       comparisons = list(c("0", "1")),
                       hide.ns = F) +
    coord_cartesian(ylim = 10^3*c(0.5e-5, 35)) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.2)),
                  labels = trans_format("log10", math_format(10^.x)),
                  breaks = 10^3*c(1e1,1e-2,1e-5)) +
    scale_x_discrete(labels = c("0" = "Survival", "1" = "Mortality")) +
    xlab('') + ylab('') + labs(title = 'Viruses')
  
  
  p = ggarrange(p1,p2,p3, ncol = 3, nrow = 1)
  return(p)
}



pd1 = myplot(df_long1 = df_long[df_long$Timepoint == 'D1',])
pd2 = myplot(df_long1 = df_long[df_long$Timepoint == 'D4',])
pd3 = myplot(df_long1 = df_long[df_long$Timepoint == 'D7',])
# !choose
pdf(paste0("Outputs/03_biomass_comparison.pdf"), width = 6, height = 2)
print(pd1)
print(pd2)
print(pd3)
dev.off()

