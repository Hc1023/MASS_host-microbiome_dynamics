rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(scales)
library(gridExtra)
library(ggforce)
library(ggpubr)
library(rstatix)
library(cowplot)

load('Inputs/1211_microbe.rdata')
load('Inputs/1211_metadata.rdata')
m = 0

myfun = function(m){
  compare_n = c('D1', 'D4', 'D7')
  
  df_long1 = df_long %>%
    filter(Timepoint %in% compare_n, Mortality28d == m) %>%
    droplevels()
  
  pathogens = rownames(data)
  data1 = data %>% select(df_long1$SampleID)
  dat = cbind(df_long1, t(data1))
  
  Fisher_results = data.frame(Pathogen = pathogens, P_val = NA)
  dat_prevalence = data.frame()
  for (i in seq_along(pathogens)) {
    bug = pathogens[i]
    dat_bug <- dat %>%
      select(Timepoint, all_of(bug)) %>%
      mutate(Present = .data[[bug]] > 0)
    tab <- table(dat_bug$Timepoint, dat_bug$Present)
    if (ncol(tab) > 1) {
      Fisher_results$P_val[i] <- fisher.test(tab)$p.value
    } else {
      Fisher_results$P_val[i] <- NA  # or some informative value like NA or 1
    }

    dat_bug = dat_bug %>% group_by(Timepoint) %>% 
      summarise(Prevalence = sum(Present)/n()*100) %>%
      mutate(Pathogen = bug)
    dat_prevalence = rbind(dat_prevalence, dat_bug)
  }
  Fisher_results %>%
    mutate(P_val = as.numeric(unlist(P_val))) %>% 
    arrange(P_val) -> Fisher_results
  Fisher_results = Fisher_results[1:18,]
  Fisher_results$P_val <- signif(Fisher_results$P_val, 3)
  levels_ = Fisher_results$Pathogen
  dat_prevalence %<>% filter(Pathogen %in% levels_)
  dat_prevalence$Pathogen = factor(dat_prevalence$Pathogen, levels = levels_)
  dat_prevalence$Timepoint = factor(dat_prevalence$Timepoint, levels = compare_n)
  
  tmp = dat_prevalence %>% group_by(Pathogen) %>% summarise(y_max = max(Prevalence)+1) %>% data.frame()
  annot_df = Fisher_results %>% mutate(x= 1) %>%
    left_join(tmp, by = 'Pathogen')
  annot_df$Pathogen = factor(annot_df$Pathogen, levels = Fisher_results$Pathogen)

  color_ = ifelse(m == 1, "#D73027", "#4575B4")
  levels(dat_prevalence$Pathogen)[levels(dat_prevalence$Pathogen) == 'HHV-4'] = 'EBV'
  levels(annot_df$Pathogen)[levels(annot_df$Pathogen) == 'HHV-4'] = 'EBV'
  p = ggplot(dat_prevalence , 
             aes(x = Timepoint, y = Prevalence)) +
    geom_bar(stat = "identity", width = 0.7, color = 'black', 
             fill = alpha(color_, 1), linewidth = 0.3) +
    # geom_line(aes(group = Pathogen)) +
    theme_bw() + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(panel.grid.minor = element_blank()) +
    geom_text(data = annot_df, 
              aes(x = x+1, y = y_max, label = paste0('p=',P_val)),
              inherit.aes = FALSE, vjust = 0, size = 3) +
    facet_wrap(~ Pathogen, scales = "free_y", ncol = 6)

  return(p)
}
p1 = myfun(1)
p0 = myfun(0)

pdf(paste0("Outputs/03_longitudinal_prevalence.pdf"), width = 8.5, height = 4.5)
print(p0)
print(p1)
dev.off()



