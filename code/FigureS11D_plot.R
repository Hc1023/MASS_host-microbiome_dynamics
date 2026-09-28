# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
library(dplyr)

library(ggplot2)

timepoints <- c("D3", "D5")

cts_levels <- paste0("CTS", 1:3)

prefix <- "Outputs"

source("code/shared_CMAISE_visit_rules.R")
e<-new.env();load("Inputs/260827_SRR_meta_model.rdata",e)
calls<-read.csv("Inputs/260830_SRR_CTS_classification.csv")
cohort<-calls %>% select(HumanID,Timepoint,CTS) %>% filter(Timepoint %in% timepoints) %>%
 left_join(e$meta_model %>% select(HumanID,Mortality28d),by="HumanID")
cohort<-bind_rows(lapply(timepoints,function(tp)cohort %>% filter(Timepoint==tp,cmaise_keep_visit(HumanID,tp))))
mortality_summary<-cohort %>% group_by(Timepoint,CTS) %>% summarise(n=n(),deaths=sum(Mortality28d==1),mortality_proportion=deaths/n,.groups="drop")
write.csv(mortality_summary,"Outputs/FigureS11D_mortality_summary.csv",row.names=FALSE)

p_heat <- ggplot(mortality_summary, aes(x = factor(CTS, levels = cts_levels), y = factor(Timepoint, levels = rev(timepoints)), 
    fill = mortality_proportion)) + geom_tile(color = "white", linewidth = 0.8) + geom_text(aes(label = ifelse(is.na(mortality_proportion), 
    "NA", sprintf("%.1f%%\n%d/%d", 100 * mortality_proportion, deaths, n))), size = 3.5) + scale_fill_gradient(low = "#f7fbff", 
    high = "#b2182b", limits = c(0, max(mortality_summary$mortality_proportion, na.rm = TRUE)), labels = scales::percent, 
    na.value = "grey90") + labs(x = NULL, y = NULL, fill = "28-day \nmortality") + 
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank())

ggsave(repo_output(paste0(prefix, "_mortality_heatmap.pdf"), "FigureS11D_plot.R"), p_heat, width = 4, 
    height = 2)

