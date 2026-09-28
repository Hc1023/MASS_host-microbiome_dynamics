source('code/shared_sources.R')
scores<-readRDS('Outputs/Figure3C_prepare_MASS_GSVA_scores.rds')
saveRDS(scores,'Outputs/Module_scores.rds')
