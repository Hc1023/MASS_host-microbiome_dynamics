# Reconstruct plotted repeat-level AUCs without fitting models or rerunning CV.
library(dplyr)
library(tidyr)
oof <- read.csv("Outputs/FigureS11E_analysis_SRR_CTS_nestedCV_oof_predictions.csv")
stopifnot(!anyDuplicated(oof[c("Timepoint","Repeat","HumanID")]),length(unique(oof$Repeat))==20)
calc_auc <- function(y,p) as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)))
performance <- oof %>% pivot_longer(c(M0_probability,M1_probability),names_to="Model",values_to="Probability") %>%
  mutate(Model=sub("_probability$","",Model)) %>% group_by(Timepoint,Repeat,Model) %>%
  summarise(AUC=calc_auc(Mortality28d,Probability),Brier=mean((Mortality28d-Probability)^2),.groups="drop")
wide <- performance %>% pivot_wider(names_from=Model,values_from=c(AUC,Brier)) %>%
  mutate(DeltaAUC=AUC_M1-AUC_M0,DeltaBrier=Brier_M1-Brier_M0)
contrasts <- wide %>% group_by(Timepoint) %>% summarise(Mean_DeltaAUC=mean(DeltaAUC),
  Repeat_SE=sd(DeltaAUC)/sqrt(n()),N_repeats=n(),.groups="drop")
write.csv(performance,"Outputs/FigureS11E_repeat_performance.csv",row.names=FALSE)
write.csv(wide %>% select(Timepoint,Repeat,DeltaAUC,DeltaBrier),"Outputs/FigureS11E_paired_differences.csv",row.names=FALSE)
write.csv(contrasts,"Outputs/FigureS11E_paired_descriptive.csv",row.names=FALSE)
