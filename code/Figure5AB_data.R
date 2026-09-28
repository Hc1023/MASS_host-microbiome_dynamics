# Derive panel data from the original saved CTS labels without reclassification.
source("code/shared_CTS_prepare.R")
states <- df_plot
flows <- meta_analysis %>% arrange(HumanID,Timepoint) %>% group_by(HumanID) %>%
  mutate(Next_timepoint=lead(Timepoint),Next_state=lead(CTS)) %>% ungroup() %>%
  filter(!is.na(Next_timepoint)) %>% count(Timepoint,CTS,Next_timepoint,Next_state,name="N")
mortality <- meta_analysis %>% filter(Timepoint %in% c("D1","D4","D7"),CTS %in% c("1","2","3")) %>%
  group_by(Timepoint,CTS) %>% summarise(N=n(),Deaths=sum(Mortality28d==1),Mortality_probability=Deaths/N,.groups="drop")
write.csv(states,"Outputs/Figure5A_reconstructed_states.csv",row.names=FALSE)
write.csv(flows,"Outputs/Figure5A_reconstructed_flows.csv",row.names=FALSE)
write.csv(mortality,"Outputs/Figure5B_reconstructed_mortality.csv",row.names=FALSE)
