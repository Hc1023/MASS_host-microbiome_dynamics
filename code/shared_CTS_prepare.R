# Use saved CTS calls and probabilities; no classifier refit.
source("code/shared_paths.R")
library(tidyverse)
library(magrittr)
load("Inputs/1211_metadata.rdata")
cts_env <- new.env()
load("Inputs/1616_meta_model.rdata", cts_env)
cts_saved <- bind_rows(lapply(c("D1", "D4", "D7"), function(tp) {
  z <- cts_env$meta_model[, c("HumanID", tp, paste0(tp, c("_CTS", "_CTS1", "_CTS2", "_CTS3")))]
  names(z) <- c("HumanID", "SampleID", "CTS", "CTS1", "CTS2", "CTS3")
  z
})) %>% filter(!is.na(SampleID))
stopifnot(!anyDuplicated(cts_saved$SampleID))
df_long <- left_join(df_long, select(cts_saved, -HumanID), by="SampleID")
meta_analysis = meta %>% filter(!is.na(D1) & !is.na(D4.death) & !is.na(D7.death))

table(meta_analysis$Mortality28d)

meta_analysis %<>% dplyr::select(HumanID, Mortality28d, D1, D4 = D4.death, D7 = D7.death) %>% pivot_longer(cols = c(D1, 
    D4, D7), names_to = "Timepoint", values_to = "SampleID")

tmp = df_long %>% dplyr::select(SampleID, CTS)

meta_analysis %<>% left_join(tmp, by = "SampleID")

meta_analysis$CTS = factor(meta_analysis$CTS, levels = c("Death", "Alive", "1", "2", "3"))

meta_analysis$CTS[meta_analysis$SampleID == "Death"] = "Death"

d28_rows <- meta_analysis %>% distinct(HumanID, Mortality28d) %>% mutate(Timepoint = "D28", SampleID = NA_character_, 
    CTS = ifelse(Mortality28d == 1, "Death", "Alive"))

meta_analysis %<>% bind_rows(d28_rows) %>% mutate(Timepoint = factor(Timepoint, levels = c("D1", "D4", 
    "D7", "D28"))) %>% arrange(HumanID, Timepoint)

df_plot = meta_analysis %>% dplyr::select(HumanID, Timepoint, State = CTS)

df_plot$State %<>% factor(., levels = c("Death", "Alive", "1", "2", "3"))


stopifnot(!anyNA(meta_analysis$CTS))
