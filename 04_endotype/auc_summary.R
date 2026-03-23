library(dplyr)
library(magrittr)
library(writexl)

# read
df1 <- read.csv("Outputs/Supplementary_data_auc_main.csv")
df2 <- read.csv("Outputs/Supplementary_data_auc_miall.csv")
df3 <- read.csv("Outputs/Supplementary_data_auc_miall_svm.csv")
df4 <- read.csv("Outputs/Supplementary_data_auc_miall_rf.csv")
df5 <- read.csv("Outputs/Supplementary_data_auc_miall_xgb.csv")
df6 <- read.csv("Outputs/Supplementary_data_auc_selected.csv")

df6$model = 'GLM'

df_all <- bind_rows(
  a_main = df1,
  b_selected = df6,
  c1_miall_glmnet = df2,
  c2_miall_svm = df3,
  c3_miall_rf = df4,
  c4_miall_xgb = df5,
  .id = "analysis"
)

# add ranking within grp for each analysis
df_all <- df_all %>%
  group_by(analysis, grp) %>%
  mutate(Rank = row_number(desc(median))) %>%
  ungroup() %>%
  arrange(analysis, Timepoint, Mortality_nd)
    

# write
write.csv(
  df_all,
  "Outputs/Supplementary_data_auc_all.csv",
  row.names = FALSE
)


df_all %>%
  filter(grepl("miall", analysis), Inputs == "CTS + microbe") %>%
  select(analysis, grp, Timepoint, Mortality_nd, median, Rank) %>%
  arrange(Timepoint, Mortality_nd)


df_all %>%
  filter(grepl("miall", analysis), Inputs == "CTS + microbe") %>%
  count(Rank)

df_all %>%
  filter(grepl("miall", analysis), Inputs == "CTS + microbe", 
         Timepoint == 'D1', Mortality_nd == 28) %>%
  count(Rank)

df_all %>%
  filter(grepl("miall", analysis), Inputs == "APACHEII", 
         Timepoint != 'D1', Mortality_nd == 28) %>%
  count(Rank)

df_all %>%
  filter(grepl("sel", analysis), Inputs == "CTS + microbe") %>%
  count(Rank)

