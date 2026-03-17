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

# put into list
dfs <- list(
  main = df1,
  miall_glmnet = df2,
  miall_svm = df3,
  miall_rf = df4,
  miall_xgb = df5,
  selected = df6
)

# arrange all
dfs <- lapply(dfs, function(x) {
  x %>% arrange(Timepoint, Mortality_nd)
})

# write to excel with multiple sheets
write_xlsx(
  dfs,
  "Outputs/Supplementary_data_auc_all.xlsx"
)