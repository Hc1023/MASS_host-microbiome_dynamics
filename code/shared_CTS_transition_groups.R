# Shared CTS grouping definition for the current S10 analysis.
df1 = meta_model %>% dplyr::select(HumanID, Mortality28d, D1_CTS, D4_CTS, D7_CTS)

df1 %<>% na.omit()

df1$grp = "Mixed"

idx = df1$D1_CTS == 3 & df1$D7_CTS != 3

df1$grp[idx] = "Deteriorating"

idx = df1$D1_CTS != 3 & df1$D7_CTS == 3

df1$grp[idx] = "Improving"

idx = df1$D1_CTS == 3 & df1$D4_CTS == 3 & df1$D7_CTS == 3

df1$grp[idx] = "Stable low-risk"

idx = df1$D1_CTS != 3 & df1$D4_CTS != 3 & df1$D7_CTS != 3

df1$grp[idx] = "Stable high-risk"

df_sum <- NULL
