library(tidyverse)
library(magrittr)
library(gtsummary)
library(gt)
library(openxlsx)

load('Inputs/1315_meta_model.rdata')
df1 = meta_model

vars = c('Age', 'Gender', 
         'SOFA_24h', 'Immunosuppression', 'MV')

df_tab = df1 %>% dplyr::select(HumanID, Mortality28d, all_of(vars))


df2 = read.csv('Inputs/df_comorbities.csv')
df2 %<>% filter(HumanID %in% df1$HumanID)
vars = c('DM','MI','COPD','HepaticImpairment','RenalDisease','Tumor','HM','ConnectiveTissueDisease')

df2_tab <- df2 %>%
  mutate(HM = as.numeric(Leukemia | Lymphoma)) %>%
  select(HumanID, all_of(vars)) %>%
  mutate(across(-HumanID, ~ as.integer(. > 0)))
df2_tab <- df2_tab %>%
  mutate(across(-HumanID, ~ factor(. , levels=c(0,1))))

df_tab %<>% left_join(df2_tab, by = 'HumanID')

load('Inputs/1211_metadata.rdata')
df_tab %<>% left_join(meta %>% dplyr::select(HumanID, PneumoniaType), by = 'HumanID')

tmp = read.csv('../../AI_MASS/MASS_metadata_database_20250607/11.Medical_history_info.csv')
df_tab %<>% left_join(tmp %>% dplyr::select(HumanID, TransplantHistory), by = 'HumanID')
df_tab$TransplantHistory <- sub("\\..*", "", df_tab$TransplantHistory)
df_tab$TransplantHistory %<>% factor()
str(df_tab)

tmp = read.csv('../../AI_MASS/MASS_metadata_database_20250607/6.Experiment_info.csv')
tmp %<>% dplyr::select(HumanID, WBC, 
                       LymphocyteCount, NeutrophilCount, 
                       PlateletCount, PCT, hsCRP)
df_tab %<>% left_join(tmp, by = 'HumanID')
df_tab %<>% left_join(meta_model %>% dplyr::select(HumanID, APACHEII_24h), by = 'HumanID')
# comorbity分类变量显示为1行
vars = c('Immunosuppression', 'MV','DM','MI','COPD','HepaticImpairment',
         'RenalDisease','Tumor','HM','ConnectiveTissueDisease',
         'TransplantHistory')
value = purrr::map(vars, ~"1") |> rlang::set_names(vars)

dat = df_tab[,-1] %>%
  tbl_summary(
    by = Mortality28d,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p})"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)   # n 不保留小数，% 保留1位
    ),
    missing = "no",
    value = value
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      all_continuous() ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = function(x) style_pvalue(x, digits = 2)
  ) %>%
  bold_labels()


df_excel <- as_tibble(dat)
colnames(df_excel) = gsub('\\*','',colnames(df_excel))
df_excel$Characteristic = gsub('_','',df_excel$Characteristic)

write.xlsx(
  df_excel,
  file = "Outputs/01_table1.xlsx",
  rowNames = FALSE
)

write.csv(df_tab, 'Outputs/Source_data_table1.csv',
          row.names = F)
