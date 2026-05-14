rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(purrr)
library(ggplot2)
library(stringr)
library(forcats)

load('Inputs/1616_meta_model.rdata')
load('Inputs/1616_microbe.rdata')
# microbe: log2 transformed mass

## logistic-----
library(jstable)       # 用于亚组分析
library(dplyr)         # 数据处理
library(forestploter)  # 绘制森林图
library(grid)          # 基础绘图支持
library(survival)      # 获取示例数据集
library(purrr)

myfun <- function(d = "D1",
                  pathogen_col = c("D1_mi2_Bacf", "D1_mi2_Viruses", 
                                   "D1_mi_HHV.4","D1_mi_HCMV"),
                  vars = c('Gender', 'Age', 'CenterGroup', 'PneumoniaTypeGroup',
                           'CCI', 'SOFA_24h', 'Immunosuppression', 'MV')
) {
  
  cts_var <- paste0(d, "_CTS")
  
  df1 <- meta_model %>%
    dplyr::select(Mortality28d, starts_with(paste0(d, "_")), all_of(vars)) %>%
    tidyr::drop_na() 
  
  df1$Mortality28d <- as.numeric(as.character(df1$Mortality28d))

  # ---- 跑模型 ----
  res_D <- map_dfr(pathogen_col, function(p) {
    
    fml <- as.formula(paste("Mortality28d ~", p))
    
    out <- tryCatch(
      {
        TableSubgroupMultiGLM(
          formula = fml,
          var_subgroups = cts_var,
          var_cov = vars,
          data = df1
        ) %>%
          mutate(pathogen = p, timepoint = d)
      },
      error = function(e) {
        tibble(
          pathogen = p,
          timepoint = d,
          error = e$message
        )
      }
    )
    
    out
  })
  return(res_D)

}

res_D1 <- myfun("D1")


#### HCMV, HHV.4: draw forest plot ####

res_all_filtered = res_D1
res_all_filtered <- res_all_filtered %>%
  mutate(
    Variable = ifelse(
      grepl("^D[147]_CTS$", Variable),
      pathogen,
      Variable
    )
  )
res_all_filtered$Variable = sub("_mi2?_", " ", res_all_filtered$Variable)
res_all_filtered$Variable = sub("HHV\\.4", "EBV", res_all_filtered$Variable)
res_all_filtered$Variable = sub("Bacf", "Bacteria/Fungi", res_all_filtered$Variable)

res_all_filtered <- res_all_filtered %>%
  mutate(
    Variable_trim = trimws(Variable),
    order_flag = case_when(
      Variable_trim == "Overall" ~ 4,
      Variable_trim %in% c("1", "2", "3") ~ as.numeric(Variable_trim),
      TRUE ~ 0
    )
  ) %>%
  group_by(timepoint, pathogen) %>%
  arrange(order_flag, .by_group = TRUE) %>%
  ungroup() %>%
  dplyr::select(-Variable_trim, -order_flag)

res = res_all_filtered[,c(1:8)]
res %<>%
  mutate(
    Variable = case_when(
      trimws(Variable) %in% c("1", "2", "3") ~ paste0("  CTS", trimws(Variable)),
      Variable == "Overall" ~ "  Total",
      TRUE ~ Variable
    )
  )
res = res[,-c(2,3)]

res <- res %>%
  mutate(across(c(OR, Lower, Upper),
                ~ suppressWarnings(as.numeric(.x))))
res$"OR (95% CI)" <- ifelse(
  is.na(res$"OR"),
  "",
  sprintf("%.2f (%.2f to %.2f)", 
          res$OR, res$Lower, res$Upper)
)
res[,c(5,6,7)][is.na(res[, c(5,6,7)])] = " "
# 添加空白列用于森林图布局
res$` ` <- paste(rep(" ", 15), collapse = " ")

colnames(res)[c(1,6)] = c("Subgroup", "P interaction")
res2 = res
res2$`P interaction` = paste0("     ", res2$`P interaction`)
# ---------------------------------------------------4.绘制森林图
# 绘制基础森林图
plot_sub <- forest(
  data = res2[, c(1, 7, 8, 5, 6)],  # 按顺序选择展示列
  lower = res2$Lower,                   # 置信区间下限
  upper = res2$Upper,                   # 置信区间上限
  est = res2$OR,                      # 点估计值（OR）
  ci_column = 3,                       # 点估计对应列位置
  ref_line = 1#,                        # 参考线位置（OR=1）
  # xlim = c(0, 3)                       # X轴范围
) %>%
  add_border(part = "header", row = 1, where = "top") %>%
  add_border(part = "header", row = 1, where = "bottom")

# 输出森林图
pdf(paste0("Outputs/04_forest.pdf"), width = 8, height = 6)
plot(plot_sub)
dev.off()
