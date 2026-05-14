rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)

load('Inputs/1616_meta_model.rdata')
sel_cts = read.csv('Outputs/Supplementary_data_selection_ctsg.csv')
sel_mi = read.csv('Outputs/Supplementary_data_selection_microbe.csv')

#### modeling ####
library(ggpubr)
library(e1071)    # SVM 模型
library(caret)    # 数据切分 + 评价
library(pROC)     # ROC AUC
library(glmnet)

pred_fun = function(var, m = 'GLM', tps = c(3,7,28)){
  set.seed(123)
  n_iter = 500
  auc_mat = matrix(nrow = n_iter, ncol = length(tps))
  colnames(auc_mat) = paste0('M', tps)
  for (i in tps) {
    print(i)
    
    formula <- as.formula(
      paste(paste0("Mortality",i,"d ~"), 
            paste(var, collapse = " + ")))
    
    
    auc_values <- numeric(n_iter)
    
    for(i_iter in 1:n_iter){
      
      # --- 随机划分 ---
      train_index <- createDataPartition(df_G1[[paste0('Mortality',i,'d')]], p = 0.7, list = FALSE)
      train <- df_G1[train_index, ]
      test  <- df_G1[-train_index, ]
      
      if(m == 'GLM'){
        # --- GLM（logistic）---
        glm_model <- glm(
          formula,
          data = train,
          family = binomial
        )
        
        # --- 预测（概率）---
        pred_prob <- predict(
          glm_model,
          newdata = test,
          type = "response"
        )
      }
      
      if(m == "GLMNET"){
        # y 必须是 0/1 numeric
        y_train <- as.integer(as.character(train[[paste0('Mortality',i,'d')]]))
        y_test  <- as.integer(as.character(test[[paste0('Mortality',i,'d')]]))  # 用于 roc（可直接用 test$Mortality28d 也行）
        
        # X matrix（去掉截距）
        X_train <- model.matrix(formula, data = train)[, -1, drop = FALSE]
        X_test  <- model.matrix(formula, data = test )[, -1, drop = FALSE]
        
        fit <- cv.glmnet(
          x = X_train,
          y = y_train,
          family = "binomial",
          alpha = 0.5
        )
        
        pred_prob <- as.numeric(
          predict(fit, newx = X_test, s = "lambda.min", type = "response")
        )
      }
      
      if (m == "RF") {
        # --- Random Forest ---
        library(randomForest)
        
        rf_model <- randomForest(
          formula,
          data = train,
          ntree = 1000,
          mtry  = floor(sqrt(ncol(train) - 1)),
          importance = TRUE
        )
        
        # --- 预测 ---
        pred_prob <- predict(
          rf_model,
          test,
          type = "prob"
        )[,"1"]
      }
      # --- AUC ---
      outcome <- paste0("Mortality", i, "d")
      auc_values[i_iter] <- auc(roc(test[[outcome]], pred_prob, quiet = TRUE))
      
    }
    auc_mat[,paste0('M',i)] = auc_values
    print(mean(auc_values))
    # cat("\n", outcome, "\n")
    # print(table(df_G1[[outcome]]))
    # print(tapply(df_G1[[var]], df_G1[[outcome]], summary))
  }
  return(auc_mat)
  
}


df_G1 = meta_model %>% dplyr::select(Mortality28d, SurvivalTimeWithin28Days, APACHEII_24h, SOFA_24h,
                                     starts_with("D1_"))
df_G1 %<>% na.omit()

for (i in c(3,7)) {
  y = df_G1$SurvivalTimeWithin28Days <= i
  df_G1[[paste0('Mortality',i,'d')]] = factor(as.numeric(y))
  # print(table(y))
}

tmp = sel_cts[,2]
names(tmp) = sel_cts[,1]
stable_vars_CTSg <- paste0('D1_CTSg_', names(rev(sort(tmp))[1:5]))
tmp = sel_mi[,2]
names(tmp) = sel_mi[,1]
stable_vars_m <- paste0('D1_mi_', names(rev(sort(tmp))[1:2]))

auc_mat1 = pred_fun(var = "APACHEII_24h")
auc_mat2 = pred_fun(var = stable_vars_CTSg)
auc_mat3 = pred_fun(var = stable_vars_m)
auc_mat4 = pred_fun(var = c(stable_vars_CTSg, stable_vars_m))

head(auc_mat1)

auc_long <- bind_rows(
  as.data.frame(auc_mat1) %>% mutate(model = "APACHEII"),
  as.data.frame(auc_mat2) %>% mutate(model = "CTS"),
  as.data.frame(auc_mat3) %>% mutate(model = "microbe"),
  as.data.frame(auc_mat4) %>% mutate(model = "CTS + microbe")
) %>%
  pivot_longer(cols = c(M3, M7, M28), names_to = "endpoint", values_to = "AUC") %>%
  mutate(
    endpoint = factor(endpoint, levels = c("M3", "M7", "M28"),
                      labels = c("Mortality 3d", "Mortality 7d", "Mortality 28d")),
    model = factor(model, levels = c("APACHEII", "CTS", "microbe", "CTS + microbe"))
  )

model_cols <- c(
  "APACHEII" = "#7F7F7F",   # baseline grey
  "CTS"         = "#1B9E77",   # green
  "microbe"            = "#7570B3",   # purple
  "CTS + microbe"     = "#D95F02"    # orange (combined model)
)

p = ggplot(auc_long, aes(x = model, y = AUC, fill = model)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(
    aes(color = model),
    width = 0.12, height = 0,
    alpha = 0.2, size = 1.2,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    fatten = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  facet_wrap(~ endpoint, nrow = 1) +
  scale_fill_manual(values = model_cols) +
  scale_color_manual(values = model_cols) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = NULL, y = "AUC") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c(
    "APACHEII" = "APACHEII",
    "CTS"  = "CTS",
    "microbe"    = "Microbes",
    "CTS + microbe"     = "CTS+Microbes"
  ))

print(p)


#### D4 and D7 ####

df_G1 = meta_model %>% dplyr::select(Mortality28d, SurvivalTimeWithin28Days, APACHEII_24h, SOFA_24h,
                                     starts_with("D4_"))
df_G1 %<>% na.omit()
table(df_G1$Mortality28d)

tmp = sel_cts[,3]
names(tmp) = sel_cts[,1]
stable_vars_CTSg <- paste0('D4_CTSg_', names(rev(sort(tmp))[1:5]))
tmp = sel_mi[,3]
names(tmp) = sel_mi[,1]
stable_vars_m <- paste0('D4_mi_', names(rev(sort(tmp))[1:2]))

auc_mat1_D4 = pred_fun(var = "APACHEII_24h", tps = 28)
auc_mat2_D4 = pred_fun(var = stable_vars_CTSg, tps = 28)
auc_mat3_D4 = pred_fun(var = stable_vars_m, tps = 28)
auc_mat4_D4 = pred_fun(var = c(stable_vars_CTSg, stable_vars_m), tps = 28)

## D7
df_G1 = meta_model %>% dplyr::select(Mortality28d, SurvivalTimeWithin28Days, APACHEII_24h, SOFA_24h,
                                     starts_with("D7_"))
df_G1 %<>% na.omit()

tmp = sel_cts[,4]
names(tmp) = sel_cts[,1]
stable_vars_CTSg <- paste0('D7_CTSg_', names(rev(sort(tmp))[1:5]))
tmp = sel_mi[,4]
names(tmp) = sel_mi[,1]
stable_vars_m <- paste0('D7_mi_', names(rev(sort(tmp))[1:2]))

auc_mat1_D7 = pred_fun(var = "APACHEII_24h", tps = 28)
auc_mat2_D7 = pred_fun(var = stable_vars_CTSg, tps = 28)
auc_mat3_D7 = pred_fun(var = stable_vars_m, tps = 28)
auc_mat4_D7 = pred_fun(var = c(stable_vars_CTSg, stable_vars_m), tps = 28)

auc_mat1_D = bind_cols(auc_mat1[,3], auc_mat1_D4, auc_mat1_D7)
colnames(auc_mat1_D) = c('D1','D4','D7')
auc_mat2_D = bind_cols(auc_mat2[,3], auc_mat2_D4, auc_mat2_D7)
colnames(auc_mat2_D) = c('D1','D4','D7')
auc_mat3_D = bind_cols(auc_mat3[,3], auc_mat3_D4, auc_mat3_D7)
colnames(auc_mat3_D) = c('D1','D4','D7')
auc_mat4_D = bind_cols(auc_mat4[,3], auc_mat4_D4, auc_mat4_D7)
colnames(auc_mat4_D) = c('D1','D4','D7')


auc_long_D <- bind_rows(
  as.data.frame(auc_mat1_D) %>% mutate(model = "APACHEII"),
  as.data.frame(auc_mat2_D) %>% mutate(model = "CTS"),
  as.data.frame(auc_mat3_D) %>% mutate(model = "microbe"),
  as.data.frame(auc_mat4_D) %>% mutate(model = "CTS + microbe")
) %>%
  pivot_longer(cols = c(D1, D4, D7), names_to = "Tp", values_to = "AUC") %>%
  mutate(
    Tp = factor(Tp, levels = c("D1", "D4", "D7")),
    model = factor(model, levels = c("APACHEII", "CTS", "microbe", "CTS + microbe"))
  )

model_cols <- c(
  "APACHEII" = "#7F7F7F",   # baseline grey
  "CTS"         = "#1B9E77",   # green
  "microbe"            = "#7570B3",   # purple
  "CTS + microbe"     = "#D95F02"    # orange (combined model)
)

p = ggplot(auc_long_D, aes(x = model, y = AUC, fill = model)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(
    aes(color = model),
    width = 0.12, height = 0,
    alpha = 0.2, size = 1.2,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    fatten = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  facet_wrap(~ Tp, nrow = 1) +
  scale_fill_manual(values = model_cols) +
  scale_color_manual(values = model_cols) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = NULL, y = "AUC") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c(
    "APACHEII" = "APACHEII",
    "CTS"  = "CTS",
    "microbe"    = "Microbes",
    "CTS + microbe"     = "CTS+Microbes"
  ))


print(p)


df_all <- bind_rows(
  auc_long %>% rename(grp = endpoint),
  auc_long_D %>% rename(grp = Tp) %>% filter(grp != 'D1')
)
df_all$grp <- as.character(df_all$grp)
df_all$grp[df_all$grp == "Mortality 28d"] <- "Mortality 28d-D1"
df_all$grp <- factor(df_all$grp, levels = as.character(unique(df_all$grp)))


p = ggplot(df_all, aes(x = model, y = AUC, fill = model)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(
    aes(color = model),
    width = 0.12, height = 0,
    alpha = 0.1, size = 1.2,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    fatten = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  facet_wrap(~ grp, nrow = 1) +
  scale_fill_manual(values = model_cols) +
  scale_color_manual(values = model_cols) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             color = "grey50",
             linewidth = 0.6) +
  # scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = NULL, y = "AUC") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = c(
    "APACHEII" = "APACHEII",
    "CTS"  = "CTS",
    "microbe"    = "Microbes",
    "CTS + microbe"     = "CTS+Microbes"
  ))

pdf(paste0("Outputs/04_auc_all_sel.pdf"), width = 7, height = 2.5)
print(p)
dev.off()

if(F){
  res1 = df_all %>%
    group_by(model, grp) %>%
    summarise(
      median = median(AUC, na.rm = TRUE),
      Q1 = quantile(AUC, 0.25, na.rm = TRUE),
      Q3 = quantile(AUC, 0.75, na.rm = TRUE),
      AUC_summary = sprintf("%.3f (%.3f-%.3f)", median, Q1, Q3),
      .groups = "drop"
    )
  colnames(res1)[1] = 'Inputs'
  res1$Timepoint = 'D1'
  res1$Timepoint[res1$grp %in% c('D4', 'D7')] = as.character(res1$grp[res1$grp %in% c('D4', 'D7')])
  
  res1$Mortality_nd = 28
  res1$Mortality_nd[grepl('3d', res1$grp)] = 3
  res1$Mortality_nd[grepl('7d', res1$grp)] = 7
  
  res1$CTSg = NA
  res1$microbe = NA
  
  tmp = sel_cts[,2]
  names(tmp) = sel_cts[,1]
  res1$CTSg[res1$Timepoint == 'D1' & grepl('CTS', res1$Inputs)] = paste(names(rev(sort(tmp))[1:5]),collapse = ",")
  
  tmp = sel_cts[,3]
  names(tmp) = sel_cts[,1]
  res1$CTSg[res1$Timepoint == 'D4' & grepl('CTS', res1$Inputs)] = paste(names(rev(sort(tmp))[1:5]),collapse = ",")
  
  tmp = sel_cts[,4]
  names(tmp) = sel_cts[,1]
  res1$CTSg[res1$Timepoint == 'D7' & grepl('CTS', res1$Inputs)] = paste(names(rev(sort(tmp))[1:5]),collapse = ",")
  
  
  tmp = sel_mi[,2]
  names(tmp) = sel_mi[,1]
  res1$microbe[res1$Timepoint == 'D1' & grepl('mi', res1$Inputs)] = paste(names(rev(sort(tmp))[1:2]),collapse = ",")
  
  tmp = sel_mi[,3]
  names(tmp) = sel_mi[,1]
  res1$microbe[res1$Timepoint == 'D4' & grepl('mi', res1$Inputs)] = paste(names(rev(sort(tmp))[1:2]),collapse = ",")
  
  tmp = sel_mi[,4]
  names(tmp) = sel_mi[,1]
  res1$microbe[res1$Timepoint == 'D7' & grepl('mi', res1$Inputs)] = paste(names(rev(sort(tmp))[1:2]),collapse = ",")
  
  
  write.csv(res1, 'Outputs/Supplementary_data_selection_auc.csv',
            row.names = F)
  
}
