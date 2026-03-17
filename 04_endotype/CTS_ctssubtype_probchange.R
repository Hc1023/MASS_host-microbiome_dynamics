rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)

load('Inputs/1315_meta_model.rdata')

meta_analysis = meta_model %>% filter(!is.na(D1) & !is.na(D4.death) & !is.na(D7.death))
meta_analysis %<>% dplyr::select(HumanID, Mortality28d, matches("CTS"))
meta_analysis %<>% dplyr::select(-matches("CTSg"))
meta_analysis %<>% dplyr::select(-matches("D14"), -matches("D21"))

cts_type = 2

df1 = meta_analysis %>% 
  dplyr::select(HumanID, Mortality28d, matches(paste0("D[147]_CTS", cts_type)), D1_CTS) %>%
  filter(D1_CTS == cts_type)

df1_long = df1 %>% pivot_longer(
  cols = matches("^D(1|4|7)_CTS[123]$"),
  names_to = c("Timepoint"),
  names_pattern = "(D[147])_CTS[123]",
  values_to = "CTS_value"
)
meta_analysis2 = df1_long
meta_analysis2$TimeDays <- as.numeric(sub('D','',meta_analysis2$Timepoint))

meta_analysis2 %<>% left_join(meta_model %>% 
                                dplyr::select(HumanID, 
                                              SurvivalTimeWithin28Days),
                              by = 'HumanID')

surv_df = meta_analysis2 %>% 
  dplyr::select(HumanID, Mortality28d, SurvivalTimeWithin28Days) %>%
  distinct() %>%
  mutate(
    Mortality28d = as.numeric(as.character(Mortality28d))  # 如果是factor: "0"/"1"
  )

ggplot(meta_analysis2,
       aes(x = Timepoint,
           y = CTS_value,
           fill = factor(Mortality28d))) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA
  ) +
  scale_fill_manual(values = alpha(c("0" = "#4575B4", "1" = "#D73027"),0.7)) +
  theme_bw()


cts = 'CTS_value'

library(nlme)
lme_fit <- lme(
  fixed  = as.formula(paste0(cts, "~ TimeDays * Mortality28d")),
  random = ~ TimeDays | HumanID,  
  data   = meta_analysis2,
  na.action = na.exclude,
  control = lmeControl(opt = "optim")
)
summary(lme_fit)

library(survival)

cox_fit <- coxph(
  Surv(SurvivalTimeWithin28Days, Mortality28d) ~ 1 + cluster(HumanID),
  data = surv_df,
  x = TRUE
)
library(JM)
jm_fit <- jointModel(
  lmeObject = lme_fit,
  survObject = cox_fit,
  timeVar = "TimeDays",
  method = "weibull-PH-aGH"   # 常用且稳；也可试 "spline-PH-aGH"
)

summary(jm_fit)

#### prediction ####
pred_grid <- expand.grid(
  TimeDays = c(1, 4, 7),
  Mortality28d = c(0, 1)
)

# JM 要求一个“真实存在过”的 ID 占位
pred_grid$HumanID <- meta_analysis2$HumanID[1]

pred_grid$Group <- ifelse(pred_grid$Mortality28d == 0,
                          "Survivor", "Non-survivor")
pred_grid$fit_jm <- as.numeric(
  predict(
    jm_fit,
    newdata = pred_grid,
    process = "Longitudinal",
    type = "Marginal"
  )
)
pred_grid$Timepoint = factor(paste0('D', pred_grid$TimeDays))

str(pred_grid)
library(ggplot2)

dodge_w <- 0.35

library(ggnewscale)
library(ggpubr)
library(scales)

cols_strong <- c("0" = "#4575B4", "1" = "#D73027")
cols_light  <- c("0" = alpha("#4575B4", 0.7),
                 "1" = alpha("#D73027", 0.7))




p1 = ggplot() +
  ## 1) boxplot（淡色边框）
  geom_boxplot(
    data = meta_analysis2,
    aes(x = Timepoint, y = .data[[cts]],
        group = interaction(Timepoint, Mortality28d),
        color = factor(Mortality28d)),
    width = 0.25,
    outlier.shape = NA,
    position = position_dodge(width = dodge_w),
    alpha = 0.1
  ) +
  scale_color_manual(
    values = cols_light,
    labels = c("0" = "Survivor", "1" = "Mortality"),
    name = 'Observed'
  ) +
  ggnewscale::new_scale_color() +
  
  ## 2) JM fitted mean（深色点线）
  geom_point(
    data = pred_grid,
    aes(x = Timepoint, y = fit_jm,
        group = factor(Mortality28d),
        color = factor(Mortality28d)),
    position = position_dodge(width = dodge_w),
    size = 3
  ) +
  geom_line(
    data = pred_grid,
    aes(x = Timepoint, y = fit_jm,
        group = factor(Mortality28d),
        color = factor(Mortality28d)),
    position = position_dodge(width = dodge_w),
    linewidth = 1
  ) +
  scale_color_manual(
    values = cols_strong,
    labels = c("0" = "Survival", "1" = "Mortality"),
    name = "JM fit"
  ) +
  
  ## 3) 每个时间点 Wilcoxon 星号（用 meta_analysis2）
  ggpubr::stat_compare_means(
    data = meta_analysis2,
    aes(x = Timepoint, y = .data[[cts]], group = Mortality28d),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 5,
    alpha = 0.7
  ) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1),
                     expand = expansion(mult = c(0.05, 0.12))) +
  theme_bw(base_size = 12) +
  labs(x = "Timepoint", y = cts, color = "Outcome")


pdf(paste0("Outputs/1329_JMsubtype2.pdf"), width = 3.5, height = 2.2)
print(p1)
dev.off()

