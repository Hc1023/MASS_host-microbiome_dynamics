rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)

load('../../AI_MASS/Utility/metadata.rdata')
metadata1 = metadata[!is.na(metadata$Mortality28d),c(1:3,8:12)]
data = read.csv('../../AI_MASS/Blood/BLgenusercc_1101.csv', row.names = 1)
mapid = read.csv('../../AI_MASS/Utility/BLmicrobe_mapIds.csv')
mapid = mapid[!duplicated(mapid$Genus),-1]
rownames(mapid) = mapid$Genus

df_long <- metadata1 %>%
  pivot_longer(
    cols = starts_with("BL_"),   # 要转换的列
    names_to = "Timepoint",           # 新列存放原列名
    values_to = "SampleID",            # 新列存放原值
    values_drop_na = TRUE             # 如果想去掉 NA 行
  )
df_long$Timepoint = sub("BL_", "", df_long$Timepoint)
df_long = data.frame(df_long)

# !choose
d = 'D1'
df_long1 = df_long[df_long$Timepoint == d,]

data1 = data[,df_long1$SampleID]
df_long1$Microbial_mass = colSums(data1)


median_ =median(df_long1$Microbial_mass)

df_long1$Group_mass = 0
df_long1$Group_mass[df_long1$Microbial_mass > median_] = 1


library(survival)
library(survminer)
table(is.na(df_long1[,7]))
df_clean <- na.omit(df_long1[, c("SurvivalTimeWithin28Days", "Mortality28d", "Group_mass")])
df_clean$Mortality28d = as.numeric(as.character(df_clean$Mortality28d))

# Create survival object
surv_obj <- Surv(time = df_clean$SurvivalTimeWithin28Days,
                 event = df_clean$Mortality28d)

# Fit model
fit <- survfit(surv_obj ~ Group_mass, data = df_clean)

# Plot



p1 = ggsurvplot(
  fit,
  data = df_clean,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.labs = c('High mass','Low mass'),
  font.legend = c(14, "plain", "black"),
  legend.title = "",
  # legend.labs = levels(df_clean$Group_mass),
  palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0, 28),
  break.time.by = 7
) 


p1

# 筛选死亡组
df_death <- df_clean[df_clean$Mortality28d == 1, ]

# 构建生存对象
surv_obj <- Surv(time = df_death$SurvivalTimeWithin28Days,
                 event = df_death$Mortality28d)

# 拟合生存曲线
fit_death <- survfit(surv_obj ~ 1, data = df_death)

p2 = ggsurvplot(
  fit_death,
  data = df_death,
  risk.table = TRUE,     # 风险表
  pval = TRUE,           # p值（对于单组无意义，可去掉）
  conf.int = TRUE,       # 置信区间
  xlab = "Time (days)",
  ylab = "Survival Probability",
  palette = c("#c43932"),  # 只用一条曲线颜色
  legend.labs = "Mortality",
  font.legend = c(14, "plain", "black"),
  legend.title = "",
  xlim = c(0, 28),
  break.time.by = 7,
  censor = FALSE         # 不显示删失点
)


pdf(paste0("Outputs/03_surv.pdf"), width = 4, height = 3)
print(p1)
print(p2)
dev.off()
