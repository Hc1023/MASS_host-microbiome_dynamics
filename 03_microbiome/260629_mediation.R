rm(list = ls())
library(tidyverse)
library(magrittr)
library(mediation)

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)

out_dir <- file.path(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260629_med"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# -------------------------
# 1) 读入 + 构建 D1/D4 wide + delta
# -------------------------
dfmm <- read.csv("Inputs/microbe_hostmodule.csv")
# 这里的microbe已经是log2(fg+1)

compare_n <- c("D1","D4")

df_sub <- dfmm %>%
  filter(Timepoint %in% compare_n) %>%
  group_by(HumanID) %>%
  filter(n_distinct(Timepoint) == length(compare_n)) %>%
  ungroup()

df_wide <- df_sub %>%
  dplyr::select(
    HumanID, Timepoint,
    up1, up2, up3, dw1,
    HCMV,
    Gender, Age, CenterGroup, PneumoniaTypeGroup,
    CCI, SOFA_24h, Immunosuppression, MV,
    Mortality28d, SurvivalTimeWithin28Days
  ) %>%
  pivot_wider(
    names_from  = Timepoint,
    values_from = c(up1, up2, up3, dw1, HCMV)
  )

df_med <- df_wide %>%
  mutate(
    d_up1  = up1_D4 - up1_D1,
    d_up2  = up2_D4 - up2_D1,
    d_up3  = up3_D4 - up3_D1,
    d_dw1  = dw1_D4 - dw1_D1,
    d_HCMV = HCMV_D4 - HCMV_D1
  )

# -------------------------
# 2) 协变量处理：factor -> dummy；连续变量 scale
# -------------------------
X <- "HCMV_D4"               # treat
Y <- "Mortality28d"          # outcome
modules <- c("up1_D4","up2_D4","up3_D4","dw1_D4")  # mediators (可替换成 d_up1... 等)

covs_cat <- c("Gender","CenterGroup","PneumoniaTypeGroup","Immunosuppression","MV")
covs_num <- c("Age","CCI","SOFA_24h")

stopifnot(all(c(X, Y, modules, covs_cat, covs_num) %in% names(df_med)))

# 结局变量确保是 0/1 numeric（如果本来就是 0/1 可省略）
df_med[[Y]] <- as.numeric(as.character(df_med[[Y]]))

# 避免 model.matrix 因 level 名字含 '-' 出问题
df_med$PneumoniaTypeGroup[df_med$PneumoniaTypeGroup == "non-CAP"] <- "non.CAP"

df_m <- df_med %>%
  mutate(across(all_of(covs_cat), as.factor))

# dummy（去掉截距列）
D <- model.matrix(
  ~ Gender + CenterGroup + PneumoniaTypeGroup + Immunosuppression + MV,
  data = df_m
)[, -1, drop = FALSE]

df_m <- bind_cols(
  df_m %>% dplyr::select(-all_of(covs_cat)),
  as.data.frame(D)
)

# scale 连续协变量（可选：也可以 scale X/M，但解释会变）
df_m[covs_num] <- lapply(df_m[covs_num], scale)

df_med <- df_m
covs_m  <- c(covs_num, colnames(D))
cov_str <- paste(covs_m, collapse = " + ")

# -------------------------
# 3) 一个小工具：安全提取系数（不存在就 NA）
# -------------------------
get_coef <- function(fit, term, col = "Estimate") {
  cf <- tryCatch(coef(summary(fit)), error = function(e) NULL)
  if (is.null(cf) || !(term %in% rownames(cf)) || !(col %in% colnames(cf))) return(NA_real_)
  unname(cf[term, col])
}

# -------------------------
# 4) 主循环：输出 a / b / c' / c + mediate 的 ACME/ADE/Total/Prop
# -------------------------
s_list <- list()

for (M in modules) {
  message("Running mediator: ", M)
  
  # (1) a: mediator model  M ~ X + cov
  med.fml <- reformulate(termlabels = c(X, covs_m), response = M)
  med.fit <- lm(med.fml, data = df_med)
  
  # (2) b, c': outcome model  Y ~ X + M + cov
  out.fml <- reformulate(termlabels = c(X, M, covs_m), response = Y)
  out.fit <- glm(out.fml, data = df_med, family = binomial())
  
  # (3) c: total effect model  Y ~ X + cov
  tot.fml <- reformulate(termlabels = c(X, covs_m), response = Y)
  tot.fit <- glm(tot.fml, data = df_med, family = binomial())
  
  # (4) mediation
  med.out <- mediate(
    model.m  = med.fit,
    model.y  = out.fit,
    treat    = X,
    mediator = M,
    boot     = TRUE,
    sims     = 2000
  )
  s <- summary(med.out)
  
  # ---- 提取路径系数 ----
  a_est <- get_coef(med.fit, term = X, col = "Estimate")
  a_se  <- get_coef(med.fit, term = X, col = "Std. Error")
  a_p   <- get_coef(med.fit, term = X, col = "Pr(>|t|)")
  
  b_est <- get_coef(out.fit, term = M, col = "Estimate")       # log-odds
  b_se  <- get_coef(out.fit, term = M, col = "Std. Error")
  b_p   <- get_coef(out.fit, term = M, col = "Pr(>|z|)")
  b_or  <- ifelse(is.na(b_est), NA_real_, exp(b_est))
  
  cprime_est <- get_coef(out.fit, term = X, col = "Estimate")  # log-odds
  cprime_se  <- get_coef(out.fit, term = X, col = "Std. Error")
  cprime_p   <- get_coef(out.fit, term = X, col = "Pr(>|z|)")
  cprime_or  <- ifelse(is.na(cprime_est), NA_real_, exp(cprime_est))
  
  c_est <- get_coef(tot.fit, term = X, col = "Estimate")       # log-odds
  c_se  <- get_coef(tot.fit, term = X, col = "Std. Error")
  c_p   <- get_coef(tot.fit, term = X, col = "Pr(>|z|)")
  c_or  <- ifelse(is.na(c_est), NA_real_, exp(c_est))
  
  # ---- mediate 的因果效应（ACME/ADE/Total/Prop）----
  s_df <- tibble(
    mediator = M,
    
    # a/b/c'/c（用于路径解释；注意 logistic 下 ACME 不是简单 a*b）
    a = a_est, a_se = a_se, a_p = a_p,
    b = b_est, b_se = b_se, b_OR = b_or, b_p = b_p,
    c_prime = cprime_est, c_prime_se = cprime_se, c_prime_OR = cprime_or, c_prime_p = cprime_p,
    c_total = c_est, c_total_se = c_se, c_total_OR = c_or, c_total_p = c_p,
    
    # 因果中介量（你之前用的）
    ACME = s$d.avg,
    ACME_low = s$d.avg.ci[1],
    ACME_high = s$d.avg.ci[2],
    ACME_p = s$d.avg.p,
    
    ADE = s$z.avg,
    ADE_low = s$z.avg.ci[1],
    ADE_high = s$z.avg.ci[2],
    ADE_p = s$z.avg.p,
    
    Total = s$tau.coef,
    Total_low = s$tau.ci[1],
    Total_high = s$tau.ci[2],
    Total_p = s$tau.p,
    
    PropMediated = s$n.avg,
    Prop_low = s$n.avg.ci[1],
    Prop_high = s$n.avg.ci[2],
    Prop_p = s$n.avg.p
  )
  
  s_list[[M]] <- s_df
}

res_med <- bind_rows(s_list)

# -------------------------
# 5) 可选：加上更好看的模块名 + 显著性标记
# -------------------------
module_map <- c(
  up1_D4 = "Inflammatory signaling",
  up2_D4 = "Phagolysosome function",
  up3_D4 = "IFN signaling",
  dw1_D4 = "Ribosome biogenesis"
)

res_med <- res_med %>%
  mutate(
    Module = factor(unname(module_map[mediator]),
                    levels = rev(unname(module_map[modules]))),
    sig = case_when(
      ACME_p < 0.001 ~ "***",
      ACME_p < 0.01  ~ "**",
      ACME_p < 0.05  ~ "*",
      TRUE           ~ ""
    ),
    ACME_ci = sprintf("%.3f (%.3f, %.3f)", ACME, ACME_low, ACME_high),
    a_ci    = sprintf("%.3f ± %.3f", a, a_se),
    b_ci    = sprintf("%.3f (OR=%.3f)", b, b_OR),
    cprime_ci = sprintf("%.3f (OR=%.3f)", c_prime, c_prime_OR),
    c_ci      = sprintf("%.3f (OR=%.3f)", c_total, c_total_OR)
  )


if(F){
  module_map <- c(
    up1_D4 = "Inflammatory signaling",
    up2_D4 = "Phagolysosome function",
    up3_D4 = "IFN-1 signaling",
    dw1_D4 = "Ribosome biogenesis"
  )
  
  res_med = read.csv(file.path(out_dir, "Supplementary_data_10a_med.csv"))
  
}

plot_df <- res_med %>%
  mutate(
    Module = factor(module_map[mediator],
                    levels = rev(module_map[modules[c(1,3,2,4)]])),
    sig = case_when(
      ACME_p < 0.001 ~ "***",
      ACME_p < 0.01  ~ "**",
      ACME_p < 0.05  ~ "*",
      TRUE           ~ ""
    ),
    p_label = ifelse(ACME_p < 0.001, "<0.001", sprintf("%.3f", ACME_p)),
    ci_label = sprintf("%.3f (%.3f, %.3f)", ACME, ACME_low, ACME_high)
  ) %>%
  mutate(p_label2 = paste0(p_label, sig))


p_acme <- ggplot(plot_df, aes(x = ACME, y = Module)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  
  geom_errorbarh(aes(xmin = ACME_low, xmax = ACME_high),
                 width = 0.18, linewidth = 0.6) +
  
  geom_point(aes(fill = ACME_p < 0.05),
             shape = 21, color = "black", size = 3.2) +
  
  # ⭐ p值标在点正上方
  geom_text(
    aes(label = paste0("p=", p_label2)),
    nudge_y = 0.35,      # 控制往上移动距离
    size = 3.2,
    hjust = 0.5
  ) +
  
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"),
                    guide = "none") +
  
  labs(x = "ACME (risk difference)", y = NULL) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  ) +
  scale_x_continuous(n.breaks = 3)

p_acme

pdf(file = file.path(out_dir, "med.pdf"), 
    width = 3.3, height = 2)
print(p_acme)
dev.off()

if(F){
  write.csv(res_med, 'Outputs/Supplementary_data_med.csv',
            row.names = F)
}

