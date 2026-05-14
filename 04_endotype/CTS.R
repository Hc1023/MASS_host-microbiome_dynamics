rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ConsensusTranscriptomicSubtype)
library(edgeR)

load('Inputs/1211_metadata.rdata')
load('Inputs/1211_transcriptome.rdata')
load('Inputs/1616_microbe.rdata')

# new_expr_data must be a data.frame or matrix
# First column = gene identifiers (optional), rownames = Ensembl IDs
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
# cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE)
cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE, log =T)

new_expr <- cpm_norm
result <- run_subtype_classifier(new_expr_data = new_expr)

## 获取每个分型的概率
rf <- result$rf_model
new_mat <- result$expression_corrected$new   # genes x samples
cts_prob <- predict(rf, t(new_mat), type = "prob")  

df_long$CTS = result[["predictions"]]$CTS
identical(df_long$SampleID, rownames(cts_prob))
df_long = cbind(df_long, cts_prob)
colnames(df_long)[6:8] = paste0('CTS', colnames(df_long)[6:8])
table(df_long$Mortality28d, df_long$CTS)

# D1/4/7 have samples or death
meta_analysis = meta %>% filter(!is.na(D1) & !is.na(D4.death) & !is.na(D7.death))
table(meta_analysis$Mortality28d)

meta_analysis %<>%
  dplyr::select(
    HumanID,
    Mortality28d,
    D1,
    D4 = D4.death,
    D7 = D7.death
  ) %>%
  pivot_longer(
    cols = c(D1, D4, D7),
    names_to = "Timepoint",
    values_to = "SampleID"
  )
tmp = df_long %>% dplyr::select(SampleID, CTS)

meta_analysis %<>% left_join(tmp, by = 'SampleID')

meta_analysis$CTS = factor(meta_analysis$CTS, 
                           levels = c("Death","Alive","1", "2", "3"))
meta_analysis$CTS[meta_analysis$SampleID == 'Death'] = 'Death'

d28_rows <- meta_analysis %>%
  distinct(HumanID, Mortality28d) %>%   # 保留每个人的结局信息
  mutate(
    Timepoint = "D28",
    SampleID  = NA_character_,
    CTS       = ifelse(Mortality28d == 1, 'Death', 'Alive')
  )

meta_analysis %<>%
  bind_rows(d28_rows) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D1","D4","D7","D28"))
  ) %>%
  arrange(HumanID, Timepoint)


df_plot = meta_analysis %>% dplyr::select(HumanID, Timepoint, State = CTS)
df_plot$State %<>% factor(., levels = c("Death","Alive","1", "2", "3"))


library(ggalluvial)

p = ggplot(df_plot,
       aes(x = Timepoint,
           stratum = State,
           alluvium = HumanID,
           fill = State,
           label = State)) +
  
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            alpha = .8, width = 1/12) +
  
  geom_stratum(width = 1/12, color = "black") +
  geom_text(stat = "stratum") +
  
  scale_x_discrete(expand = c(.1, .1)) +
  
  labs(title = "CTS Trajectory Alluvial Plot",
       x = "Timepoint",
       y = "Number of Patients",
       fill = "State") +
  
  theme_classic() +
  scale_fill_manual(values = c(
    "Death" = "#d73027",   
    "1"     = "#FC8D62",   
    "2"     = "#66C2A5",   
    "3"     = "#8DA0CB",   
    "Alive" = "#4575b4"    
  )) 


if(F){
  pdf(paste0("Outputs/04_CTS.pdf"), width = 5, height = 3)
  print(p)
  dev.off()
}


#### State mortality at each timepoint ####

df_mort_prop <- meta_analysis %>%
  filter(Timepoint %in% c("D1", "D4", "D7"),
         CTS %in% c("1", "2", "3")) %>%
  mutate(
    CTS = factor(CTS, levels = c("1", "2", "3")),
    Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7"))
  )

mort_summary <- df_mort_prop %>%
  group_by(Timepoint, CTS) %>%
  summarise(
    n_total = n(),
    n_dead  = sum(Mortality28d == 1),
    prop_dead = n_dead / n_total,
    .groups = "drop"
  )

mort_summary

# 转成 matrix：行 = Timepoint，列 = CTS
mat_prob <- mort_summary %>%
  dplyr::select(Timepoint, CTS, prop_dead) %>%
  pivot_wider(
    names_from  = CTS,
    values_from = prop_dead
  ) %>%
  column_to_rownames("Timepoint") %>%
  as.matrix()


library(ComplexHeatmap)
library(circlize)
library(grid)
library(scales)

col_fun <- colorRamp2(
  c(0, 0.3, 0.6),
  c("#F7FBFF", "#FDAE61", "#B2182B")
)


p = Heatmap(
  mat_prob,
  name = "P(death)",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title = "Timepoint",
  column_title = "CTS",
  
  row_names_side = "left",
  column_names_side = "top",
  
  rect_gp = gpar(col = "white", lwd = 1),
  
  # 在格子里写百分比
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(
      percent(mat_prob[i, j], accuracy = 0.1),
      x = x,
      y = y,
      gp = gpar(
        fontsize = 10,
        col = "black"
      )
    )
  },
  
  heatmap_legend_param = list(
    at = c(0, 0.3, 0.6),
    labels = percent(c(0, 0.3, 0.6)),
    title = "Mortality\nprobability"
  )
)

pdf(paste0("Outputs/04_CTS_mort.pdf"), width = 3, height = 2.2)
print(p)
dev.off()

#### State transition probability ####

df_next <- df_plot %>%
  group_by(HumanID) %>%
  mutate(
    Next_Timepoint = lead(Timepoint),
    Next_State     = lead(State)
  ) %>%
  ungroup()

df_trans <- df_next %>%
  filter(
    Timepoint %in% c("D1","D4","D7"),
    State != 'Death',
    !is.na(Next_State)
  )

library(broom)

prob_death_ci <- df_trans %>%
  group_by(Timepoint, State) %>%
  summarise(
    tidy(
      prop.test(
        x = sum(Next_State == "Death"),
        n = n()
      )
    ),
    .groups = "drop"
  ) %>%
  transmute(
    Timepoint,
    State,
    prob_death = estimate,
    conf.low,
    conf.high
  )

## visualization

# 确保顺序正确
prob_death_ci2 <- prob_death_ci %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D1", "D4", "D7")),
    State     = factor(State, levels = c("1", "2", "3"))
  )

# 转成 matrix：行 = Timepoint，列 = CTS
mat_prob <- prob_death_ci2 %>%
  dplyr::select(Timepoint, State, prob_death) %>%
  pivot_wider(
    names_from  = State,
    values_from = prob_death
  ) %>%
  column_to_rownames("Timepoint") %>%
  as.matrix()



col_fun <- colorRamp2(
  c(0, 0.3, 0.6),
  c("#F7FBFF", "#FDAE61", "#B2182B")
)

library(ComplexHeatmap)
library(circlize)
library(grid)
library(scales)

p1 = Heatmap(
  mat_prob,
  name = "P(death)",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_title = "Timepoint",
  column_title = "CTS",
  
  row_names_side = "left",
  column_names_side = "top",
  
  rect_gp = gpar(col = "white", lwd = 1),

  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(
      percent(mat_prob[i, j], accuracy = 0.1),
      x = x,
      y = y,
      gp = gpar(
        fontsize = 10,
        col = "black"
      )
    )
  },
  
  heatmap_legend_param = list(
    at = c(0, 0.3, 0.6),
    labels = percent(c(0, 0.3, 0.6)),
    title = "Next-timepoint\nmortality"
  )
)

p2 =ggplot(prob_death_ci,
       aes(x = Timepoint,
           y = prob_death,
           color = State,
           group = State)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.1,
    linewidth = 0.6
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    x = "Timepoint",
    y = "Next-timepoint mortality",
    color = "CTS"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values = c(
    "1"     = "#FC8D62",   
    "2"     = "#66C2A5",   
    "3"     = "#8DA0CB" ) 
  )

p2

pdf(paste0("Outputs/04_CTS_mort_next.pdf"), width = 3, height = 2)
print(p1)
print(p2)
dev.off()


#### CTS 1/2/3 probability longitudinal trend in mortality/survival ####
meta_analysis2 = meta_analysis %>% filter(SampleID != 'Death', !is.na(SampleID))
tmp = df_long %>% dplyr::select(SampleID, CTS1, CTS2, CTS3)
meta_analysis2 %<>% left_join(tmp, by = 'SampleID')
meta_analysis2 %<>% droplevels()
str(meta_analysis2$Timepoint)
head(meta_analysis2)

# 构造数值型时间变量（Joint model 必须）
meta_analysis2$TimeDays <- as.numeric(sub('D','',meta_analysis2$Timepoint))
# 存活时间
meta_analysis2 %<>% left_join(meta %>% 
                         dplyr::select(HumanID, 
                                       SurvivalTimeWithin28Days),
                       by = 'HumanID')

surv_df = meta_analysis2 %>% 
  dplyr::select(HumanID, Mortality28d, SurvivalTimeWithin28Days) %>%
  distinct() %>%
  mutate(
    Mortality28d = as.numeric(as.character(Mortality28d))  # 如果是factor: "0"/"1"
  )

## 纵向子模型 lme

cts = 'CTS1'

get_jm = function(cts){
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
  
  return(jm_fit)
}
jm_fit_cts1 = get_jm(cts = 'CTS1')
jm_fit_cts2 = get_jm(cts = 'CTS2')
jm_fit_cts3 = get_jm(cts = 'CTS3')


get_pred = function(cts, jm_fit, convergence){
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
  
  if(convergence == 1){
    p2 = ggplot() +
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
    return(p2)
  }
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
  return(p1)

}
cts= 'CTS1'
jm_fit = jm_fit_cts1

p1 = get_pred(cts= 'CTS1', jm_fit = jm_fit_cts1, convergence = 0)
# cts2 is not convergent, likely because its early death
p2 = get_pred(cts= 'CTS2', jm_fit = jm_fit_cts2, convergence = 1)
p3 = get_pred(cts= 'CTS3', jm_fit = jm_fit_cts3, convergence = 0)

pdf(paste0("Outputs/04_CTS_JM.pdf"), width = 3.5, height = 2.2)
print(p1)
print(p2)
print(p3)
dev.off()


if(F){
  library(openxlsx)
  library(dplyr)
  
  library(dplyr)
  library(openxlsx)
  
  extract_jm_one_sheet <- function(jm_fit, model_name) {
    txt <- capture.output(summary(jm_fit))
    
    ## helper to parse coefficient block
    parse_block <- function(lines, section_name) {
      lines <- trimws(lines)
      lines <- lines[lines != ""]
      lines <- lines[!grepl("^Value", lines)]
      
      out <- lapply(lines, function(x) {
        parts <- strsplit(x, "\\s+")[[1]]
        n <- length(parts)
        if (n < 5) return(NULL)
        data.frame(
          Model = model_name,
          Section = section_name,
          Term = paste(parts[1:(n-4)], collapse = " "),
          Estimate = parts[n-3],
          Std_Error = parts[n-2],
          z_value = parts[n-1],
          p_value = parts[n],
          stringsAsFactors = FALSE
        )
      })
      bind_rows(out)
    }
    
    ## longitudinal block
    i_long <- grep("^Longitudinal Process$", txt)
    i_event <- grep("^Event Process$", txt)
    long_tab <- parse_block(txt[(i_long + 1):(i_event - 1)], "Longitudinal")
    
    ## event block
    i_scale <- grep("^Scale:", txt)
    event_tab <- parse_block(txt[(i_event + 1):(i_scale - 1)], "Event")
    
    ## variance components
    i_var <- grep("^Variance Components:", txt)
    i_coef <- grep("^Coefficients:", txt)
    var_lines <- txt[(i_var + 1):(i_coef - 1)]
    var_lines <- trimws(var_lines)
    var_lines <- var_lines[var_lines != ""]
    var_lines <- var_lines[!grepl("^StdDev", var_lines)]
    
    var_tab <- lapply(var_lines, function(x) {
      parts <- strsplit(x, "\\s+")[[1]]
      if (length(parts) == 2) {
        data.frame(
          Model = model_name,
          Section = "Variance",
          Term = parts[1],
          Estimate = parts[2],
          Std_Error = NA,
          z_value = NA,
          p_value = NA,
          stringsAsFactors = FALSE
        )
      } else if (length(parts) >= 3) {
        data.frame(
          Model = model_name,
          Section = "Variance",
          Term = parts[1],
          Estimate = parts[2],
          Std_Error = parts[3],
          z_value = NA,
          p_value = NA,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }) %>% bind_rows()
    
    ## fit statistics
    i_fit <- grep("^   log.Lik", txt)
    fit_vals <- strsplit(trimws(txt[i_fit + 1]), "\\s+")[[1]]
    fit_tab <- data.frame(
      Model = model_name,
      Section = "Model fit",
      Term = c("logLik", "AIC", "BIC"),
      Estimate = fit_vals,
      Std_Error = NA,
      z_value = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    )
    
    bind_rows(long_tab, event_tab, var_tab, fit_tab)
  }
  
  cts1_sheet <- extract_jm_one_sheet(jm_fit_cts1, "CTS1")
  cts1_sheet[nrow(cts1_sheet)+1,1:3] = c('CTS1', 'Convergence', T)
  cts2_sheet <- extract_jm_one_sheet(jm_fit_cts2, "CTS2")
  cts2_sheet[nrow(cts2_sheet)+1,1:3] = c('CTS2', 'Convergence', F)
  cts3_sheet <- extract_jm_one_sheet(jm_fit_cts3, "CTS3")
  cts3_sheet[nrow(cts3_sheet)+1,1:3] = c('CTS3', 'Convergence', T)
  
  cts_sheets = bind_rows(cts1_sheet, cts2_sheet, cts3_sheet)
  cts_sheets[is.na(cts_sheets)] <- ""
  write.csv(cts_sheets, 'Outputs/Supplementary_data_cts.csv', row.names = F)
}


