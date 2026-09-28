source("code/shared_CTS_prepare.R")
library(JM)

meta_analysis2 = meta_analysis %>% filter(SampleID != "Death", !is.na(SampleID))

tmp = df_long %>% dplyr::select(SampleID, CTS1, CTS2, CTS3)

meta_analysis2 %<>% left_join(tmp, by = "SampleID")

meta_analysis2 %<>% droplevels()

str(meta_analysis2$Timepoint)

write.csv(meta_analysis2,"Outputs/Figure5C_model_rows.csv",row.names=FALSE)

meta_analysis2$TimeDays <- as.numeric(sub("D", "", meta_analysis2$Timepoint)) - 1

meta_analysis2 %<>% left_join(meta %>% dplyr::select(HumanID, SurvivalTimeWithin28Days), by = "HumanID")

surv_df = meta_analysis2 %>% dplyr::select(HumanID, Mortality28d, SurvivalTimeWithin28Days) %>% distinct() %>% 
    mutate(Mortality28d = as.numeric(as.character(Mortality28d)))

cts = "CTS1"

get_jm = function(cts) {
    library(nlme)
    lme_fit <- lme(fixed = as.formula(paste0(cts, "~ TimeDays * Mortality28d")), random = ~TimeDays | 
        HumanID, data = meta_analysis2, na.action = na.exclude, control = lmeControl(opt = "optim"))
    summary(lme_fit)
    library(survival)
    cox_fit <- coxph(Surv(SurvivalTimeWithin28Days, Mortality28d) ~ 1 + cluster(HumanID), data = surv_df, 
        x = TRUE)
    library(JM)
    # Same model and quadrature; CTS1 needed a larger optimizer iteration limit.
    jm_fit <- jointModel(lmeObject = lme_fit, survObject = cox_fit, timeVar = "TimeDays", method = "weibull-PH-aGH",
        control = list(iter.qN = if (cts == "CTS1") 2000L else 350L))
    return(jm_fit)
}

if ("--refit" %in% commandArgs(TRUE)) {
    jm_fit_cts1 <- get_jm("CTS1")
    jm_fit_cts3 <- get_jm("CTS3")
    saveRDS(list(CTS1=jm_fit_cts1,CTS3=jm_fit_cts3),"Outputs/Figure5C_joint_models.rds")
} else {
    saved_models <- readRDS("Outputs/Figure5C_joint_models.rds")
    jm_fit_cts1 <- saved_models$CTS1
    jm_fit_cts3 <- saved_models$CTS3
}
stopifnot(jm_fit_cts1$convergence == 0, jm_fit_cts3$convergence == 0)

get_pred = function(cts, jm_fit, convergence) {
    pred_grid <- expand.grid(TimeDays = c(0, 3, 6), Mortality28d = c(0, 1))
    pred_grid$HumanID <- meta_analysis2$HumanID[1]
    pred_grid$Group <- ifelse(pred_grid$Mortality28d == 0, "Survivor", "Non-survivor")
    pred_grid$fit_jm <- as.numeric(predict(jm_fit, newdata = pred_grid, process = "Longitudinal", type = "Marginal"))
    pred_grid$Timepoint = factor(paste0("D", pred_grid$TimeDays + 1))
    str(pred_grid)
    library(ggplot2)
    dodge_w <- 0.35
    library(ggnewscale)
    library(ggpubr)
    library(scales)
    cols_strong <- c(`0` = "#4575B4", `1` = "#D73027")
    cols_light <- c(`0` = alpha("#4575B4", 0.7), `1` = alpha("#D73027", 0.7))
    if (convergence == 1) {
        p2 = ggplot() + geom_boxplot(data = meta_analysis2, aes(x = Timepoint, y = .data[[cts]], group = interaction(Timepoint, 
            Mortality28d), color = factor(Mortality28d)), width = 0.25, outlier.shape = NA, position = position_dodge(width = dodge_w), 
            alpha = 0.1) + scale_color_manual(values = cols_light, labels = c(`0` = "Survivor", `1` = "Mortality"), 
            name = "Observed") + ggpubr::stat_compare_means(data = meta_analysis2, aes(x = Timepoint, 
            y = .data[[cts]], group = Mortality28d), method = "wilcox.test", label = "p.signif", hide.ns = TRUE, 
            tip.length = 0.01, size = 5, alpha = 0.7) + scale_y_continuous(limits = c(0, 1), breaks = c(0, 
            0.5, 1), expand = expansion(mult = c(0.05, 0.12))) + theme_bw(base_size = 12) + labs(x = "Timepoint", 
            y = cts, color = "Outcome")
        return(p2)
    }
    p1 = ggplot() + geom_boxplot(data = meta_analysis2, aes(x = Timepoint, y = .data[[cts]], group = interaction(Timepoint, 
        Mortality28d), color = factor(Mortality28d)), width = 0.25, outlier.shape = NA, position = position_dodge(width = dodge_w), 
        alpha = 0.1) + scale_color_manual(values = cols_light, labels = c(`0` = "Survivor", `1` = "Mortality"), 
        name = "Observed") + ggnewscale::new_scale_color() + geom_point(data = pred_grid, aes(x = Timepoint, 
        y = fit_jm, group = factor(Mortality28d), color = factor(Mortality28d)), position = position_dodge(width = dodge_w), 
        size = 3) + geom_line(data = pred_grid, aes(x = Timepoint, y = fit_jm, group = factor(Mortality28d), 
        color = factor(Mortality28d)), position = position_dodge(width = dodge_w), linewidth = 1) + scale_color_manual(values = cols_strong, 
        labels = c(`0` = "Survival", `1` = "Mortality"), name = "JM fit") + ggpubr::stat_compare_means(data = meta_analysis2, 
        aes(x = Timepoint, y = .data[[cts]], group = Mortality28d), method = "wilcox.test", label = "p.signif", 
        hide.ns = TRUE, tip.length = 0.01, size = 5, alpha = 0.7) + scale_y_continuous(limits = c(0, 
        1), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.05, 0.12))) + theme_bw(base_size = 12) + 
        labs(x = "Timepoint", y = cts, color = "Outcome")
    return(p1)
}

cts = "CTS1"

jm_fit = jm_fit_cts1

p1 = get_pred(cts = "CTS1", jm_fit = jm_fit_cts1, convergence = 0)


p3 = get_pred(cts = "CTS3", jm_fit = jm_fit_cts3, convergence = 0)

pdf(repo_output(paste0("Outputs/04_CTS_JM.pdf"), "Figure5C_analysis.R"), width = 3.5, height = 2.2)

print(p1)


print(p3)

dev.off()
