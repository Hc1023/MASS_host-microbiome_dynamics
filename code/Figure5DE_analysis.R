# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(magrittr)

library(data.table)

library(purrr)

library(ggplot2)

library(stringr)

library(forcats)

library(grid)

library(survival)

project_dir <- "."

invisible(NULL)

out_dir <- repo_path(".")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

load("Inputs/1616_meta_model.rdata")

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

df_sum <- df1 %>% filter(!is.na(grp), !is.na(Mortality28d)) %>% group_by(grp) %>% summarise(n = n(), 
    death_rate = mean(Mortality28d == 1), .groups = "drop")

df_sum <- df_sum %>% mutate(death_pct = death_rate * 100, n_label = paste0("n = ", n))

df_sum$grp = factor(df_sum$grp, levels = c("Stable high-risk", "Stable low-risk", "Deteriorating", "Improving", 
    "Mixed"))

p = ggplot(df_sum, aes(x = grp, y = death_pct)) + geom_col(fill = "#D55E00", color = "black", width = 0.6, 
    alpha = 0.85) + geom_text(aes(label = n_label), hjust = -0.1, size = 3, color = "black") + scale_y_continuous(name = "Mortality (%)", 
    limits = c(0, max(df_sum$death_pct) * 1.3), expand = expansion(mult = c(0, 0.05))) + theme_bw() + 
    theme(axis.title.y = element_text(size = 11), panel.grid.minor = element_blank()) + labs(x = NULL) + 
    coord_flip()

tmp = meta_model %>% dplyr::select(HumanID, matches("D[147]_mi"))

df1_mi = df1 %>% left_join(tmp, by = "HumanID")

delta_long <- df1_mi %>% pivot_longer(cols = matches("^D(1|4|7)_(mi2?|mi)_"), names_to = c("Day", "feature"), 
    names_pattern = "^D(1|4|7)_(.+)$", values_to = "value") %>% mutate(Day = paste0("D", Day)) %>% pivot_wider(names_from = Day, 
    values_from = value) %>% mutate(delta = D7 - D1)

delta_long$feature <- gsub("^(mi2|mi)_", "", delta_long$feature)

delta_long %<>% filter(grp != "Mixed")

library(rstatix)

res_delta <- delta_long %>% filter(!is.na(D1), !is.na(D7)) %>% group_by(grp, feature) %>% summarise(n = n(), 
    median_delta = median(delta, na.rm = TRUE), mean_delta = mean(delta, na.rm = TRUE), p_value = tryCatch(wilcox.test(D7, 
        D1, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_), .groups = "drop") %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH"))

mean_delta_mat <- res_delta %>% dplyr::select(feature, grp, mean_delta) %>% tidyr::pivot_wider(names_from = grp, 
    values_from = mean_delta)

pval_mat <- res_delta %>% dplyr::select(feature, grp, p_value) %>% tidyr::pivot_wider(names_from = grp, 
    values_from = p_value)

pval_mat_2 = pval_mat %>% filter(!is.na(Deteriorating)) %>% filter(Deteriorating < 1)

pval_mat_2 %<>% arrange(Deteriorating)

bubble_df <- pval_mat_2 %>% pivot_longer(cols = -feature, names_to = "grp", values_to = "p") %>% left_join(mean_delta_mat %>% 
    pivot_longer(cols = -feature, names_to = "grp", values_to = "delta"), by = c("feature", "grp")) %>% 
    mutate(neg_log10_p = -log10(p))

bubble_df$feature = factor(bubble_df$feature, levels = rev(pval_mat_2$feature))

bubble_df <- bubble_df %>% mutate(sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", 
    p < 0.1 ~ "<U+00B7>", TRUE ~ ""))

levels(bubble_df$feature)[levels(bubble_df$feature) == "HHV.4"] <- "EBV"

levels(bubble_df$feature)[levels(bubble_df$feature) == "Influenza.A"] <- "Influenza A"

levels(bubble_df$feature)[levels(bubble_df$feature) == "Bacf"] <- "Bacteria/Fungi"

feature_x = "Klebsiella"

feature_x = "Total"

delta_long$feature[delta_long$feature == "HHV.4"] = "EBV"

delta_long$feature = gsub("\\.", " ", delta_long$feature)

delta_long$feature[delta_long$feature == "SARS CoV 2"] = "SARS-CoV-2"

delta_long$feature[delta_long$feature == "HAdV B"] = "HAdV-B"

delta_long$feature[delta_long$feature == "HHV 6B"] = "HHV-6B"

delta_long$feature[delta_long$feature == "Bacf"] = "Bacteria/Fungi"

draw_feature_x = function(feature_x) {
    group_cols <- c(Deteriorating = "#8DA0CB", Improving = "#D95F02")
    dfx = delta_long %>% filter(feature == feature_x, grp %in% names(group_cols)) %>% mutate(grp = factor(grp, 
        levels = names(group_cols)))
    df_long <- dfx %>% pivot_longer(cols = c(D1, D4, D7), names_to = "Day", values_to = "value") %>% 
        mutate(Day = factor(Day, levels = c("D1", "D4", "D7")), Day_num = as.numeric(Day))
    stat_line <- df_long %>% group_by(Day, grp) %>% summarise(m = mean(value, na.rm = TRUE), se = sd(value, 
        na.rm = TRUE)/sqrt(sum(!is.na(value))), .groups = "drop") %>% mutate(lower = m - se, upper = m + 
        se)
    y_range <- range(c(stat_line$lower, stat_line$upper), na.rm = TRUE)
    y_span <- diff(y_range)
    if (!is.finite(y_span) || y_span == 0) 
        y_span <- 1
    stat_grp <- dfx %>% filter(!is.na(D1), !is.na(D7)) %>% group_by(grp) %>% summarise(n = n(), log2FC = mean(D7 - 
        D1, na.rm = TRUE), p = tryCatch(wilcox.test(D7, D1, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_), 
        .groups = "drop") %>% mutate(sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", 
        p < 0.1 ~ "<U+00B7>", TRUE ~ "ns"), label = paste0("Change = ", sprintf("%.2f", log2FC), " ", 
        sig), x.position = 2, y.position = case_when(grp == "Deteriorating" ~ y_range[2] - y_span * 0.15, 
        grp == "Improving" ~ y_range[1] + y_span * 0.05, TRUE ~ NA_real_))
    p = ggplot(df_long, aes(x = Day_num, y = value, color = grp, group = grp)) + stat_summary(fun.data = mean_se, 
        geom = "errorbar", width = 0.15, linewidth = 0.6) + stat_summary(fun = mean, geom = "line", linewidth = 0.8) + 
        stat_summary(fun = mean, geom = "point", size = 2.5) + theme_bw() + labs(x = NULL, y = "log2(mass+1)", 
        color = "Trajectory group") + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
        ggtitle(feature_x) + geom_text(data = stat_grp, aes(x = x.position, y = y.position, label = label, 
        color = grp), inherit.aes = FALSE, size = 3.2, show.legend = FALSE) + scale_color_manual(values = group_cols) + 
        scale_x_continuous(breaks = 1:3, labels = c("D1", "D4", "D7"), limits = c(0.9, 3.1)) + scale_y_continuous(expand = expansion(mult = c(0.05, 
        0.1)))
    p
    return(p)
}

features_to_plot <- c("Total", "Bacteria/Fungi", "Viruses")

plots <- setNames(map(features_to_plot, draw_feature_x), features_to_plot)

draw_feature_x("Total")

draw_feature_x("Bacteria/Fungi")

draw_feature_x("Viruses")

group_cols <- c(Deteriorating = "#8DA0CB", Improving = "#D95F02")

legend_df <- tibble(Day = rep(c(1, 2), 2), value = rep(c(1, 1), 2), grp = factor(rep(names(group_cols), 
    each = 2), levels = names(group_cols)))

legend_plot <- ggplot(legend_df, aes(x = Day, y = value, color = grp, group = grp)) + geom_line(linewidth = 0.8, 
    show.legend = TRUE) + geom_point(size = 2.5, show.legend = TRUE) + scale_color_manual(values = group_cols, 
    labels = c("Deteriorating (CTS3 -> 1/2)", "Improving (CTS1/2 -> 3)")) + labs(color = NULL) + theme_void() + 
    theme(legend.position = "top", legend.text = element_text(size = 9), legend.key.width = unit(1.1, 
        "cm")) + guides(color = guide_legend(override.aes = list(linewidth = 0.8, size = 2.5)))

legend_grob <- ggplotGrob(legend_plot)

legend_grob <- cowplot::get_legend(legend_plot)

ggsave(filename = repo_output(repo_path(out_dir, "trajectory_two_groups_legend.pdf"), "Figure5DE_analysis.R"), 
    plot = legend_grob, width = 4.8, height = 0.35)

walk2(plots, names(plots), ~ggsave(filename = repo_output(repo_path(out_dir, paste0("trajectory_two_groups_", 
    str_replace_all(.y, "[/ ]+", "_"), ".pdf")), "Figure5DE_analysis.R"), plot = .x, width = 2.3, 
    height = 2.5))

surv_df <- df1 %>% filter(grp %in% c("Deteriorating", "Improving")) %>% dplyr::select(HumanID, grp) %>% 
    left_join(meta_model %>% dplyr::select(HumanID, Mortality28d, SurvivalTimeWithin28Days), by = "HumanID") %>% 
    mutate(grp = factor(grp, levels = names(group_cols)), event = as.integer(as.character(Mortality28d) == 
        "1"), time = as.numeric(SurvivalTimeWithin28Days)) %>% filter(!is.na(grp), !is.na(event), !is.na(time))

surv_fit <- survfit(Surv(time, event) ~ grp, data = surv_df)

surv_diff <- survdiff(Surv(time, event) ~ grp, data = surv_df)

logrank_p <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)

surv_sum <- summary(surv_fit)

surv_plot_df <- tibble(time = surv_sum$time, surv = surv_sum$surv, n.risk = surv_sum$n.risk, n.event = surv_sum$n.event, 
    n.censor = surv_sum$n.censor, grp = str_remove(surv_sum$strata, "^grp=")) %>% mutate(grp = factor(grp, 
    levels = names(group_cols))) %>% bind_rows(tibble(time = 0, surv = 1, n.risk = as.integer(table(surv_df$grp)[names(group_cols)]), 
    n.event = 0, n.censor = 0, grp = factor(names(group_cols), levels = names(group_cols)))) %>% arrange(grp, 
    time)

censor_df <- surv_plot_df %>% filter(n.censor > 0)

surv_stats <- surv_df %>% group_by(grp) %>% summarise(n = n(), events = sum(event == 1), event_rate = events/n, 
    .groups = "drop") %>% mutate(logrank_p = logrank_p)

write.csv(surv_stats, file = repo_output(repo_path(out_dir, "survival_two_groups_summary.csv"), "Figure5DE_analysis.R"), 
    row.names = FALSE)

surv_labels <- c(Deteriorating = paste0("Deteriorating (CTS3 -> 1/2), n = ", surv_stats$n[match("Deteriorating", 
    surv_stats$grp)]), Improving = paste0("Improving (CTS1/2 -> 3), n = ", surv_stats$n[match("Improving", 
    surv_stats$grp)]))

p_surv <- ggplot(surv_plot_df, aes(x = time, y = surv, color = grp, group = grp)) + geom_step(linewidth = 0.8) + 
    geom_point(data = censor_df, shape = 3, size = 1.8, stroke = 0.6, show.legend = FALSE) + annotate("text", 
    x = 1, y = 0.12, hjust = 0, size = 3.2, label = paste0("Log-rank p = ", signif(logrank_p, 3))) + 
    scale_color_manual(values = group_cols, labels = surv_labels) + scale_x_continuous(name = "Days after enrollment", 
    limits = c(0, 28), breaks = c(0, 7, 14, 21, 28)) + scale_y_continuous(name = "Survival probability", 
    limits = c(0, 1), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 
        0.04))) + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(), 
    legend.text = element_text(size = 8), legend.margin = margin(0, 0, 0, 0)) + guides(color = guide_legend(nrow = 1, 
    byrow = TRUE, override.aes = list(linewidth = 0.8)))

ggsave(filename = repo_output(repo_path(out_dir, "survival_two_groups_KM.pdf"), "Figure5DE_analysis.R"), 
    plot = p_surv, width = 5.8, height = 2.6)

ggsave(filename = repo_output(repo_path(out_dir, "survival_two_groups_KM2.pdf"), "Figure5DE_analysis.R"), 
    plot = p_surv, width = 5, height = 2)

d1_surv_df <- meta_model %>% dplyr::select(HumanID, Mortality28d, SurvivalTimeWithin28Days, D1_CTS) %>% 
    mutate(grp = case_when(D1_CTS %in% c(1, 2) ~ "D1 CTS1/2", D1_CTS == 3 ~ "D1 CTS3", TRUE ~ NA_character_), 
        grp = factor(grp, levels = c("D1 CTS1/2", "D1 CTS3")), event = as.integer(as.character(Mortality28d) == 
            "1"), time = as.numeric(SurvivalTimeWithin28Days)) %>% filter(!is.na(grp), !is.na(event), 
    !is.na(time))

d1_cols <- c(`D1 CTS1/2` = "#D95F02", `D1 CTS3` = "#8DA0CB")

d1_surv_fit <- survfit(Surv(time, event) ~ grp, data = d1_surv_df)

d1_surv_diff <- survdiff(Surv(time, event) ~ grp, data = d1_surv_df)

d1_logrank_p <- pchisq(d1_surv_diff$chisq, df = length(d1_surv_diff$n) - 1, lower.tail = FALSE)

d1_surv_sum <- summary(d1_surv_fit)

d1_surv_plot_df <- tibble(time = d1_surv_sum$time, surv = d1_surv_sum$surv, n.risk = d1_surv_sum$n.risk, 
    n.event = d1_surv_sum$n.event, n.censor = d1_surv_sum$n.censor, grp = str_remove(d1_surv_sum$strata, 
        "^grp=")) %>% mutate(grp = factor(grp, levels = names(d1_cols))) %>% bind_rows(tibble(time = 0, 
    surv = 1, n.risk = as.integer(table(d1_surv_df$grp)[names(d1_cols)]), n.event = 0, n.censor = 0, 
    grp = factor(names(d1_cols), levels = names(d1_cols)))) %>% arrange(grp, time)

d1_censor_df <- d1_surv_plot_df %>% filter(n.censor > 0)

d1_surv_stats <- d1_surv_df %>% group_by(grp) %>% summarise(n = n(), events = sum(event == 1), event_rate = events/n, 
    events_within_7d = sum(event == 1 & time <= 7), .groups = "drop") %>% mutate(logrank_p = d1_logrank_p)

write.csv(d1_surv_stats, file = repo_output(repo_path(out_dir, "survival_D1_CTS_summary.csv"), "Figure5DE_analysis.R"), 
    row.names = FALSE)

d1_surv_labels <- c(`D1 CTS1/2` = paste0("D1 CTS1/2, n = ", d1_surv_stats$n[match("D1 CTS1/2", d1_surv_stats$grp)]), 
    `D1 CTS3` = paste0("D1 CTS3, n = ", d1_surv_stats$n[match("D1 CTS3", d1_surv_stats$grp)]))

p_d1_surv <- ggplot(d1_surv_plot_df, aes(x = time, y = surv, color = grp, group = grp)) + geom_step(linewidth = 0.8) + 
    geom_point(data = d1_censor_df, shape = 3, size = 1.8, stroke = 0.6, show.legend = FALSE) + annotate("text", 
    x = 1, y = 0.12, hjust = 0, size = 3.2, label = paste0("Log-rank p = ", signif(d1_logrank_p, 3))) + 
    scale_color_manual(values = d1_cols, labels = d1_surv_labels) + scale_x_continuous(name = "Days after enrollment", 
    limits = c(0, 28), breaks = c(0, 7, 14, 21, 28)) + scale_y_continuous(name = "Survival probability", 
    limits = c(0, 1), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 
        0.04))) + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(), 
    legend.text = element_text(size = 8), legend.margin = margin(0, 0, 0, 0)) + guides(color = guide_legend(nrow = 1, 
    byrow = TRUE, override.aes = list(linewidth = 0.8)))

ggsave(filename = repo_output(repo_path(out_dir, "survival_D1_CTS_KM.pdf"), "Figure5DE_analysis.R"), 
    plot = p_d1_surv, width = 5.8, height = 2.6)

write.csv(d1_surv_plot_df,"Outputs/Figure5D_baseline_KM.csv",row.names=FALSE)
write.csv(surv_plot_df,"Outputs/Figure5D_transition_KM.csv",row.names=FALSE)
for (nm in c("d1_surv_fit","surv_fit")) {
    ss <- summary(get(nm),times=c(0,7,14,21,28),extend=TRUE)
    write.csv(data.frame(Time=ss$time,Group=ss$strata,N_risk=ss$n.risk,Survival=ss$surv),
              paste0("Outputs/Figure5D_",nm,"_risk.csv"),row.names=FALSE)
}
trajectory_data <- delta_long %>% filter(feature %in% features_to_plot,
    grp %in% c("Deteriorating","Improving")) %>%
    dplyr::select(HumanID,grp,feature,D1,D4,D7) %>%
    pivot_longer(c(D1,D4,D7),names_to="Timepoint",values_to="Log2_biomass_plus1")
write.csv(trajectory_data,"Outputs/Figure5E_observations.csv",row.names=FALSE)
trajectory_summary <- trajectory_data %>% group_by(grp,feature,Timepoint) %>%
    summarise(N=sum(!is.na(Log2_biomass_plus1)),Mean=mean(Log2_biomass_plus1,na.rm=TRUE),
      SE=sd(Log2_biomass_plus1,na.rm=TRUE)/sqrt(N),.groups="drop")
write.csv(trajectory_summary,"Outputs/Figure5E_summary.csv",row.names=FALSE)
trajectory_tests <- delta_long %>% filter(feature %in% features_to_plot,
    grp %in% c("Deteriorating","Improving")) %>% group_by(grp,feature) %>%
    summarise(N=sum(complete.cases(D1,D7)),Mean_change=mean(D7-D1,na.rm=TRUE),
      P_value=wilcox.test(D7,D1,paired=TRUE,exact=FALSE)$p.value,.groups="drop")
write.csv(trajectory_tests,"Outputs/Figure5E_tests.csv",row.names=FALSE)
