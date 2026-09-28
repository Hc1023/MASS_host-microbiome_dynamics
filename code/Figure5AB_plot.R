source("code/shared_CTS_prepare.R")

library(ggalluvial)

p = ggplot(df_plot, aes(x = Timepoint, stratum = State, alluvium = HumanID, fill = State, label = State)) + 
    geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.8, width = 1/12) + geom_stratum(width = 1/12, 
    color = "black") + geom_text(stat = "stratum") + scale_x_discrete(expand = c(0.1, 0.1)) + labs(title = "CTS Trajectory Alluvial Plot", 
    x = "Timepoint", y = "Number of Patients", fill = "State") + theme_classic() + scale_fill_manual(values = c(Death = "#d73027", 
    `1` = "#FC8D62", `2` = "#66C2A5", `3` = "#8DA0CB", Alive = "#4575b4"))

if (TRUE) {
    pdf(repo_output(paste0("Outputs/1312_CTS.pdf"), "Figure5AB_plot.R"), width = 5, height = 3)
    print(p)
    dev.off()
}

df_mort_prop <- meta_analysis %>% filter(Timepoint %in% c("D1", "D4", "D7"), CTS %in% c("1", "2", "3")) %>% 
    mutate(CTS = factor(CTS, levels = c("1", "2", "3")), Timepoint = factor(Timepoint, levels = c("D1", 
        "D4", "D7")))

mort_summary <- df_mort_prop %>% group_by(Timepoint, CTS) %>% summarise(n_total = n(), n_dead = sum(Mortality28d == 
    1), prop_dead = n_dead/n_total, .groups = "drop")

mort_summary

mat_prob <- mort_summary %>% dplyr::select(Timepoint, CTS, prop_dead) %>% pivot_wider(names_from = CTS, 
    values_from = prop_dead) %>% column_to_rownames("Timepoint") %>% as.matrix()

col_fun <- circlize::colorRamp2(c(0, 0.3, 0.6), c("#F7FBFF", "#FDAE61", "#B2182B"))

library(ComplexHeatmap)

library(circlize)

library(grid)

library(scales)

p = Heatmap(mat_prob, name = "P(death)", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
    row_title = "Timepoint", column_title = "CTS", row_names_side = "left", column_names_side = "top", 
    rect_gp = gpar(col = "white", lwd = 1), cell_fun = function(j, i, x, y, w, h, fill) {
        grid.text(percent(mat_prob[i, j], accuracy = 0.1), x = x, y = y, gp = gpar(fontsize = 10, col = "black"))
    }, heatmap_legend_param = list(at = c(0, 0.3, 0.6), labels = percent(c(0, 0.3, 0.6)), title = "Mortality\nprobability"))

pdf(repo_output(paste0("Outputs/1312_CTS_mort.pdf"), "Figure5AB_plot.R"), width = 3, height = 2.2)

print(p)

dev.off()


write.csv(df_plot,"Outputs/Figure5A_states.csv",row.names=FALSE)
write.csv(mort_summary,"Outputs/Figure5B_mortality.csv",row.names=FALSE)
