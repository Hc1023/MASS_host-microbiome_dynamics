# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(magrittr)

library(data.table)

library(ggpubr)

library(gridExtra)

load(repo_input("Inputs/1616_microbe.rdata"))

load(repo_input("Inputs/1211_metadata.rdata"))

{
    genus_sum2 = data
    group_vec <- mapid_vec[rownames(genus_sum2)]
    prev <- rowSums(genus_sum2 > 0)
    bac_df <- genus_sum2[group_vec == "Bacteria", , drop = FALSE]
    fun_df <- genus_sum2[group_vec == "Fungi", , drop = FALSE]
    virus_df <- genus_sum2[group_vec == "Viruses", , drop = FALSE]
    virus_ordered <- virus_df[order(prev[rownames(virus_df)], decreasing = TRUE), ]
    bac_ordered <- bac_df[order(prev[rownames(bac_df)], decreasing = TRUE), ]
    virus_top10 <- head(virus_ordered, 10)
    bac_top10 <- head(bac_ordered, 10)
    genus_sum3 <- rbind(bac_top10, fun_df, virus_top10)
    genus_sum = genus_sum3
    pathogens = rownames(genus_sum3)
}

d = "D1"

myplot = function(d) {
    df_long1 = df_long[df_long$Timepoint == d, ]
    genus_sum3 = genus_sum[, df_long1$SampleID]
    genus_long <- genus_sum3 %>% tibble::rownames_to_column(var = "Pathogen") %>% pivot_longer(cols = -Pathogen, 
        names_to = "SampleID", values_to = "Abundance")
    genus_long %<>% filter(Abundance > 0)
    genus_long$Pathogen %<>% factor(., levels = pathogens)
    genus_long %<>% left_join(df_long %>% select(Mortality28d, SampleID), by = "SampleID")
    genus_long$Abundance = log2(genus_long$Abundance + 1)
    genus_long %<>% mutate(Abundance = ifelse(Mortality28d == "1", -Abundance, Abundance))
    results <- genus_long %>% group_by(Pathogen) %>% summarise(n0 = sum(Mortality28d == 0), n1 = sum(Mortality28d == 
        1), P_value = ifelse(length(unique(Mortality28d)) < 2 | min(n0, n1) < 3, NA, wilcox.test(abs(Abundance) ~ 
        Mortality28d, exact = FALSE)$p.value)) %>% as.data.frame()
    results <- results %>% mutate(sig_label = case_when(P_value < 0.001 ~ "***", P_value < 0.01 ~ "**", 
        P_value < 0.05 ~ "*", P_value < 0.1 ~ ".", TRUE ~ ""))
    max_y <- genus_long %>% group_by(Pathogen) %>% summarise(max_y = max(Abundance))
    annot_df <- max_y %>% left_join(results %>% select(Pathogen, sig_label), by = "Pathogen")
    annot_df$max_y = ifelse(annot_df$sig_label %in% c("*", "**", "***"), annot_df$max_y * 0.9, annot_df$max_y * 
        1.1)
    df_longplot = genus_long %>% filter(Pathogen %in% pathogens[1:11])
    p1 = ggplot(df_longplot, aes(x = Pathogen, y = Abundance, fill = Mortality28d)) + geom_boxplot(position = position_identity(), 
        width = 0.7, color = "black", linewidth = 0.3, outlier.alpha = 0.2) + scale_fill_manual(values = c(`0` = "#4575B4", 
        `1` = "#D73027")) + labs(x = "", y = "", fill = "") + theme_bw() + labs(title = d) + scale_y_continuous(labels = function(x) abs(x), 
        limits = c(-6, 6)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "None", 
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", 
            colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA), plot.margin = margin(l = 0.5, 
            r = 0.1, t = 0.1, unit = "cm")) + geom_text(data = annot_df[1:11, ], aes(x = Pathogen, y = max_y + 
        0.01, label = sig_label), inherit.aes = FALSE, vjust = 0, size = 5)
    df_longplot = genus_long %>% filter(Pathogen %in% pathogens[12:21])
    levels(df_longplot$Pathogen)[levels(df_longplot$Pathogen) == "HHV-4"] = "EBV"
    levels(annot_df$Pathogen)[levels(annot_df$Pathogen) == "HHV-4"] = "EBV"
    p2 = ggplot(df_longplot, aes(x = Pathogen, y = Abundance, fill = Mortality28d)) + geom_boxplot(position = position_identity(), 
        width = 0.7, color = "black", linewidth = 0.3, outlier.alpha = 0.2) + scale_fill_manual(values = c(`0` = "#4575B4", 
        `1` = "#D73027")) + labs(x = "", y = "", fill = "") + theme_bw() + labs(title = d) + scale_y_continuous(labels = function(x) abs(x), 
        limits = c(-15, 15)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "None", 
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA), legend.background = element_rect(fill = "transparent", 
            colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA), plot.margin = margin(l = 0.5, 
            r = 0.1, t = 0.1, unit = "cm")) + geom_text(data = annot_df[12:21, ], aes(x = Pathogen, y = max_y + 
        0.01, label = sig_label), inherit.aes = FALSE, vjust = 0, size = 5)
    return(list(p1, p2))
}

plist1 = myplot(d = "D1")

plist2 = myplot(d = "D4")

plist3 = myplot(d = "D7")

p1 = plist1[[1]]

p2 = plist1[[2]]

p3 = plist2[[1]]

p4 = plist2[[2]]

p5 = plist3[[1]]

p6 = plist3[[2]]

pdf(repo_output(paste0("Outputs/03_microbe_mortality_association mass.pdf"), "FigureS4B_analysis.R"), 
    width = 2.6, height = 2.4)

print(p1)

print(p2)

print(p3)

print(p4)

print(p5)

print(p6)

dev.off()

