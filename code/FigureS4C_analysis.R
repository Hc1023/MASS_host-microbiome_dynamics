# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(data.table)

library(maaslin3)

library(ggplot2)

library(magrittr)

load(repo_input("Inputs/1616_microbe.rdata"))

load(repo_input("Inputs/1211_metadata.rdata"))

taxa_table = data.frame(t(data[pathogens, ]))

taxa_table = log2(taxa_table + 1)

colnames(taxa_table) = pathogens

d = 1

df_long1 = df_long %>% filter(Timepoint == paste0("D", d)) %>% left_join(meta[, -2], by = "HumanID")

df_long1$PneumoniaType = ifelse(df_long1$PneumoniaType == "1.CAP", "CAP", "Others") %>% as.factor()

df_long1 %<>% mutate(Center = ifelse(Center %in% names(which(table(Center) < 10)), "Others", Center))

vars = c("Gender", "Age", "CenterGroup", "PneumoniaTypeGroup", "CCI", "SOFA_24h", "Immunosuppression", 
    "MV")

var_str = paste(vars, collapse = "+")

df_long1$Gender %<>% as.factor()

rownames(df_long1) = df_long1$SampleID

taxa1 = taxa_table[df_long1$SampleID, ]

fit_out <- maaslin3(input_data = taxa1, input_metadata = df_long1, output = "Outputs", formula = paste("~ Mortality28d +", 
    var_str), normalization = "NONE", transform = "NONE", augment = TRUE, standardize = F, max_significance = 0.1, 
    median_comparison_abundance = F, median_comparison_prevalence = FALSE, cores = 1, warn_prevalence = F)

all_mas = read.delim2("Outputs/03_mas_output/all_results.tsv")

all_mas %<>% filter(metadata == "Mortality28d")

all_mas$coef %<>% as.numeric()

all_mas$null_hypothesis %<>% as.numeric()

all_mas$xmin = all_mas$coef - as.numeric(all_mas$stderr)

all_mas$xmax = all_mas$coef + as.numeric(all_mas$stderr)

all_mas$pval_individual %<>% as.numeric()

head(all_mas)

median_df = data.frame(model = c("Prevalence", "Abundance"), null = c(all_mas$null_hypothesis[all_mas$model == 
    "prevalence"][1], all_mas$null_hypothesis[all_mas$model == "abundance"][1]))

all_mas$feature %<>% factor(., levels = rev(pathogens))

mytheme = theme_bw() + theme(axis.text = element_text(color = "black", size = 13), axis.ticks.length = unit(0.2, 
    "lines"), axis.ticks = element_line(size = 0.5, color = "black"), title = element_text(color = "black", 
    size = 13), panel.border = element_rect(linewidth = 1), legend.text = element_text(size = 13), legend.title = element_text(size = 13), 
    legend.position = "right", plot.margin = margin(l = 0.1, r = 0.1, t = 1, b = 0.1, unit = "cm"), legend.spacing.y = unit(2, 
        "pt"))

levels(all_mas$feature)[levels(all_mas$feature) == "HHV-4"] = "EBV"

pp = ggplot(all_mas, aes(x = coef, y = feature)) + ggplot2::guides(linetype = ggplot2::guide_legend(title = "Null hypothesis", 
    order = 1), ) + ggplot2::geom_vline(data = median_df, ggplot2::aes(xintercept = null, linetype = model), 
    color = "darkgray", size = 0.3) + ggplot2::scale_linetype_manual(values = c(Prevalence = "dashed", 
    Abundance = "solid")) + geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0.2, linewidth = 0.4) + 
    geom_point(aes(shape = model), size = 3.2) + ggplot2::scale_shape_manual(name = "Association", values = c(21, 
    24)) + ggplot2::guides(shape = ggplot2::guide_legend(order = 2), ) + ggplot2::labs(x = expression(paste(beta, 
    " coefficient")), y = "") + geom_point(all_mas[all_mas$model == "abundance", ], mapping = aes(shape = model, 
    fill = pval_individual), size = 3.2, shape = 21, stroke = 0.3) + ggplot2::scale_fill_gradient(low = "#8B008B", 
    high = "white", limits = c(1e-05, 1), breaks = c(1e-05, 0.05, 1), labels = c(1e-05, 0.05, 1), transform = scales::pseudo_log_trans(sigma = 0.001), 
    name = "Abundance P") + ggnewscale::new_scale_fill() + geom_point(all_mas[all_mas$model == "prevalence", 
    ], mapping = aes(shape = model, fill = pval_individual), size = 3.2, shape = 24, stroke = 0.3) + 
    ggplot2::scale_fill_gradient(low = "#008B8B", high = "white", limits = c(1e-05, 1), breaks = c(1e-05, 
        0.05, 1), labels = c(1e-05, 0.05, 1), transform = scales::pseudo_log_trans(sigma = 0.001), name = "Prevalence P") + 
    mytheme + theme(strip.text = element_text(size = 13), panel.grid = element_line(linewidth = 0.3))

pp

pdf(file = repo_output("Outputs/03_maaslin.pdf", "FigureS4C_analysis.R"), width = 6, 
    height = 5.5)

print(pp)

dev.off()

tmp = all_mas %>% filter(pval_individual < 0.05)

