# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(ComplexHeatmap))

suppressPackageStartupMessages(library(circlize))

suppressPackageStartupMessages(library(grid))

project_dir <- "."

invisible(NULL)

out_dir <- repo_path(".")

repo_dir(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- "Inputs/metadata_table_source_1007.csv"

min_positive_n <- 10

exclude_microbes <- c("Lachnospiraceae", "Aspergillus")

load(repo_input("Inputs/1211_metadata.rdata"))

load(repo_input("Inputs/1616_microbe.rdata"))

metadata <- read_csv(repo_input(metadata_file), show_col_types = FALSE)

clinical_history_vars <- c("Mortality28d", "Age", "SOFA_24h", "APACHEII_24h", "WBC", "LymphocyteCount", 
    "NeutrophilCount", "PlateletCount", "PCT", "hsCRP", "PaO2_FiO2", "Immunosuppression", "MV", "DM", 
    "MI", "COPD", "HepaticImpairment", "RenalDisease", "Tumor", "HM", "ConnectiveTissueDisease", "TransplantHistory")

clinical_history_labels <- c(Mortality28d = "28-day mortality", Age = "Age", SOFA_24h = "SOFA 24h", APACHEII_24h = "APACHE-II 24h", 
    WBC = "WBC", LymphocyteCount = "Lymphocyte count", NeutrophilCount = "Neutrophil count", PlateletCount = "Platelet count", 
    PCT = "PCT", hsCRP = "hsCRP", PaO2_FiO2 = "PaO2/FiO2", Immunosuppression = "Immunosuppression", MV = "Mechanical ventilation", 
    DM = "Diabetes mellitus", MI = "Myocardial infarction", COPD = "Chronic pulmonary disease", HepaticImpairment = "Liver disease", 
    RenalDisease = "Chronic kidney disease", Tumor = "Solid tumour", HM = "Haematologic malignancy", 
    ConnectiveTissueDisease = "Connective tissue disease", TransplantHistory = "Transplantation")

p_signif <- function(p) {
    case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

spearman_pair <- function(dat, x, y) {
    tmp <- dat %>% select(all_of(c(x, y))) %>% drop_na()
    if (nrow(tmp) < 3 || n_distinct(tmp[[x]]) < 2 || n_distinct(tmp[[y]]) < 2) {
        return(tibble(var_x = x, var_y = y, rho = NA_real_, p = NA_real_))
    }
    test <- suppressWarnings(cor.test(tmp[[x]], tmp[[y]], method = "spearman", exact = FALSE))
    tibble(var_x = x, var_y = y, rho = unname(test$estimate), p = test$p.value)
}

d1_map <- df_long %>% filter(Timepoint == "D1") %>% distinct(HumanID, SampleID)

d1_samples <- intersect(d1_map$SampleID, colnames(data))

d1_map <- d1_map %>% filter(SampleID %in% d1_samples)

d1_data <- data[, d1_map$SampleID, drop = FALSE]

keep_features <- rownames(d1_data)[rowSums(d1_data > 0, na.rm = TRUE) >= min_positive_n]

keep_features <- setdiff(keep_features, exclude_microbes)

d1_data <- d1_data[keep_features, , drop = FALSE]

microbe_labels <- rownames(d1_data)

microbe_labels[microbe_labels == "HHV-4"] <- "EBV"

microbe_vars <- paste0("mi_", make.names(microbe_labels, unique = TRUE))

microbe_label_map <- setNames(microbe_labels, microbe_vars)

d1_biomass <- as_tibble(t(log2(d1_data + 1)), .name_repair = "minimal") %>% setNames(microbe_vars) %>% 
    mutate(SampleID = colnames(d1_data), .before = 1) %>% inner_join(d1_map, by = "SampleID") %>% select(HumanID, 
    SampleID, all_of(microbe_vars))

analysis_df <- metadata %>% select(HumanID, all_of(clinical_history_vars)) %>% inner_join(d1_biomass, 
    by = "HumanID") %>% mutate(across(all_of(clinical_history_vars), as.numeric))

metadata_cor <- expand_grid(microbe = microbe_vars, marker = clinical_history_vars) %>% pmap_dfr(~spearman_pair(analysis_df, 
    ..1, ..2)) %>% mutate(p.signif = p_signif(p), microbe_label = recode(var_x, !!!microbe_label_map), 
    marker_label = recode(var_y, !!!clinical_history_labels), microbe_label = factor(microbe_label, levels = rev(microbe_label_map)), 
    marker_label = factor(marker_label, levels = unname(clinical_history_labels)))

metadata_sig_microbes <- metadata_cor %>% group_by(var_x) %>% summarise(any_significant = any(!is.na(p) & 
    p < 0.05), .groups = "drop") %>% filter(any_significant) %>% pull(var_x)

metadata_cor_plot <- metadata_cor %>% filter(var_x %in% metadata_sig_microbes) %>% mutate(microbe_label = factor(microbe_label, 
    levels = rev(microbe_label_map[metadata_sig_microbes])))

metadata_rho_mat <- metadata_cor_plot %>% select(microbe_label, marker_label, rho) %>% pivot_wider(names_from = marker_label, 
    values_from = rho) %>% column_to_rownames("microbe_label") %>% as.matrix()

metadata_sig_mat <- metadata_cor_plot %>% select(microbe_label, marker_label, p.signif) %>% pivot_wider(names_from = marker_label, 
    values_from = p.signif) %>% column_to_rownames("microbe_label") %>% as.matrix()

history_labels <- unname(clinical_history_labels[c("Immunosuppression", "MV", "DM", "MI", "COPD", "HepaticImpairment", 
    "RenalDisease", "Tumor", "HM", "ConnectiveTissueDisease", "TransplantHistory")])

metadata_col_group <- if_else(colnames(metadata_rho_mat) %in% history_labels, "Medical history", "Clinical marker") %>% 
    factor(levels = c("Clinical marker", "Medical history"))

metadata_top_anno <- HeatmapAnnotation(Group = metadata_col_group, col = list(Group = c(`Clinical marker` = "#9ECAE1", 
    `Medical history` = "#FDD0A2")), annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(3, 
    "mm"))

metadata_heatmap <- Heatmap(metadata_rho_mat, name = "Spearman\nrho", col = colorRamp2(c(-0.35, 0, 0.35), 
    c("#3B6FB6", "white", "#C4413A")), rect_gp = gpar(col = "white", lwd = 0.7), na_col = "grey95", cluster_rows = TRUE, 
    cluster_columns = TRUE, top_annotation = metadata_top_anno, column_split = metadata_col_group, cluster_column_slices = FALSE, 
    row_names_side = "left", row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), 
    column_names_rot = 45, heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), 
        labels_gp = gpar(fontsize = 8)), cell_fun = function(j, i, x, y, width, height, fill) {
        star <- metadata_sig_mat[i, j]
        if (!is.na(star) && star != "") {
            grid.text(star, x, y, gp = gpar(fontsize = 8, fontface = "bold"))
        }
    })

pdf(repo_output(repo_path(out_dir, "D1_microbe_clinical_history_spearman_heatmap.pdf"), "Figure2C_analysis.R"), 
    width = 10, height = 4.4)

draw(metadata_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right", column_title = "D1 microbe biomass versus clinical and medical history variables\nRaw Spearman P: * < 0.05, ** < 0.01, *** < 0.001")

dev.off()


write.csv(metadata_cor, "Outputs/Figure2C_all_tests.csv", row.names=FALSE)
write.csv(metadata_cor_plot, "Outputs/Figure2C_plot.csv", row.names=FALSE)
