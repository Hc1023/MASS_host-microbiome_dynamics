# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(glmnet)
    library(pROC)
    library(tidyr)
})

project_dir <- "."

input_file <- repo_path(project_dir, "Inputs/260823_meta_model.rdata")

step2_cohort_file <- repo_path(project_dir, "Inputs/260825_D4_prepared_data.rds")

legacy_outer_file <- repo_path(project_dir, "Inputs/260825_D4_outer_fold_assignments_used.csv")

legacy_inner_file <- repo_path(project_dir, "Inputs/260825_D4_inner_fold_assignments.csv")

output_dir <- repo_path(project_dir, "Outputs")

repo_dir(output_dir, recursive = TRUE, showWarnings = FALSE)

base_seed <- 260901L

n_repeats <- 20L

n_outer_folds <- 5L

n_inner_folds <- 5L

alpha_value <- 0.5

stable_frequency_threshold <- 0.5

input_env <- new.env(parent = emptyenv())

load(repo_input(input_file), envir = input_env)

stopifnot(all(c("meta_model", "model_prepare_info") %in% ls(input_env)))

if (!repo_exists(step2_cohort_file)) stop("Missing Step2 paired-cohort source: ", step2_cohort_file)

meta_model <- input_env$meta_model

model_prepare_info <- input_env$model_prepare_info

stopifnot(!anyDuplicated(meta_model$HumanID))

step2_prepared <- readRDS(repo_input(step2_cohort_file))

paired_ids <- as.character(step2_prepared$main$HumanID)

if (length(paired_ids) != 236L || anyDuplicated(paired_ids)) {
    stop("The Step2 paired-cohort source does not contain 236 unique patients")
}

cts_genes <- as.character(model_prepare_info$cts_genes_detected)

trajectory_genes <- as.character(model_prepare_info$module_genes_detected)

cts_vars <- paste0("D4_CTSg_", cts_genes)

trajectory_vars <- paste0("D4_DYNg_", trajectory_genes)

module_vars <- c("D4_mod_up1", "D4_mod_up2", "D4_mod_up3", "D4_mod_dw1")

required_vars <- c("HumanID", "D4", "Mortality28d", "SOFA_96h", module_vars, cts_vars, trajectory_vars)

missing_vars <- setdiff(required_vars, names(meta_model))

if (length(missing_vars)) stop("Missing required variables: ", paste(missing_vars, collapse = ", "))

if (length(cts_genes) != 18L || length(trajectory_genes) != 371L || anyDuplicated(cts_genes) || anyDuplicated(trajectory_genes) || 
    length(intersect(cts_genes, trajectory_genes))) {
    stop("Unexpected CTS/trajectory feature definitions")
}

analysis_cohort <- meta_model %>% filter(HumanID %in% paired_ids) %>% select(all_of(required_vars)) %>% 
    mutate(Mortality28d = as.integer(as.character(Mortality28d)), across(-c(HumanID, D4, Mortality28d), 
        as.numeric)) %>% arrange(HumanID)

if (nrow(analysis_cohort) != 236L || !setequal(analysis_cohort$HumanID, paired_ids) || anyNA(analysis_cohort) || 
    !identical(sort(unique(analysis_cohort$Mortality28d)), c(0L, 1L)) || any(!is.finite(as.matrix(analysis_cohort[, 
    setdiff(names(analysis_cohort), c("HumanID", "D4"))])))) {
    stop("The common Day-4 cohort failed expected n/outcome/completeness checks")
}

feature_metadata <- bind_rows(tibble(Feature = cts_vars, Gene = cts_genes, Gene_source = "CTS"), tibble(Feature = trajectory_vars, 
    Gene = trajectory_genes, Gene_source = "trajectory"))

make_stratified_folds <- function(y, k, seed) {
    set.seed(seed)
    fold <- integer(length(y))
    for (event in sort(unique(y))) {
        index <- sample(which(y == event))
        fold[index] <- rep(seq_len(k), length.out = length(index))
    }
    if (any(tabulate(fold, nbins = k) == 0L)) 
        stop("An empty fold was created")
    fold
}

if (!repo_exists(legacy_outer_file) || !repo_exists(legacy_inner_file)) {
    stop("Missing legacy frozen fold assignments")
}

outer_folds <- read.csv(repo_input(legacy_outer_file), check.names = FALSE) %>% transmute(HumanID = as.character(HumanID), 
    Repeat = as.integer(Repeat), Outer_fold = as.integer(Fold))

inner_folds <- read.csv(repo_input(legacy_inner_file), check.names = FALSE) %>% transmute(HumanID = as.character(HumanID), 
    Repeat = as.integer(Repeat), Outer_fold = as.integer(OuterFold), Inner_fold = as.integer(InnerFold))

if (nrow(outer_folds) != nrow(analysis_cohort) * n_repeats || any(count(outer_folds, Repeat, HumanID)$n != 
    1L) || n_distinct(outer_folds$Repeat) != n_repeats || any((outer_folds %>% count(Repeat, Outer_fold) %>% 
    count(Repeat))$n != n_outer_folds) || any(count(inner_folds, Repeat, Outer_fold, HumanID)$n != 1L) || 
    any(inner_folds$Inner_fold < 1L | inner_folds$Inner_fold > n_inner_folds) || !setequal(outer_folds$HumanID, 
    analysis_cohort$HumanID)) {
    stop("Frozen fold assignments failed QC")
}

write.csv(outer_folds, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_outer_folds.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(inner_folds, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_inner_folds.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

calc_auc <- function(y, probability) {
    as.numeric(pROC::auc(pROC::roc(response = y, predictor = probability, levels = c(0, 1), direction = "<", 
        quiet = TRUE)))
}

fit_outer_elastic <- function(training_data, test_data, gene_features, repeat_id, outer_fold_id) {
    feature_order <- c("SOFA_96h", gene_features)
    fold_table <- inner_folds %>% filter(Repeat == repeat_id, Outer_fold == outer_fold_id)
    inner_foldid <- fold_table$Inner_fold[match(training_data$HumanID, fold_table$HumanID)]
    if (anyNA(inner_foldid) || n_distinct(inner_foldid) != n_inner_folds) {
        stop("Failed to map frozen inner folds")
    }
    fit <- cv.glmnet(x = as.matrix(training_data[, feature_order, drop = FALSE]), y = training_data$Mortality28d, 
        family = "binomial", alpha = alpha_value, foldid = inner_foldid, nfolds = n_inner_folds, type.measure = "deviance", 
        standardize = TRUE, intercept = TRUE, penalty.factor = c(0, rep(1, length(gene_features))))
    coefficient <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
    selected <- setdiff(names(coefficient)[coefficient != 0], c("(Intercept)", "SOFA_96h"))
    list(probability = as.numeric(predict(fit, newx = as.matrix(test_data[, feature_order, drop = FALSE]), 
        s = "lambda.min", type = "response")), lambda = fit$lambda.min, selected = selected, coefficients = coefficient[selected])
}

model_definitions <- tibble(Model = c("M0", "M1", "M2", "M3"), Display = c("SOFA96", "SOFA96 + 18 D4 CTS genes", 
    "SOFA96 + 371 D4 trajectory genes", "SOFA96 + 18 D4 CTS genes + 371 D4 trajectory genes"), Algorithm = c("Logistic regression", 
    rep("Elastic net (alpha=0.5; inner 5-fold lambda.min)", 3)))

oof_list <- vector("list", n_repeats * n_outer_folds)

selection_list <- vector("list", n_repeats * n_outer_folds * 3L)

outer_model_list <- vector("list", n_repeats * n_outer_folds * 3L)

oof_index <- 1L

selection_index <- 1L

for (repeat_id in seq_len(n_repeats)) {
    for (outer_fold_id in seq_len(n_outer_folds)) {
        test_ids <- outer_folds$HumanID[outer_folds$Repeat == repeat_id & outer_folds$Outer_fold == outer_fold_id]
        training_data <- analysis_cohort %>% filter(!HumanID %in% test_ids)
        test_data <- analysis_cohort %>% filter(HumanID %in% test_ids)
        m0_fit <- glm(Mortality28d ~ SOFA_96h, data = training_data, family = binomial())
        elastic_fits <- list(M1 = fit_outer_elastic(training_data, test_data, cts_vars, repeat_id, outer_fold_id), 
            M2 = fit_outer_elastic(training_data, test_data, trajectory_vars, repeat_id, outer_fold_id), 
            M3 = fit_outer_elastic(training_data, test_data, c(cts_vars, trajectory_vars), repeat_id, 
                outer_fold_id))
        oof_list[[oof_index]] <- tibble(HumanID = test_data$HumanID, Mortality28d = test_data$Mortality28d, 
            Repeat = repeat_id, Outer_fold = outer_fold_id, M0_probability = as.numeric(predict(m0_fit, 
                test_data, type = "response")), M1_probability = elastic_fits$M1$probability, M2_probability = elastic_fits$M2$probability, 
            M3_probability = elastic_fits$M3$probability)
        oof_index <- oof_index + 1L
        fit_results <- elastic_fits
        for (model in c("M1", "M2", "M3")) {
            fit_result <- fit_results[[model]]
            selected_metadata <- feature_metadata %>% filter(Feature %in% fit_result$selected) %>% mutate(Repeat = repeat_id, 
                Outer_fold = outer_fold_id, Model = model, Lambda_min = fit_result$lambda, Coefficient = unname(fit_result$coefficients[Feature])) %>% 
                select(Repeat, Outer_fold, Model, Lambda_min, Feature, Gene, Gene_source, Coefficient)
            selection_list[[selection_index]] <- selected_metadata
            outer_model_list[[selection_index]] <- tibble(Repeat = repeat_id, Outer_fold = outer_fold_id, 
                Model = model, Lambda_min = fit_result$lambda, Selected_gene_n = length(fit_result$selected), 
                Selected_CTS_n = sum(fit_result$selected %in% cts_vars), Selected_trajectory_n = sum(fit_result$selected %in% 
                  trajectory_vars), Selected_genes = paste(feature_metadata$Gene[match(fit_result$selected, 
                  feature_metadata$Feature)], collapse = ";"))
            selection_index <- selection_index + 1L
        }
    }
}

oof_predictions <- bind_rows(oof_list) %>% arrange(Repeat, Outer_fold, HumanID)

gene_selection <- bind_rows(selection_list) %>% arrange(Model, Repeat, Outer_fold, Gene_source, Gene)

outer_model_summary <- bind_rows(outer_model_list) %>% arrange(Model, Repeat, Outer_fold)

probability_cols <- paste0(model_definitions$Model, "_probability")

if (nrow(oof_predictions) != nrow(analysis_cohort) * n_repeats || any(count(oof_predictions, Repeat, 
    HumanID)$n != 1L) || any(!is.finite(as.matrix(oof_predictions[, probability_cols]))) || nrow(outer_model_summary) != 
    3L * n_repeats * n_outer_folds || any(count(outer_model_summary, Model, Repeat, Outer_fold)$n != 
    1L)) {
    stop("OOF predictions or outer-model audit failed QC")
}

repeat_performance <- lapply(seq_len(n_repeats), function(repeat_id) {
    repeat_data <- oof_predictions %>% filter(Repeat == repeat_id)
    bind_rows(lapply(model_definitions$Model, function(model) {
        probability <- repeat_data[[paste0(model, "_probability")]]
        tibble(Repeat = repeat_id, Model = model, AUC = calc_auc(repeat_data$Mortality28d, probability), 
            Brier = mean((repeat_data$Mortality28d - probability)^2))
    }))
}) %>% bind_rows()

performance <- repeat_performance %>% group_by(Model) %>% summarise(N = nrow(analysis_cohort), Deaths = sum(analysis_cohort$Mortality28d == 
    1L), Survivors = sum(analysis_cohort$Mortality28d == 0L), Repeats = n(), Outer_folds_per_repeat = n_outer_folds, 
    Mean_CV_AUC = mean(AUC), AUC_Q25 = quantile(AUC, 0.25), AUC_Q75 = quantile(AUC, 0.75), AUC_IQR = IQR(AUC), 
    Mean_CV_Brier = mean(Brier), Brier_Q25 = quantile(Brier, 0.25), Brier_Q75 = quantile(Brier, 0.75), 
    Brier_IQR = IQR(Brier), .groups = "drop") %>% left_join(model_definitions, by = "Model") %>% select(Model, 
    Display, Algorithm, everything())

contrast_definitions <- tibble(Contrast = c("M1 - M0", "M2 - M0", "M3 - M0", "M2 - M1", "M3 - M2", "M3 - M1"), 
    New_model = c("M1", "M2", "M3", "M2", "M3", "M3"), Reference_model = c("M0", "M0", "M0", "M1", "M2", 
        "M1"))

contrast_distributions <- lapply(seq_len(nrow(contrast_definitions)), function(i) {
    definition <- contrast_definitions[i, ]
    wide <- repeat_performance %>% filter(Model %in% c(definition$New_model, definition$Reference_model)) %>% 
        pivot_wider(names_from = Model, values_from = c(AUC, Brier))
    tibble(Repeat = wide$Repeat, Contrast = definition$Contrast, New_model = definition$New_model, Reference_model = definition$Reference_model, 
        DeltaAUC = wide[[paste0("AUC_", definition$New_model)]] - wide[[paste0("AUC_", definition$Reference_model)]], 
        DeltaBrier = wide[[paste0("Brier_", definition$New_model)]] - wide[[paste0("Brier_", definition$Reference_model)]])
}) %>% bind_rows()

contrasts <- contrast_distributions %>% group_by(Contrast, New_model, Reference_model) %>% summarise(Repeats = n(), 
    Mean_DeltaAUC = mean(DeltaAUC), DeltaAUC_Q25 = quantile(DeltaAUC, 0.25), DeltaAUC_Q75 = quantile(DeltaAUC, 
        0.75), DeltaAUC_IQR = IQR(DeltaAUC), Mean_DeltaBrier = mean(DeltaBrier), DeltaBrier_Q25 = quantile(DeltaBrier, 
        0.25), DeltaBrier_Q75 = quantile(DeltaBrier, 0.75), DeltaBrier_IQR = IQR(DeltaBrier), .groups = "drop") %>% 
    mutate(Contrast = factor(Contrast, levels = contrast_definitions$Contrast)) %>% arrange(Contrast) %>% 
    mutate(Contrast = as.character(Contrast))

candidate_grid <- bind_rows(crossing(Model = "M1", feature_metadata %>% filter(Gene_source == "CTS")), 
    crossing(Model = "M2", feature_metadata %>% filter(Gene_source == "trajectory")), crossing(Model = "M3", 
        feature_metadata))

gene_frequency <- candidate_grid %>% left_join(gene_selection %>% count(Model, Feature, name = "Selection_count"), 
    by = c("Model", "Feature")) %>% mutate(Selection_count = replace_na(Selection_count, 0L), Outer_model_n = n_repeats * 
    n_outer_folds, Selection_frequency = Selection_count/Outer_model_n, Stable_50pct = Selection_frequency >= 
    stable_frequency_threshold) %>% arrange(Model, desc(Selection_frequency), Gene_source, Gene)

full_inner_fold <- make_stratified_folds(analysis_cohort$Mortality28d, n_inner_folds, base_seed + 999999L)

fit_full_elastic <- function(model, gene_features) {
    feature_order <- c("SOFA_96h", gene_features)
    x <- as.matrix(analysis_cohort[, feature_order, drop = FALSE])
    fit <- cv.glmnet(x = x, y = analysis_cohort$Mortality28d, family = "binomial", alpha = alpha_value, 
        foldid = full_inner_fold, nfolds = n_inner_folds, type.measure = "deviance", standardize = TRUE, 
        intercept = TRUE, penalty.factor = c(0, rep(1, length(gene_features))))
    coefficient <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
    selected <- setdiff(names(coefficient)[coefficient != 0], c("(Intercept)", "SOFA_96h"))
    coefficient_table <- tibble(Model = model, Feature = names(coefficient), Coefficient = unname(coefficient), 
        Lambda_min = fit$lambda.min) %>% filter(Feature %in% c("(Intercept)", "SOFA_96h", selected)) %>% 
        left_join(feature_metadata, by = "Feature") %>% mutate(Gene = case_when(Feature == "(Intercept)" ~ 
        "(Intercept)", Feature == "SOFA_96h" ~ "SOFA96", TRUE ~ Gene), Gene_source = case_when(Feature == 
        "(Intercept)" ~ "Intercept", Feature == "SOFA_96h" ~ "Clinical", TRUE ~ Gene_source), Penalized = !Feature %in% 
        c("(Intercept)", "SOFA_96h"))
    list(model = model, fit = fit, lambda_min = fit$lambda.min, selected_features = selected, coefficients = coefficient, 
        coefficient_table = coefficient_table, full_inner_fold = full_inner_fold, feature_order = feature_order, 
        penalty_factor = c(0, rep(1, length(gene_features))))
}

full_fits <- list(M1 = fit_full_elastic("M1", cts_vars), M2 = fit_full_elastic("M2", trajectory_vars), 
    M3 = fit_full_elastic("M3", c(cts_vars, trajectory_vars)))

final_coefficients <- bind_rows(lapply(full_fits, `[[`, "coefficient_table")) %>% arrange(Model, factor(Gene_source, 
    c("Intercept", "Clinical", "CTS", "trajectory")), Gene)

final_genes <- final_coefficients %>% filter(Gene_source %in% c("CTS", "trajectory")) %>% left_join(gene_frequency %>% 
    select(Model, Feature, Selection_count, Outer_model_n, Selection_frequency, Stable_50pct), by = c("Model", 
    "Feature")) %>% select(Model, Gene, Gene_source, Feature, Coefficient, Lambda_min, Selection_count, 
    Outer_model_n, Selection_frequency, Stable_50pct)

frozen_bundle <- list(freeze_id = "260901_Step3_D4_signature", source_input = input_file, cohort_definition = c("Step2 paired D1+D4 transcriptome cohort", 
    "SOFA_96h", "Mortality28d"), cohort = analysis_cohort %>% select(HumanID, Mortality28d, SOFA_96h), 
    feature_metadata = feature_metadata, cv_specification = list(fold_source = c(outer = legacy_outer_file, 
        inner = legacy_inner_file), full_data_fit_seed = base_seed + 999999L, repeats = n_repeats, outer_folds = n_outer_folds, 
        inner_folds = n_inner_folds, stratified = TRUE, alpha = alpha_value, lambda = "lambda.min", SOFA_unpenalized = TRUE, 
        genes_penalized = TRUE, pooled_OOF_per_repeat = TRUE), outer_folds = outer_folds, inner_folds = inner_folds, 
    full_data_inner_fold = tibble(HumanID = analysis_cohort$HumanID, Inner_fold = full_inner_fold), full_fits = full_fits, 
    trajectory_signature = full_fits$M2, final_genes = final_genes, final_coefficients = final_coefficients)

write.csv(performance, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_performance.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(contrasts, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_contrasts.csv"), "Figure6BC_FigureS11B_analysis.R"), 
    row.names = FALSE)

write.csv(oof_predictions, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_oof_predictions.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(gene_selection, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_gene_selection.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(gene_frequency, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_gene_frequency.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(final_genes, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_final_genes.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(final_coefficients, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_final_coefficients.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(repeat_performance, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_repeat_performance.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(contrast_distributions, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_contrast_distributions.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

write.csv(outer_model_summary, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_outer_model_summary.csv"), 
    "Figure6BC_FigureS11B_analysis.R"), row.names = FALSE)

saveRDS(frozen_bundle, repo_output(repo_path(output_dir, "260901_Step3_D4_signature_frozen_model.rds"), 
    "Figure6BC_FigureS11B_analysis.R"), compress = "xz")

get_perf <- function(model, field) performance[[field]][performance$Model == model]

get_contrast <- function(contrast, field) contrasts[[field]][contrasts$Contrast == contrast]

qtext <- function(x) sprintf("%.0f (IQR %.0f-%.0f)", median(x), quantile(x, 0.25), quantile(x, 0.75))

gene_list <- function(model) {
    genes <- final_genes$Gene[final_genes$Model == model]
    if (length(genes)) 
        paste(genes, collapse = ", ")
    else "none"
}

stable_list <- function(model, source = NULL) {
    x <- gene_frequency %>% filter(Model == model, Stable_50pct)
    if (!is.null(source)) 
        x <- x %>% filter(Gene_source == source)
    if (!nrow(x)) 
        return("none at the prespecified >=50% threshold")
    paste(sprintf("%s (%.0f%%)", x$Gene, 100 * x$Selection_frequency), collapse = ", ")
}

count_lines <- vapply(c("M1", "M2", "M3"), function(model) {
    x <- outer_model_summary %>% filter(Model == model)
    if (model == "M3") {
        sprintf("- M3 selected %s total genes: %s CTS and %s trajectory genes.", qtext(x$Selected_gene_n), 
            qtext(x$Selected_CTS_n), qtext(x$Selected_trajectory_n))
    }
    else {
        sprintf("- %s selected %s genes.", model, qtext(x$Selected_gene_n))
    }
}, character(1))

performance_lines <- vapply(model_definitions$Model, function(model) {
    sprintf("- %s: mean AUC %.3f (IQR %.3f-%.3f); mean Brier %.3f (IQR %.3f-%.3f).", model, get_perf(model, 
        "Mean_CV_AUC"), get_perf(model, "AUC_Q25"), get_perf(model, "AUC_Q75"), get_perf(model, "Mean_CV_Brier"), 
        get_perf(model, "Brier_Q25"), get_perf(model, "Brier_Q75"))
}, character(1))

contrast_lines <- vapply(contrast_definitions$Contrast, function(contrast) {
    sprintf("- %s: mean DeltaAUC %+.3f (IQR %.3f to %.3f); mean DeltaBrier %+.3f (IQR %.3f to %.3f).", 
        contrast, get_contrast(contrast, "Mean_DeltaAUC"), get_contrast(contrast, "DeltaAUC_Q25"), get_contrast(contrast, 
            "DeltaAUC_Q75"), get_contrast(contrast, "Mean_DeltaBrier"), get_contrast(contrast, "DeltaBrier_Q25"), 
        get_contrast(contrast, "DeltaBrier_Q75"))
}, character(1))

m3_m2_auc <- get_contrast("M3 - M2", "Mean_DeltaAUC")

m3_stable_cts <- gene_frequency %>% filter(Model == "M3", Gene_source == "CTS", Stable_50pct) %>% nrow()

trajectory_close <- abs(m3_m2_auc) < 0.01

why_m3 <- if (trajectory_close) {
    "M2 was chosen because its trajectory-focused feature space achieved discrimination within 0.01 AUC of M3 while avoiding expansion into the combined CTS-plus-trajectory space."
} else {
    "M2 was retained as the trajectory-focused signature to preserve the prespecified biological trajectory construct; its performance difference from M3 is reported without claiming equivalence."
}

m4_stability_text <- if (m3_stable_cts == 0L) {
    "No CTS gene in M3 reached the prespecified >=50% outer-model selection-frequency threshold."
} else {
    sprintf("M3 contained %d CTS gene(s) selected in >=50%% of outer models: %s.", m3_stable_cts, stable_list("M3", 
        "CTS"))
}

