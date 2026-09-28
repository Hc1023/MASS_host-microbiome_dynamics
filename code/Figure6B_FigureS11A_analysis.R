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

legacy_dir <- repo_path(project_dir, "Outputs")

output_dir <- repo_path(project_dir, "Outputs")

repo_dir(output_dir, recursive = TRUE, showWarnings = FALSE)

n_repeats <- 20L

n_outer_folds <- 5L

n_inner_folds <- 5L

load(repo_input(input_file))
source("code/shared_sources.R")
meta_model <- module_overlay(meta_model)

stopifnot(exists("meta_model"), !anyDuplicated(meta_model$HumanID))

cts_genes <- c("SLC4A1", "SERPINB1", "FECH", "ACER3", "NLRC4", "BTN3A3", "SNX3", "GADD45A", "CA1", "PGD", 
    "STOM", "TDRD9", "HK3", "EPB42", "BPGM", "GLRX5", "UBE2H", "METTL9")

d1_cts_vars <- paste0("D1_CTSg_", cts_genes)

d4_cts_vars <- paste0("D4_CTSg_", cts_genes)

module_vars <- c("D4_mod_up1", "D4_mod_up2", "D4_mod_up3", "D4_mod_dw1")

viral_vars <- c("D4_mi_HHV.4", "D4_mi_HCMV")

required_vars <- c("HumanID", "Mortality28d", "SOFA_96h", d1_cts_vars, d4_cts_vars, module_vars, viral_vars)

missing_vars <- setdiff(required_vars, names(meta_model))

if (length(missing_vars) > 0L) {
    stop("Required variables are absent: ", paste(missing_vars, collapse = ", "))
}

if (length(d1_cts_vars) != 18L || length(d4_cts_vars) != 18L || anyDuplicated(d1_cts_vars) || anyDuplicated(d4_cts_vars) || 
    length(module_vars) != 4L || anyDuplicated(module_vars)) {
    stop("Fixed CTS/module feature definitions are invalid")
}

legacy_prepared_file <- repo_path(legacy_dir, "Inputs/260825_D4_prepared_data.rds")

legacy_outer_file <- repo_path(legacy_dir, "Inputs/260825_D4_outer_fold_assignments_used.csv")

legacy_inner_file <- repo_path(legacy_dir, "Inputs/260825_D4_inner_fold_assignments.csv")

stopifnot(repo_exists(legacy_prepared_file), repo_exists(legacy_outer_file), repo_exists(legacy_inner_file))

legacy_prepared <- readRDS(repo_input(legacy_prepared_file))

outer_folds <- read.csv(repo_input(legacy_outer_file), check.names = FALSE)

inner_folds <- read.csv(repo_input(legacy_inner_file), check.names = FALSE)

paired_ids <- legacy_prepared$main$HumanID

analysis_cohort <- meta_model %>% filter(HumanID %in% paired_ids) %>% select(all_of(required_vars)) %>% 
    mutate(Mortality28d = as.integer(as.character(Mortality28d)), across(-c(HumanID, Mortality28d), as.numeric)) %>% 
    arrange(HumanID)

if (!setequal(analysis_cohort$HumanID, paired_ids) || nrow(analysis_cohort) != 236L) {
    stop("The reused paired cohort is not the expected 236 patients")
}

if (anyNA(analysis_cohort) || !identical(sort(unique(analysis_cohort$Mortality28d)), c(0L, 1L)) || any(!is.finite(as.matrix(analysis_cohort[, 
    -1])))) {
    stop("The paired cohort has missing/non-finite data or an invalid outcome")
}

if (!setequal(outer_folds$HumanID, analysis_cohort$HumanID) || n_distinct(outer_folds$Repeat) != n_repeats || 
    any(count(outer_folds, Repeat, HumanID)$n != 1L) || any(count(outer_folds, Repeat)$n != nrow(analysis_cohort)) || 
    any((outer_folds %>% group_by(Repeat) %>% summarise(K = n_distinct(Fold), .groups = "drop"))$K != 
        n_outer_folds)) {
    stop("Legacy outer-fold assignments fail reuse checks")
}

outer_fold_sizes <- outer_folds %>% count(Repeat, Fold)

if (nrow(outer_fold_sizes) != n_repeats * n_outer_folds || any(outer_fold_sizes$n == 0L)) stop("Legacy outer folds are incomplete")

expected_inner <- lapply(seq_len(n_repeats), function(repeat_id) {
    lapply(seq_len(n_outer_folds), function(fold_id) {
        training_ids <- setdiff(analysis_cohort$HumanID, outer_folds$HumanID[outer_folds$Repeat == repeat_id & 
            outer_folds$Fold == fold_id])
        tibble(Repeat = repeat_id, OuterFold = fold_id, HumanID = training_ids)
    }) %>% bind_rows()
}) %>% bind_rows()

if (!setequal(paste(expected_inner$Repeat, expected_inner$OuterFold, expected_inner$HumanID), paste(inner_folds$Repeat, 
    inner_folds$OuterFold, inner_folds$HumanID)) || any(count(inner_folds, Repeat, OuterFold, HumanID)$n != 
    1L) || any(inner_folds$InnerFold < 1L | inner_folds$InnerFold > n_inner_folds)) {
    stop("Legacy inner-fold assignments do not match the outer training sets")
}

inner_fold_sizes <- inner_folds %>% count(Repeat, OuterFold, InnerFold)

if (nrow(inner_fold_sizes) != n_repeats * n_outer_folds * n_inner_folds || any(inner_fold_sizes$n == 
    0L)) stop("Legacy inner folds are incomplete")

write.csv(outer_folds, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_outer_folds_reused.csv"), 
    "Figure6B_FigureS11A_analysis.R"), row.names = FALSE)

write.csv(inner_folds, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_inner_folds_reused.csv"), 
    "Figure6B_FigureS11A_analysis.R"), row.names = FALSE)

calc_auc <- function(y, probability) {
    as.numeric(pROC::auc(pROC::roc(response = y, predictor = probability, levels = c(0, 1), direction = "<", 
        quiet = TRUE)))
}

fit_nested_glmnet <- function(train, test, feature_vars, unpenalized_vars, repeat_id, outer_fold_id) {
    fold_table <- inner_folds %>% filter(Repeat == repeat_id, OuterFold == outer_fold_id)
    inner_foldid <- fold_table$InnerFold[match(train$HumanID, fold_table$HumanID)]
    if (anyNA(inner_foldid) || n_distinct(inner_foldid) != n_inner_folds) {
        stop("Failed to map reused inner folds")
    }
    penalty <- ifelse(feature_vars %in% unpenalized_vars, 0, 1)
    cv_fit <- cv.glmnet(x = as.matrix(train[, feature_vars, drop = FALSE]), y = train$Mortality28d, family = "binomial", 
        alpha = 0.5, type.measure = "deviance", nfolds = n_inner_folds, foldid = inner_foldid, penalty.factor = penalty, 
        standardize = TRUE, intercept = TRUE)
    list(probability = as.numeric(predict(cv_fit, newx = as.matrix(test[, feature_vars, drop = FALSE]), 
        s = "lambda.min", type = "response")), lambda_min = cv_fit$lambda.min)
}

model_definitions <- tibble(Model = paste0("M", 0:5), Display = c("SOFA96", "SOFA96 + D1 CTS genes", 
    "SOFA96 + D4 CTS genes", "SOFA96 + D4 CTS genes + 4 D4 trajectory modules", "SOFA96 + D4 EBV/HCMV", 
    "SOFA96 + D4 EBV/HCMV + D4 CTS genes"), Feature_space = c("SOFA_96h", "SOFA_96h + 18 D1 CTS genes", 
    "SOFA_96h + 18 D4 CTS genes", "SOFA_96h + 18 D4 CTS genes + 4 D4 trajectory modules", "SOFA_96h + D4 EBV + D4 HCMV", 
    "SOFA_96h + D4 EBV + D4 HCMV + 18 D4 CTS genes"), Algorithm = c("Logistic regression", rep("Elastic net (alpha=0.5; reused inner 5-fold lambda.min)", 
    3), "Logistic regression", "Elastic net (alpha=0.5; reused inner 5-fold lambda.min)"), Analysis_role = c("Contemporaneous clinical baseline", 
    "Retained admission CTS state", "Updated Day-4 CTS state", "Day-4 trajectory-module increment", "Day-4 EBV/HCMV beyond contemporaneous SOFA96", 
    "Day-4 CTS genes beyond Day-4 EBV/HCMV and SOFA96"))

# Refit M3 only; all other out-of-fold predictions are immutable archived results.
oof_predictions <- read.csv("../MASS_mortality-main/Outputs/260831_Step2_D4_D1/260831_Step2_D4_D1_oof_predictions.csv")
old_oof <- oof_predictions
for(repeat_id in seq_len(n_repeats)) for(fold_id in seq_len(n_outer_folds)) {
 test_ids <- outer_folds$HumanID[outer_folds$Repeat==repeat_id & outer_folds$Fold==fold_id]
 train <- analysis_cohort %>% filter(!HumanID %in% test_ids)
 test <- analysis_cohort %>% filter(HumanID %in% test_ids)
 m3 <- fit_nested_glmnet(train,test,c("SOFA_96h",module_vars,d4_cts_vars),c("SOFA_96h",module_vars),repeat_id,fold_id)
 ix <- match(paste(repeat_id, test$HumanID),paste(oof_predictions$Repeat,oof_predictions$HumanID))
 stopifnot(!anyNA(ix))
 oof_predictions$M3_probability[ix] <- m3$probability
 oof_predictions$M3_lambda_min[ix] <- m3$lambda_min
}
stopifnot(identical(old_oof[,setdiff(names(old_oof),c("M3_probability","M3_lambda_min"))],oof_predictions[,setdiff(names(old_oof),c("M3_probability","M3_lambda_min"))]))
probability_cols <- paste0(model_definitions$Model, "_probability")

if (nrow(oof_predictions) != nrow(analysis_cohort) * n_repeats || any(count(oof_predictions, Repeat, 
    HumanID)$n != 1L) || any(!is.finite(as.matrix(oof_predictions[, probability_cols])))) {
    stop("Repeated outer held-out predictions failed QC")
}

repeat_performance <- lapply(seq_len(n_repeats), function(repeat_id) {
    repeat_data <- oof_predictions %>% filter(Repeat == repeat_id)
    bind_rows(lapply(model_definitions$Model, function(model) {
        probability <- repeat_data[[paste0(model, "_probability")]]
        tibble(Repeat = repeat_id, Model = model, AUC = calc_auc(repeat_data$Mortality28d, probability), 
            Brier = mean((repeat_data$Mortality28d - probability)^2))
    }))
}) %>% bind_rows()

cohort_n <- nrow(analysis_cohort)

cohort_deaths <- sum(analysis_cohort$Mortality28d == 1L)

cohort_survivors <- sum(analysis_cohort$Mortality28d == 0L)

performance <- repeat_performance %>% group_by(Model) %>% summarise(N = cohort_n, Deaths = cohort_deaths, 
    Survivors = cohort_survivors, Repeats = n(), Outer_folds_per_repeat = n_outer_folds, Mean_CV_AUC = mean(AUC), 
    AUC_Q25 = quantile(AUC, 0.25), AUC_Q75 = quantile(AUC, 0.75), AUC_IQR = IQR(AUC), SD_CV_AUC = sd(AUC), 
    Mean_CV_Brier = mean(Brier), Brier_Q25 = quantile(Brier, 0.25), Brier_Q75 = quantile(Brier, 0.75), 
    Brier_IQR = IQR(Brier), SD_CV_Brier = sd(Brier), .groups = "drop") %>% left_join(model_definitions, 
    by = "Model") %>% select(Model, Display, Analysis_role, Feature_space, Algorithm, everything())

contrast_definitions <- tibble(Contrast = c("M1 - M0", "M2 - M1", "M3 - M2", "M3 - M0", "M4 - M0", "M5 - M4", 
    "M5 - M2", "M5 - M3"), New_model = c("M1", "M2", "M3", "M3", "M4", "M5", "M5", "M5"), Reference_model = c("M0", 
    "M1", "M2", "M0", "M0", "M4", "M2", "M3"), Comparison_role = c("Admission CTS state beyond SOFA96", 
    "Updated D4 CTS versus retained D1 CTS", "D4 trajectory modules beyond updated D4 CTS", "Full D4 host reassessment versus SOFA96", 
    "D4 EBV/HCMV beyond contemporaneous SOFA96", "D4 CTS genes beyond D4 EBV/HCMV and SOFA96", "D4 EBV/HCMV beyond D4 CTS genes and SOFA96", 
    "Viral-plus-CTS path versus trajectory-module path"))

contrast_distributions <- lapply(seq_len(nrow(contrast_definitions)), function(i) {
    definition <- contrast_definitions[i, ]
    wide <- repeat_performance %>% filter(Model %in% c(definition$New_model, definition$Reference_model)) %>% 
        pivot_wider(names_from = Model, values_from = c(AUC, Brier))
    tibble(Repeat = wide$Repeat, Contrast = definition$Contrast, New_model = definition$New_model, Reference_model = definition$Reference_model, 
        Comparison_role = definition$Comparison_role, DeltaAUC = wide[[paste0("AUC_", definition$New_model)]] - 
            wide[[paste0("AUC_", definition$Reference_model)]], DeltaBrier = wide[[paste0("Brier_", definition$New_model)]] - 
            wide[[paste0("Brier_", definition$Reference_model)]])
}) %>% bind_rows()

mean_ci <- function(x) {
    estimate <- mean(x)
    half_width <- qt(0.975, df = length(x) - 1L) * sd(x)/sqrt(length(x))
    c(lower = estimate - half_width, upper = estimate + half_width)
}

contrasts <- contrast_distributions %>% group_by(Contrast, New_model, Reference_model, Comparison_role) %>% 
    summarise(Repeats = n(), Mean_paired_DeltaAUC = mean(DeltaAUC), DeltaAUC_Q25 = quantile(DeltaAUC, 
        0.25), DeltaAUC_Q75 = quantile(DeltaAUC, 0.75), DeltaAUC_IQR = IQR(DeltaAUC), DeltaAUC_95CI_lower = mean_ci(DeltaAUC)["lower"], 
        DeltaAUC_95CI_upper = mean_ci(DeltaAUC)["upper"], Paired_AUC_ttest_P = t.test(DeltaAUC, mu = 0)$p.value, 
        Mean_paired_DeltaBrier = mean(DeltaBrier), DeltaBrier_Q25 = quantile(DeltaBrier, 0.25), DeltaBrier_Q75 = quantile(DeltaBrier, 
            0.75), DeltaBrier_IQR = IQR(DeltaBrier), DeltaBrier_95CI_lower = mean_ci(DeltaBrier)["lower"], 
        DeltaBrier_95CI_upper = mean_ci(DeltaBrier)["upper"], Paired_Brier_ttest_P = t.test(DeltaBrier, 
            mu = 0)$p.value, DeltaBrier_definition = "new model minus reference model; negative favors new model", 
        .groups = "drop") %>% mutate(Contrast = factor(Contrast, levels = contrast_definitions$Contrast)) %>% 
    arrange(Contrast) %>% mutate(Contrast = as.character(Contrast))

write.csv(performance, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_performance.csv"), "Figure6B_FigureS11A_analysis.R"), 
    row.names = FALSE)

write.csv(contrasts, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_contrasts.csv"), "Figure6B_FigureS11A_analysis.R"), 
    row.names = FALSE)

write.csv(contrast_distributions, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_contrast_distributions.csv"), 
    "Figure6B_FigureS11A_analysis.R"), row.names = FALSE)

write.csv(repeat_performance, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_repeat_performance.csv"), 
    "Figure6B_FigureS11A_analysis.R"), row.names = FALSE)

write.csv(oof_predictions, repo_output(repo_path(output_dir, "260831_Step2_D4_D1_oof_predictions.csv"), 
    "Figure6B_FigureS11A_analysis.R"), row.names = FALSE)

