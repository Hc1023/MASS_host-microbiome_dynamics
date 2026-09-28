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

output_dir <- repo_path(project_dir, "Outputs")

repo_dir(output_dir, recursive = TRUE, showWarnings = FALSE)

base_seed <- 260829L

n_repeats <- 20L

n_outer_folds <- 5L

n_inner_folds <- 5L

load(repo_input(input_file))

stopifnot(exists("meta_model"), !anyDuplicated(meta_model$HumanID))

cts_genes <- c("SLC4A1", "SERPINB1", "FECH", "ACER3", "NLRC4", "BTN3A3", "SNX3", "GADD45A", "CA1", "PGD", 
    "STOM", "TDRD9", "HK3", "EPB42", "BPGM", "GLRX5", "UBE2H", "METTL9")

cts_vars <- paste0("D1_CTSg_", cts_genes)

ebv_var <- "D1_mi_HHV.4"

hcmv_var <- "D1_mi_HCMV"

viral_vars <- c(ebv_var, hcmv_var)

required_vars <- c("HumanID", "Mortality28d", "SOFA_24h", cts_vars, viral_vars)

missing_vars <- setdiff(required_vars, names(meta_model))

if (length(missing_vars) > 0L) {
    stop("Required variables are absent: ", paste(missing_vars, collapse = ", "))
}

if (length(cts_vars) != 18L || anyDuplicated(cts_vars)) {
    stop("The fixed D1 CTS feature space must contain exactly 18 unique genes")
}

analysis_cohort <- meta_model %>% select(all_of(required_vars)) %>% filter(if_all(everything(), ~!is.na(.x))) %>% 
    mutate(Mortality28d = as.integer(as.character(Mortality28d)), SOFA_24h = as.numeric(SOFA_24h), across(all_of(c(cts_vars, 
        viral_vars)), as.numeric))

if (!identical(sort(unique(analysis_cohort$Mortality28d)), c(0L, 1L))) {
    stop("Mortality28d must be coded 0/1")
}

if (any(!is.finite(as.matrix(analysis_cohort[, -1])))) {
    stop("The complete-case analysis cohort contains non-finite values")
}

if (any(vapply(analysis_cohort[cts_vars], sd, numeric(1)) == 0)) {
    stop("At least one D1 CTS gene has zero variance")
}

cohort_n <- nrow(analysis_cohort)

cohort_deaths <- sum(analysis_cohort$Mortality28d == 1L)

cohort_survivors <- sum(analysis_cohort$Mortality28d == 0L)

make_stratified_folds <- function(y, k, seed) {
    set.seed(seed)
    fold <- integer(length(y))
    for (event in sort(unique(y))) {
        idx <- sample(which(y == event))
        fold[idx] <- rep(seq_len(k), length.out = length(idx))
    }
    if (any(tabulate(fold, nbins = k) == 0L)) 
        stop("An empty fold was created")
    fold
}

calc_auc <- function(y, probability) {
    as.numeric(pROC::auc(pROC::roc(response = y, predictor = probability, levels = c(0, 1), direction = "<", 
        quiet = TRUE)))
}

fit_nested_glmnet <- function(train, test, feature_vars, unpenalized_vars, inner_seed) {
    x_train <- as.matrix(train[, feature_vars, drop = FALSE])
    x_test <- as.matrix(test[, feature_vars, drop = FALSE])
    y_train <- train$Mortality28d
    penalty <- ifelse(feature_vars %in% unpenalized_vars, 0, 1)
    inner_foldid <- make_stratified_folds(y_train, n_inner_folds, inner_seed)
    cv_fit <- cv.glmnet(x = x_train, y = y_train, family = "binomial", alpha = 0.5, type.measure = "deviance", 
        nfolds = n_inner_folds, foldid = inner_foldid, penalty.factor = penalty, standardize = TRUE, 
        intercept = TRUE)
    coefficient <- as.matrix(coef(cv_fit, s = "lambda.min"))[-1, 1]
    list(probability = as.numeric(predict(cv_fit, newx = x_test, s = "lambda.min", type = "response")), 
        lambda = cv_fit$lambda.min, selected_CTS_n = sum(coefficient[cts_vars] != 0))
}

model_definitions <- tibble(Model = c("M0", "M1", "M2", "M3"), Display = c("SOFA24", "SOFA24 + EBV/HCMV", 
    "SOFA24 + CTS genes", "SOFA24 + CTS genes + EBV/HCMV"), Feature_space = c("SOFA_24h", "SOFA_24h + D1 EBV + D1 HCMV", 
    "SOFA_24h + 18 D1 CTS genes", "SOFA_24h + 18 D1 CTS genes + D1 EBV + D1 HCMV"), Algorithm = c("Logistic regression", 
    "Logistic regression", "Elastic net (alpha=0.5; inner 5-fold lambda.min)", "Elastic net (alpha=0.5; inner 5-fold lambda.min)"), 
    Analysis_role = c("Primary baseline", "Primary viral model", "Supplementary CTS comparator", "Supplementary combined comparator"))

oof_predictions <- lapply(seq_len(n_repeats), function(repeat_id) {
    outer_fold <- make_stratified_folds(analysis_cohort$Mortality28d, n_outer_folds, base_seed + repeat_id * 
        1000L)
    repeat_data <- analysis_cohort %>% mutate(Outer_fold = outer_fold)
    lapply(seq_len(n_outer_folds), function(fold_id) {
        train <- repeat_data %>% filter(Outer_fold != fold_id)
        test <- repeat_data %>% filter(Outer_fold == fold_id)
        d0_fit <- glm(Mortality28d ~ SOFA_24h, data = train, family = binomial())
        d0_probability <- as.numeric(predict(d0_fit, test, type = "response"))
        d1_fit <- glm(Mortality28d ~ SOFA_24h + D1_mi_HHV.4 + D1_mi_HCMV, data = train, family = binomial())
        d1_probability <- as.numeric(predict(d1_fit, test, type = "response"))
        d2 <- fit_nested_glmnet(train, test, c("SOFA_24h", cts_vars), "SOFA_24h", base_seed + repeat_id * 
            100L + fold_id)
        d3 <- fit_nested_glmnet(train, test, c("SOFA_24h", cts_vars, viral_vars), c("SOFA_24h", viral_vars), 
            base_seed + 100000L + repeat_id * 100L + fold_id)
        tibble(HumanID = test$HumanID, Mortality28d = test$Mortality28d, Repeat = repeat_id, Outer_fold = fold_id, 
            M0_probability = d0_probability, M1_probability = d1_probability, M2_probability = d2$probability, 
            M3_probability = d3$probability, M2_lambda_min = d2$lambda, M3_lambda_min = d3$lambda, M2_selected_CTS_n = d2$selected_CTS_n, 
            M3_selected_CTS_n = d3$selected_CTS_n)
    }) %>% bind_rows()
}) %>% bind_rows() %>% arrange(Repeat, Outer_fold, HumanID)

if (nrow(oof_predictions) != cohort_n * n_repeats || any(count(oof_predictions, Repeat, HumanID)$n != 
    1L) || any(!is.finite(as.matrix(oof_predictions[, grep("_probability$", names(oof_predictions))])))) {
    stop("Repeated outer held-out predictions failed completeness/uniqueness checks")
}

repeat_performance <- lapply(seq_len(n_repeats), function(repeat_id) {
    repeat_predictions <- oof_predictions %>% filter(Repeat == repeat_id)
    bind_rows(lapply(model_definitions$Model, function(model) {
        probability <- repeat_predictions[[paste0(model, "_probability")]]
        tibble(Repeat = repeat_id, Model = model, AUC = calc_auc(repeat_predictions$Mortality28d, probability), 
            Brier = mean((repeat_predictions$Mortality28d - probability)^2))
    }))
}) %>% bind_rows()

performance <- repeat_performance %>% group_by(Model) %>% summarise(N = cohort_n, Deaths = cohort_deaths, 
    Survivors = cohort_survivors, Repeats = n(), Outer_folds_per_repeat = n_outer_folds, Mean_CV_AUC = mean(AUC), 
    AUC_Q25 = quantile(AUC, 0.25), AUC_Q75 = quantile(AUC, 0.75), AUC_IQR = IQR(AUC), SD_CV_AUC = sd(AUC), 
    Mean_CV_Brier = mean(Brier), Brier_Q25 = quantile(Brier, 0.25), Brier_Q75 = quantile(Brier, 0.75), 
    Brier_IQR = IQR(Brier), SD_CV_Brier = sd(Brier), .groups = "drop") %>% left_join(model_definitions, 
    by = "Model") %>% select(Model, Display, Analysis_role, Feature_space, Algorithm, everything())

contrast_definitions <- tibble(Contrast = c("M1 - M0", "M2 - M0", "M3 - M1", "M3 - M2"), New_model = c("M1", 
    "M2", "M3", "M3"), Reference_model = c("M0", "M0", "M1", "M2"), Comparison_role = c("Primary: EBV/HCMV beyond SOFA24", 
    "Supplementary: CTS genes beyond SOFA24", "Primary: CTS genes beyond SOFA24 and EBV/HCMV", "Supplementary: EBV/HCMV beyond SOFA24 and CTS genes"))

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
    standard_error <- sd(x)/sqrt(length(x))
    critical <- qt(0.975, df = length(x) - 1L)
    c(mean = estimate, lower = estimate - critical * standard_error, upper = estimate + critical * standard_error)
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

write.csv(performance, repo_output(repo_path(output_dir, "260829_Step1_D1_performance.csv"), "Figure6A_analysis.R"), 
    row.names = FALSE)

write.csv(contrasts, repo_output(repo_path(output_dir, "260829_Step1_D1_contrasts.csv"), "Figure6A_analysis.R"), 
    row.names = FALSE)

write.csv(oof_predictions, repo_output(repo_path(output_dir, "260829_Step1_D1_oof_predictions.csv"), 
    "Figure6A_analysis.R"), row.names = FALSE)

get_performance <- function(model, variable) {
    performance[[variable]][performance$Model == model]
}

get_contrast <- function(contrast, variable) {
    contrasts[[variable]][contrasts$Contrast == contrast]
}

primary_delta_auc <- get_contrast("M1 - M0", "Mean_paired_DeltaAUC")

primary_delta_brier <- get_contrast("M1 - M0", "Mean_paired_DeltaBrier")

primary_ci_lower <- get_contrast("M1 - M0", "DeltaAUC_95CI_lower")

primary_ci_upper <- get_contrast("M1 - M0", "DeltaAUC_95CI_upper")

primary_ci_crosses_zero <- primary_ci_lower <= 0 && primary_ci_upper >= 0

if (primary_delta_auc > 0 && primary_ci_crosses_zero) {
    viral_interpretation <- paste0("At ICU admission, addition of EBV and HCMV to SOFA-24h modestly improved ", 
        "discrimination for 28-day mortality, although the paired confidence interval included zero.")
} else if (primary_delta_auc > 0) {
    viral_interpretation <- paste0("At ICU admission, addition of EBV and HCMV to SOFA-24h improved ", 
        "discrimination for 28-day mortality.")
} else {
    viral_interpretation <- paste0("At ICU admission, addition of EBV and HCMV to SOFA-24h did not improve ", 
        "discrimination for 28-day mortality.")
}

cts_delta_auc <- get_contrast("M3 - M1", "Mean_paired_DeltaAUC")

cts_ci_lower <- get_contrast("M3 - M1", "DeltaAUC_95CI_lower")

cts_ci_upper <- get_contrast("M3 - M1", "DeltaAUC_95CI_upper")

if (cts_delta_auc <= 0.01 || (cts_ci_lower <= 0 && cts_ci_upper >= 0)) {
    cts_interpretation <- paste0("Addition of CTS-associated genes provided little further discrimination ", 
        "beyond the clinical and viral features.")
} else {
    cts_interpretation <- paste0("Addition of CTS-associated genes further improved discrimination beyond ", 
        "the clinical and viral features.")
}

results_paragraph <- paste0(viral_interpretation, " Mean repeat-level outer held-out AUC increased from ", 
    sprintf("%.3f", get_performance("M0", "Mean_CV_AUC")), " for SOFA-24h alone to ", sprintf("%.3f", 
        get_performance("M1", "Mean_CV_AUC")), " after addition of EBV and HCMV (paired Delta AUC ", 
    sprintf("%+.3f", primary_delta_auc), ", 95% CI ", sprintf("%.3f", primary_ci_lower), " to ", sprintf("%.3f", 
        primary_ci_upper), "; paired Delta Brier ", sprintf("%+.3f", primary_delta_brier), ") (Figure 6A). ", 
    cts_interpretation, " The combined model had an AUC of ", sprintf("%.3f", get_performance("M3", "Mean_CV_AUC")), 
    " (paired Delta AUC versus SOFA-24h plus EBV/HCMV ", sprintf("%+.3f", cts_delta_auc), ") (Figure 6A and Supplementary Figure).")

performance_lines <- vapply(model_definitions$Model, function(model) {
    sprintf("- %s (%s): mean AUC %.3f (IQR %.3f-%.3f); mean Brier %.3f (IQR %.3f-%.3f).", model, model_definitions$Display[model_definitions$Model == 
        model], get_performance(model, "Mean_CV_AUC"), get_performance(model, "AUC_Q25"), get_performance(model, 
        "AUC_Q75"), get_performance(model, "Mean_CV_Brier"), get_performance(model, "Brier_Q25"), get_performance(model, 
        "Brier_Q75"))
}, character(1))

contrast_lines <- vapply(contrast_definitions$Contrast, function(contrast) {
    sprintf("- %s: paired Delta AUC %+.3f (IQR %.3f to %.3f; 95%% CI %.3f to %.3f); paired Delta Brier %+.3f (IQR %.3f to %.3f).", 
        contrast, get_contrast(contrast, "Mean_paired_DeltaAUC"), get_contrast(contrast, "DeltaAUC_Q25"), 
        get_contrast(contrast, "DeltaAUC_Q75"), get_contrast(contrast, "DeltaAUC_95CI_lower"), get_contrast(contrast, 
            "DeltaAUC_95CI_upper"), get_contrast(contrast, "Mean_paired_DeltaBrier"), get_contrast(contrast, 
            "DeltaBrier_Q25"), get_contrast(contrast, "DeltaBrier_Q75"))
}, character(1))

summary_lines <- c("# Figure 6A: Day-1 admission prediction analysis", "", "## Common Day-1 cohort", 
    "", sprintf("- N: %d", cohort_n), sprintf("- Deaths by 28 days: %d", cohort_deaths), sprintf("- Survivors at 28 days: %d", 
        cohort_survivors), "- Required data: D1 transcriptome, Mortality28d, SOFA_24h, all 18 D1 CTS genes, D1 EBV, and D1 HCMV.", 
    "- EBV and HCMV are D1_mi_HHV.4 and D1_mi_HCMV, already transformed as log2(microbial mass + 1).", 
    "", "## Cross-validation design", "", "- 20 repeats x stratified 5-fold outer CV; all four models use identical outer folds.", 
    "- Each repeat-level AUC and Brier score combines all five outer held-out folds.", "- M2/M3 use stratified 5-fold inner CV within every outer training fold to select lambda.min.", 
    "- M2: SOFA unpenalized and 18 CTS genes penalized.", "- M3: SOFA, EBV, and HCMV unpenalized and 18 CTS genes penalized.", 
    "- Delta Brier is new minus reference; negative values favor the new model.", "", "## Model performance", 
    "", performance_lines, "", "## Paired contrasts", "", contrast_lines, "", "## Main-text interpretation", 
    "", results_paragraph, "", "## Placement", "", "- Main Figure 6A: M0-M3 mean AUC and IQR.", "- Supplementary Figure: paired Delta AUC and Delta Brier distributions.", 
    "- CTS-associated models are supplementary comparators.")

writeLines(summary_lines, repo_output(repo_path(output_dir, "260829_Step1_D1_summary.md"), "Figure6A_analysis.R"), 
    useBytes = TRUE)

message(sprintf(paste0("Completed repeated nested CV: n=%d, deaths=%d, survivors=%d; ", "mean AUC M0/M1/M2/M3 = %.3f/%.3f/%.3f/%.3f; M1-M0 Delta AUC = %+.3f"), 
    cohort_n, cohort_deaths, cohort_survivors, get_performance("M0", "Mean_CV_AUC"), get_performance("M1", 
        "Mean_CV_AUC"), get_performance("M2", "Mean_CV_AUC"), get_performance("M3", "Mean_CV_AUC"), primary_delta_auc))

