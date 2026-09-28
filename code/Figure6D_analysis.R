# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
source("code/shared_CMAISE_visit_rules.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(pROC)
    library(glmnet)
})

project_dir <- "."

input_file <- repo_path(project_dir, "Inputs/260827_SRR_meta_model.rdata")

frozen_file <- repo_path(project_dir, "Inputs/260901_Step3_D4_signature_frozen_model.rds")

output_dir <- repo_path(project_dir, "Outputs")

repo_dir(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "260901_Step4_CMAISE_validation"

model_levels <- c("M1", "M2", "M3")

if (!repo_exists(input_file) || !repo_exists(frozen_file)) stop("Missing CMAISE input or Step3 frozen model")

input_env <- new.env(parent = emptyenv())

load(repo_input(input_file), envir = input_env)

stopifnot(all(c("meta_model", "model_prepare_info") %in% ls(input_env)))

cmaise <- as_tibble(input_env$meta_model)

cmaise_info <- input_env$model_prepare_info

frozen <- readRDS(repo_input(frozen_file))

if (nrow(cmaise) != 402L || anyDuplicated(cmaise$HumanID)) stop("Expected 402 unique CMAISE patients")

if (!identical(cmaise_info$expression_scale, "TMM-normalized log2-CPM (prior.count = 1)")) stop("Unexpected expression scale")

if (!all(model_levels %in% names(frozen$full_fits))) stop("Step3 frozen bundle lacks M1-M3")

frozen_models <- frozen$full_fits[model_levels]

coefficient_tables <- bind_rows(lapply(frozen_models, `[[`, "coefficient_table"))

selected_counts <- coefficient_tables %>% filter(Gene_source %in% c("CTS", "trajectory")) %>% count(Model, 
    name = "Selected_gene_n")

if (!identical(selected_counts$Model, model_levels) || any(selected_counts$Selected_gene_n < 1L)) stop("Invalid frozen gene sets")

availability <- cmaise_info$expression_availability

stopifnot(all(c("HumanID", "has_expression_D3", "has_expression_D5") %in% names(availability)))

cmaise <- cmaise %>% left_join(availability, by = "HumanID")

if (anyNA(cmaise$has_expression_D3) || anyNA(cmaise$has_expression_D5)) stop("Expression availability mapping failed")

map_feature <- function(feature, timepoint) {
    if (feature == "SOFA_96h") 
        return(paste0("SOFA_", timepoint))
    if (grepl("^D4_(CTSg|DYNg)_", feature)) 
        return(sub("^D4_", paste0(timepoint, "_"), feature))
    stop("Cannot map frozen feature: ", feature)
}

timepoints <- c("D3", "D5")

# Preserve the original common-cohort eligibility after removing the independent sensitivity fit.
all_terms <- unique(c(
  coefficient_tables$Feature[coefficient_tables$Feature != "(Intercept)"],
  frozen$feature_metadata$Feature[frozen$feature_metadata$Gene_source == "CTS"]
))

mapped_features <- crossing(Timepoint = timepoints, Frozen_feature = all_terms) %>% rowwise() %>% mutate(CMAISE_feature = map_feature(Frozen_feature, 
    Timepoint)) %>% ungroup()

missing_features <- setdiff(mapped_features$CMAISE_feature, names(cmaise))

if (length(missing_features)) stop("CMAISE lacks frozen feature(s): ", paste(missing_features, collapse = ", "))

audit_one <- function(timepoint) {
    required <- mapped_features %>% filter(Timepoint == timepoint) %>% pull(CMAISE_feature)
    expr_flag <- cmaise[[paste0("has_expression_", timepoint)]]
    masks <- list(rep(TRUE, nrow(cmaise)), !is.na(cmaise$Mortality28d), !is.na(cmaise$Mortality28d) & 
        expr_flag, !is.na(cmaise$Mortality28d) & expr_flag & !is.na(cmaise[[paste0("SOFA_", timepoint)]]), 
        !is.na(cmaise$Mortality28d) & expr_flag & rowSums(is.na(cmaise[, required, drop = FALSE])) == 
            0L)
    steps <- c("CMAISE patients", "Mortality28d available", paste0(timepoint, " transcriptome available"), 
        paste0(timepoint, " contemporaneous SOFA available"), paste0(timepoint, " common M1-M3 model cohort"))
    Map(function(step, keep) tibble(Timepoint = timepoint, Audit_step = step, N = sum(keep), Deaths = sum(as.integer(as.character(cmaise$Mortality28d[keep])) == 
        1L, na.rm = TRUE), Survivors = sum(as.integer(as.character(cmaise$Mortality28d[keep])) == 0L, 
        na.rm = TRUE)), steps, masks) %>% bind_rows() %>% mutate(Removed_from_previous = c(NA_integer_, 
        -diff(N)))
}

data_audit <- bind_rows(lapply(timepoints, audit_one))

feature_coverage <- lapply(timepoints, function(timepoint) bind_rows(lapply(model_levels, function(model) {
    terms <- coefficient_tables %>% filter(Model == model, Gene_source %in% c("CTS", "trajectory"))
    mapped <- vapply(terms$Feature, map_feature, character(1), timepoint = timepoint)
    tibble(Timepoint = timepoint, Model = model, Frozen_gene_n = nrow(terms), CTS_gene_n = sum(terms$Gene_source == 
        "CTS"), Trajectory_gene_n = sum(terms$Gene_source == "trajectory"), Present_gene_n = sum(mapped %in% 
        names(cmaise)), Present_gene_coverage_pct = 100 * Present_gene_n/Frozen_gene_n, Complete_gene_patient_n = sum(cmaise[[paste0("has_expression_", 
        timepoint)]] & rowSums(is.na(cmaise[, mapped, drop = FALSE])) == 0L))
}))) %>% bind_rows()

if (any(feature_coverage$Present_gene_coverage_pct != 100)) stop("Incomplete frozen-feature coverage")

score_model <- function(data, timepoint, model) {
    fit_object <- frozen_models[[model]]
    terms <- fit_object$coefficient_table
    intercept <- terms$Coefficient[terms$Feature == "(Intercept)"]
    terms <- terms %>% filter(Feature != "(Intercept)") %>% mutate(CMAISE_feature = vapply(Feature, map_feature, 
        character(1), timepoint = timepoint))
    lp <- rep(intercept, nrow(data))
    for (i in seq_len(nrow(terms))) lp <- lp + terms$Coefficient[i] * data[[terms$CMAISE_feature[i]]]
    if (inherits(fit_object$fit, "glm")) {
        newdata <- data.frame(matrix(nrow = nrow(data), ncol = 0))
        for (feature in fit_object$feature_order) newdata[[feature]] <- data[[map_feature(feature, timepoint)]]
        native_lp <- as.numeric(predict(fit_object$fit, newdata = newdata, type = "link"))
    }
    else {
        newx <- matrix(NA_real_, nrow(data), length(fit_object$feature_order), dimnames = list(NULL, 
            fit_object$feature_order))
        for (feature in fit_object$feature_order) newx[, feature] <- data[[map_feature(feature, timepoint)]]
        native_lp <- as.numeric(predict(fit_object$fit, newx = newx, s = fit_object$lambda_min, type = "link"))
    }
    qc <- max(abs(lp - native_lp))
    if (!is.finite(qc) || qc > 1e-10) 
        stop("Score reconstruction failed for ", model, " ", timepoint)
    list(lp = lp, probability = plogis(lp), qc = qc)
}

score_timepoint <- function(timepoint) {
    required <- mapped_features %>% filter(Timepoint == timepoint) %>% pull(CMAISE_feature)
    data <- cmaise %>% filter(.data[[paste0("has_expression_", timepoint)]], !is.na(Mortality28d), if_all(all_of(required), 
        ~!is.na(.x))) %>% filter(cmaise_keep_visit(HumanID,timepoint)) %>% arrange(HumanID)
    out <- data.frame(HumanID = data$HumanID, Timepoint = timepoint, Mortality28d = as.integer(as.character(data$Mortality28d)), 
        SOFA_contemporaneous = data[[paste0("SOFA_", timepoint)]])
    for (model in model_levels) {
        score <- score_model(data, timepoint, model)
        out[[paste0("LP_", model)]] <- score$lp
        out[[paste0("Probability_", model)]] <- score$probability
        out[[paste0("QC_", model)]] <- score$qc
    }
    out
}

predictions <- bind_rows(lapply(timepoints, score_timepoint))

roc_metrics <- function(y, score) {
    roc <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
    list(roc = roc, auc = as.numeric(pROC::auc(roc)), lower = ci[1], upper = ci[3])
}

performance <- lapply(timepoints, function(timepoint) {
    dat <- predictions %>% filter(Timepoint == timepoint)
    y <- dat$Mortality28d
    sofa <- roc_metrics(y, dat$SOFA_contemporaneous)
    baseline <- tibble(Timepoint = timepoint, Model = "M0", Model_label = "Contemporaneous SOFA", N = nrow(dat), 
        Deaths = sum(y), Survivors = sum(y == 0), AUC = sofa$auc, AUC_CI_lower = sofa$lower, AUC_CI_upper = sofa$upper, 
        Brier = NA_real_, Calibration_intercept = NA_real_, Calibration_slope = NA_real_)
    molecular <- bind_rows(lapply(model_levels, function(model) {
        probability <- dat[[paste0("Probability_", model)]]
        lp <- dat[[paste0("LP_", model)]]
        roc <- roc_metrics(y, probability)
        tibble(Timepoint = timepoint, Model = model, Model_label = c(M1 = "CTS elastic-net model", M2 = "Trajectory model", 
            M3 = "Combined model")[[model]], N = nrow(dat), Deaths = sum(y), Survivors = sum(y == 0), 
            AUC = roc$auc, AUC_CI_lower = roc$lower, AUC_CI_upper = roc$upper, Brier = mean((y - probability)^2), 
            Calibration_intercept = unname(coef(glm(y ~ 1, family = binomial(), offset = lp))[1]), Calibration_slope = unname(coef(glm(y ~ 
                lp, family = binomial()))[2]))
    }))
    bind_rows(baseline, molecular)
}) %>% bind_rows()

comparison_pairs <- tribble(~New_model, ~Reference_model, "M1", "M0", "M2", "M0", "M3", "M0", "M2", "M1", 
    "M3", "M2", "M3", "M1")

comparisons <- bind_rows(lapply(timepoints, function(timepoint) {
    dat <- predictions %>% filter(Timepoint == timepoint)
    bind_rows(Map(function(new_model, reference_model) {
        score <- function(model) if (model == "M0") 
            dat$SOFA_contemporaneous
        else dat[[paste0("Probability_", model)]]
        rn <- pROC::roc(dat$Mortality28d, score(new_model), levels = c(0, 1), direction = "<", quiet = TRUE)
        rr <- pROC::roc(dat$Mortality28d, score(reference_model), levels = c(0, 1), direction = "<", 
            quiet = TRUE)
        test <- pROC::roc.test(rn, rr, paired = TRUE, method = "delong", conf.int = TRUE)
        auc_new <- as.numeric(pROC::auc(rn))
        auc_reference <- as.numeric(pROC::auc(rr))
        tibble(Timepoint = timepoint, Contrast = paste(new_model, "-", reference_model), New_model = new_model, 
            Reference_model = reference_model, N = nrow(dat), Deaths = sum(dat$Mortality28d), AUC_new = auc_new, 
            AUC_reference = auc_reference, DeltaAUC = auc_new - auc_reference, DeltaAUC_CI_lower = unname(test$conf.int[1]), 
            DeltaAUC_CI_upper = unname(test$conf.int[2]), DeLong_p = test$p.value)
    }, comparison_pairs$New_model, comparison_pairs$Reference_model))
}))

write.csv(data_audit, repo_output(repo_path(output_dir, paste0(prefix, "_data_audit.csv")), "Figure6D_analysis.R"), 
    row.names = FALSE)

write.csv(feature_coverage, repo_output(repo_path(output_dir, paste0(prefix, "_feature_coverage.csv")), 
    "Figure6D_analysis.R"), row.names = FALSE)

write.csv(predictions, repo_output(repo_path(output_dir, paste0(prefix, "_predictions.csv")), "Figure6D_analysis.R"), 
    row.names = FALSE)

write.csv(performance, repo_output(repo_path(output_dir, paste0(prefix, "_performance.csv")), "Figure6D_analysis.R"), 
    row.names = FALSE)

write.csv(comparisons, repo_output(repo_path(output_dir, paste0(prefix, "_comparisons.csv")), "Figure6D_analysis.R"), 
    row.names = FALSE)

write.csv(coefficient_tables, repo_output(repo_path(output_dir, paste0(prefix, "_frozen_coefficients_applied.csv")), 
    "Figure6D_analysis.R"), row.names = FALSE)

performance_lines <- apply(performance, 1, function(x) if (x[["Model"]] == "M0") sprintf("- %s %s: AUC %.3f (95%% CI %.3f-%.3f); n=%s, deaths=%s.", 
    x[["Timepoint"]], x[["Model"]], as.numeric(x[["AUC"]]), as.numeric(x[["AUC_CI_lower"]]), as.numeric(x[["AUC_CI_upper"]]), 
    x[["N"]], x[["Deaths"]]) else sprintf("- %s %s: AUC %.3f (95%% CI %.3f-%.3f); Brier %.3f; calibration intercept %.2f; calibration slope %.2f; n=%s, deaths=%s.", 
    x[["Timepoint"]], x[["Model"]], as.numeric(x[["AUC"]]), as.numeric(x[["AUC_CI_lower"]]), as.numeric(x[["AUC_CI_upper"]]), 
    as.numeric(x[["Brier"]]), as.numeric(x[["Calibration_intercept"]]), as.numeric(x[["Calibration_slope"]]), 
    x[["N"]], x[["Deaths"]]))

comparison_lines <- apply(comparisons, 1, function(x) sprintf("- %s %s: DeltaAUC %+.3f (95%% CI %.3f to %.3f); DeLong p=%.4f.", 
    x[["Timepoint"]], x[["Contrast"]], as.numeric(x[["DeltaAUC"]]), as.numeric(x[["DeltaAUC_CI_lower"]]), 
    as.numeric(x[["DeltaAUC_CI_upper"]]), as.numeric(x[["DeLong_p"]])))


write.csv(bind_rows(lapply(timepoints,function(tp){z<-predictions[predictions$Timepoint==tp,];data.frame(Timepoint=tp,Audit_step="Eligible measurement cohort",N=nrow(z),Deaths=sum(z$Mortality28d),Survivors=sum(z$Mortality28d==0))})),"Outputs/CMAISE_final_cohort_counts.csv",row.names=FALSE)
