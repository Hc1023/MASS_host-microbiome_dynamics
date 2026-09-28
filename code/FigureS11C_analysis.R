# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(pROC)
})

project_dir <- "."

input_file <- repo_path(project_dir, "Inputs/260827_SRR_meta_model.rdata")

step4_dir <- repo_path(project_dir, "Outputs")

prediction_file <- "Outputs/Figure6D_analysis_Step4_CMAISE_validation_predictions.csv"

output_dir <- repo_path(project_dir, "Outputs")

repo_dir(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "260901_Step4_2_CMAISE_subgroups"

input_env <- new.env(parent = emptyenv())

load(repo_input(input_file), envir = input_env)

cmaise <- as_tibble(input_env$meta_model)

predictions <- read.csv(repo_input(prediction_file), check.names = FALSE)

if (!"lung_infection" %in% names(cmaise) || !is.logical(cmaise$lung_infection) || anyNA(cmaise$lung_infection)) {
    stop("CMAISE lung_infection must be a complete logical variable")
}

predictions <- predictions %>% left_join(cmaise %>% select(HumanID, lung_infection), by = "HumanID")

if (anyNA(predictions$lung_infection) || any(count(predictions, Timepoint, HumanID)$n != 1L)) {
    stop("Step4 prediction-to-subgroup mapping failed")
}

timepoints <- c("D3", "D5")

model_levels <- paste0("M", 0:3)

subgroup_definitions <- c(`Lung infection` = TRUE, `Non-lung infection` = FALSE)

roc_metrics <- function(y, score) {
    roc <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
    list(roc = roc, auc = as.numeric(pROC::auc(roc)), lower = ci[1], upper = ci[3])
}

get_score <- function(data, model) {
    if (model == "M0") 
        data$SOFA_contemporaneous
    else data[[paste0("Probability_", model)]]
}

subgroup_performance <- bind_rows(lapply(timepoints, function(timepoint) {
    bind_rows(Map(function(subgroup, flag) {
        dat <- predictions %>% filter(Timepoint == timepoint, lung_infection == flag)
        y <- dat$Mortality28d
        if (n_distinct(y) != 2L) 
            stop("Both outcomes required in ", timepoint, " ", subgroup)
        bind_rows(lapply(model_levels, function(model) {
            score <- get_score(dat, model)
            roc <- roc_metrics(y, score)
            if (model == "M0") {
                brier <- intercept <- slope <- NA_real_
            }
            else {
                lp <- dat[[paste0("LP_", model)]]
                brier <- mean((y - score)^2)
                intercept <- unname(coef(glm(y ~ 1, family = binomial(), offset = lp))[1])
                slope <- unname(coef(glm(y ~ lp, family = binomial()))[2])
            }
            tibble(Timepoint = timepoint, Subgroup = subgroup, Model = model, N = nrow(dat), Deaths = sum(y), 
                Survivors = sum(y == 0), AUC = roc$auc, AUC_CI_lower = roc$lower, AUC_CI_upper = roc$upper, 
                Brier = brier, Calibration_intercept = intercept, Calibration_slope = slope)
        }))
    }, names(subgroup_definitions), subgroup_definitions))
}))

comparison_pairs <- tribble(~New_model, ~Reference_model, "M1", "M0", "M2", "M0", "M3", "M0", "M2", "M1", 
    "M3", "M2", "M3", "M1")

subgroup_comparisons <- bind_rows(lapply(timepoints, function(timepoint) {
    bind_rows(Map(function(subgroup, flag) {
        dat <- predictions %>% filter(Timepoint == timepoint, lung_infection == flag)
        bind_rows(Map(function(new_model, reference_model) {
            rn <- pROC::roc(dat$Mortality28d, get_score(dat, new_model), levels = c(0, 1), direction = "<", 
                quiet = TRUE)
            rr <- pROC::roc(dat$Mortality28d, get_score(dat, reference_model), levels = c(0, 1), direction = "<", 
                quiet = TRUE)
            test <- pROC::roc.test(rn, rr, paired = TRUE, method = "delong", conf.int = TRUE)
            auc_new <- as.numeric(pROC::auc(rn))
            auc_reference <- as.numeric(pROC::auc(rr))
            tibble(Timepoint = timepoint, Subgroup = subgroup, Contrast = paste(new_model, "-", reference_model), 
                New_model = new_model, Reference_model = reference_model, N = nrow(dat), Deaths = sum(dat$Mortality28d), 
                AUC_new = auc_new, AUC_reference = auc_reference, DeltaAUC = auc_new - auc_reference, 
                DeltaAUC_CI_lower = unname(test$conf.int[1]), DeltaAUC_CI_upper = unname(test$conf.int[2]), 
                DeLong_p = test$p.value)
        }, comparison_pairs$New_model, comparison_pairs$Reference_model))
    }, names(subgroup_definitions), subgroup_definitions))
}))

write.csv(subgroup_performance, repo_output(repo_path(output_dir, paste0(prefix, "_performance.csv")), 
    "FigureS11C_analysis.R"), row.names = FALSE)

write.csv(subgroup_comparisons, repo_output(repo_path(output_dir, paste0(prefix, "_comparisons.csv")), 
    "FigureS11C_analysis.R"), row.names = FALSE)

write.csv(predictions %>% select(HumanID, Timepoint, Mortality28d, lung_infection), repo_output(repo_path(output_dir, 
    paste0(prefix, "_cohort.csv")), "FigureS11C_analysis.R"), row.names = FALSE)

