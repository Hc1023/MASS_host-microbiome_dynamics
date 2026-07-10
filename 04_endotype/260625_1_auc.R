rm(list = ls())

library(tidyverse)
library(e1071)
library(caret)
library(pROC)
library(glmnet)
library(xgboost)

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_host-microbiome_dynamics"
setwd(project_dir)

out_dir <- file.path(
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_auc_sel"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load("Inputs/1616_meta_model.rdata")

n_boot <- as.integer(Sys.getenv("N_BOOT", "20"))
auc_iter <- as.integer(Sys.getenv("AUC_ITER", as.character(n_boot)))
run_auc_models <- Sys.getenv("RUN_AUC_MODELS", "TRUE") == "TRUE"

forecast_tasks <- tibble::tribble(
  ~task_id,  ~task_label,  ~day,  ~horizon_days,
  "D1_M3",  "D1 3-day",   "D1",  3,
  "D1_M7",  "D1 7-day",   "D1",  7,
  "D1_M28", "D1 28-day",  "D1",  28,
  "D4_M28", "D4 28-day",  "D4",  28,
  "D7_M28", "D7 28-day",  "D7",  28
)

feature_spaces <- tibble::tribble(
  ~feature_space,      ~feature_space_label,
  "APACHEII",          "APACHE II",
  "CTS",               "CTS",
  "microbe",           "Microbes",
  "CTS + microbe",     "CTS + Microbes"
)

auc_algorithms <- tibble::tribble(
  ~algorithm,      ~algorithm_label,
  "randomForest", "Random Forest",
  "xgboost",      "XGBoost",
  "svm",          "SVM-RBF",
  "glmnet",       "Elastic Net"
)

get_task_data <- function(day, horizon_days) {
  df <- meta_model %>%
    dplyr::select(
      Mortality28d,
      SurvivalTimeWithin28Days,
      APACHEII_24h,
      starts_with(paste0(day, "_CTSg")),
      starts_with(paste0(day, "_mi"))
    ) %>%
    na.omit()

  if (horizon_days == 28) {
    outcome <- as.integer(as.character(df$Mortality28d))
  } else {
    outcome <- as.integer(df$SurvivalTimeWithin28Days <= horizon_days)
  }

  df$Outcome <- factor(
    ifelse(outcome == 1, "Yes", "No"),
    levels = c("No", "Yes")
  )
  df$sample_id <- rownames(df)
  df
}

get_feature_cols <- function(df, day, feature_space) {
  cts_cols <- grep(paste0("^", day, "_CTSg"), names(df), value = TRUE)
  microbe_cols <- grep(paste0("^", day, "_mi"), names(df), value = TRUE)

  if (feature_space == "APACHEII") {
    return("APACHEII_24h")
  }
  if (feature_space == "CTS") {
    return(cts_cols)
  }
  if (feature_space == "microbe") {
    return(microbe_cols)
  }
  if (feature_space == "CTS + microbe") {
    return(c(cts_cols, microbe_cols))
  }
  stop("Unknown feature space: ", feature_space)
}

calc_auc <- function(y, prob) {
  if (anyNA(prob) || length(unique(y)) < 2 || length(unique(prob)) < 2) {
    return(NA_real_)
  }
  as.numeric(pROC::auc(
    pROC::roc(
      response = y,
      predictor = prob,
      levels = c("No", "Yes"),
      direction = "<",
      quiet = TRUE
    )
  ))
}

scale_train_test <- function(train_x, test_x) {
  center <- vapply(train_x, mean, numeric(1), na.rm = TRUE)
  scale_val <- vapply(train_x, sd, numeric(1), na.rm = TRUE)
  scale_val[!is.finite(scale_val) | scale_val == 0] <- 1

  train_scaled <- sweep(as.matrix(train_x), 2, center, "-")
  train_scaled <- sweep(train_scaled, 2, scale_val, "/")
  test_scaled <- sweep(as.matrix(test_x), 2, center, "-")
  test_scaled <- sweep(test_scaled, 2, scale_val, "/")

  list(train = train_scaled, test = test_scaled)
}

make_inner_folds <- function(y, k = 3) {
  createFolds(y, k = k, returnTrain = FALSE)
}

fit_auc_model <- function(train_x, train_y, test_x, algorithm) {
  nzv <- caret::nearZeroVar(train_x)
  if (length(nzv) > 0) {
    train_x <- train_x[, -nzv, drop = FALSE]
    test_x <- test_x[, colnames(train_x), drop = FALSE]
  }

  if (ncol(train_x) == 0) {
    stop("All features removed by near-zero variance filtering")
  }

  folds <- make_inner_folds(train_y, k = 3)
  y_num <- as.integer(train_y == "Yes")

  if (algorithm == "glmnet") {
    scaled <- scale_train_test(train_x, test_x)

    if (ncol(scaled$train) == 1) {
      train_dat <- data.frame(y = train_y, x = scaled$train[, 1])
      test_dat <- data.frame(x = scaled$test[, 1])
      fit <- glm(y ~ x, data = train_dat, family = binomial)
      pred_prob <- as.numeric(predict(
        fit,
        newdata = test_dat,
        type = "response"
      ))
      return(list(
        prob = pred_prob,
        best_tune = "single_feature_logistic_fallback",
        n_features_used = ncol(train_x)
      ))
    }

    alpha_grid <- c(0, 0.5, 1)
    cv_scores <- sapply(alpha_grid, function(alpha) {
      fit <- cv.glmnet(
        x = scaled$train,
        y = y_num,
        family = "binomial",
        alpha = alpha,
        type.measure = "auc",
        nfolds = 3
      )
      max(fit$cvm, na.rm = TRUE)
    })
    best_alpha <- alpha_grid[which.max(cv_scores)]
    fit <- cv.glmnet(
      x = scaled$train,
      y = y_num,
      family = "binomial",
      alpha = best_alpha,
      type.measure = "auc",
      nfolds = 3
    )
    pred_prob <- as.numeric(predict(
      fit,
      newx = scaled$test,
      s = "lambda.min",
      type = "response"
    ))
    return(list(
      prob = pred_prob,
      best_tune = paste0("alpha=", best_alpha, "; lambda=", signif(fit$lambda.min, 4)),
      n_features_used = ncol(train_x)
    ))
  }

  if (algorithm == "xgboost") {
    grid <- expand.grid(
      nrounds = c(50, 100),
      max_depth = c(2, 3),
      eta = c(0.03, 0.1),
      subsample = 0.8,
      colsample_bytree = 0.8
    )
    dtrain <- xgb.DMatrix(data = as.matrix(train_x), label = y_num)
    fold_idx <- lapply(folds, function(idx) idx)
    cv_scores <- rep(NA_real_, nrow(grid))

    for (i in seq_len(nrow(grid))) {
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = grid$max_depth[i],
        eta = grid$eta[i],
        subsample = grid$subsample[i],
        colsample_bytree = grid$colsample_bytree[i]
      )
      cv <- xgb.cv(
        params = params,
        data = dtrain,
        nrounds = grid$nrounds[i],
        folds = fold_idx,
        verbose = 0
      )
      cv_scores[i] <- max(cv$evaluation_log$test_auc_mean, na.rm = TRUE)
    }

    best <- grid[which.max(cv_scores), , drop = FALSE]
    fit <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = best$max_depth,
        eta = best$eta,
        subsample = best$subsample,
        colsample_bytree = best$colsample_bytree
      ),
      data = dtrain,
      nrounds = best$nrounds,
      verbose = 0
    )
    pred_prob <- as.numeric(predict(fit, newdata = as.matrix(test_x)))
    return(list(
      prob = pred_prob,
      best_tune = paste(names(best), best[1, ], sep = "=", collapse = "; "),
      n_features_used = ncol(train_x)
    ))
  }

  if (algorithm == "randomForest") {
    p <- ncol(train_x)
    grid <- expand.grid(
      mtry = sort(unique(pmax(1, c(floor(sqrt(p)), floor(p / 3), p)))),
      min.node.size = c(1, 5)
    )
    cv_scores <- rep(NA_real_, nrow(grid))

    for (i in seq_len(nrow(grid))) {
      fold_scores <- sapply(folds, function(valid_idx) {
        fit <- ranger::ranger(
          x = train_x[-valid_idx, , drop = FALSE],
          y = train_y[-valid_idx],
          probability = TRUE,
          num.trees = 500,
          mtry = grid$mtry[i],
          min.node.size = grid$min.node.size[i],
          classification = TRUE,
          seed = 123
        )
        prob <- predict(fit, data = train_x[valid_idx, , drop = FALSE])$predictions[, "Yes"]
        calc_auc(train_y[valid_idx], prob)
      })
      cv_scores[i] <- mean(fold_scores, na.rm = TRUE)
    }

    best <- grid[which.max(cv_scores), , drop = FALSE]
    fit <- ranger::ranger(
      x = train_x,
      y = train_y,
      probability = TRUE,
      num.trees = 500,
      mtry = best$mtry,
      min.node.size = best$min.node.size,
      classification = TRUE,
      seed = 123
    )
    pred_prob <- predict(fit, data = test_x)$predictions[, "Yes"]
    return(list(
      prob = pred_prob,
      best_tune = paste(names(best), best[1, ], sep = "=", collapse = "; "),
      n_features_used = ncol(train_x)
    ))
  }

  if (algorithm == "svm") {
    scaled <- scale_train_test(train_x, test_x)
    sigma <- tryCatch(
      as.numeric(kernlab::sigest(scaled$train, frac = 1)),
      error = function(e) c(0.01, 0.05, 0.1)
    )
    gamma_grid <- sort(unique(pmax(sigma, .Machine$double.eps)))
    grid <- expand.grid(gamma = gamma_grid, cost = c(0.25, 1, 4))
    cv_scores <- rep(NA_real_, nrow(grid))

    for (i in seq_len(nrow(grid))) {
      fold_scores <- sapply(folds, function(valid_idx) {
        fit <- svm(
          x = scaled$train[-valid_idx, , drop = FALSE],
          y = train_y[-valid_idx],
          kernel = "radial",
          gamma = grid$gamma[i],
          cost = grid$cost[i],
          probability = TRUE,
          scale = FALSE
        )
        pred <- predict(fit, scaled$train[valid_idx, , drop = FALSE], probability = TRUE)
        prob <- attr(pred, "probabilities")[, "Yes"]
        calc_auc(train_y[valid_idx], prob)
      })
      cv_scores[i] <- mean(fold_scores, na.rm = TRUE)
    }

    best <- grid[which.max(cv_scores), , drop = FALSE]
    fit <- svm(
      x = scaled$train,
      y = train_y,
      kernel = "radial",
      gamma = best$gamma,
      cost = best$cost,
      probability = TRUE,
      scale = FALSE
    )
    pred <- predict(fit, scaled$test, probability = TRUE)
    pred_prob <- attr(pred, "probabilities")[, "Yes"]
    return(list(
      prob = pred_prob,
      best_tune = paste(names(best), best[1, ], sep = "=", collapse = "; "),
      n_features_used = ncol(train_x)
    ))
  }

  stop("Unknown algorithm: ", algorithm)
}

run_auc_benchmark <- function(B = auc_iter, seed = 123) {
  set.seed(seed)
  auc_rows <- list()
  pred_rows <- list()
  row_id <- 1
  pred_id <- 1

  for (task_idx in seq_len(nrow(forecast_tasks))) {
    task <- forecast_tasks[task_idx, ]
    task_df <- get_task_data(task$day, task$horizon_days)

    for (iter in seq_len(B)) {
      split_idx <- createDataPartition(task_df$Outcome, p = 0.7, list = FALSE)
      train_df <- task_df[split_idx, , drop = FALSE]
      test_df <- task_df[-split_idx, , drop = FALSE]

      for (space_idx in seq_len(nrow(feature_spaces))) {
        feature_space <- feature_spaces$feature_space[space_idx]
        feature_cols <- get_feature_cols(task_df, task$day, feature_space)

        train_x <- train_df[, feature_cols, drop = FALSE]
        test_x <- test_df[, feature_cols, drop = FALSE]

        for (alg_idx in seq_len(nrow(auc_algorithms))) {
          algorithm <- auc_algorithms$algorithm[alg_idx]
          message(
            "AUC task=", task$task_id,
            " iter=", iter,
            " feature=", feature_space,
            " algorithm=", algorithm
          )

          fit_res <- tryCatch(
            fit_auc_model(train_x, train_df$Outcome, test_x, algorithm),
            error = function(e) list(
              prob = rep(NA_real_, nrow(test_df)),
              best_tune = NA_character_,
              n_features_used = NA_integer_,
              error = conditionMessage(e)
            )
          )

          if (is.null(fit_res$error)) {
            auc_val <- as.numeric(pROC::auc(
              pROC::roc(
                response = test_df$Outcome,
                predictor = fit_res$prob,
                levels = c("No", "Yes"),
                direction = "<",
                quiet = TRUE
              )
            ))
            err_msg <- NA_character_
          } else {
            auc_val <- NA_real_
            err_msg <- fit_res$error
          }

          auc_rows[[row_id]] <- data.frame(
            iteration = iter,
            task_id = task$task_id,
            task_label = task$task_label,
            day = task$day,
            horizon_days = task$horizon_days,
            algorithm = algorithm,
            algorithm_label = auc_algorithms$algorithm_label[alg_idx],
            feature_space = feature_space,
            feature_space_label = feature_spaces$feature_space_label[space_idx],
            auc = auc_val,
            best_tune = fit_res$best_tune,
            n_total = nrow(task_df),
            n_train = nrow(train_df),
            n_test = nrow(test_df),
            event_rate_train = mean(train_df$Outcome == "Yes"),
            event_rate_test = mean(test_df$Outcome == "Yes"),
            n_features_started = length(feature_cols),
            n_features_used = fit_res$n_features_used,
            error = err_msg,
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1

          pred_rows[[pred_id]] <- data.frame(
            iteration = iter,
            task_id = task$task_id,
            task_label = task$task_label,
            sample_id = test_df$sample_id,
            observed = test_df$Outcome,
            pred_prob = fit_res$prob,
            algorithm = algorithm,
            algorithm_label = auc_algorithms$algorithm_label[alg_idx],
            feature_space = feature_space,
            feature_space_label = feature_spaces$feature_space_label[space_idx],
            stringsAsFactors = FALSE
          )
          pred_id <- pred_id + 1
        }
      }
    }
  }

  auc_long <- bind_rows(auc_rows)
  pred_long <- bind_rows(pred_rows)

  auc_summary <- auc_long %>%
    group_by(task_id, task_label, day, horizon_days, algorithm, algorithm_label,
             feature_space, feature_space_label) %>%
    summarise(
      n_iter = sum(!is.na(auc)),
      median_auc = median(auc, na.rm = TRUE),
      mean_auc = mean(auc, na.rm = TRUE),
      q25_auc = quantile(auc, 0.25, na.rm = TRUE),
      q75_auc = quantile(auc, 0.75, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      n_total = first(n_total),
      n_train = median(n_train),
      n_test = median(n_test),
      event_rate_train = median(event_rate_train),
      event_rate_test = median(event_rate_test),
      n_features_started = first(n_features_started),
      n_features_used = median(n_features_used, na.rm = TRUE),
      n_errors = sum(!is.na(error)),
      .groups = "drop"
    ) %>%
    mutate(
      auc_label = sprintf("%.3f (%.3f-%.3f)", median_auc, q25_auc, q75_auc)
    )

  auc_feature_summary <- auc_long %>%
    group_by(task_id, task_label, day, horizon_days,
             feature_space, feature_space_label) %>%
    summarise(
      n_values = sum(!is.na(auc)),
      median_auc = median(auc, na.rm = TRUE),
      mean_auc = mean(auc, na.rm = TRUE),
      q25_auc = quantile(auc, 0.25, na.rm = TRUE),
      q75_auc = quantile(auc, 0.75, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      n_total = first(n_total),
      n_errors = sum(!is.na(error)),
      .groups = "drop"
    ) %>%
    mutate(
      auc_label = sprintf("%.3f (%.3f-%.3f)", median_auc, q25_auc, q75_auc)
    )

  task_levels <- forecast_tasks$task_id
  feature_levels <- feature_spaces$feature_space
  algorithm_levels <- auc_algorithms$algorithm

  auc_plotting_table <- auc_summary %>%
    mutate(
      task_order = match(task_id, task_levels),
      feature_space_order = match(feature_space, feature_levels),
      algorithm_order = match(algorithm, algorithm_levels)
    ) %>%
    arrange(task_order, feature_space_order, algorithm_order)

  list(
    auc_long = auc_long,
    pred_long = pred_long,
    auc_summary = auc_summary,
    auc_feature_summary = auc_feature_summary,
    auc_plotting_table = auc_plotting_table
  )
}

if (run_auc_models) {
  auc_outputs <- run_auc_benchmark(B = auc_iter)
  write.csv(
    auc_outputs$auc_long,
    file.path(out_dir, "auc_validation_long.csv"),
    row.names = FALSE
  )
  write.csv(
    auc_outputs$auc_summary,
    file.path(out_dir, "auc_validation_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    auc_outputs$pred_long,
    file.path(out_dir, "auc_predictions_long.csv"),
    row.names = FALSE
  )
  write.csv(
    auc_outputs$auc_feature_summary,
    file.path(out_dir, "auc_feature_space_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    auc_outputs$auc_plotting_table,
    file.path(out_dir, "auc_plotting_table.csv"),
    row.names = FALSE
  )
  saveRDS(
    auc_outputs,
    file.path(out_dir, "auc_modeling_outputs.rds")
  )
}
