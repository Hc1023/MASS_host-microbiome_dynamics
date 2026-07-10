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
  "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main/Outputs/260625_feature_sel"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load("Inputs/1616_meta_model.rdata")

fi_boot <- as.integer(Sys.getenv("FI_BOOT", "100"))
fi_n_perm <- as.integer(Sys.getenv("FI_N_PERM", "1"))
fi_top_prop <- as.numeric(Sys.getenv("FI_TOP_PROP", "0.25"))
fi_seed <- as.integer(Sys.getenv("FI_SEED", "123"))
glmnet_alpha <- as.numeric(Sys.getenv("FI_GLMNET_ALPHA", "0.5"))
svm_cost <- as.numeric(Sys.getenv("FI_SVM_COST", "1"))
xgb_nrounds <- as.integer(Sys.getenv("FI_XGB_NROUNDS", "100"))
rf_num_trees <- as.integer(Sys.getenv("FI_RF_NUM_TREES", "500"))

parse_env_list <- function(var, default) {
  value <- Sys.getenv(var, default)
  strsplit(value, ",", fixed = TRUE)[[1]] %>%
    trimws() %>%
    discard(~.x == "")
}

selected_feature_spaces <- parse_env_list("FI_FEATURE_SPACES", "CTS + microbe")
selected_algorithms <- parse_env_list(
  "FI_ALGORITHMS",
  "randomForest,xgboost,svm,glmnet"
)

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
) %>%
  filter(feature_space %in% selected_feature_spaces)

importance_algorithms <- tibble::tribble(
  ~algorithm,      ~algorithm_label,
  "randomForest", "Random Forest",
  "xgboost",      "XGBoost",
  "svm",          "SVM-RBF",
  "glmnet",       "Elastic Net"
) %>%
  filter(algorithm %in% selected_algorithms)

clean_feature_names <- function(x, d) {
  case_when(
    x == "APACHEII_24h" ~ "APACHE II",
    str_detect(x, paste0("^", d, "_CTSg_")) ~
      sub(paste0("^", d, "_CTSg_"), "", x),
    str_detect(x, paste0("^", d, "_mi2?_")) ~
      sub(paste0("^", d, "_mi2?_"), "", x),
    TRUE ~ x
  ) %>%
    str_replace_all("\\.", "-") %>%
    recode(
      "Influenza-A" = "Influenza A",
      "HHV-4" = "EBV",
      "Bacf" = "Bacteria/Fungi"
    )
}

get_feature_type <- function(x, d) {
  case_when(
    x == "APACHEII_24h" ~ "APACHEII",
    str_detect(x, paste0("^", d, "_CTSg")) ~ "CTSg",
    str_detect(x, paste0("^", d, "_mi")) ~ "Microbe",
    TRUE ~ "Other"
  )
}

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

scale_train <- function(train_x) {
  center <- vapply(train_x, mean, numeric(1), na.rm = TRUE)
  scale_val <- vapply(train_x, sd, numeric(1), na.rm = TRUE)
  scale_val[!is.finite(scale_val) | scale_val == 0] <- 1
  train_scaled <- sweep(as.matrix(train_x), 2, center, "-")
  train_scaled <- sweep(train_scaled, 2, scale_val, "/")
  list(train = train_scaled, center = center, scale = scale_val)
}

scale_new <- function(new_x, center, scale_val) {
  new_scaled <- sweep(as.matrix(new_x), 2, center, "-")
  sweep(new_scaled, 2, scale_val, "/")
}

extract_yes_prob <- function(pred, positive_class = "Yes") {
  prob <- attr(pred, "probabilities")
  if (is.null(prob)) {
    return(rep(NA_real_, length(pred)))
  }
  if (positive_class %in% colnames(prob)) {
    return(prob[, positive_class])
  }
  prob[, ncol(prob)]
}

fit_importance_model <- function(train_x, train_y, algorithm) {
  nzv <- caret::nearZeroVar(train_x)
  if (length(nzv) > 0) {
    train_x <- train_x[, -nzv, drop = FALSE]
  }
  if (ncol(train_x) == 0) {
    stop("All features removed by near-zero variance filtering")
  }

  feature_names <- colnames(train_x)
  y_num <- as.integer(train_y == "Yes")

  if (algorithm == "glmnet") {
    scaled <- scale_train(train_x)
    fit <- cv.glmnet(
      x = scaled$train,
      y = y_num,
      family = "binomial",
      alpha = glmnet_alpha,
      type.measure = "auc",
      nfolds = 3
    )
    predict_fun <- function(new_x) {
      new_scaled <- scale_new(
        new_x[, feature_names, drop = FALSE],
        scaled$center,
        scaled$scale
      )
      as.numeric(predict(
        fit,
        newx = new_scaled,
        s = "lambda.min",
        type = "response"
      ))
    }
    return(list(
      predict = predict_fun,
      features_used = feature_names,
      tune = paste0("alpha=", glmnet_alpha, "; lambda=", signif(fit$lambda.min, 4))
    ))
  }

  if (algorithm == "xgboost") {
    dtrain <- xgb.DMatrix(data = as.matrix(train_x), label = y_num)
    fit <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = 3,
        eta = 0.1,
        subsample = 0.8,
        colsample_bytree = 0.8
      ),
      data = dtrain,
      nrounds = xgb_nrounds,
      verbose = 0
    )
    predict_fun <- function(new_x) {
      as.numeric(predict(fit, newdata = as.matrix(new_x[, feature_names, drop = FALSE])))
    }
    return(list(
      predict = predict_fun,
      features_used = feature_names,
      tune = paste0("nrounds=", xgb_nrounds, "; max_depth=3; eta=0.1")
    ))
  }

  if (algorithm == "randomForest") {
    p <- ncol(train_x)
    mtry <- max(1, floor(sqrt(p)))
    fit <- ranger::ranger(
      x = train_x,
      y = train_y,
      probability = TRUE,
      num.trees = rf_num_trees,
      mtry = mtry,
      min.node.size = 5,
      classification = TRUE,
      seed = fi_seed
    )
    predict_fun <- function(new_x) {
      predict(fit, data = new_x[, feature_names, drop = FALSE])$predictions[, "Yes"]
    }
    return(list(
      predict = predict_fun,
      features_used = feature_names,
      tune = paste0("num.trees=", rf_num_trees, "; mtry=", mtry, "; min.node.size=5")
    ))
  }

  if (algorithm == "svm") {
    scaled <- scale_train(train_x)
    sigma <- tryCatch(
      as.numeric(kernlab::sigest(scaled$train, frac = 1)),
      error = function(e) c(0.01, 0.05, 0.1)
    )
    gamma <- sort(unique(pmax(sigma, .Machine$double.eps)))[2]
    if (is.na(gamma)) {
      gamma <- sort(unique(pmax(sigma, .Machine$double.eps)))[1]
    }
    fit <- svm(
      x = scaled$train,
      y = train_y,
      kernel = "radial",
      gamma = gamma,
      cost = svm_cost,
      probability = TRUE,
      scale = FALSE
    )
    predict_fun <- function(new_x) {
      new_scaled <- scale_new(
        new_x[, feature_names, drop = FALSE],
        scaled$center,
        scaled$scale
      )
      pred <- predict(fit, new_scaled, probability = TRUE)
      extract_yes_prob(pred)
    }
    return(list(
      predict = predict_fun,
      features_used = feature_names,
      tune = paste0("kernel=radial; gamma=", signif(gamma, 4), "; cost=", svm_cost)
    ))
  }

  stop("Unknown algorithm: ", algorithm)
}

make_bootstrap_split <- function(y, max_tries = 100) {
  class_idx <- split(seq_along(y), y)
  all_idx <- seq_along(y)

  for (i in seq_len(max_tries)) {
    train_idx <- unlist(lapply(class_idx, function(idx) {
      sample(idx, length(idx), replace = TRUE)
    }), use.names = FALSE)
    test_idx <- setdiff(all_idx, unique(train_idx))

    if (length(test_idx) > 0 && length(unique(y[test_idx])) == 2) {
      return(list(train = train_idx, test = test_idx))
    }
  }

  split_idx <- caret::createDataPartition(y, p = 0.7, list = FALSE)
  list(train = as.integer(split_idx), test = setdiff(all_idx, as.integer(split_idx)))
}

permute_feature_auc <- function(test_x, test_y, fit_res, feature, n_perm = 1) {
  scores <- rep(NA_real_, n_perm)
  for (i in seq_len(n_perm)) {
    perm_x <- test_x
    perm_x[[feature]] <- sample(perm_x[[feature]])
    scores[i] <- calc_auc(test_y, fit_res$predict(perm_x))
  }
  mean(scores, na.rm = TRUE)
}

run_feature_importance <- function(B = fi_boot, seed = fi_seed) {
  set.seed(seed)
  importance_rows <- list()
  model_rows <- list()
  row_id <- 1
  model_id <- 1

  for (task_idx in seq_len(nrow(forecast_tasks))) {
    task <- forecast_tasks[task_idx, ]
    task_df <- get_task_data(task$day, task$horizon_days)

    for (iter in seq_len(B)) {
      split <- make_bootstrap_split(task_df$Outcome)
      train_df <- task_df[split$train, , drop = FALSE]
      test_df <- task_df[split$test, , drop = FALSE]

      for (space_idx in seq_len(nrow(feature_spaces))) {
        feature_space <- feature_spaces$feature_space[space_idx]
        feature_cols <- get_feature_cols(task_df, task$day, feature_space)
        train_x <- train_df[, feature_cols, drop = FALSE]
        test_x <- test_df[, feature_cols, drop = FALSE]

        for (alg_idx in seq_len(nrow(importance_algorithms))) {
          algorithm <- importance_algorithms$algorithm[alg_idx]
          message(
            "Feature importance task=", task$task_id,
            " iter=", iter,
            " feature=", feature_space,
            " algorithm=", algorithm
          )

          fit_res <- tryCatch(
            fit_importance_model(train_x, train_df$Outcome, algorithm),
            error = function(e) list(error = conditionMessage(e))
          )

          if (!is.null(fit_res$error)) {
            model_rows[[model_id]] <- data.frame(
              iteration = iter,
              task_id = task$task_id,
              task_label = task$task_label,
              day = task$day,
              horizon_days = task$horizon_days,
              algorithm = algorithm,
              algorithm_label = importance_algorithms$algorithm_label[alg_idx],
              feature_space = feature_space,
              feature_space_label = feature_spaces$feature_space_label[space_idx],
              baseline_auc = NA_real_,
              n_total = nrow(task_df),
              n_train = nrow(train_df),
              n_test = nrow(test_df),
              n_features_started = length(feature_cols),
              n_features_used = NA_integer_,
              tune = NA_character_,
              error = fit_res$error,
              stringsAsFactors = FALSE
            )
            model_id <- model_id + 1
            next
          }

          baseline_prob <- fit_res$predict(test_x)
          baseline_auc <- calc_auc(test_df$Outcome, baseline_prob)

          model_rows[[model_id]] <- data.frame(
            iteration = iter,
            task_id = task$task_id,
            task_label = task$task_label,
            day = task$day,
            horizon_days = task$horizon_days,
            algorithm = algorithm,
            algorithm_label = importance_algorithms$algorithm_label[alg_idx],
            feature_space = feature_space,
            feature_space_label = feature_spaces$feature_space_label[space_idx],
            baseline_auc = baseline_auc,
            n_total = nrow(task_df),
            n_train = nrow(train_df),
            n_test = nrow(test_df),
            n_features_started = length(feature_cols),
            n_features_used = length(fit_res$features_used),
            tune = fit_res$tune,
            error = NA_character_,
            stringsAsFactors = FALSE
          )
          model_id <- model_id + 1

          for (feature in fit_res$features_used) {
            permuted_auc <- permute_feature_auc(
              test_x = test_x,
              test_y = test_df$Outcome,
              fit_res = fit_res,
              feature = feature,
              n_perm = fi_n_perm
            )
            importance_rows[[row_id]] <- data.frame(
              iteration = iter,
              task_id = task$task_id,
              task_label = task$task_label,
              day = task$day,
              horizon_days = task$horizon_days,
              algorithm = algorithm,
              algorithm_label = importance_algorithms$algorithm_label[alg_idx],
              feature_space = feature_space,
              feature_space_label = feature_spaces$feature_space_label[space_idx],
              feature = feature,
              feature_label = clean_feature_names(feature, task$day),
              feature_type = get_feature_type(feature, task$day),
              baseline_auc = baseline_auc,
              permuted_auc = permuted_auc,
              importance_auc_drop = baseline_auc - permuted_auc,
              n_perm = fi_n_perm,
              stringsAsFactors = FALSE
            )
            row_id <- row_id + 1
          }
        }
      }
    }
  }

  importance_long <- bind_rows(importance_rows) %>%
    group_by(iteration, task_id, algorithm, feature_space) %>%
    mutate(
      rank_in_iteration = rank(-importance_auc_drop, ties.method = "first"),
      top_n_cutoff = max(1, ceiling(n() * fi_top_prop)),
      top_feature = rank_in_iteration <= top_n_cutoff & importance_auc_drop > 0
    ) %>%
    ungroup()

  model_qc <- bind_rows(model_rows)

  importance_summary <- importance_long %>%
    group_by(
      task_id, task_label, day, horizon_days,
      algorithm, algorithm_label,
      feature_space, feature_space_label,
      feature, feature_label, feature_type
    ) %>%
    summarise(
      n_iter = sum(!is.na(importance_auc_drop)),
      median_importance = median(importance_auc_drop, na.rm = TRUE),
      mean_importance = mean(importance_auc_drop, na.rm = TRUE),
      q25_importance = quantile(importance_auc_drop, 0.25, na.rm = TRUE),
      q75_importance = quantile(importance_auc_drop, 0.75, na.rm = TRUE),
      positive_importance_frequency = mean(importance_auc_drop > 0, na.rm = TRUE),
      top_feature_frequency = mean(top_feature, na.rm = TRUE),
      median_baseline_auc = median(baseline_auc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(task_id, feature_space, algorithm, desc(median_importance), feature_label)

  list(
    importance_long = importance_long,
    importance_summary = importance_summary,
    model_qc = model_qc
  )
}

outputs <- run_feature_importance()

write.csv(
  outputs$importance_long,
  file.path(out_dir, "feature_importance_permutation_long.csv"),
  row.names = FALSE
)
write.csv(
  outputs$importance_summary,
  file.path(out_dir, "feature_importance_permutation_summary.csv"),
  row.names = FALSE
)
write.csv(
  outputs$model_qc,
  file.path(out_dir, "feature_importance_model_qc.csv"),
  row.names = FALSE
)
saveRDS(
  outputs,
  file.path(out_dir, "feature_importance_permutation_outputs.rds")
)

message("Feature importance outputs written to: ", out_dir)
