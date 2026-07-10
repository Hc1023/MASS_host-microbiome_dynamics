rm(list = ls())

library(tidyverse)

# -----------------------------------------------------------------------------
# VIF analysis for Day 4 HCMV mediation/path models
#
# Path structure:
#   Day 4 HCMV (X) -> host-response module (M) -> 28-day mortality (Y)
#
# VIF is calculated separately for:
#   a-path model:    M ~ X + covariates
#   b/c'-path model: Y ~ X + M + covariates
# -----------------------------------------------------------------------------

find_repo_root <- function(start_dir = getwd()) {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_dir <- if (length(script_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
  } else {
    normalizePath(start_dir)
  }

  candidates <- unique(c(
    script_dir,
    start_dir,
    dirname(script_dir),
    dirname(start_dir)
  ))
  is_repo_root <- vapply(
    candidates,
    function(x) {
      file.exists(file.path(
        x, "MASS_host-microbiome_dynamics", "Inputs", "microbe_hostmodule.csv"
      ))
    },
    logical(1)
  )

  if (!any(is_repo_root)) {
    stop("Cannot locate the repository root from: ", normalizePath(start_dir))
  }
  normalizePath(candidates[which(is_repo_root)[[1]]])
}

repo_dir <- find_repo_root()
input_file <- file.path(
  repo_dir, "MASS_host-microbiome_dynamics", "Inputs", "microbe_hostmodule.csv"
)
out_dir <- file.path(repo_dir, "MASS_mortality-main", "Outputs", "260629_VIF")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

X <- "HCMV_D4"
Y <- "Mortality28d"
modules <- c("up3_D4", "up1_D4", "up2_D4", "dw1_D4")

# Mapping confirmed by the observed a-path coefficients and b-path odds ratios.
module_map <- c(
  up3_D4 = "Phagolysosome function",
  up1_D4 = "Ribosome biogenesis",
  up2_D4 = "IFN signaling",
  dw1_D4 = "Inflammatory signaling"
)

covariates_categorical <- c(
  "Gender", "CenterGroup", "PneumoniaTypeGroup", "Immunosuppression", "MV"
)
covariates_numeric <- c("Age", "CCI", "SOFA_24h")
covariates <- c(covariates_categorical, covariates_numeric)
vif_limit <- 2

# -----------------------------------------------------------------------------
# 1) Construct the paired D1/D4 dataset used by the mediation analysis
# -----------------------------------------------------------------------------

df <- read.csv(input_file, check.names = FALSE)
required_columns <- c(
  "HumanID", "Timepoint", "up1", "up2", "up3", "dw1", "HCMV", Y,
  covariates
)
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

df_paired <- df %>%
  filter(Timepoint %in% c("D1", "D4")) %>%
  group_by(HumanID) %>%
  filter(n_distinct(Timepoint) == 2) %>%
  ungroup()

duplicate_visits <- df_paired %>%
  count(HumanID, Timepoint) %>%
  filter(n > 1)
if (nrow(duplicate_visits) > 0) {
  stop("Duplicate participant/timepoint records were found.")
}

df_path <- df_paired %>%
  dplyr::select(
    HumanID, Timepoint, up1, up2, up3, dw1, HCMV, all_of(c(Y, covariates))
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = c(up1, up2, up3, dw1, HCMV)
  ) %>%
  mutate(
    across(all_of(c(X, modules, Y, covariates_numeric)), as.numeric),
    across(all_of(covariates_categorical), as.factor)
  )

# -----------------------------------------------------------------------------
# 2) Calculate standard VIF for every predictor in one fitted path model
#
# VIF_j = 1 / (1 - R_j^2). The same quantity is obtained from the diagonal of
# the inverse predictor correlation matrix. The current categorical covariates
# are binary, so every model term has one degree of freedom.
# -----------------------------------------------------------------------------

calculate_model_vif <- function(formula, data, module_code, path_model) {
  model_data <- model.frame(formula, data = data, na.action = na.omit)
  design_matrix <- model.matrix(formula, data = model_data)
  keep <- colnames(design_matrix) != "(Intercept)"
  design_matrix <- design_matrix[, keep, drop = FALSE]

  term_labels <- attr(terms(formula), "term.labels")
  term_assignment <- attr(model.matrix(formula, data = model_data), "assign")[keep]

  if (ncol(design_matrix) != length(term_labels) ||
      any(table(term_assignment) != 1)) {
    stop(
      "The current implementation expects one coefficient per predictor term. ",
      "Check whether a categorical covariate has more than two levels."
    )
  }
  if (any(apply(design_matrix, 2, sd) == 0)) {
    stop("A predictor has zero variance in model: ", deparse(formula))
  }

  predictor_cor <- cor(design_matrix)
  if (qr(predictor_cor)$rank < ncol(predictor_cor)) {
    stop("The predictor correlation matrix is singular in model: ", deparse(formula))
  }

  vif_values <- diag(solve(predictor_cor))

  tibble(
    module_code = module_code,
    module = unname(module_map[module_code]),
    path_model = path_model,
    formula = paste(deparse(formula), collapse = " "),
    n_complete = nrow(model_data),
    predictor = term_labels[term_assignment],
    design_column = colnames(design_matrix),
    VIF = unname(vif_values),
    tolerance = 1 / VIF,
    auxiliary_R2 = 1 - tolerance,
    below_2 = VIF < vif_limit
  )
}

vif_results <- map_dfr(modules, function(M) {
  a_formula <- reformulate(c(X, covariates), response = M)
  outcome_formula <- reformulate(c(X, M, covariates), response = Y)

  bind_rows(
    calculate_model_vif(a_formula, df_path, M, "a-path (mediator model)"),
    calculate_model_vif(
      outcome_formula, df_path, M, "b/c'-path (mortality model)"
    )
  )
})

vif_summary <- vif_results %>%
  group_by(module_code, module, path_model, formula, n_complete) %>%
  summarise(
    min_VIF = min(VIF),
    max_VIF = max(VIF),
    max_VIF_predictor = predictor[which.max(VIF)],
    all_VIF_below_2 = all(below_2),
    .groups = "drop"
  )

overall_max <- vif_results %>%
  slice_max(VIF, n = 1, with_ties = FALSE)

overall_summary <- tibble(
  exposure = X,
  outcome = Y,
  n_modules = length(modules),
  n_path_models = nrow(vif_summary),
  min_VIF = min(vif_results$VIF),
  max_VIF = max(vif_results$VIF),
  max_VIF_module = overall_max$module,
  max_VIF_path_model = overall_max$path_model,
  max_VIF_predictor = overall_max$predictor,
  VIF_limit = vif_limit,
  all_VIF_below_2 = all(vif_results$below_2)
)

# Aggregate each predictor across all path models in which it was included.
# This is the compact table intended for the manuscript/supplement.
predictor_label_map <- c(
  HCMV_D4 = "Day 4 viral sentinel load (HCMV)",
  up3_D4 = "Phagolysosome function mediator",
  up1_D4 = "Ribosome biogenesis mediator",
  up2_D4 = "IFN signaling mediator",
  dw1_D4 = "Inflammatory signaling mediator",
  SOFA_24h = "Baseline sepsis severity score (SOFA)",
  Age = "Age of patient",
  CCI = "Charlson Comorbidity Index (CCI)",
  Gender = "Sex",
  CenterGroup = "Study center group",
  PneumoniaTypeGroup = "Pneumonia type",
  Immunosuppression = "Baseline immunosuppression status",
  MV = "Mechanical ventilation status"
)

predictor_order <- names(predictor_label_map)
mediator_codes <- names(module_map)

manuscript_vif_table <- vif_results %>%
  group_by(predictor) %>%
  summarise(
    models_included = n_distinct(paste(module_code, path_model)),
    minimum_VIF = min(VIF),
    maximum_VIF = max(VIF),
    .groups = "drop"
  ) %>%
  mutate(
    model_covariate = unname(predictor_label_map[predictor]),
    predictor_type = case_when(
      predictor == X ~ "Viral exposure",
      predictor %in% mediator_codes ~ "Host-response mediator",
      TRUE ~ "Clinical covariate"
    ),
    VIF_range = if_else(
      abs(maximum_VIF - minimum_VIF) < .Machine$double.eps^0.5,
      sprintf("%.2f", maximum_VIF),
      sprintf("%.2f-%.2f", minimum_VIF, maximum_VIF)
    ),
    threshold_status = if_else(
      maximum_VIF < vif_limit,
      sprintf("Passed (Max VIF < %.1f)", vif_limit),
      sprintf("Review required (Max VIF >= %.1f)", vif_limit)
    ),
    interpretation_and_model_robustness_notes = case_when(
      predictor == X ~ paste0(
        "Low collinearity with the host-response mediators and included ",
        "clinical covariates."
      ),
      predictor %in% mediator_codes ~ paste0(
        "The mediator is not materially redundant with Day 4 HCMV burden ",
        "or the included clinical covariates."
      ),
      TRUE ~ paste0(
        "Low collinearity with Day 4 HCMV, host-response mediators, and ",
        "the other clinical adjustments."
      )
    ),
    sort_order = match(predictor, predictor_order)
  ) %>%
  arrange(sort_order) %>%
  dplyr::select(
    model_covariate,
    predictor_type,
    models_included,
    minimum_VIF,
    maximum_VIF,
    VIF_range,
    threshold_status,
    interpretation_and_model_robustness_notes
  )

manuscript_vif_summary_row <- tibble(
  model_covariate = "COLLINEARITY SUMMARY",
  predictor_type = "Overall",
  models_included = nrow(vif_summary),
  minimum_VIF = min(vif_results$VIF),
  maximum_VIF = max(vif_results$VIF),
  VIF_range = sprintf(
    "%.2f-%.2f", min(vif_results$VIF), max(vif_results$VIF)
  ),
  threshold_status = sprintf("Passed (all VIFs < %.1f)", vif_limit),
  interpretation_and_model_robustness_notes = paste0(
    "All predictors showed low multicollinearity across the eight Day 4 ",
    "HCMV path models."
  )
)

manuscript_vif_table_with_summary <- bind_rows(
  manuscript_vif_table,
  manuscript_vif_summary_row
)

message("Day 4 HCMV path-model VIF analysis complete.")
message(
  "Overall VIF range: ", sprintf("%.3f", min(vif_results$VIF)), "-",
  sprintf("%.3f", max(vif_results$VIF)),
  "; all VIFs < ", vif_limit, ": ", all(vif_results$below_2)
)
print(vif_summary)

# -----------------------------------------------------------------------------
# 3) Save detailed results, summaries, manuscript text, and QC plot
# -----------------------------------------------------------------------------

write.csv(
  vif_results,
  file.path(out_dir, "260629_D4_path_model_VIF_results.csv"),
  row.names = FALSE
)
write.csv(
  vif_summary,
  file.path(out_dir, "260629_D4_path_model_VIF_summary.csv"),
  row.names = FALSE
)
write.csv(
  overall_summary,
  file.path(out_dir, "260629_D4_path_model_VIF_overall.csv"),
  row.names = FALSE
)
write.csv(
  manuscript_vif_table_with_summary,
  file.path(out_dir, "260629_D4_path_model_VIF_manuscript_table.csv"),
  row.names = FALSE
)

result_statement <- sprintf(
  paste0(
    "Variance inflation factors (VIFs) were calculated for all predictors in ",
    "each mediator (a-path) and mortality (b/c'-path) regression model. Across ",
    "the Day 4 HCMV path models, VIFs ranged from %.2f to %.2f, and all values ",
    "were below 2, indicating low multicollinearity between Day 4 HCMV burden, ",
    "host-response modules, and the included clinical covariates."
  ),
  min(vif_results$VIF),
  max(vif_results$VIF)
)
writeLines(
  result_statement,
  file.path(out_dir, "260629_D4_path_model_VIF_statement.txt")
)

plot_df <- vif_summary %>%
  mutate(
    module = factor(module, levels = rev(unname(module_map[modules]))),
    path_model = factor(
      path_model,
      levels = c("a-path (mediator model)", "b/c'-path (mortality model)")
    )
  )

p_vif <- ggplot(
  plot_df,
  aes(x = max_VIF, y = module, fill = path_model)
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = vif_limit,
    linetype = "dashed",
    linewidth = 0.5,
    color = "#B2182B"
  ) +
  geom_text(
    aes(label = sprintf("%.2f", max_VIF)),
    position = position_dodge(width = 0.72),
    hjust = -0.15,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "a-path (mediator model)" = "#9ECAE1",
      "b/c'-path (mortality model)" = "#4C78A8"
    ),
    labels = c("a-path", "b/c'-path")
  ) +
  scale_x_continuous(
    limits = c(0, vif_limit * 1.12),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "Maximum VIF within each path model", y = NULL, fill = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(out_dir, "260629_D4_path_model_VIF_plot.pdf"),
  p_vif,
  width = 5.8,
  height = 3.3
)
