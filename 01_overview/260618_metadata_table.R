rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

project_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome"
dynamics_dir <- file.path(project_dir, "MASS_host-microbiome_dynamics")
mortality_dir <- file.path(project_dir, "MASS_mortality-main")
metadata_dir <- "/Users/huangsisi/workspace/MASS/AI_MASS/MASS_metadata_database_20260614"

core_path <- file.path(
  mortality_dir,
  "Outputs/260616_metadata/step1_patient_core_1007.csv"
)

out_dir <- file.path(mortality_dir, "Outputs/260618_metadata_table")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_meta_csv <- function(path) {
  x <- read_csv(
    path,
    show_col_types = FALSE,
    locale = locale(encoding = "UTF-8"),
    col_types = cols(.default = col_character())
  )
  names(x) <- str_remove(names(x), "^\\ufeff")
  x
}

num_prefix <- function(x) {
  suppressWarnings(as.integer(str_extract(as.character(x), "^\\d+")))
}

positive01 <- function(x) {
  case_when(
    is.na(x) ~ NA_integer_,
    num_prefix(x) == 0L ~ 0L,
    num_prefix(x) >= 1L ~ 1L,
    str_detect(as.character(x), "有|是|阳性") ~ 1L,
    str_detect(as.character(x), "无|否|阴性") ~ 0L,
    TRUE ~ NA_integer_
  )
}

yes_no_factor <- function(x) {
  factor(as.character(x), levels = c("0", "1"))
}

fmt_n_pct <- function(x, level = "1") {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  n <- sum(x == level)
  sprintf("%d (%.1f)", n, 100 * n / length(x))
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  fmt <- paste0("%.", digits, "f")
  paste0(
    sprintf(fmt, median(x)),
    " (",
    sprintf(fmt, quantile(x, 0.25)),
    ", ",
    sprintf(fmt, quantile(x, 0.75)),
    ")"
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  if (p < 0.1) return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

p_continuous <- function(dat, var) {
  x <- dat[[var]]
  g <- dat$Mortality28d
  keep <- !is.na(x) & !is.na(g)
  if (length(unique(g[keep])) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(as.numeric(x[keep]) ~ g[keep])$p.value)
}

p_continuous_by <- function(dat, var, by) {
  x <- dat[[var]]
  g <- dat[[by]]
  keep <- !is.na(x) & !is.na(g)
  if (length(unique(g[keep])) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(as.numeric(x[keep]) ~ g[keep])$p.value)
}

p_categorical <- function(dat, var) {
  x <- dat[[var]]
  g <- dat$Mortality28d
  keep <- !is.na(x) & !is.na(g)
  tab <- table(x[keep], g[keep])
  if (nrow(tab) < 1 || ncol(tab) < 2) return(NA_real_)
  fisher.test(tab)$p.value
}

p_categorical_by <- function(dat, var, by) {
  x <- dat[[var]]
  g <- dat[[by]]
  keep <- !is.na(x) & !is.na(g)
  tab <- table(x[keep], g[keep])
  if (nrow(tab) < 1 || ncol(tab) < 2) return(NA_real_)
  fisher.test(tab)$p.value
}

add_section <- function(label) {
  tibble(
    Characteristics = label,
    Total = NA_character_,
    Survival = NA_character_,
    Mortality = NA_character_,
    `P value` = NA_character_
  )
}

add_continuous_row <- function(dat, label, var, digits = 1) {
  tibble(
    Characteristics = label,
    Total = fmt_median_iqr(dat[[var]], digits),
    Survival = fmt_median_iqr(dat[[var]][dat$Mortality28d == "0"], digits),
    Mortality = fmt_median_iqr(dat[[var]][dat$Mortality28d == "1"], digits),
    `P value` = fmt_p(p_continuous(dat, var))
  )
}

add_binary_row <- function(dat, label, var) {
  tibble(
    Characteristics = label,
    Total = fmt_n_pct(dat[[var]], "1"),
    Survival = fmt_n_pct(dat[[var]][dat$Mortality28d == "0"], "1"),
    Mortality = fmt_n_pct(dat[[var]][dat$Mortality28d == "1"], "1"),
    `P value` = fmt_p(p_categorical(dat, var))
  )
}

add_continuous_row_by <- function(dat, label, var, by, group1, group2, digits = 1) {
  tibble(
    Characteristics = label,
    Total = fmt_median_iqr(dat[[var]], digits),
    Group1 = fmt_median_iqr(dat[[var]][dat[[by]] == group1], digits),
    Group2 = fmt_median_iqr(dat[[var]][dat[[by]] == group2], digits),
    `P value` = fmt_p(p_continuous_by(dat, var, by))
  )
}

add_binary_row_by <- function(dat, label, var, by, group1, group2) {
  tibble(
    Characteristics = label,
    Total = fmt_n_pct(dat[[var]], "1"),
    Group1 = fmt_n_pct(dat[[var]][dat[[by]] == group1], "1"),
    Group2 = fmt_n_pct(dat[[var]][dat[[by]] == group2], "1"),
    `P value` = fmt_p(p_categorical_by(dat, var, by))
  )
}

add_categorical_block <- function(dat, section_label, var, levels, labels) {
  dat <- dat %>%
    mutate(.category = factor(.data[[var]], levels = levels, labels = labels))

  p_value <- fmt_p(p_categorical(dat, ".category"))

  rows <- purrr::map2_dfr(levels(dat$.category), levels(dat$.category), function(level, label) {
    if (!any(dat$.category == level, na.rm = TRUE)) return(NULL)
    tibble(
      Characteristics = paste0("    ", label),
      Total = fmt_n_pct(dat$.category, level),
      Survival = fmt_n_pct(dat$.category[dat$Mortality28d == "0"], level),
      Mortality = fmt_n_pct(dat$.category[dat$Mortality28d == "1"], level),
      `P value` = NA_character_
    )
  })

  bind_rows(
    tibble(
      Characteristics = section_label,
      Total = NA_character_,
      Survival = NA_character_,
      Mortality = NA_character_,
      `P value` = p_value
    ),
    rows
  )
}

add_categorical_block_by <- function(dat, section_label, var, levels, labels, by, group1, group2) {
  dat <- dat %>%
    mutate(.category = factor(.data[[var]], levels = levels, labels = labels))

  p_value <- fmt_p(p_categorical_by(dat, ".category", by))

  rows <- purrr::map2_dfr(levels(dat$.category), levels(dat$.category), function(level, label) {
    if (!any(dat$.category == level, na.rm = TRUE)) return(NULL)
    tibble(
      Characteristics = paste0("    ", label),
      Total = fmt_n_pct(dat$.category, level),
      Group1 = fmt_n_pct(dat$.category[dat[[by]] == group1], level),
      Group2 = fmt_n_pct(dat$.category[dat[[by]] == group2], level),
      `P value` = NA_character_
    )
  })

  bind_rows(
    tibble(
      Characteristics = section_label,
      Total = NA_character_,
      Group1 = NA_character_,
      Group2 = NA_character_,
      `P value` = p_value
    ),
    rows
  )
}

make_mortality_table_excel <- function(dat) {
  dat <- dat %>%
    mutate(across(where(is.factor), droplevels))

  out <- bind_rows(
    add_continuous_row(dat, "Age, years, median (IQR)", "Age", 1),
    add_categorical_block(
      dat,
      "Gender, n (%)",
      "Gender",
      levels = c("F", "M"),
      labels = c("Female", "Male")
    ),
    add_categorical_block(
      dat,
      "Type of pneumonia, n (%)",
      "PneumoniaType",
      levels = c("1.CAP", "2.HAP", "3.VAP", "4.SAP"),
      labels = c("CAP", "HAP", "VAP", "SAP")
    ),
    add_section("Comorbities, n (%)"),
    add_binary_row(dat, "    Diabetes mellitus", "DM"),
    add_binary_row(dat, "    Myocardial infarction", "MI"),
    add_binary_row(dat, "    Chronic pulmonary disease", "COPD"),
    add_binary_row(dat, "    Liver disease", "HepaticImpairment"),
    add_binary_row(dat, "    Chronic kidney disease", "RenalDisease"),
    add_binary_row(dat, "    Solid tumour", "Tumor"),
    add_binary_row(dat, "    Haematologic malignancy", "HM"),
    add_binary_row(dat, "    Connective tissue disease", "ConnectiveTissueDisease"),
    add_binary_row(dat, "    Transplantation", "TransplantHistory"),
    add_section("Laboratory indicators, median (IQR)"),
    add_continuous_row(dat, "    White blood cell (10^9/L)", "WBC", 1),
    add_continuous_row(dat, "    Lymphocyte (10^9/L)", "LymphocyteCount", 1),
    add_continuous_row(dat, "    Neutrophil (10^9/L)", "NeutrophilCount", 1),
    add_continuous_row(dat, "    Platelet (10^9/L)", "PlateletCount", 1),
    add_continuous_row(dat, "    Procalcitonin (ng/ml)", "PCT", 1),
    add_continuous_row(dat, "    C reactive protein (mg/L)", "hsCRP", 1),
    add_continuous_row(dat, "    PaO2/FiO2 (mmHg)", "PaO2_FiO2", 1),
    add_binary_row(dat, "Immunosuppression, n (%)", "Immunosuppression"),
    add_continuous_row(dat, "APACHE II score, median (IQR)", "APACHEII_24h", 1),
    add_continuous_row(dat, "SOFA score, median (IQR)", "SOFA_24h", 0),
    add_binary_row(dat, "Mechanical ventilation, n (%)", "MV")
  )

  names(out) <- c(
    "Characteristics",
    paste0("Total (n = ", nrow(dat), ")"),
    paste0("Survival (n = ", sum(dat$Mortality28d == "0", na.rm = TRUE), ")"),
    paste0("Mortality (n = ", sum(dat$Mortality28d == "1", na.rm = TRUE), ")"),
    "P value"
  )
  out
}

make_cap_hapvap_table_excel <- function(dat) {
  dat <- dat %>%
    mutate(
      PneumoniaGroup = case_when(
        PneumoniaType == "1.CAP" ~ "CAP",
        PneumoniaType %in% c("2.HAP", "3.VAP") ~ "HAP/VAP",
        TRUE ~ NA_character_
      ),
      PneumoniaGroup = factor(PneumoniaGroup, levels = c("CAP", "HAP/VAP"))
    ) %>%
    filter(!is.na(PneumoniaGroup)) %>%
    mutate(across(where(is.factor), droplevels))

  by <- "PneumoniaGroup"
  group1 <- "CAP"
  group2 <- "HAP/VAP"

  out <- bind_rows(
    add_continuous_row_by(dat, "Age, years, median (IQR)", "Age", by, group1, group2, 1),
    add_categorical_block_by(
      dat,
      "Gender, n (%)",
      "Gender",
      levels = c("F", "M"),
      labels = c("Female", "Male"),
      by = by,
      group1 = group1,
      group2 = group2
    ),
    add_binary_row_by(dat, "Mortality, n (%)", "Mortality28d", by, group1, group2),
    add_section("Comorbities, n (%)") %>% rename(Group1 = Survival, Group2 = Mortality),
    add_binary_row_by(dat, "    Diabetes mellitus", "DM", by, group1, group2),
    add_binary_row_by(dat, "    Myocardial infarction", "MI", by, group1, group2),
    add_binary_row_by(dat, "    Chronic pulmonary disease", "COPD", by, group1, group2),
    add_binary_row_by(dat, "    Liver disease", "HepaticImpairment", by, group1, group2),
    add_binary_row_by(dat, "    Chronic kidney disease", "RenalDisease", by, group1, group2),
    add_binary_row_by(dat, "    Solid tumour", "Tumor", by, group1, group2),
    add_binary_row_by(dat, "    Haematologic malignancy", "HM", by, group1, group2),
    add_binary_row_by(dat, "    Connective tissue disease", "ConnectiveTissueDisease", by, group1, group2),
    add_binary_row_by(dat, "    Transplantation", "TransplantHistory", by, group1, group2),
    add_section("Laboratory indicators, median (IQR)") %>% rename(Group1 = Survival, Group2 = Mortality),
    add_continuous_row_by(dat, "    White blood cell (10^9/L)", "WBC", by, group1, group2, 1),
    add_continuous_row_by(dat, "    Lymphocyte (10^9/L)", "LymphocyteCount", by, group1, group2, 1),
    add_continuous_row_by(dat, "    Neutrophil (10^9/L)", "NeutrophilCount", by, group1, group2, 1),
    add_continuous_row_by(dat, "    Platelet (10^9/L)", "PlateletCount", by, group1, group2, 1),
    add_continuous_row_by(dat, "    Procalcitonin (ng/ml)", "PCT", by, group1, group2, 1),
    add_continuous_row_by(dat, "    C reactive protein (mg/L)", "hsCRP", by, group1, group2, 1),
    add_continuous_row_by(dat, "    PaO2/FiO2 (mmHg)", "PaO2_FiO2", by, group1, group2, 1),
    add_binary_row_by(dat, "Immunosuppression, n (%)", "Immunosuppression", by, group1, group2),
    add_continuous_row_by(dat, "APACHE II score, median (IQR)", "APACHEII_24h", by, group1, group2, 1),
    add_continuous_row_by(dat, "SOFA score, median (IQR)", "SOFA_24h", by, group1, group2, 0),
    add_binary_row_by(dat, "Mechanical ventilation, n (%)", "MV", by, group1, group2)
  )

  names(out) <- c(
    "Characteristics",
    paste0("Total (n = ", nrow(dat), ")"),
    paste0("CAP (n = ", sum(dat$PneumoniaGroup == group1, na.rm = TRUE), ")"),
    paste0("HAP/VAP (n = ", sum(dat$PneumoniaGroup == group2, na.rm = TRUE), ")"),
    "P value"
  )
  out
}


# ---- Step 1. Read source files ----

metadata_core <- read_meta_csv(core_path)
stopifnot(nrow(metadata_core) == 1007)
stopifnot(n_distinct(metadata_core$HumanID) == nrow(metadata_core))

medical_history <- read_meta_csv(file.path(metadata_dir, "11.Medical_history_info.csv"))
experiment <- read_meta_csv(file.path(metadata_dir, "6.Experiment_info.csv"))
sofa_apa <- read_meta_csv(file.path(metadata_dir, "7.Sofa_apa_info.csv"))

load(file.path(dynamics_dir, "Inputs/1211_metadata.rdata")) # creates df_long and meta
study_ids <- df_long %>%
  distinct(HumanID) %>%
  pull(HumanID)
stopifnot(length(study_ids) == 417)

# ---- Step 2. Prepare variables used by the old table1.R ----

comorbidity_tbl <- medical_history %>%
  transmute(
    HumanID,
    DM = positive01(Diabetes),
    MI = positive01(MyocardialInfarction),
    COPD = positive01(ChronicLungDiseaseOrAsthma),
    HepaticImpairment = positive01(HepaticImpairment),
    RenalDisease = positive01(ModerateToSevereRenalImpairment),
    Tumor = positive01(SolidTumor),
    HM = as.integer(
      replace_na(positive01(Leukemia), 0L) == 1L |
        replace_na(positive01(Lymphoma), 0L) == 1L
    ),
    ConnectiveTissueDisease = positive01(RheumaticOrConnectiveTissueDisease),
    TransplantHistory = positive01(TransplantHistory)
  ) %>%
  distinct(HumanID, .keep_all = TRUE)

experiment_tbl <- experiment %>%
  transmute(
    HumanID,
    WBC = as.numeric(WBC),
    LymphocyteCount = as.numeric(LymphocyteCount),
    NeutrophilCount = as.numeric(NeutrophilCount),
    PlateletCount = as.numeric(PlateletCount),
    PCT = as.numeric(PCT),
    hsCRP = as.numeric(hsCRP)
  ) %>%
  distinct(HumanID, .keep_all = TRUE)

oxygen_tbl <- sofa_apa %>%
  transmute(
    HumanID,
    PaO2_FiO2 = as.numeric(PaO2_24h) / (as.numeric(FiO2_24h) / 100)
  ) %>%
  mutate(PaO2_FiO2 = if_else(is.finite(PaO2_FiO2), PaO2_FiO2, NA_real_)) %>%
  distinct(HumanID, .keep_all = TRUE)

metadata_table_source <- metadata_core %>%
  transmute(
    HumanID,
    StudyCohort417 = as.integer(HumanID %in% study_ids),
    cohort = if_else(HumanID %in% study_ids, "study_417", "whole_only"),
    Mortality28d = DeathWithin28DaysAfterEnrollment,
    Age = as.numeric(Age),
    Gender,
    SOFA_24h = as.numeric(SOFA_24h),
    Immunosuppression,
    MV = InvasiveMechanicalVentilation,
    PneumoniaType,
    APACHEII_24h = as.numeric(APACHEII_24h)
  ) %>%
  left_join(comorbidity_tbl, by = "HumanID") %>%
  left_join(experiment_tbl, by = "HumanID") %>%
  left_join(oxygen_tbl, by = "HumanID") %>%
  mutate(
    Mortality28d = yes_no_factor(Mortality28d),
    Gender = factor(Gender),
    PneumoniaType = factor(PneumoniaType),
    across(
      c(
        Immunosuppression,
        MV,
        DM,
        MI,
        COPD,
        HepaticImpairment,
        RenalDisease,
        Tumor,
        HM,
        ConnectiveTissueDisease,
        TransplantHistory
      ),
      yes_no_factor
    )
  )

stopifnot(nrow(metadata_table_source) == 1007)
stopifnot(sum(metadata_table_source$StudyCohort417 == 1L) == 417)

# ---- Step 3. Generate gtsummary tables using only old table1.R variables ----

table_vars <- c(
  "Age",
  "Gender",
  "SOFA_24h",
  "Immunosuppression",
  "MV",
  "DM",
  "MI",
  "COPD",
  "HepaticImpairment",
  "RenalDisease",
  "Tumor",
  "HM",
  "ConnectiveTissueDisease",
  "PneumoniaType",
  "TransplantHistory",
  "WBC",
  "LymphocyteCount",
  "NeutrophilCount",
  "PlateletCount",
  "PCT",
  "hsCRP",
  "APACHEII_24h"
)

binary_vars <- c(
  "Immunosuppression",
  "MV",
  "DM",
  "MI",
  "COPD",
  "HepaticImpairment",
  "RenalDisease",
  "Tumor",
  "HM",
  "ConnectiveTissueDisease",
  "TransplantHistory"
)

binary_value <- purrr::map(binary_vars, ~ "1") |>
  rlang::set_names(binary_vars)

whole_source <- metadata_table_source
study_source <- metadata_table_source %>%
  filter(StudyCohort417 == 1L) %>%
  mutate(across(where(is.factor), droplevels))

whole_excel <- make_mortality_table_excel(whole_source)
study_excel <- make_mortality_table_excel(study_source)
whole_cap_hapvap_excel <- make_cap_hapvap_table_excel(whole_source)
study_cap_hapvap_excel <- make_cap_hapvap_table_excel(study_source)

write_csv(metadata_table_source, file.path(out_dir, "metadata_table_source_1007.csv"))

write.xlsx(
  list(
    "Whole cohort 1007" = whole_excel,
    "Study cohort 417" = study_excel,
    "Whole CAP vs HAPVAP" = whole_cap_hapvap_excel,
    "Study CAP vs HAPVAP" = study_cap_hapvap_excel
  ),
  file = file.path(out_dir, "260618_characteristic_table.xlsx"),
  overwrite = TRUE
)
