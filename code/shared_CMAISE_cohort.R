# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
source("code/shared_CMAISE_visit_rules.R")
prepare_cmaise_cohort <- function(root = ".") {
    score_file <- "Outputs/CMAISE_joint_models/scores_long.csv"
    scores <- read.csv(repo_input(score_file), check.names = FALSE, stringsAsFactors = FALSE)
    rules <- cmaise_visit_rules()
    eligible <- rules[rules$Eligible, c("HumanID", "Timepoint")]
    stopifnot(all(paste(scores$HumanID, scores$Timepoint) %in% paste(eligible$HumanID, eligible$Timepoint)))
    module_labels <- c(M1 = "PRR/TLR/NF-kB signaling", M2 = "Phagolysosomal/autophagy", M3 = "Type I IFN response", 
        M4 = "Ribosome biogenesis")
    old_to_id <- c(`PRR/TLR/TNF/NF-kB signaling` = "M1", `Phagolysosome/autophagy` = "M2", `Type I interferon response` = "M3", 
        `Ribosome biogenesis` = "M4")
    scores$Timepoint <- factor(scores$Timepoint, levels = c("D1", "D3", "D5"))
    scores$Module_ID <- factor(scores$Module_ID, levels = names(module_labels))
    scores$Outcome <- factor(scores$Outcome, levels = c("Survival", "Mortality"))
    stopifnot(!anyDuplicated(scores[c("SampleID", "Module_ID")]))
    list(scores = scores, audit = scores, score_file = score_file, clinical_file = score_file, old_to_id = old_to_id, 
        module_labels = module_labels)
}

assert_cmaise_same_cohort <- function(a, b) {
    cols <- c("HumanID", "SampleID", "Module_ID", "Timepoint", "Subgroup", "Score", "Mortality28d", "Outcome", 
        "event_time", "TimeDays", "TimeSinceD1")
    normalize <- function(d) d %>% dplyr::select(all_of(cols)) %>% mutate(across(where(is.factor), as.character)) %>% 
        arrange(Module_ID, HumanID, TimeDays) %>% as.data.frame()
    stopifnot(isTRUE(all.equal(normalize(a), normalize(b), tolerance = 1e-12, check.attributes = FALSE)))
    invisible(TRUE)
}

