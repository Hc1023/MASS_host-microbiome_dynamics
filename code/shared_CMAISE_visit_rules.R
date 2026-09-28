# Measurement-level eligibility. Nominal Day 1 is elapsed day 0.
# Original audit is preserved; raw collection timestamps are not available here.
cmaise_visit_rules <- function() {
 a <- read.csv('Inputs/CMAISE_original_visit_audit.csv',stringsAsFactors=FALSE)
 cols <- c('HumanID','SampleID','Timepoint','TimeDays','TimeSinceD1','hospital_days','status','hospital_death','Mortality28d','event_time')
 a <- unique(a[cols]);stopifnot(!anyDuplicated(a[c('HumanID','Timepoint')]))
 a$Original_hospital_days <- a$hospital_days
 a$hospital_days[a$HumanID %in% c('JinHua_53','SSR_1')] <- 5
 a$event_time <- pmin(a$hospital_days,28)
 a$Timing_precision <- 'Common origin confirmed by investigator; corrected hospital duration; exact timestamps unavailable'
 stopifnot(!anyNA(a$Mortality28d), all(a$Mortality28d == as.integer(a$hospital_death==1 & a$hospital_days<=28)))
 a$Eligible <- !a$SampleID %in% c('SSR_144_d3','SSR_183_d3')
 a$Reason <- ifelse(a$Eligible,'Included: available binary endpoint and specified visit', 'Excluded Day-3 measurement after recorded event time')
 # Only these measurements are excluded; all other visits remain independently eligible.
 stopifnot(identical(sort(a$SampleID[!a$Eligible]),c('SSR_144_d3','SSR_183_d3')))
 write.csv(a,'Outputs/CMAISE_visit_eligibility.csv',row.names=FALSE)
 write.csv(a[!a$Eligible,],'Outputs/CMAISE_measurement_exclusion_log.csv',row.names=FALSE)
 a
}
cmaise_keep_visit <- function(ids,tp) {
 a <- cmaise_visit_rules();idx <- match(paste(ids,tp),paste(a$HumanID,a$Timepoint))
 !is.na(idx) & a$Eligible[idx]
}
