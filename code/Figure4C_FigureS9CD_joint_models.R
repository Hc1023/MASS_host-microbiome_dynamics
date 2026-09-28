# Researcher-confirmed common time origin and 5-day discharge corrections.
# Original snapshots and v3 results remain unchanged.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(nlme);library(survival);library(JM);library(ggplot2)})
set.seed(20260917)
dir.create("Outputs/CMAISE_joint_models",recursive=TRUE,showWarnings=FALSE)
out <- "Outputs/CMAISE_joint_models"
a <- read.csv("Inputs/CMAISE_original_visit_audit.csv",stringsAsFactors=FALSE)
a <- distinct(a,HumanID,SampleID,Timepoint,TimeDays,TimeSinceD1,hospital_days,hospital_death,Mortality28d,event_time)
a$Original_hospital_days <- a$hospital_days
a$hospital_days[a$HumanID %in% c("JinHua_53","SSR_1")] <- 5
log <- filter(a,hospital_days!=Original_hospital_days) %>% mutate(Reason="Researcher correction: hospital duration 5 days for both participants, all visits; latest instruction supersedes 3/5/5")
write.csv(log,file.path(out,"duration_correction_log.csv"),row.names=FALSE)
a$Followup <- pmin(a$hospital_days,28)
a$Event <- as.integer(a$hospital_death==1 & a$hospital_days<=28)
a$Eligible <- !a$SampleID %in% c("SSR_144_d3","SSR_183_d3")
a$After_followup <- a$TimeSinceD1 > a$Followup
stopifnot(all(a$Event==a$Mortality28d),all(a$Followup>0),!any(a$After_followup[a$Eligible]))
write.csv(a,file.path(out,"visit_audit.csv"),row.names=FALSE)
d <- read.csv("Inputs/260911_CMAISE_scores_long.csv",stringsAsFactors=FALSE)
i <- match(d$SampleID,a$SampleID);stopifnot(!anyNA(i),all(a$Eligible[i]))
d$hospital_days <- a$hospital_days[i];d$event_time <- a$Followup[i]
d$Timepoint <- factor(d$Timepoint,levels=c("D1","D3","D5"))
d$Mortality28d <- factor(d$Mortality28d,levels=0:1)
stopifnot(all(d$TimeSinceD1==d$TimeDays-1),all(d$TimeSinceD1<=d$event_time))
write.csv(d,file.path(out,"scores_long.csv"),row.names=FALSE)
labels <- c(M1="PRR/TLR/NF-kB signaling",M2="Phagolysosomal/autophagy",M3="Type I IFN response",M4="Ribosome biogenesis")
alltests <- list();allpred <- list();counts <- list()
for(cohort in c("overall","lung","nonlung")){
 dat <- if(cohort=="overall")d else filter(d,Subgroup==if(cohort=="lung")"lung infection" else "non-lung infection")
 fits <- list()
 for(id in names(labels)){
  fpath <- file.path(out,paste0(cohort,"_",id,"_fit.rds"))
  if(file.exists(fpath)){cached<-readRDS(fpath);expected<-filter(dat,Module_ID==id) %>% arrange(HumanID,TimeSinceD1);stopifnot(isTRUE(all.equal(cached$input,expected,check.attributes=FALSE)));fits[[id]]<-cached;next}
  x <- filter(dat,Module_ID==id) %>% arrange(HumanID,TimeSinceD1)
  message("Fitting ",cohort," ",id," n=",nrow(x))
  warn <- character()
  ans <- withCallingHandlers({
   lf <- try(lme(Score~TimeSinceD1*Mortality28d,random=~TimeSinceD1|HumanID,data=x,method="REML",control=lmeControl(opt="optim",maxIter=200,msMaxIter=200)),silent=TRUE)
   rs <- "intercept and slope"
   if(inherits(lf,"try-error")){
    rs <- "intercept only (initial slope fit failed)"
    lf <- lme(Score~TimeSinceD1*Mortality28d,random=~1|HumanID,data=x,method="REML",control=lmeControl(opt="optim",maxIter=200,msMaxIter=200))
   }
   sd <- distinct(x,HumanID,event_time,Mortality28d) %>% mutate(event=as.integer(as.character(Mortality28d)))
   stopifnot(nrow(sd)==n_distinct(x$HumanID))
   cf <- coxph(Surv(event_time,event)~1,data=sd,x=TRUE)
   jf <- jointModel(lf,cf,timeVar="TimeSinceD1",method="weibull-PH-aGH")
   retry <- FALSE
   if(!is.null(jf$convergence)&&jf$convergence!=0){
    retry <- TRUE
    jf <- jointModel(lf,cf,timeVar="TimeSinceD1",method="weibull-PH-aGH",control=list(iter.EM=200,iter.qN=2000))
   }
   lt <- summary(jf)$`CoefTable-Long`;et <- summary(jf)$`CoefTable-Event`
   ir <- grep("TimeSinceD1:Mortality28d",rownames(lt));ar <- grep("Assoct|association",rownames(et),ignore.case=TRUE)
   stopifnot(length(ir)==1,length(ar)==1)
   pg <- expand_grid(TimeDays=c(1,3,5),Mortality28d=factor(0:1,levels=0:1)) %>% mutate(TimeSinceD1=TimeDays-1,HumanID=x$HumanID[1],Score=0)
   pr <- as.data.frame(predict(jf,newdata=pg,type="Marginal",interval="confidence",returnData=TRUE))
   tab <- data.frame(Cohort=cohort,Module_ID=id,N=nrow(sd),Measurements=nrow(x),Events=sum(sd$event),Random_structure=rs,Convergence=if(is.null(jf$convergence))NA else jf$convergence,Retry=retry,Slope=lt[ir,"Value"],SE=lt[ir,"Std.Err"],P=lt[ir,"p-value"],Alpha=et[ar,"Value"],Alpha_SE=et[ar,"Std.Err"],Alpha_P=et[ar,"p-value"])
   list(jm=jf,longitudinal=lf,survival=cf,table=tab,predictions=pr,input=x)
  },warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")})
  ans$warnings <- warn;saveRDS(ans,fpath);fits[[id]]<-ans
 }
 tt <- bind_rows(lapply(fits,`[[`,"table")) %>% mutate(FDR=p.adjust(P,"BH"),CI_low=Slope-1.96*SE,CI_high=Slope+1.96*SE,Alpha_CI_low=Alpha-1.96*Alpha_SE,Alpha_CI_high=Alpha+1.96*Alpha_SE)
 pp <- bind_rows(lapply(names(fits),function(id)mutate(fits[[id]]$predictions,Module_ID=id,Cohort=cohort)))
 write.csv(tt,file.path(out,paste0(cohort,"_estimates.csv")),row.names=FALSE)
 write.csv(pp,file.path(out,paste0(cohort,"_predictions.csv")),row.names=FALSE)
 writeLines(unlist(lapply(names(fits),function(id)c(id,fits[[id]]$warnings))),file.path(out,paste0(cohort,"_warnings.txt")))
 alltests[[cohort]]<-tt;allpred[[cohort]]<-pp
 counts[[cohort]]<-dat %>% filter(Module_ID=="M1") %>% count(Timepoint,Mortality28d) %>% mutate(Cohort=cohort)
}
write.csv(bind_rows(alltests),file.path(out,"all_estimates.csv"),row.names=FALSE)
write.csv(bind_rows(allpred),file.path(out,"all_predictions.csv"),row.names=FALSE)
write.csv(bind_rows(counts),file.path(out,"sample_counts.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines(c("Common origin confirmed by researcher.","D1/D3/D5 elapsed time: 0/2/4 days.","JinHua_53 and SSR_1 discharge duration corrected to 5 days at every visit.","Event: recorded in-hospital death within 28 days.","Follow-up: min(hospital duration,28), with live discharge censored at discharge.","Discharge censoring requires a noninformative-censoring assumption conditional on model structure; no post-discharge survival assumed.","Outcome in longitudinal submodel: retrospective descriptive sensitivity analysis, not prospective prediction.","BH across four slope interactions separately within overall/lung/nonlung."),file.path(out,"analysis_specification.txt"))
