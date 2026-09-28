suppressPackageStartupMessages(library(JM))
out<-"Outputs/CMAISE_joint_models";paths<-list.files(out,pattern="_fit.rds$",full.names=TRUE)
z<-lapply(paths,function(p){
 x<-readRDS(p);j<-x$jm;d<-x$input
 sd<-unique(d[c("HumanID","event_time","Mortality28d")])
 stopifnot(j$n==nrow(sd),j$N==nrow(d),sum(j$y$d)==sum(as.character(sd$Mortality28d)=="1"),max(abs(exp(j$y$logT)-sd$event_time))<1e-10,all(is.finite(vcov(j))),min(eigen(vcov(j),symmetric=TRUE,only.values=TRUE)$values)>0,j$convergence==0)
 # JM::initial.surv constructs a final start-stop segment at each last observation.
 # A measurement at the discharge instant produces a zero-length segment in
 # that initialization-only Cox fit, not a removed final likelihood observation.
 last<-d[!duplicated(d$HumanID,fromLast=TRUE),]
 data.frame(Model=basename(p),N=j$n,Measurements=j$N,Events=sum(j$y$d),Convergence=j$convergence,Minimum_covariance_eigenvalue=min(eigen(vcov(j),symmetric=TRUE,only.values=TRUE)$values),Initialization_zero_length_intervals=sum(last$TimeSinceD1==last$event_time),All_measurements_retained=TRUE,Warnings=paste(unique(x$warnings),collapse="; "))
})
write.csv(do.call(rbind,z),file.path(out,"model_acceptance.csv"),row.names=FALSE)
cat("All 12 models: final likelihood includes all participants, events and measurements; finite positive-definite covariance; convergence zero.\n")
