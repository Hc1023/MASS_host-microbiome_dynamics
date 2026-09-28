# Export the necessary refit snapshot. Nonconverged estimates are not publication results.
library(JM)
models <- readRDS("Outputs/Figure5C_joint_models.rds")
rows <- lapply(names(models),function(cts) {
  fit <- models[[cts]]
  d <- expand.grid(TimeDays=c(0,3,6),Mortality28d=c(0,1))
  d$HumanID <- fit$data$HumanID[1]
  d$CTS <- cts
  d$Timepoint <- paste0("D",d$TimeDays+1)
  d$Convergence_code <- fit$convergence
  d$Prediction <- d$SE <- d$CI_low <- d$CI_high <- NA_real_
  if (fit$convergence==0) {
    d$Prediction <- as.numeric(predict(fit,newdata=d,process="Longitudinal",type="Marginal"))
    X <- model.matrix(~TimeDays*Mortality28d,d)
    idx <- paste0("Y.",names(fit$coefficients$betas))
    V <- vcov(fit)[idx,idx,drop=FALSE]
    d$SE <- sqrt(rowSums((X%*%V)*X))
    d$CI_low <- d$Prediction-qnorm(.975)*d$SE
    d$CI_high <- d$Prediction+qnorm(.975)*d$SE
  }
  d$HumanID <- NULL
  d
})
write.csv(do.call(rbind,rows),"Outputs/Figure5C_reconstructed_predictions.csv",row.names=FALSE,na="")
