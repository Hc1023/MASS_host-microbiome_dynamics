# Original Table1 categorical test: two-sided Fisher exact test.
e <- new.env();load("Inputs/1211_metadata.rdata",e)
d <- e$meta
d$Pneumonia <- ifelse(d$PneumoniaType=="1.CAP","CAP",
  ifelse(d$PneumoniaType %in% c("2.HAP","3.VAP"),"NP",NA_character_))
stopifnot(nrow(d)==417,!anyNA(d$Pneumonia))
tab <- table(factor(d$Pneumonia,levels=c("CAP","NP")),d$Mortality28d)
stopifnot(all(tab==matrix(c(202,63,124,28),nrow=2)))
result <- data.frame(Pneumonia=rownames(tab),Total_N=rowSums(tab),Survival_N=tab[,"0"],Mortality_N=tab[,"1"])
result$Total_percent <- 100*result$Total_N/nrow(d)
result$Survival_percent <- 100*result$Survival_N/sum(tab[,"0"])
result$Mortality_percent <- 100*result$Mortality_N/sum(tab[,"1"])
result$Fisher_P <- fisher.test(tab)$p.value
write.csv(result,"Outputs/Table1_pneumonia.csv",row.names=FALSE)
