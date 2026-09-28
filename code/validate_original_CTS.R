# Audit original stored labels without modifying them or retraining a classifier.
source('code/shared_CTS_prepare.R')
p <- as.matrix(cts_saved[,c('CTS1','CTS2','CTS3')])
label <- as.integer(as.character(cts_saved$CTS))
valid <- complete.cases(p) & !is.na(label)
audit <- cts_saved
audit$Label_probability_mismatch <- NA
audit$Label_probability_mismatch[valid] <- p[cbind(which(valid),label[valid])] < apply(p[valid,,drop=FALSE],1,max)-1e-12
write.csv(audit,'Outputs/CTS_original_label_probability_audit.csv',row.names=FALSE)
e<-new.env();load('Inputs/1616_meta_model.rdata',e)
a<-e$meta_model
stopifnot(sum(a$D1_CTS %in% c(1,2))==266,sum(a$D1_CTS==3,na.rm=TRUE)==110)
states<-read.csv('Outputs/Figure5A_reconstructed_states.csv')
stopifnot(all(table(states$Timepoint)==214),!anyDuplicated(states[c('HumanID','Timepoint')]))
flows<-read.csv('Outputs/Figure5A_reconstructed_flows.csv')
stopifnot(all(tapply(flows$N,flows$Timepoint,sum)==214))
z<-read.csv('Outputs/FigureS11D_mortality_summary.csv')
stopifnot(identical(as.integer(tapply(z$n,z$Timepoint,sum)),c(347L,290L)))
write.csv(data.frame(Check=c('Original D1 group sizes 266/110','Alluvial state and flow conservation','CMAISE D3/D5 counts 347/290'),Pass=TRUE),'Outputs/CTS_original_validation.csv',row.names=FALSE)
cat('Original labels retained. Unresolved label/probability mismatches:',sum(audit$Label_probability_mismatch,na.rm=TRUE),'\n')
