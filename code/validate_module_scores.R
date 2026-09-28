library(edgeR);library(limma);library(GSVA)
load('Inputs/1211_metadata.rdata');load('Inputs/1211_transcriptome.rdata')
m<-readRDS('Outputs/Figure3C_prepare_MASS_GSVA_scores.rds');saveRDS(m,'Outputs/Module_scores.rds');old<-as.matrix(read.csv('Inputs/Supplementary_data_6_module_gsva.csv',row.names=1,check.names=FALSE));old<-old[rownames(m),colnames(m)]
a<-read.csv('Inputs/MASS_GSVA_scores_long_with_metadata.csv');b<-read.csv('Outputs/Figure3C_prepare_MASS_GSVA_scores_long_with_metadata.csv');ix<-match(paste(b$SampleID,b$Module),paste(a$SampleID,a$Module));stopifnot(!anyNA(ix));main_diff<-max(abs(a$Score[ix]-b$Score));stopifnot(main_diff<1e-12)
write.csv(data.frame(Comparison=c('Rebuilt versus saved main trajectory','Rebuilt versus legacy subgroup'),Max_abs_difference=c(main_diff,max(abs(m-old)))),'Outputs/Module_score_comparison.csv',row.names=FALSE)
write.csv(m,'Outputs/Module_scores.csv')
dge<-calcNormFactors(DGEList(as.matrix(counts[,colnames(m)])))
e<-voom(dge,plot=FALSE)$E;sym<-gene_attr[rownames(e),'SYMBOL'];ok<-!is.na(sym)&!grepl('^ENSG',sym);e<-e[ok,];rownames(e)<-sym[ok];e<-e[rownames(e)%in%names(table(rownames(e)))[table(rownames(e))==1],]
g<-read.csv('Inputs/selected_genes_by_module_wide.csv');g<-setNames(lapply(g,function(x)unique(trimws(x[!is.na(x)&trimws(x)!='']))),rownames(m))
v<-gsva(gsvaParam(e,g,kcdf='Gaussian',minSize=10,maxSize=500),verbose=FALSE);saveRDS(v,'../MASS_host-microbiome_dynamics_archive/2026-09-17_directory_cleanup/Module_legacy_reconstruction.rds')
write.csv(data.frame(Module=rownames(m),Legacy_reproduction_max_abs=apply(abs(v-old),1,max),Main_vs_legacy_max_abs=apply(abs(m-old),1,max),Main_vs_legacy_median_abs=apply(abs(m-old),1,median)),'Outputs/Module_preprocessing_audit.csv',row.names=FALSE)
