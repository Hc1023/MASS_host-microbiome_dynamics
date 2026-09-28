suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
out<-"Outputs/CMAISE_joint_models";dir.create("Outputs/plot_components",recursive=TRUE,showWarnings=FALSE)
for(cohort in c("overall","lung","nonlung")){
 p<-read.csv(file.path(out,paste0(cohort,"_predictions.csv")))
 t<-read.csv(file.path(out,paste0(cohort,"_estimates.csv")))
 p$Outcome<-factor(p$Mortality28d,0:1,c("No event","Death <=28 d"))
 t$label<-ifelse(t$FDR<.001,"JM interaction FDR <0.001",sprintf("JM interaction FDR %.3f",t$FDR))
 g<-ggplot(p,aes(TimeDays,pred,color=Outcome,fill=Outcome,group=Outcome))+
 geom_ribbon(aes(ymin=low,ymax=upp),alpha=.16,color=NA)+geom_line(linewidth=.65)+geom_point(size=1.7)+
 geom_label(data=t,aes(x=1,y=Inf,label=label),inherit.aes=FALSE,hjust=0,vjust=1.15,size=2.3,fontface="bold",linewidth=0)+
 facet_wrap(~Module_ID,ncol=4,scales="free_y")+scale_x_continuous(breaks=c(1,3,5),labels=c("D1","D3","D5"))+
 scale_color_manual(values=c("#4575B4","#D73027"),name=NULL)+scale_fill_manual(values=c("#4575B4","#D73027"),name=NULL)+
 scale_y_continuous(expand=expansion(mult=c(.08,.28)))+labs(x=NULL,y="JM fitted GSVA score")+
 theme_bw(base_size=9)+theme(legend.position="top",panel.grid.minor=element_blank(),strip.text=element_blank(),strip.background=element_blank(),plot.margin=margin(8,8,6,8))
 ggsave(paste0("Outputs/plot_components/",cohort,"_JM.pdf"),g,width=12,height=2.5)
}
