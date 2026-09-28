# Run explicitly selected current analysis steps from the repository root.
steps <- c('Figure3C_prepare.R','validate_module_scores.R',
 'prepare_analysis_sources.R','Figure5AB_data.R','Figure5AB_plot.R','validate_original_CTS.R',
 'Figure5C_analysis.R','Figure5C_data.R','Figure5DE_analysis.R',
 'FigureS10ABC_analysis.R','FigureS11D_plot.R','Figure3DE_FigureS5C_analysis.R',
 'Figure3DE_FigureS5C_plot.R','FigureS6_FigureS7_FigureS8_analysis.R',
 'Figure6B_FigureS11A_analysis.R','Figure6BC_FigureS11B_plot.R',
 'Figure4C_FigureS9CD_trajectories.R','FigureS9AB_plot.R','Figure4C_FigureS9CD_joint_models.R',
 'validate_CMAISE_models.R','Figure4C_FigureS9CD_joint_plot.R')

args <- commandArgs(TRUE)
if(!length(args) || "--list" %in% args){cat(paste(steps,collapse="\n"),"\n");quit(status=0)}
selected <- if("--all" %in% args) steps else args
stopifnot(all(selected %in% steps))
for(s in selected){
 message("Running ",s)
 extra <- if(s=="Figure5C_analysis.R")"--refit" else character()
 status<-system2(file.path(R.home("bin"),"Rscript"),c(file.path("code",s),extra),env="LC_ALL=en_US.UTF-8")
 if(status!=0)stop("Failed: ",s)
}
