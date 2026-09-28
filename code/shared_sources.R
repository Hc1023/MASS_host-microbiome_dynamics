# Module scores shared across analyses. CTS labels are read from the original input snapshot.
module_overlay <- function(x) {
 z <- readRDS('Outputs/Module_scores.rds')
 for(tp in c('D1','D4','D7')) for(m in rownames(z))
  x[[paste0(tp,'_mod_',m)]] <- as.numeric(z[m,match(x[[tp]],colnames(z))])
 x
}
