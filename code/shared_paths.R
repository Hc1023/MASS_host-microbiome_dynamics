# Shared filesystem routing. Input snapshots live once under Inputs/.
repo_path <- function(...) {
  parts <- list(...)
  p <- do.call(base::file.path, parts)
  vapply(p, function(x) {
    if (grepl('(^|/)Inputs/', x)) return(file.path('Inputs', basename(x)))
    if (grepl('(^|/)code/', x)) return(file.path('code', basename(x)))
    if (grepl('\\.[[:alnum:]]+$', basename(x))) {
      if (file.exists(file.path('Inputs', basename(x)))) return(file.path('Inputs', basename(x)))
      return(file.path(if(grepl('\\.(pdf|png|svg|jpg)$',x,ignore.case=TRUE)) 'Outputs/plot_components' else 'Outputs', basename(x)))
    }
    if (grepl('^[0-9]{6}_', basename(x))) return(file.path('Outputs', basename(x)))
    if (grepl('(^|/)Outputs($|/)', x)) return('Outputs')
    if (grepl('(^|/)Figures($|/)', x)) return('Figures')
    if (grepl('(^|/)Inputs($|/)', x)) return('Inputs')
    '.'
  }, character(1), USE.NAMES=FALSE)
}
repo_input <- function(path) {
  if (!is.character(path)) return(path)
  vapply(path, function(p) {
    if (file.exists(p)) return(p)
    b <- basename(p)
    for (d in c('Inputs','Outputs','code','Figures','Outputs/plot_components')) {
      candidate <- file.path(d,b)
      if(file.exists(candidate)) return(candidate)
    }
    file.path('Inputs', b)
  }, character(1), USE.NAMES=FALSE)
}
repo_exists <- function(...) file.exists(repo_input(unlist(list(...), use.names=FALSE)))
repo_output <- function(path, script) {
  b <- basename(path)
  b <- gsub('[^A-Za-z0-9_.+ -]', '_', b)
  dir <- if(grepl('\\.(pdf|png|svg|jpe?g)$',b,ignore.case=TRUE)) 'Outputs/plot_components' else 'Outputs'
  b <- sub("^[0-9]{4,6}_", "", b)
  if (!grepl("^Figure[0-9S]", b)) b <- paste0(tools::file_path_sans_ext(basename(script)), "_", b)

  file.path(dir,b)
}
repo_dir <- function(path, ...) {
  d <- if(grepl('Figures',path)) 'Outputs/plot_components' else if(grepl('Inputs',path)) 'Inputs' else 'Outputs'
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
}
