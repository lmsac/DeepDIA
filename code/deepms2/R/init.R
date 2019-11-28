#' Peptide MS2 Spectrum Prediction
#' v1

local({
  if(!exists('getCurrentScript')) {
    getCurrentScript = function() {
      lastScriptSourced = tail(unlist(lapply(sys.frames(), function(env) env$ofile)), 1)
      
      if (is.null(lastScriptSourced)) {
        # No script sourced, checking invocation through Rscript
        cmdArgs = commandArgs(trailingOnly = FALSE)
        needle = "--file="
        match = grep(needle, cmdArgs)
        if (length(match) > 0) {
          return(normalizePath(sub(needle, "", cmdArgs[match]), winslash = .Platform$file.sep, mustWork = TRUE))
        }
      } else {
        # 'source'd via R console
        return(normalizePath(lastScriptSourced, winslash = .Platform$file.sep, mustWork = TRUE))
      }
    }
  }
  
  currentScript = getCurrentScript()
  
  tempwd = getwd()
  setwd(dirname(currentScript))
  
  message(paste0('Loading script: ', currentScript))
  
  SOURCE_PATHS = list.files(
    path = 'functions',
    pattern = '\\.R',
    recursive = TRUE,
    full.names = TRUE
  )
  SOURCE_PATHS = SOURCE_PATHS[SOURCE_PATHS != paste0('./', basename(currentScript))]
  lapply(SOURCE_PATHS, source)
  
  setwd(tempwd)
})