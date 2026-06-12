#!/usr/bin/env Rscript
# Extended Data Fig. 10 wrapper

get_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) return(normalizePath(frame$ofile, winslash = "/", mustWork = FALSE))
      ""
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) return(dirname(tail(frame_files, 1)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  normalizePath(".")
}

script_dir <- get_script_dir()
source(file.path(script_dir, "ed10_usa_focal_country_impl.R"))
