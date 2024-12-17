library(parallel)

#' Set the Number of Cores for Parallel Processing
#'
#' This function sets the number of cores to be used for parallel processing.
#' If the requested number exceeds the available cores, it will default to the 16 cores.
#'
#' @param n Integer specifying the number of cores. Default is 16.
#' @export
set_cores <- function(n = 16) {
  if (n > detectCores()) {
    warning("Number of cores exceeds available cores. Setting to maximum available.")
    n <- detectCores()
  }
  options(scuRMST.cores = n)
}

#' Get the Current Number of Cores for Parallel Processing
#'
#' This function retrieves the current number of cores set for parallel processing.
#'
#' @return Integer indicating the number of cores.
#' @export
get_cores <- function() {
  return(getOption("scuRMST.cores", 16))
}
