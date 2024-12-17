# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html


library(testthat)
library(devtools)

# packages <- c("survRM2", "mvtnorm", "ggplot2", "MASS", "tidyr", "survival", "nph", "tidyverse",
#               "foreach", "doParallel", "cowplot", "IRdisplay", "rlang", "simtrial", "ggrepel")
# install_and_load <- function(package) {
#   if (!require(package, character.only = TRUE)) {
#     install.packages(package)
#     library(package, character.only = TRUE)
#   }
# }

# invisible(lapply(packages, install_and_load))

# n_cores <- detectCores()
# cluster <- makeCluster(min(16, n_cores))  # Use the minimum of 32 or the number of available cores
# registerDoParallel(cluster)
# invisible(clusterEvalQ(cluster, {
#   library(survRM2)
#   library(mvtnorm)
#   library(cubature)
#   library(survival)
#   library(nph)
#   library(simtrial)
#   library(foreach)
# }))


# test
test_that('Test some functions',{
  result <- tryCatch({ 
set_cores(16)
 set.seed(2024)
lambda_H1 <- 0.9
 HR <- 1.7
 lambda_H0 <- lambda_H1 * HR
 sim_size <- 2000
 r <- 60
 cen_time <- 1
 alpha <- 0.1
 power <- 0.85
 result_scu <- c()
 # Give a searching range of total sample size N:
n <- ceiling(100 / 2) 
     acc_time <- 100/r
     opt_rmst <- get_opt_des(n = n, sim_size = sim_size, 
                 acc_time = acc_time, cen_time = cen_time, int_step = 4, 
                 lambda_H0 = lambda_H0, lambda_H1 = lambda_H1, H1_type = 'PH', 
                 alpha = alpha, power = power, method = 'Complex') 
  

  }, error = function(e) {
    stop("An error occurred: ", e$message)
  })
  
  # If no error occurs, the test will pass
  expect_true(TRUE)

})
