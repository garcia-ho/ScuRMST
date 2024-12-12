# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html


library(testthat)
library(devtools)
#devtools::load_all()

packages <- c("survRM2", "mvtnorm", "ggplot2", "MASS", "tidyr", "survival", "nph", "tidyverse",
              "foreach", "doParallel", "cowplot", "IRdisplay", "rlang", "simtrial", "ggrepel")
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

invisible(lapply(packages, install_and_load))

n_cores <- detectCores()
cluster <- makeCluster(min(16, n_cores))  # Use the minimum of 32 or the number of available cores
registerDoParallel(cluster)
invisible(clusterEvalQ(cluster, {
  library(survRM2)
  library(mvtnorm)
  library(cubature)
  library(survival)
  library(nph)
  library(simtrial)
  library(foreach)
}))

# Export necessary functions and variables to the cluster
clusterExport(cluster, "expo_gen_2stages")

# test
test_that('Testing expo_gen_2stages and RMST_sim_test',{
  result <- tryCatch({

  median_con <- 10 # month
  lambda_H0 <- log(2)/median_con
  lambda_H1 <- lambda_H0 * 0.67
  sim_size <- 1000 
  acc_time <- 24
  cen_time <- 12
  tau <- 24
  n <- 100  
  set.seed(2024)
  data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                           cen_time = cen_time,arm = 0, interim = 0)[ , c(4,5,1)]  
  data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H1, dist = 'exp', 
                            cen_time = cen_time,arm = 1, interim = 0)[ , c(4,5,1)] 
  simple_rmst_2 <- RMST_sim_test(data_C = data_C, data_E = data_E_H1, sim_size = sim_size, tau = tau,
                            n = n, alpha = 0.05 ,sided = 'two_sided')
  print(simple_rmst_2$test_result$rejection)
  
  }, error = function(e) {
    stop("An error occurred: ", e$message)
  })
  
  # If no error occurs, the test will pass
  expect_true(TRUE)

})

# Stop the cluster after tests
stopCluster(cluster)