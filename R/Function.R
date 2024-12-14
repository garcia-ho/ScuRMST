#—————————————————————————————————————————————————————
# All functions used in this project are stored here!! 
#—————————————————————————————————————————————————————

# Import packages
#' @import survRM2
#' @import mvtnorm
#' @import ggplot2
#' @import MASS
#' @import truncnorm
#' @import tidyr
#' @import survival
#' @import nph
#' @import tidyverse
#' @import foreach
#' @import doParallel
#' @import cowplot
#' @import IRdisplay
#' @import rlang
#' @import simtrial
#' @import ggrepel
#' @import viridis
#' @importFrom stats dnorm filter integrate pnorm qnorm quantile rexp runif
#' @importFrom dplyr mutate
NULL
#‘
#’
utils::globalVariables(c("k", "variable", "value", "linetype_group"))

#' expo_gen_2stages
#'
#' This function generate exponential dist survival data for 2 stages test.
#' The censoring distribution in interim period is different from the whole trail 
#' If 'interim' is a number, it returns an (N * 5) array : 
#' \code{[arm, obs_time_int, event_int, obs_time_fin, event_fin]}
#' If 'interim' is a list or c() meaning a series of interim timepoints, 
#' the third dimension of result array return the survival result of ith interim (N * 5) 

#' @param  N Number of patients
#' @param  dist 'exp' for exponential, 'pcw_exp' for piecewise exponential 
#' @param  acc_time Accrual time period with constant rate
#' @param  lambda The parameter for exponential distribution
#' @param  cen_time Extra (minimum) censoring period after accrual period
#' @param  arm Group label(0,1)
#' @param  interim Interim time point
#' @param  HR1 For piecewise exponential only. The hazard ratio before change_time
#' @param  HR2 For piecewise exponential only. The hazard ratio after change_time 
#' @param  change_time For piecewise exponential only. The time when hazard ratio changes from HR1 to HR2  
#' 
#' @examples
#' set.seed(2024)
#' sim_size <- 5000
#' N <- 112
#' n <- ceiling(N / 2)  # total sample size per arm
#' r <- 60
#' acc_time <- N / r
#' cen_time <- 1
#' lambda_H0 <- 0.9 * 1.5
#' lambda_H1 <- 0.9
#' int_step <- 4
#'
#' int_factor <- seq(0.4, 0.7, by = int_step / N)  
#' # Each time interim sample size increase by 4
#' interim_list <- int_factor * acc_time
#'
#' # parameter interim can be a list
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, 
#'                      dist = 'exp', cen_time = cen_time, arm = 0, interim = interim_list) 
#' print(data_C[ , ,1])
#'
#' # parameter interim is a number
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, 
#'                    dist = 'exp', cen_time = cen_time, arm = 0, interim = interim_list[1]) 
#' print(data_C)
#' 
#' @export
expo_gen_2stages <- function(N,dist,acc_time,cen_time,lambda,HR1,HR2,arm,interim,change_time)
{
  if (dist == 'exp') {
    survival_time_all <- rexp(N, rate = lambda)
  } else if (dist == 'pcw_exp') { #piecewise exponential with one change time point
    survival_time_all <- rpwexp(n = N,fail_rate = data.frame(rate = c(lambda * HR1, lambda * HR2), 
                                            duration = c(change_time)))
}

arm_all <- rep(arm, N)
entry_time_all <- runif(N, min = 0, max = acc_time)
censor_time_fin <- runif(N, min = cen_time, max = acc_time + cen_time)
all_data <- cbind(survival_time_all, entry_time_all, censor_time_fin, arm_all)

# _______________________Define function for event counting________________________________
# filter the sample die before interim
cal_event_int <- function(row) { 
  survival_time_all <- row[1]
  censor_time_fin <- row[3]
  entry_time_all <- row[2]
  if (entry_time_all >= interim_val) { # Not in the interim
    return(c(0,0))
  } 
  else if (entry_time_all + min(survival_time_all, censor_time_fin) < interim_val) {
    return(c(min(survival_time_all, censor_time_fin), 1))
  } 
  else if (entry_time_all + min(survival_time_all, censor_time_fin) > interim_val &&
           entry_time_all < interim_val) {
    return(c(interim_val - entry_time_all, 0))
  }
}
# calculate the event in final study
cal_event_fin <- function(row) { 
  survival_time_all <- row[1]
  entry_time_all <- row[2]
  censor_time_fin <- row[3]
  obs_time_int <- row[5]
  event_int <- row[6]
if(entry_time_all >= interim_val) {  
  return(c(min(survival_time_all,censor_time_fin), 
          as.integer(survival_time_all <= censor_time_fin)))
}
else if (event_int == 0 && obs_time_int != 0) { 
  return(c(min(survival_time_all, censor_time_fin), 
          as.integer(survival_time_all <= censor_time_fin)))
}
else {
  return(c(obs_time_int,event_int))
}
}
#_________________________________________________________________

if (length(interim) == 1) {
    interim_val <- interim
    obs_event_int <- apply(all_data, 1, cal_event_int)
    sur_data_int <- cbind(all_data,t(obs_event_int))
    obs_event_fin <- apply(sur_data_int, 1, cal_event_fin) 
    all_data <- cbind(sur_data_int,t(obs_event_fin)) 
    return(all_data[ ,c(4,5,6,7,8)])
# In all_data 4th column is arm, 5th obs_time_int, 
# 6th event_int, 7th obs_time_fin, 8th event_fin
  }

else 
  {  # interim is a c() it return different interim event 
    result <- array(NA, dim = c(N , 5, length(interim)))
    for (i in 1 : length(interim))
    {
      interim_val <- interim[i]
      obs_event_int <- apply(all_data, 1, cal_event_int)
      sur_data_int <- cbind(all_data,t(obs_event_int))
      obs_event_fin <- apply(sur_data_int, 1, cal_event_fin) 
      result[, , i] <- cbind(sur_data_int,t(obs_event_fin))[ ,c(4,5,6,7,8)]
    }
   return(result)
  }
}





#' RMST_sim_cal
#' 
#' This function estimate the RMST values of each arm based on generated survival data.
#' It returns 2 * sim_size matrix. First row control group, second one experiment group.
#' There is no guaranteed that the maximum simulated survival time is larger than tau. 
#' Please refer to the 'tau' in survRM2 package: rmst2 function
#' The output of this function is array to accelerate the following grid search code.

#' @param n Sample size in each arm
#' @param data_E Survival data of experiment group generated by expo_gen_2stages
#' @param data_C Survival data of control group generated by expo_gen_2stages
#' @param tau Cut-off time for RMST
#' @param sim_size Simulation times
#' 
#' @examples
#' sim_size <- 5000 
#' N <- 100
#' r <- 60
#' acc_time <- N / r
#' cen_time <- 1
#' lambda_H1 <- 0.9
#' HR <- 1.7
#' lambda_H0 <- 0.9 * 1.7
#' change_time <- 1
#' interim <- 0.6 * acc_time
#' n <- ceiling(N / 2)
#' alpha <- 0.05
#' tau_f <- 2.5
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                            lambda = lambda_H0, dist = 'exp', 
#'                            cen_time = cen_time,arm = 0, interim = interim)    
#' data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                            lambda = lambda_H0, dist = 'exp', 
#'                            cen_time = cen_time,arm = 1, interim = interim)
#' rmst_h0_int <- RMST_sim_cal(n = n, data_E = data_E_H0[ , c(2,3,1)], 
#'                            data_C = data_C[ , c(2,3,1)],
#'                            tau = tau_f,sim_size = sim_size)
#' print(rmst_h0_int)
#' @export

RMST_sim_cal <- function(n,data_E,data_C,tau,sim_size)
{
    sim_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'survRM2') %dopar% {
    pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
    # the simulation data is not guaranteed to be larger than tau
    if (tau < min(max(data_E[((k-1)*n+1):(k*n),1]),  max(data_C[((k-1)*n+1):(k*n),1]))) {  
        rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], tau = tau) # No need for adjustment
      }  else {
        rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3])
        # The rmst2 function will automatically adjust tau for us if it's not specified
      }
    c(rmst_result$RMST.arm0$rmst[1],rmst_result$RMST.arm1$rmst[1])
    }
    
    return(sim_result)
}





#' RMST_sim_test
#' 
#' Different from RMST_sim_cal, It return a dataframe of p-value and rejection times of single-stage RMST test .
#' It also counts the times of tau adjustment (adjusted tau is the minimax survival time of two groups)
#' This function can be used to compare our rejection method with classical RMST difference test 
#' 
#' @param data_E Survival data of experiment group generated by expo_gen_2stages
#' @param data_C Survival data of control group generated by expo_gen_2stages
#' @param sim_size Simulation times
#' @param n Sample size in each arm
#' @param alpha Stated type I error level
#' @param sided "two_sided" for two-sided test, "greater" for one-sided superiority test
#' @param tau Prespecified cut-off time for RMST(used for final stage)
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{test_result}{A data frame with the following columns:
#'     \describe{
#'       \item{rejection}{The proportion of simulations where the RMST test p-value is less than or equal to alpha.}
#'       \item{tau adjustment}{The proportion of simulations where tau adjustment occurred.}
#'     }
#'   }
#'   \item{p_value}{A numeric vector of the p-values for each simulation.}
#' }
#' 
#' @examples
#'median_con <- 10 # month
#' lambda_H0 <- log(2)/median_con
#' lambda_H1 <- lambda_H0 * 0.67
#' sim_size <- 5000 
#' acc_time <- 24
#' cen_time <- 12
#' tau <- 24
#' n <- 100  
#' set.seed(2024)
#' 
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                    lambda = lambda_H0, dist = 'exp', 
#'                    cen_time = cen_time,arm = 0, interim = 0)[ , c(4,5,1)]    
#' data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                    lambda = lambda_H0, dist = 'exp', 
#'                    cen_time = cen_time,arm = 1, interim = 0)[ , c(4,5,1)]                    
#' simple_rmst <- RMST_sim_test(data_C = data_C, data_E = data_E_H0, 
#'                    sim_size = sim_size, tau = tau, 
#'                    n = n, alpha = 0.05 ,sided = 'greater')
#' 
#' print(simple_rmst$test_result)
#' print(simple_rmst$p_value)
#' @export
#' 
RMST_sim_test <- function(n, data_E, data_C, tau, sim_size, alpha, sided)
{
    tau_adj_count <- 0
    sim_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'survRM2') %dopar% {
          test_res <- 0
          p <- 0
          pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
          
        if (tau < min(max(data_E[((k-1)*n+1):(k*n),1]),  max(data_C[((k-1)*n+1):(k*n),1]))) {
            rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], tau = tau, alpha = alpha)
            # the simulation data is not guaranteed to be larger than tau
          }  else {
            rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], alpha = alpha)
            tau_adj_count <- tau_adj_count + 1   # tau is adjusted automatically
          }
        if (sided == 'two_sided') {
          p <- rmst_result$unadjusted.result[1,4]  
          # The p value of RMST difference test. It's a two sided test in the package
          } else if (sided == 'greater') {
          diff <- rmst_result$unadjusted.result[1,1]
          std <- (rmst_result$unadjusted.result[1,1] - 
                  rmst_result$unadjusted.result[1,2]) / qnorm(1 - alpha/2)
          p <- 1 - pnorm(diff / std)
          }
        if ( p <= alpha ) {
            test_res <- 1
          } else {
            test_res <- 0
          }
        c(test_res,tau_adj_count, p )  # The last element of the 2nd row is tau_adj_count
        }

    return(list(test_result = data.frame('rejection' = sum(sim_result[1, ]) / sim_size, 
                                      'tau adjustment' = sim_result[2, sim_size] / sim_size),
                  p_value = sim_result[3, ])) 
}





#' log_rank_sim 
#' 
#' For log-rank test simulation. Single-stage p-value or two-stage statistics. 
#' Return the simulated type I error, empirical mean and variance of Z-statistics of log-rank test
#' @param data_E Survival data of experiment group generated by expo_gen_2stages
#' @param data_C Survival data of control group generated by expo_gen_2stages
#' @param sim_size simulation times
#' @param n Sample size in each arm
#' @param alpha Stated type I error level
#' @param sided "two_sided" for two-sided test, "greater" for one-sided superiority test
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{rejection}{The proportion of simulations where the log-rank test statistic is less than or equal to alpha.}
#'   \item{z_stats}{A numeric vector of the z statistics (W/sigma) for each simulation.}
#'   \item{var_w}{The variance of statistics W for correlation calculation.}
#' }
#' 
#' @examples
#' sim_size <- 5000 
#' N <- 100
#' r <- 60
#' acc_time <- N / r
#' cen_time <- 1
#' lambda_H1 <- 0.9
#' HR <- 1.7
#' lambda_H0 <- 0.9 * 1.7
#' change_time <- 1
#' interim <- 0.6 * acc_time
#' n <- ceiling(N / 2)
#' alpha <- 0.05
#' 
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                        lambda = lambda_H0, dist = 'exp', 
#'                        cen_time = cen_time,arm = 0, interim = interim)    
#' data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                        lambda = lambda_H0, dist = 'exp', 
#'                        cen_time = cen_time,arm = 1, interim = interim)
#' lr_h0_int <- log_rank_sim(data_C = data_C[ , c(2,3,1)], data_E = data_E_H0[ , c(2,3,1)], 
#'                         sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')                                
#' print(lr_h0_int$rejection)
#' print(lr_h0_int$z_stats)
#' print(lr_h0_int$var_w)
#' 
#' @export
#' 
log_rank_sim <- function(data_C, data_E, sim_size, n, alpha, sided)
{
  logrank_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'nph') %dopar% {
    pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
    if (sided == 'greater') {
        result <- logrank.test(pre_data[,1],pre_data[,2],pre_data[,3],alternative = "greater")
      } 
    else if (sided == 'two_sided') {
        result <- logrank.test(pre_data[,1],pre_data[,2],pre_data[,3],alternative = "two.sided")
      }
    p <- result$test$p
    z_stats <- result$test$z  # return the z statistics for two stages trials
    var_w <- sum(result$D$w^2 * result$D$var)
    c(p, z_stats, var_w)
  }
  return(list(rejection = sum(logrank_result[1, ] <= alpha) / sim_size, # rejection times
              z_stats = logrank_result[2, ], # the z statistics W/sigma of every simulation
              var_w = logrank_result[3, ]) )  # variance of W for correlation 
}





#' theo_RMST
#' Calculate the theorecital RMST and RSDST of a explicit survival function(exponential only)
#' @param lambda The lambda parameter of exponential distribution
#' @param dist Only 'exp' now, require further development for other distributions
#' @param tau Cut-off time for RMST
#' 
theo_RMST <- function(lambda, dist, tau) 
{
  if (dist == 'exp') 
    {
    surv_fun <- function(t) {
             exp(-lambda * t)
        }
    RMST <- integrate(surv_fun, lower = 0, upper = tau)$value
    }
    return(RMST)
}





#' PET_norm
#' 
#' It returns \code{c(PET0, PET1)}, where \code{PET0 = prob(E-C>m1 & E>t1|H0)} using bivariate normal density.
#' Input the estimated RMST mean and variance of each group (control and experiment) 
#' and the critical value \code{(m1,t1)}
#' @param mu_c Mean of the control group.
#' @param var_c Variance of the control group.
#' @param mu_e Mean of the experimental group.
#' @param var_e Variance of the experimental group.
#' @param m1 Threshold for RMST difference E-C at interim analysis.
#' @param t1 Time point for Experimental group RMST E at interim analysis.
#' @return Probability of Early Termination (PET).
#' 

PET_norm <- function(mu_c,var_c,mu_e,var_e,m1,t1)
{
    mu_h0 <- c(0, mu_c)
    sigma_h0 <- matrix(c(2*var_c, var_c, var_c, var_c), nrow = 2)
    mu_h1 <- c(mu_e-mu_c, mu_e)
    sigma_h1 <- matrix(c(var_e+var_c, var_e, var_e, var_c), nrow = 2)
    upper <- c(Inf, Inf)
    lower <- c(m1,t1)
    p_rj_h0 <- 1 - pmvnorm(lower, upper, mean = mu_h0, sigma = sigma_h0) # p(E-C>m1 & E>t1|H0)
    p_rj_h1 <- 1 - pmvnorm(lower, upper, mean = mu_h1, sigma = sigma_h1)
    return(c(p_rj_h0,p_rj_h1))
}





#' mu_cov_mc 
# 
#' This function uses Monte Carlo simulation to calculate the variance-covariance matrix of \code{[E1-C1, E1, E2-C2, E2]}. 
#' The calculation formula refers to Lu (2021) sequential trials.
#' The input \code{rmst_int} should be the interim RMST data of two groups (generated by \code{RMST_sim_cal}).
#' \code{rmst_fin} is the final RMST data of two groups (generated by \code{RMST_sim_cal}).
#'
#' @param rmst_int A matrix of interim RMST data generated by \code{RMST_sim_cal} function.
#' @param rmst_fin A matrix of final stage RMST data generated by \code{RMST_sim_cal} function.
#' @param sim_size n integer representing the number of Monte Carlo simulation iterations.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{mu}{A numeric vector of the means \code{[mean(E1-C1), mean(E1), mean(E2-C2), mean(E2)]}.}
#'   \item{sigma}{A 4x4 variance-covariance matrix of \code{[E1-C1, E1, E2-C2, E2]}.}
#' }
#'
#' @examples
#' # Example usage of mu_cov_mc
#' sim_size <- 5000 
#' N <- 100
#' r <- 60
#' acc_time <- N / r
#' cen_time <- 1
#' lambda_H1 <- 0.9
#' HR <- 1.7
#' lambda_H0 <- 0.9 * 1.7
#' change_time <- 1
#' interim <- 0.6 * acc_time
#' n <- ceiling(N / 2)
#' alpha <- 0.05
#' tau_f <- 2.5
#' 
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                            lambda = lambda_H0, dist = 'exp', 
#'                            cen_time = cen_time,arm = 0, interim = interim)    
#' data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                              lambda = lambda_H0, dist = 'exp', 
#'                              cen_time = cen_time,arm = 1, interim = interim)
#' rmst_h0_int <- RMST_sim_cal(n = n, data_E = data_E_H0[ , c(2,3,1)], 
#'                            data_C = data_C[ , c(2,3,1)],
#'                            tau = interim, sim_size = sim_size)
#' rmst_h0_fin <- RMST_sim_cal(n = n, data_E = data_E_H0[ , c(4,5,1)], 
#'                            data_C = data_C[ , c(4,5,1)],
#'                            tau = tau_f, sim_size = sim_size)
#' mu_cov_h0 <- mu_cov_mc(rmst_int = rmst_h0_int, rmst_fin = rmst_h0_fin, sim_size = sim_size)
#' 
#' print(mu_cov_h0$mu)
#' print(mu_cov_h0$sigma)
#'
#' @export

mu_cov_mc <- function(rmst_int, rmst_fin, sim_size){
    diff_C <- rbind(rmst_int[1, ] - mean(rmst_int[1, ]), rmst_fin[1, ] - mean(rmst_fin[1, ]))
    cov_C <- matrix(0, nrow = 2, ncol = 2)
    for (i in 1:sim_size) 
      {
        product <- diff_C[, i] %*% t(diff_C[, i])  
        cov_C <- cov_C + product  
      }
    cov_C <- cov_C / sim_size   # [ Var(C1)  Cov(C1, C2)
                                #  Cov(C1, C2)  Var(C2) ]

    diff_E <- rbind(rmst_int[2, ] - mean(rmst_int[2, ]), rmst_fin[2, ] - mean(rmst_fin[2, ]))
    cov_E <- matrix(0, nrow = 2, ncol = 2)
    for (i in 1:sim_size) 
      {
        product <- diff_E[,i] %*% t(diff_E[,i])  
        cov_E <- cov_E + product  
      }
    cov_E <- cov_E / sim_size   # [ Var(E1)  Cov(E1, E2)
                                #  Cov(E1, E2)  Var(E2) ]

    var_cov_all <- matrix(
            c(cov_E[1,1]+cov_C[1,1], cov_E[1,1], cov_C[1,2]+cov_E[1,2], cov_E[1,2],
              cov_E[1,1], cov_E[1,1], cov_E[1,2], cov_E[1,2],
              cov_C[1,2]+cov_E[1,2], cov_E[1,2], cov_E[2,2]+cov_C[2,2], cov_E[2,2],
              cov_E[1,2], cov_E[1,2], cov_E[2,2], cov_E[2,2]), nrow = 4, ncol = 4)
    mu <- c(mean(rmst_int[2,] - rmst_int[1,]), mean(rmst_int[2,]),
            mean(rmst_fin[2,] - rmst_fin[1,]), mean(rmst_fin[2,]))

    return(list(mu = mu,
                sigma = var_cov_all))
  }






#' find_m_logrank
#' 
#' This function performs a grid search to find the critical values (m1, m2) for the log-rank test in a two-stage design.
#' It makes use of probability \code{Prob(W1/sigma1 > m1 & W2/sigma2 > m2 | H_0) = alpha} and
#' the normality of Z-statistics (Kwak and Jung(2017))
#' (W/sigma) are estimated by simulation. The z statistics can be obtained using the \code{log_rank_sim} function 
#' and accessed via its \code{$z_stats} output.
#' When \code{power} is not given, the function finds the critical values (m1, m2) that control the overall type I error
#' while maximizing the power. When \code{power} is given, it searches for valid critical values such that the empirical
#' power is greater than the specified \code{power}. If \code{power} is given, \code{int_n} and \code{fin_n} are required
#' for the expected sample size calculation.
#'
#' @param logrank_data A matrix of log-rank test z statistics for each simulation.
#' @param corr_h0 A numeric value representing the correlation of two stages, calculated as \code{sqrt(var(W1) / var(W))}.
#' @param search_times An integer specifying the number of grid search iterations.
#' @param int_n An optional integer specifying the interim sample size (required if \code{power} is given).
#' @param fin_n An optional integer specifying the total sample size in both arms (required if \code{power} is given).
#' @param alpha A numeric value specifying the significance level for hypothesis testing.
#' @param sim_size An integer specifying the number of simulations to perform, say 10000.
#' @param power An optional numeric value specifying the desired power for the test.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{m1}{Critical value for the first stage test.}
#'   \item{m2}{Critical value for the second stage test.}
#'   \item{PET0}{Probability of early termination under H0.}
#'   \item{PET1}{Probability of early termination under H1.}
#'   \item{alpha}{Empirical overall type I error rate.}
#'   \item{power}{Empirical power of the test.}
#'   \item{PET}{Probability of early termination.}
#'   \item{EN0}{Expected sample size under H0.}
#'   \item{EN1}{Expected sample size under H1.}
#'   \item{EN}{Overall expected sample size.}
#' }
#' If no valid critical values are found, the function returns a data frame with: 
#' \code{m1 = 0, m2 = 0, PET0 = 0, PET1 = 0, alpha = 0, power = 0, PET = 0, EN0 = NA, EN1 = NA, EN = NA}.
#'
#' @examples
#' sim_size <- 5000 
#' N <- 100
#' r <- 60
#' acc_time <- N / r
#' cen_time <- 1
#' lambda_H1 <- 0.9
#' HR <- 1.7
#' lambda_H0 <- 0.9 * 1.7
#' change_time <- 1
#' interim <- 0.6 * acc_time
#' n <- ceiling(N / 2)
#' alpha <- 0.05
#' tau_f <- 2.5
#' 
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                        lambda = lambda_H0, dist = 'exp', 
#'                        cen_time = cen_time,arm = 0, interim = interim)    
#' data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                        lambda = lambda_H0, dist = 'exp', 
#'                        cen_time = cen_time,arm = 1, interim = interim)
#' data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                        lambda = lambda_H1, dist = 'exp', 
#'                        cen_time = cen_time,arm = 1, interim = interim)
#' lr_h0_int <- log_rank_sim(data_C = data_C[ , c(2,3,1)], data_E = data_E_H0[ , c(2,3,1)], 
#'                         sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
#' lr_h1_int <- log_rank_sim(data_C = data_C[ , c(2,3,1)], data_E = data_E_H1[ , c(2,3,1)], 
#'                         sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
#' lr_h0_fin <- log_rank_sim(data_C = data_C[ , c(4,5,1)], data_E = data_E_H0[ , c(4,5,1)], 
#'                         sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
#' lr_h1_fin <- log_rank_sim(data_C = data_C[ , c(4,5,1)], data_E = data_E_H1[ , c(4,5,1)], 
#'                         sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
#' # Get W/sigma
#' z_stats_h1_int <- lr_h1_int$z_stats
#' z_stats_h1_fin <- lr_h1_fin$z_stats
#' z_stats_h0_int <- lr_h0_int$z_stats
#' z_stats_h0_fin <- lr_h0_fin$z_stats
#' logrank_data <- rbind(z_stats_h0_int, z_stats_h1_int, z_stats_h0_fin, z_stats_h1_fin) 
#' corr_h0 <- sqrt(mean(lr_h0_int$var_w) / mean(lr_h0_fin$var_w)) 
#' 
#' best_lr <- find_m_logrank(logrank_data = logrank_data, search_times = 200, corr_h0 = corr_h0,
#'                            alpha = 0.05, sim_size = sim_size)
#'
#' @export
#' 
find_m_logrank <- function(logrank_data, corr_h0, search_times, int_n = NULL, 
                          fin_n = NULL, alpha, sim_size, power = NULL) 
{
  z_stats_h0_int <-  logrank_data[1, ]
  z_stats_h1_int <-  logrank_data[2, ]
  z_stats_h0_fin <-  logrank_data[3, ]
  z_stats_h1_fin <-  logrank_data[4, ]

  # Cov matrix of (W1, W | H0)
  sigma_h0 <- matrix(c(1, corr_h0, corr_h0, 1), nrow = 2)
  mean_h0 <- c(0, 0)

  #Function to  solve m2 in P( W1/sigma1 > m1 & W/sigma > m2 | H0) = alpha given m1
  norm_2d <- function(m2, m1, mean, sigma, alpha) 
    {
      prob <- pmvnorm(lower = c(m1, m2), 
                    upper = rep(Inf, 2), 
                    mean = mean, 
                    sigma = sigma)
      return (prob - alpha)
    }

  cal_proc <- function(lr_int, lr_fin, m1, m2) {
                       sum((lr_int > m1) & (lr_fin > m2))/ sim_size
                      }
  cal_pet <- function(lr_int, m1){
                      sum(lr_int <= m1) / sim_size
                      }
  ub_m1 <- quantile(z_stats_h1_int, 0.7)
  lb_m1 <- quantile(z_stats_h0_int, 0.3)
  ub_m2 <- quantile(z_stats_h1_fin, 0.7)
  lb_m2 <- quantile(z_stats_h0_fin, 0.3)
  m1_values <- seq(lb_m1, ub_m1, by = (ub_m1 - lb_m1) / search_times) 
  m2_values <- seq(lb_m2, ub_m2, by = (ub_m2 - lb_m2) / search_times)
  combinations <- expand.grid(m1 = m1_values, m2 = m2_values)
  combinations$PET0 <- sapply(1:nrow(combinations), function(i) {
                              cal_pet(z_stats_h0_int, combinations$m1[i])})
  combinations$PET1 <- sapply(1:nrow(combinations), function(i) {
                              cal_pet(z_stats_h1_int, combinations$m1[i]) }) 
  combinations$alpha <- sapply(1:nrow(combinations), function(i) {
                              cal_proc(z_stats_h0_int, z_stats_h0_fin, 
                                    combinations$m1[i], combinations$m2[i]) })
  combinations$power <- sapply(1:nrow(combinations), function(i) {
                              cal_proc(z_stats_h1_int, z_stats_h1_fin, 
                                    combinations$m1[i], combinations$m2[i]) })
  if (is.null(power)) { # find the most powerful one
      fil_combs <- combinations[abs(combinations$alpha - alpha) < 0.05 * alpha &
                                combinations$alpha < alpha, ]
      crit_val_res <- fil_combs[which.max(fil_combs$power), ]
            }
  else { 
      crit_val_res <- combinations[abs(combinations$alpha - alpha) < 0.05 * alpha & 
                                (combinations$alpha < alpha) & (combinations$power > power)&
                                (abs(combinations$power - power) < 0.05 * power) , ]
            }

  if (is.null(crit_val_res) || dim(crit_val_res)[1] == 0) {   
      return(data.frame(m1 = 0, m2 = 0, PET0 = 0, PET1 = 0, alpha = 0, power = 0, 
                        PET = 0, EN0 = NA, EN1 = NA, EN = NA))
      }
  crit_val_res <- data.frame(crit_val_res)
  colnames(crit_val_res) <- c('m1', 'm2', 'PET0', 'PET1', 'alpha', 'power')
  
  if(is.null(power)) # Power is not specified, return the most powerfule result
    {
      powerful_m1 <- crit_val_res[which(crit_val_res$power == max(crit_val_res$power)), ]

      if(is.null(dim(powerful_m1))){ #unique solution
        return(powerful_m1)
      }
      else {  # find the smallest m1 if multiply solution exist
       return(powerful_m1[which(powerful_m1$m1 == max(powerful_m1$m1)), ])
      }
    }

  else  # When power is given, find the min(E(N)) design under (alpha, power) constraint
  {
    best_res <- crit_val_res[which(crit_val_res$power >= power), ]
    if (dim(best_res)[1] == 0){ #no valid result
      return(data.frame(m1 = 0, m2 = 0, PET0 = 0, PET1 = 0, alpha = 0, power = 0, 
                        PET = 0, EN0 = NA, EN1 = NA, EN = NA))
    }
    best_res$PET <- rowMeans(best_res[, c('PET0','PET1')])
    best_res$EN0 <- best_res$PET0 * int_n + (1 - best_res$PET0) * fin_n
    best_res$EN1 <- best_res$PET1 * int_n + (1 - best_res$PET1) * fin_n
    best_res$EN <- rowMeans(best_res[, c('EN0', 'EN1')])
  }
    return(best_res)
 }




#' adp_grid_src
#' 
#' This is the function for adaptive grid search of RMST two-stage design with one searching parameters gamma.
#' The grid search make use of the conditional distribution (Normal) of E | D, which is illustrated in paper.
#' When \code{power} is not given, the function finds the critical values that control the overall type I error \code{<=alpha}
#' while maximizing the empirical power. When \code{power} is given, it searches for valid critical values such that the empirical
#' power is greater than the specified \code{power}. If \code{power} is given, \code{int_n} and \code{fin_n} are required
#' 
#' @param rmst_data Combined RMST data generated by the \code{RMST_sim_cal} function.
#' @param mu_cov_h0 Variance-covariance matrix of \code{[E1-C1, E1, E2-C2, E2]} under the null hypothesis (H0), obtained from the \code{mu_cov_mc} function.
#' @param mu_cov_h1 Variance-covariance matrix of \code{[E1-C1, E1, E2-C2, E2]} under the alternative hypothesis (H1), obtained from the \code{mu_cov_mc} function.
#' @param int_n Interim sample size (n).
#' @param fin_n Total sample size (N).
#' @param sim_size Number of simulations to perform in the grid search.
#' @param method Method used for the RMST design. "Simple" indicates a Simple RMST with the rejection region \code{(D1 > m1 & D2 > m2)}, where \code{D = E - C} is the RMST difference between groups. "Complex" indicates a Sculpted RMST with the rejection region \code{(D1 > m1 & E1 > q1 & D2 > m2 & E2 > q2)}, where E is the RMST in the experimental group.
#' @param alpha Stated type I error level.
#' @param power Desired power of the test (optional).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{m1}{Critical value for RMST difference at the first stage}
#'   \item{m2}{Critical value for RMST difference at the second stage}
#'   \item{q1}{Critical value for Experiment group RMST at the first stage } 
#'   \item{q2}{Critical value for Experiment group RMST at the second stage } 
#'   \item{PET0}{Probability of early termination under H0.}
#'   \item{PET1}{Probability of early termination under H1.}
#'   \item{alpha}{Empirical overall type I error rate.}
#'   \item{power}{Empirical power of the test.}
#'   \item{PET}{Probability of early termination \code{(PET0+PET1)/2}.}
#'   \item{EN0}{Expected sample size under H0.}
#'   \item{EN1}{Expected sample size under H1.}
#'   \item{EN}{Average expected sample size.}
#' }
#' If no valid critical values are found, the function returns a data frame with \code{m1 = 0, m2 = 0, q1 = 0, q2 = 0, PET0 = 0, PET1 = 0, alpha = 0, power = 0, PET = 0, EN0 = NA, EN1 = NA, EN = NA}.
#'
#' @examples
#' # Example usage of adp_grid_src function
#' sim_size <- 5000 
#' N <- 100
#' r <- 60
#' acc_time <- N / r
#' cen_time <- 1
#' lambda_H1 <- 0.9
#' HR <- 1.7
#' lambda_H0 <- 0.9 * 1.7
#' change_time <- 1
#' interim <- 0.6 * acc_time
#' n <- ceiling(N / 2)
#' alpha <- 0.05
#' tau_f <- 2.5
#' 
#' # Generate data for control arm (C) and experimental arm (E) under H0 and H1
#' data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                            lambda = lambda_H0, dist = 'exp', 
#'                            cen_time = cen_time, arm = 0, interim = interim)
#' data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                            lambda = lambda_H0, dist = 'exp',
#'                            cen_time = cen_time, arm = 1, interim = interim)
#' data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, 
#'                            lambda = lambda_H1, dist = 'exp', 
#'                            cen_time = cen_time, arm = 1, interim = interim)
#'
#' # Calculate RMST data for interim and final stages under H0 and H1
#' rmst_h0_int <- RMST_sim_cal(n = n, data_E = data_E_H0[ , c(2,3,1)], 
#'                data_C = data_C[ , c(2,3,1)], tau = interim, sim_size = sim_size)
#' rmst_h0_fin <- RMST_sim_cal(n = n, data_E = data_E_H0[ , c(4,5,1)], 
#'                data_C = data_C[ , c(4,5,1)], tau = tau_f, sim_size = sim_size)
#' rmst_h1_int <- RMST_sim_cal(n = n, data_E = data_E_H1[ , c(2,3,1)], 
#'                data_C = data_C[ , c(2,3,1)], tau = interim, sim_size = sim_size)
#' rmst_h1_fin <- RMST_sim_cal(n = n, data_E = data_E_H1[ , c(4,5,1)], 
#'                data_C = data_C[ , c(4,5,1)], tau = tau_f, sim_size = sim_size)
#'
#' # Combine RMST data
#' rmst_data <- rbind(rmst_h0_int, rmst_h1_int, rmst_h0_fin, rmst_h1_fin)
#'
#' # Calculate variance-covariance matrices under H0 and H1
#' mu_cov_h0 <- mu_cov_mc(rmst_int = rmst_h0_int, rmst_fin = rmst_h0_fin, sim_size = sim_size)
#' mu_cov_h1 <- mu_cov_mc(rmst_int = rmst_h1_int, rmst_fin = rmst_h1_fin, sim_size = sim_size)
#'
#' # Perform adaptive grid search for the best RMST design
#' best_RMST <- adp_grid_src(rmst_data = rmst_data, mu_cov_h0 = mu_cov_h0, 
#'                mu_cov_h1 = mu_cov_h1, int_n = interim * r, fin_n = 2 * n, 
#'                alpha = 0.05, sim_size = sim_size, method = 'Complex')
#' 
#' @export
#' 
adp_grid_src <- function(rmst_data, mu_cov_h0, mu_cov_h1, int_n, fin_n, 
                        sim_size, method, alpha, power = NULL) 
  {
      # Interim
      mu1 <- mu_cov_h1$mu[c(1,2)]
      sigma1 <- mu_cov_h1$sigma[1:2, 1:2]
      # Final
      mu2 <- mu_cov_h1$mu[c(3,4)]
      sigma2 <- mu_cov_h1$sigma[3:4, 3:4]

      rmst_h0_int <- rmst_data[c(1,2) , ]
      rmst_h1_int <- rmst_data[c(3,4) , ]
      rmst_h0_fin <- rmst_data[c(5,6) , ]
      rmst_h1_fin <- rmst_data[c(7,8) , ]

  cal_q <- function(m, tar_prob, mu, sigma)  # conditional normal dist
      {
        mu_D <- mu[1]
        mu_E <- mu[2]
        sigma_D <- sqrt(sigma[1, 1])
        sigma_E <- sqrt(sigma[2, 2])
        rho <- sigma[1, 2] / (sigma_D * sigma_E) #corr
        #truncated normal
        alpha <- (m - mu_D) / sigma_D  
        # Mean and variance of the truncated normal distribution D | D > m
        mean_D_given_D_gt_m <- mu_D + sigma_D * dnorm(alpha) / (1 - pnorm(alpha))
        var_D_given_D_gt_m <- sigma_D^2 * (1 - (alpha * dnorm(alpha) / (1 - pnorm(alpha))) - 
                                          (dnorm(alpha) / (1 - pnorm(alpha)))^2)
        # Mean of E given D > m
        mean_E_given_D_gt_m <- mu_E + rho * (sigma_E / sigma_D) * 
                              (mean_D_given_D_gt_m - mu_D)
        # Variance of E given D > m
        var_E_given_D_gt_m <- (1 - rho^2) * sigma_E^2 + 
                              (rho * sigma_E / sigma_D)^2 * var_D_given_D_gt_m
        # Calculate q such that P(E > q | D > m) = p
        q <- qnorm(tar_prob, mean = mean_E_given_D_gt_m, sd = sqrt(var_E_given_D_gt_m), 
                    lower.tail = FALSE)
        if (is.nan(q)) {
          return(NA)
          } 
        else {
          return(q)
          }
      }

  safe_cal_q <- function(m, tar_prob, mu, sigma) {
      suppressWarnings(cal_q(m, tar_prob, mu, sigma))
    }

  cal_proc <- function(rmst_int, rmst_fin, m1, q1, m2, q2) {
                       sum((rmst_int[2, ] - rmst_int[1, ] > m1) & (rmst_int[2, ] > q1) &
                            (rmst_fin[2, ] - rmst_fin[1, ] > m2) & (rmst_fin[2, ] > q2)) / sim_size
                      }
  cal_pet <- function(rmst_int, m1, q1){
                      sum((rmst_int[2, ] - rmst_int[1, ] < m1) | 
                          (rmst_int[2, ] < q1)) / sim_size
                      }

  ub_m1 <- quantile(rmst_h1_int[2,] - rmst_h1_int[1, ], 0.7)
  lb_m1 <- quantile(rmst_h0_int[2,] - rmst_h0_int[1, ], 0.3)
  ub_m2 <- quantile(rmst_h1_fin[2,] - rmst_h1_fin[1, ], 0.7)
  lb_m2 <- quantile(rmst_h0_fin[2,] - rmst_h0_fin[1, ], 0.3)
  m1_values <- seq(lb_m1, ub_m1, by = (ub_m1 - lb_m1) / 100) 
  m2_values <- seq(lb_m2, ub_m2, by = (ub_m2 - lb_m2) / 100)

  # D1>m1, D2>m2
  if(method == 'Simple')
  {
    combinations <- expand.grid(m1 = m1_values, m2 = m2_values)
    combinations$q1 <- -Inf
    combinations$q2 <- -Inf
    combinations$gamma <- 0
    combinations$PET0 <- sapply(1:nrow(combinations), function(i) {
                              cal_pet(rmst_h0_int, combinations$m1[i], -Inf)})
    combinations$PET1 <- sapply(1:nrow(combinations), function(i) {
                              cal_pet(rmst_h1_int, combinations$m1[i], -Inf) }) 
    combinations$alpha <- sapply(1:nrow(combinations), function(i) {
                              cal_proc(rmst_h0_int, rmst_h0_fin, combinations$m1[i], -Inf, 
                                      combinations$m2[i], -Inf) })
    combinations$power <- sapply(1:nrow(combinations), function(i) {
                              cal_proc(rmst_h1_int, rmst_h1_fin, combinations$m1[i], -Inf, 
                                      combinations$m2[i], -Inf) }) 
    if (is.null(power)) { # find the most powerful one
      fil_combs <- combinations[abs(combinations$alpha - alpha) < 0.05 * alpha &
                                combinations$alpha < alpha, ]
      crit_val_res <- fil_combs[which.max(fil_combs$power), ]
            }
    else { 
      crit_val_res <- combinations[abs(combinations$alpha - alpha) < 0.05 * alpha & 
                                combinations$alpha < alpha & (combinations$power > power), ]
            }
  }

  # D1>m1, E1>q1, D2>m2, E2>q2
  if (method == 'Complex') 
  {
    crit_val_res <- foreach(gamma = seq(0, 0.1, by = 0.005), .combine = 'rbind') %dopar%
      {
        tar_prob_int <- exp(-gamma * (int_n / fin_n)) 
        tar_prob_fin <- exp(-gamma * (fin_n / fin_n)) 
        
        #interim
        q1_values <- sapply(m1_values, safe_cal_q, tar_prob = tar_prob_int, mu = mu1, sigma = sigma1)
        mq1 <- data.frame(m1_values = m1_values, q1_values = q1_values)
        mq1 <- data.frame(mq1[!is.na(mq1[,2]),])
          
        # #final
        q2_values <- sapply(m2_values, safe_cal_q, tar_prob = tar_prob_fin, mu = mu2, sigma = sigma2)
        mq2 <- data.frame(m2_values = m2_values, q2_values = q2_values)
        mq2 <- data.frame(mq2[!is.na(mq2[,2]),])

        combinations <- expand.grid(m1 = mq1$m1_values, m2 = mq2$m2_values)
        combinations$q1 <- mq1$q1_values[match(combinations$m1, mq1$m1_values)]
        combinations$q2 <- mq2$q2_values[match(combinations$m2, mq2$m2_values)]
        combinations$gamma <- gamma
        combinations$PET0 <- sapply(1:nrow(combinations), function(i) {
                                    cal_pet(rmst_h0_int, combinations$m1[i], combinations$q1[i])})
        combinations$PET1 <- sapply(1:nrow(combinations), function(i) {
                                    cal_pet(rmst_h1_int, combinations$m1[i], combinations$q1[i]) }) 
        combinations$alpha <- sapply(1:nrow(combinations), function(i) {
                                    cal_proc(rmst_h0_int, rmst_h0_fin, combinations$m1[i], combinations$q1[i], 
                                      combinations$m2[i], combinations$q2[i]) })
        combinations$power <- sapply(1:nrow(combinations), function(i) {
                                    cal_proc(rmst_h1_int, rmst_h1_fin, combinations$m1[i], combinations$q1[i], 
                                      combinations$m2[i], combinations$q2[i]) }) 
          
        if (is.null(power)) { # find the most powerful one
              fil_combs <- combinations[abs(combinations$alpha - alpha) < 0.05 * alpha &
                                        combinations$alpha < alpha, ]
              best_gamma <- fil_combs[which.max(fil_combs$power), ]
            }
        else { 
              best_gamma <- combinations[abs(combinations$alpha - alpha) < 0.05 * alpha & 
                                        combinations$alpha < alpha & (combinations$power > power), ]
            }
        best_gamma
      }   
  }

#filter and output
    if (is.null(power))  
      {
        if (dim(crit_val_res)[1] == 0 ) {   # Return NULL when something goes wrong
            return(data.frame(m1 = 0, q1 = 0, m2 = 0, q2 = 0, gamma = 0, 
                              PET0 = 0, PET1 = 0, alpha = 0, power = 0))
          }
        else {
            best_res <- crit_val_res[ which(crit_val_res[, 'power'] == max(crit_val_res[, 'power'])), ]
            return(data.frame(best_res))
          } 
      }

    else # find the min E(N)|H0 critical values
      {   
        if (dim(crit_val_res)[1] == 0) 
          {   # Return NULL when something goes wrong
            return(data.frame(m1 = 0, q1 = 0, m2 = 0, q2 = 0, gamma = 0, PET0 = 0, PET1 = 0, 
            alpha = 0, power = 0, PET = 0, EN0 = NA, EN1 = NA, EN = NA))
          }
          # calculate E(N)
          crit_val_res$PET <- (crit_val_res$PET0 + crit_val_res$PET1) / 2
          crit_val_res$EN0 <- crit_val_res$PET0  * int_n + (1 - crit_val_res$PET0 ) * fin_n
          crit_val_res$EN1 <- crit_val_res$PET1  * int_n + (1 - crit_val_res$PET1 ) * fin_n
          crit_val_res$EN <- (crit_val_res$EN0 + crit_val_res$EN1) / 2
          return(crit_val_res)
      }
  }






#' compare_line_plot
#' Used to draw line plot for comparing three methods (2 RMST and log-rank test) under different scenario
#' Input is the dataframe of 3m_comparison output. ** The order of variable is fixed **
#' @param data Please make sure the data structure is the same as the example on Github ipynb file
#' @param var_name The variable name shown on X-axis of plots
#' 
compare_line_plot <- function(data, var_name) 
  { 
    options(repr.plot.width = 20, repr.plot.height = 8)

    color_palette <- c("ScuRMST_power" = "darkred", "ScuRMST_alpha" = "darkred", 
                      "LR_power" = "lightgreen", "LR_alpha" = "lightgreen",
                      "SimRMST_power" = "blue", "SimRMST_alpha" = "blue",
                      "ScuRMST_PET0" = "darkred", "ScuRMST_PET1" = "darkred", 
                      "LR_PET0" = "lightgreen", "LR_PET1" = "lightgreen",
                      "SimRMST_PET0" = "blue", "SimRMST_PET1" = "blue")

    a_power_delta <- data.frame(data[, c(1,2,3,4,5,6,7)])
    colnames(a_power_delta) <- c(var_name,'LR_alpha','SimRMST_alpha','ScuRMST_alpha',
                          'LR_power', 'SimRMST_power','ScuRMST_power')
    a_power_long <- a_power_delta %>%
        pivot_longer(cols = -!!sym(var_name), names_to = "variable", values_to = "value")%>%
        mutate(linetype_group = ifelse(variable %in% 
            c("LR_alpha", "SimRMST_alpha", "ScuRMST_alpha"), "Alpha", "Power"))
    a_power_long <- a_power_long %>% filter(value != 0)   # 0 means could not find critical values

    plot1 <- ggplot(a_power_long, aes(x = !!sym(var_name), y = value, 
        color = variable, linetype = linetype_group)) +
    geom_point(size = 3) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, 
                aes(group = interaction(variable, linetype_group)), 
                linewidth = 1) + 
    #geom_vline(xintercept = 0.7, color = "red", linetype = "dashed", size = 1) +
    scale_linetype_manual(values = c("Alpha" = "solid", "Power" = "dotted")) +
    labs( linetype = "Line Type", color = "Variable",
          title = "Line Plot with Different Line Types") +
    scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), limits = c(0, 1)) +
    scale_color_manual(values = color_palette) +
    labs(x = var_name, y = "Value", color = "Variable",
      title = 'Type I error and Power') +
    theme_minimal(base_size = 18) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    guides(linetype = guide_legend(override.aes = list(color = "black")),
        color = guide_legend(override.aes = list(linetype = "solid")))

#_________next plot_________
    pet_delta <- data.frame(data[, c(1,8,9,10,11,12,13)])
    colnames(pet_delta) <- c(var_name,'LR_PET0', 'SimRMST_PET0', 'ScuRMST_PET0',
                      'LR_PET1', 'SimRMST_PET1', 'ScuRMST_PET1')
    pet_long <- pet_delta %>%
        pivot_longer(cols = -!!sym(var_name), names_to = "variable", values_to = "value")%>%
        mutate(linetype_group = ifelse(variable %in% 
              c("LR_PET0", "SimRMST_PET0", "ScuRMST_PET0"), "PET0", "PET1"))
    pet_long <- pet_long %>% filter(value != 0) 

    plot2 <- ggplot(pet_long, aes(x = !!sym(var_name), y = value, 
        color = variable, linetype = linetype_group)) +
        geom_point(size = 3) +
        geom_smooth(method = "loess", formula = y ~ x, se = FALSE, 
                aes(group = interaction(variable, linetype_group)), 
                linewidth = 1) + 
        #geom_vline(xintercept = 0.67, color = "red", linetype = "dashed", size = 1) +
        scale_linetype_manual(values = c("PET0" = "solid", "PET1" = "dotted")) +
        labs( linetype = "Line Type", color = "Variable",
          title = "Line Plot with Different Line Types") +

        scale_y_continuous(limits = c(0, 1)) +
        scale_color_manual(values = color_palette) +
        labs(x = var_name, y = "Value", color = "Variable",
        title = 'PET0 and PET1') +
        theme_minimal(base_size = 18) + 
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(fill = "white", color = NA),
              plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
        guides(linetype = guide_legend(override.aes = list(color = "black")),
              color = guide_legend(override.aes = list(linetype = "solid")))

    plot_grid(plot1, plot2, ncol = 2)

  }