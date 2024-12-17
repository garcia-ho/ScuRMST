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


#' Get Optimal Design
#'
#' This function search for the optimal design for a study, returning the interim sample size, total sample size
#' and other measurements of the design minimizing the expected sample size under H0: EN0.
#' For any overall sample size N, the function only search for interim sample size ranging from 0.3*N to 0.7*N.
#'
#' @param n Overall sample size for each group. This should maintain consistency in inputs.
#' @param sim_size Number of simulations to perform.
#' @param acc_time Accrual time for the study.
#' @param cen_time Time until censoring occurs.
#' @param int_step The difference in interim sample size between two searching steps.
#'                 For example, `int_step = 4` means increasing `int_n` by 4 (2 for each group) each time.
#' @param method Design method to use. Options are:
#'               - "logrank": Two-stage log-rank test
#'               - "Simple": Simple RMST difference test
#'               - "Complex": RMST with sculpted critical region.
#' @param lambda_H0 Rate parameter(exponential) under the null hypothesis.
#' @param lambda_H1 Rate parameter(exponential) of experimental group under the alternative hypothesis.
#' @param H1_type Distribution of H1:
#'                - "PH": Proportional hazards
#'                _ "NPH": Non-proportional hazards (piecewise exponential for experimental group)
#' @param HR1 For piecewise exponential only. The hazard ratio before change_time
#' @param HR2 For piecewise exponential only. The hazard ratio after change_time 
#' @param change_time For piecewise exponential only. The time when hazard ratio changes from HR1 to HR2  
#' @param tau Cut-off time for RMST at the final stage. Default is NULL.
#' @param alpha Significance level for the test.
#' @param power Desired power of the test. If NULL, returns the most powerful result.
#'
#' @return A dataframe containing the optimal interim sample size and other relevant design measurements.
#' 
#' @examples 
#' 
#' set.seed(2024)
#' lambda_H1 <- 0.9
#' HR <- 1.7
#' lambda_H0 <- lambda_H1 * HR
#' sim_size <- 5000
#' r <- 60
#' cen_time <- 1
#' alpha <- 0.1
#' power <- 0.85
#' result_scu <- c()
#' # Give a searching range of total sample size N:
#' for (N in seq(from = 74, to = 76, 2))
#' {
#'     n <- ceiling(N / 2) 
#'     acc_time <- N/r
#'     opt_rmst <- get_opt_des(n = n, sim_size = sim_size, 
#'                 acc_time = acc_time, cen_time = cen_time, int_step = 4, 
#'                 lambda_H0 = lambda_H0, lambda_H1 = lambda_H1, H1_type = 'PH', 
#'                 alpha = alpha, power = power, method = 'Complex') 
#'     res <- our_rmst[, c('PET0','PET1','alpha','power',
#'                         'PET','EN0','EN1','EN','interim_n')]
#'     res$N <- c(N)
#'     result_scu <- rbind(result_scu, res)
#'     saveRDS(result_scu, "tem_result_scu.rds")
#' }
#' 
#' # All result in result_scu
#' @export
#' 
get_opt_des <- function(n, sim_size, acc_time, cen_time, int_step, method, lambda_H0, 
                        lambda_H1, H1_type, HR1, HR2, change_time, tau = NULL, alpha, power = NULL) 
{
    N <- 2 * n #overall sample size of two groups
    r <- N / acc_time
    int_factor <- seq(0.3, 0.7, by = int_step / N)  # Each time interim sample size increase by 6
    interim_list <- int_factor * acc_time

    data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                            cen_time = cen_time, arm = 0, interim = interim_list)    
    data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                            cen_time = cen_time, arm = 1, interim = interim_list) 
    if (H1_type == 'PH'){
    data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H1, dist = 'exp', 
                            cen_time = cen_time, arm = 1, interim = interim_list)
    }
    else if (H1_type == 'NPH'){
    data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, 
                            dist = 'pcw_exp', cen_time = cen_time, HR1 = HR1, HR2 = HR2, 
                            change_time = change_time, arm = 1, interim = interim_list)
    }

    all_result <- data.frame()

    if (is.null(tau)){ # the preset cut off tau for final stage
        tau_f <- acc_time + cen_time
    }
    else{
        tau_f <- tau
    }

    for (i in 1 : length(interim_list))
    {  
        interim <- interim_list[i]  # interim is the length of interim period
        if (method != 'logrank')
        {
            rmst_h0_int <- RMST_sim_cal(n = n,data_E = data_E_H0[ , c(2,3,1), i], 
                                data_C = data_C[ , c(2,3,1), i],tau = interim, sim_size = sim_size)
            rmst_h1_int <- RMST_sim_cal(n = n,data_E = data_E_H1[ , c(2,3,1), i], 
                                data_C = data_C[ , c(2,3,1), i],tau = interim, sim_size = sim_size)
            rmst_h0_fin <- RMST_sim_cal(n = n,data_E = data_E_H0[ , c(4,5,1), i], 
                                data_C = data_C[ , c(4,5,1), i], tau = tau_f, sim_size = sim_size)
            rmst_h1_fin <- RMST_sim_cal(n = n,data_E = data_E_H1[ , c(4,5,1), i], 
                                data_C = data_C[ , c(4,5,1), i], tau = tau_f, sim_size = sim_size)
            rmst_data <- rbind(rmst_h0_int, rmst_h1_int, rmst_h0_fin, rmst_h1_fin)

            mu_cov_h0 <- mu_cov_mc(rmst_int = rmst_h0_int, rmst_fin = rmst_h0_fin, sim_size = sim_size)
            mu_cov_h1 <- mu_cov_mc(rmst_int = rmst_h1_int, rmst_fin = rmst_h1_fin, sim_size = sim_size)

            best_our <- adp_grid_src(rmst_data = rmst_data, mu_cov_h0 = mu_cov_h0, mu_cov_h1 = mu_cov_h1, 
                            int_n = interim * r, fin_n = N, sim_size = sim_size, method = method,
                            alpha = alpha, power = power)
        }

        else if (method == 'logrank')   # search the min(E(N)) using log rank test
        {
            try_lr <- tryCatch({
                lr_h0_int <- log_rank_sim(data_C = data_C[, c(2, 3, 1), i], 
                                   data_E = data_E_H0[, c(2, 3, 1), i], 
                                   sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
                lr_h1_int <- log_rank_sim(data_C = data_C[, c(2, 3, 1), i], 
                                   data_E = data_E_H1[, c(2, 3, 1), i], 
                                   sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
                lr_h0_fin <- log_rank_sim(data_C = data_C[, c(4, 5, 1), i], 
                                   data_E = data_E_H0[, c(4, 5, 1), i], 
                                   sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
                lr_h1_fin <- log_rank_sim(data_C = data_C[, c(4, 5, 1), i], 
                                   data_E = data_E_H1[, c(4, 5, 1), i], 
                                   sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
        list(lr_h0_int = lr_h0_int, lr_h1_int = lr_h1_int, lr_h0_fin = lr_h0_fin, lr_h1_fin = lr_h1_fin)
    }, error = function(e) {
        return(NULL)
    })

    if (is.null(try_lr)) {
        next 
    }
        lr_h0_int <- try_lr$lr_h0_int
        lr_h1_int <- try_lr$lr_h1_int
        lr_h0_fin <- try_lr$lr_h0_fin
        lr_h1_fin <- try_lr$lr_h1_fin

        if (is.null(lr_h0_int) || is.null(lr_h1_int) || is.null(lr_h0_fin) || is.null(lr_h1_fin)) {
            next
        }
        # Get W/sigma
        z_stats_h1_int <- lr_h1_int$z_stats
        z_stats_h1_fin <- lr_h1_fin$z_stats
        z_stats_h0_int <- lr_h0_int$z_stats
        z_stats_h0_fin <- lr_h0_fin$z_stats
        logrank_data <- rbind(z_stats_h0_int, z_stats_h1_int, z_stats_h0_fin, z_stats_h1_fin) 

        # corr(W1, W | H0)
        corr_h0 <- sqrt(mean(lr_h0_int$var_w) / mean(lr_h0_fin$var_w))         
        best_our <- find_m_logrank(logrank_data = logrank_data, sim_size = sim_size, corr_h0 = corr_h0,
                            search_times = 150, alpha = alpha, power = power, int_n = interim * r, fin_n = N)
        }
        best_our$interim_n <- ceiling(interim * r)
        all_result <- rbind(all_result, best_our)
    }

    all_result <- na.omit(all_result)    # drop NA 
    if (dim(all_result)[1] == 0) {
        if(method == 'logrank'){
            return(data.frame(m1 = 0, m2 = 0, PET0 = 0, PET1 = 0, alpha = 0, 
                        power = 0, PET = 0, EN0 = NA, EN1 = NA, EN = NA, interim_n = NA))
        }
        else{
             return(data.frame(m1 = 0, q1 = 0, m2 = 0, q2 = 0, gamma = 0, 
                        PET0 = 0, PET1 = 0, alpha = 0, power = 0, 
                        PET = 0, EN0 = NA, EN1 = NA, EN = NA, interim_n = NA))
        }
    }
    if (is.null(power)){  #return the most powerful result
        all_result <- all_result[which(all_result$power == max(all_result$power, na.rm = TRUE)), ]
        return(all_result[1, ])
    }
    else{
        all_result <- all_result[which(all_result$EN0 == min(all_result$EN0, na.rm = TRUE)), ]
        return(all_result[1, ])
        # return the result with minimal EN. only return the first row when multiple results
    }
   
}