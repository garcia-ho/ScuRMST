% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_opt_des.R
\name{get_opt_des}
\alias{get_opt_des}
\title{Get Optimal Design}
\usage{
get_opt_des(
  n,
  sim_size,
  acc_time,
  cen_time,
  int_step,
  method,
  lambda_H0,
  lambda_H1,
  H1_type,
  HR1,
  HR2,
  change_time,
  tau = NULL,
  alpha,
  power = NULL
)
}
\arguments{
\item{n}{Overall sample size for each group. This should maintain consistency in inputs.}

\item{sim_size}{Number of simulations to perform.}

\item{acc_time}{Accrual time for the study.}

\item{cen_time}{Time until censoring occurs.}

\item{int_step}{The difference in interim sample size between two searching steps.
For example, \code{int_step = 4} means increasing \code{int_n} by 4 (2 for each group) each time.}

\item{method}{Design method to use. Options are:
- "logrank": Two-stage log-rank test
- "Simple": Simple RMST difference test
- "Complex": RMST with sculpted critical region.}

\item{lambda_H0}{Rate parameter(exponential) under the null hypothesis.}

\item{lambda_H1}{Rate parameter(exponential) of experimental group under the alternative hypothesis.}

\item{H1_type}{Distribution of H1:
- "PH": Proportional hazards
_ "NPH": Non-proportional hazards (piecewise exponential for experimental group)}

\item{HR1}{For piecewise exponential only. The hazard ratio before change_time}

\item{HR2}{For piecewise exponential only. The hazard ratio after change_time}

\item{change_time}{For piecewise exponential only. The time when hazard ratio changes from HR1 to HR2}

\item{tau}{Cut-off time for RMST at the final stage. Default is NULL.}

\item{alpha}{Significance level for the test.}

\item{power}{Desired power of the test. If NULL, returns the most powerful result.}
}
\value{
A dataframe containing the optimal interim sample size and other relevant design measurements.
}
\description{
This function search for the optimal design for a study, returning the interim sample size, total sample size
and other measurements of the design minimizing the expected sample size under H0: EN0.
For any overall sample size N, the function only search for interim sample size ranging from 0.3\emph{N to 0.7}N.
}
\examples{

set.seed(2024)
lambda_H1 <- 0.9
HR <- 1.7
lambda_H0 <- lambda_H1 * HR
sim_size <- 5000
r <- 60
cen_time <- 1
alpha <- 0.1
power <- 0.85
result_scu <- c()
# Give a searching range of total sample size N:
for (N in seq(from = 74, to = 76, 2))
{
    n <- ceiling(N / 2) 
    acc_time <- N/r
    opt_rmst <- get_opt_des(n = n, sim_size = sim_size, 
                acc_time = acc_time, cen_time = cen_time, int_step = 4, 
                lambda_H0 = lambda_H0, lambda_H1 = lambda_H1, H1_type = 'PH', 
                alpha = alpha, power = power, method = 'Complex') 
    res <- our_rmst[, c('PET0','PET1','alpha','power',
                        'PET','EN0','EN1','EN','interim_n')]
    res$N <- c(N)
    result_scu <- rbind(result_scu, res)
    saveRDS(result_scu, "tem_result_scu.rds")
}

# All result in result_scu
}
