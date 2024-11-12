#' 100 points simulated from a Gaussian distribution.
#'
#' Generated using rnorm(100).
#'
#' @format A vector of 100 variables
"gaussian"

#' 10000 points simulated from a constant velocity linear-Gaussian dynamic model (process noise standard deviation = 1).
#'
#' @format A 10000*2 matrix.
"constant_velocity_x"

#' 10000 noisy observations (Gaussian noise, standard deviation 1) of 10000 points simulated from a constant velocity linear-Gaussian dynamic model (process noise standard deviation = 1).
#'
#' @format A 10000*1 matrix.
"constant_velocity_y"

#' 100000 points simulated from a two-dimensional posterior using MCMC, in the format required by the ilike.output and ggsmc packages.
#'
#' @format A 200000 row data frame.
"toy_model_rwm_ilike"

#' 100000 points simulated from a two-dimensional posterior using MCMC, in the format required by the ggmcmc package.
#'
#' @format A 200000 row data frame.
"toy_model_rwm_ggmcmc"
