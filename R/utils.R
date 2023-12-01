#' Gets a seed for random number generation by using the processor time stamp (calling rdtsc in C).
#'
#' @return An integer.
#' @export
rdtsc_seed <- function()
{
  return(ilike_rdtsc())
}


log_sum_exp <- function(x){
  xmax = which.max(x)
  if (sum(xmax)==0)
  {
    -Inf
  }
  else
  {
    result = log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
    if (is.nan(result))
    {
      -Inf
    }
    else
    {
      result
    }
  }
}


stratified_resample_u<-function(log_weights,u)
{
  norm_log_weights = log_weights - log_sum_exp(log_weights)
  P = length(norm_log_weights)
  W = exp(norm_log_weights)
  cw = cumsum(W / sum(W))
  n = length(W)
  u = u/n
  v = u + (0:(n-1))/n

  j = 1
  indices = matrix(0,n)
  for (i in 1:n)
  {
    while(cw[j] < v[i])
    {
      j = j+1
    }
    indices[i] = j
  }

  return(indices)
}


stratified_resample<-function(log_weights)
{
  n = length(log_weights)
  u = stats::runif(n,0,1)/n
  return(stratified_resample_u(log_weights,u))
}


# can take unnormalised weights
ess <- function(log_weights)
{
  the_ess = exp(2*log_sum_exp(log_weights) - log_sum_exp(2*log_weights))

  if (is.nan(the_ess))
  {
    return(0)
  }
  else
  {
    return(the_ess)
  }
}


cess <- function(normalised_log_weights, log_incremental_weights)
{
  if (all(log_incremental_weights==-Inf))
  {
    return(0)
  }
  else if (all(log_incremental_weights==Inf))
  {
    return(length(log_incremental_weights))
  }
  else
  {
    return(exp(log(length(log_incremental_weights)) + 2*(log_sum_exp(normalised_log_weights + log_incremental_weights)) - log_sum_exp(normalised_log_weights + 2*log_incremental_weights)))
  }
}
