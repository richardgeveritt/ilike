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
  P = length(log_weights)
  W = exp(log_weights)
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

  indices
}
