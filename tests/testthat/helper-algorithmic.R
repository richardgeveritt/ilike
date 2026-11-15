# Outer model whose likelihood is computed by `method` (kf, pf or enkf); the
# fixture names the inner SSM by a path relative to itself.
algorithmic_model_file <- function(method)
{
  system.file("extdata", paste("test_algorithmic_", method, ".ilike", sep = ""), package = "ilike")
}

# Exact evidence and posterior mean of the transition coefficient a ~ U(0,1)
algorithmic_truth <- function()
{
  f = function(a) sapply(a, function(aa) exp(kf_log_likelihood(ssm_y, a = aa)))
  Z = integrate(f, 0, 1)$value
  c(log_evidence = log(Z), mean = integrate(function(a) a*f(a), 0, 1)$value/Z)
}

# Exact mean and variance of the final state x_T given the data, averaging the
# Kalman filter's final filtering distribution over the posterior of a ~ U(0,1)
algorithmic_final_state_truth <- function()
{
  f = function(a) sapply(a, function(aa) exp(kf_log_likelihood(ssm_y, a = aa)))
  m = function(a) sapply(a, function(aa) kf_final_filter(ssm_y, a = aa)["mean"])
  v = function(a) sapply(a, function(aa) kf_final_filter(ssm_y, a = aa)["var"])
  Z = integrate(f, 0, 1)$value
  mean = integrate(function(a) m(a)*f(a), 0, 1)$value/Z
  second = integrate(function(a) (v(a) + m(a)^2)*f(a), 0, 1)$value/Z
  c(mean = mean, var = second - mean^2)
}

# weighted mean and variance of the nested state x at the final time
nested_final_state_moments <- function(nested)
{
  last = nested[nested$Target == max(nested$Target), ]
  w = exp(last$LogWeight - max(last$LogWeight)); w = w/sum(w)
  m = sum(w*last$Value)
  c(mean = m, var = sum(w*(last$Value - m)^2))
}
