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
