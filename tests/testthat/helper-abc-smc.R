# Run SMC with MCMC moves on an sbi model with the ABC tolerance as the SMC
# sequence variable; returns the log evidence and the weighted mean and
# variance of theta at the final target.
run_abc_smc <- function(model, algorithm, results_name, model_parameter_list = list(), seed = 1)
{
  recipe = ilike::compile(system.file("extdata", c(model, algorithm), package = "ilike"),
                          model_parameter_list = model_parameter_list)
  z = ilike::SMC_with_MCMC(recipe, results_name = results_name, number_of_particles = 2000,
                           smc_iterations_to_store = 1000, model_parameter_list = model_parameter_list,
                           seed = seed)
  o = ilike::load_smc_output(results_name)
  final = max(o$Target)
  last = o[o$Target == final & o$Iteration == max(o$Iteration[o$Target == final]), ]
  w = exp(last$LogWeight - max(last$LogWeight)); w = w/sum(w)
  m = sum(w*last$Value)
  list(log_evidence = z, mean = m, var = sum(w*(last$Value - m)^2),
       tolerances = unique(o$TargetParameters))
}

# Exact log evidence, posterior mean and variance with prior N(1,1) for an ABC
# likelihood lik(theta), by quadrature
abc_truth <- function(lik)
{
  f = function(t) dnorm(t, 1, 1)*lik(t)
  Z = integrate(f, -Inf, Inf)$value
  m = integrate(function(t) t*f(t), -Inf, Inf)$value/Z
  c(log_evidence = log(Z), mean = m, var = integrate(function(t) (t - m)^2*f(t), -Inf, Inf)$value/Z)
}

# Gaussian-kernel ABC likelihood for observed summary 1 with s ~ N(theta, v)
gaussian_abc_likelihood <- function(epsilon, v = 1)
{
  function(t) dnorm(1, t, sqrt(v + epsilon^2))
}

expect_abc_smc <- function(result, truth)
{
  expect_equal(result$log_evidence, unname(truth["log_evidence"]), tolerance = 0.03)
  expect_equal(result$mean, unname(truth["mean"]), tolerance = 0.1)
  expect_equal(result$var, unname(truth["var"]), tolerance = 0.1)
}
