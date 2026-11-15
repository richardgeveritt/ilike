# log density of MVN(m, S) at x, without needing the mvtnorm package
log_dmvnorm <- function(x, m, S)
{
  r = x - m
  -0.5*(length(x)*log(2*pi) + log(det(S)) + sum(r*solve(S, r)))
}

# posterior mean and variance of theta for prior density `prior` on (0, Inf)
# and one observation 1 with likelihood N(theta, 1), by quadrature
positive_posterior_moments <- function(prior)
{
  f = function(t) prior(t)*dnorm(1, t, 1)
  Z = integrate(f, 0, Inf)$value
  m = integrate(function(t) t*f(t), 0, Inf)$value/Z
  v = integrate(function(t) (t - m)^2*f(t), 0, Inf)$value/Z
  c(mean = m, var = v)
}

run_is_fixture <- function(files, results_name, number_of_importance_points = 5000, seed = 1)
{
  recipe = ilike::compile(system.file("extdata", files, package = "ilike"))
  ilike::IS(recipe, results_name = results_name,
            number_of_importance_points = number_of_importance_points, seed = seed)
}
