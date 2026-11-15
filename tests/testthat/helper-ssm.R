run_kf_fixture <- function(fixture, results_name)
{
  filename = system.file("extdata", fixture, package = "ilike")
  ilike::KF(ilike::compile(filename), results_name = results_name)
}

# Exact Kalman filter log-likelihood for the scalar model in
# test_kf_linear_gaussian.ilike (x_0 ~ N(m0,v0), x_t = a x_{t-1} + N(0,q),
# y_t = x_t + N(0,r)), with y_0 observed at the first index.
kf_log_likelihood <- function(y, a = 0.9, q = 0.5, r = 1, m0 = 0, v0 = 1)
{
  m = m0; v = v0; ll = 0
  for (t in seq_along(y))
  {
    if (t > 1) { m = a*m; v = a*a*v + q }
    s = v + r
    ll = ll + dnorm(y[t], m, sqrt(s), log = TRUE)
    k = v/s; m = m + k*(y[t] - m); v = (1 - k)*v
  }
  ll
}

run_pf_fixture <- function(fixture, results_name, number_of_particles = 5000, seed = 1)
{
  filename = system.file("extdata", fixture, package = "ilike")
  ilike::PF(ilike::compile(filename), results_name = results_name,
            number_of_particles = number_of_particles, seed = seed)
}

# Exact final filtering mean and variance for the same model
kf_final_filter <- function(y, a = 0.9, q = 0.5, r = 1, m0 = 0, v0 = 1)
{
  m = m0; v = v0
  for (t in seq_along(y))
  {
    if (t > 1) { m = a*m; v = a*a*v + q }
    s = v + r; k = v/s; m = m + k*(y[t] - m); v = (1 - k)*v
  }
  c(mean = m, var = v)
}

ssm_y <- c(0.5, 1.0, 0.2, -0.3, 0.8)
