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

# Exact Kalman filter log-likelihood for test_enkf_bivariate.ilike
kf_bivariate_log_likelihood <- function()
{
  A = matrix(c(0.9, -0.1, 0.2, 0.8), 2); Q = matrix(c(0.5, 0.2, 0.2, 0.4), 2)
  H = matrix(c(1, 0.3, 0.5, 1), 2); R = matrix(c(1, 0.6, 0.6, 1), 2)
  y = matrix(c(-1.3703, 0.1127, 1.5672, -7e-04, -0.1496, 0.1355, 0.4484, 1.6447, -0.7621, 1.3297,
               0.7706, 0.6335, 1.1424, -0.9968, 0.1188, -1.3231, -0.7974, 1.1012, -0.6907, 1.4045), 10, 2)
  m = c(0, 0); P = diag(2); ll = 0
  for (t in 1:10)
  {
    if (t > 1) { m = A %*% m; P = A %*% P %*% t(A) + Q }
    S = H %*% P %*% t(H) + R; e = y[t, ] - H %*% m
    ll = ll - 0.5*(2*log(2*pi) + log(det(S)) + sum(e*solve(S, e)))
    K = P %*% t(H) %*% solve(S); m = m + K %*% e; P = (diag(2) - K %*% H) %*% P
  }
  ll
}
