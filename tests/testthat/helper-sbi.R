# Likelihood of theta, for one observation 1, implied by each simulation-based
# method (simulator N(theta,1); epsilon = 0.5 for ABC)
sbi_likelihood_function <- function(kind)
{
  switch(kind,
         sl = function(t) dnorm(1, t, 1),
         gaussian_abc = function(t) dnorm(1, t, sqrt(1 + 0.5^2)),
         # uniform kernel on |x - 1| < 0.5, normalised by its width 1
         function(t) (pnorm(1.5 - t) - pnorm(0.5 - t))/1)
}

# log evidence, posterior mean and variance with prior N(1, 1), by quadrature
sbi_truth <- function(kind)
{
  lik = sbi_likelihood_function(kind)
  f = function(t) dnorm(t, 1, 1)*lik(t)
  Z = integrate(f, -Inf, Inf)$value
  m = integrate(function(t) t*f(t), -Inf, Inf)$value/Z
  v = integrate(function(t) (t - m)^2*f(t), -Inf, Inf)$value/Z
  c(log_evidence = log(Z), mean = m, var = v)
}
