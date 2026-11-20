# Converting Stan programs to the ilike language (dev/plan.md 3.4, C230): the
# translated models' log densities against R's, and refusals.

stan_compile = function(program, data)
{
  f = tempfile(fileext = ".ilike")
  ilike::convert_stan_model(program, data, output = f)
  ilike::compile(f)
}

test_that("eight schools is translated, with the constraint on tau giving a half-Cauchy prior", {
  program = '
data { int<lower=0> J; array[J] real y; array[J] real<lower=0> sigma; }
parameters { real mu; real<lower=0> tau; vector[J] eta; }
transformed parameters { vector[J] theta = mu + tau * eta; }
model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 5);
  eta ~ std_normal();
  y ~ normal(theta, sigma);
}
generated quantities { real z = mu; }'
  d = list(J = 8, y = c(28, 8, -3, 7, -1, 1, 18, 12), sigma = c(15, 10, 16, 11, 9, 11, 10, 18))
  expect_warning(text <- ilike::convert_stan_model(program, d), "generated quantities")
  expect_match(text, "tau ~ halfcauchy(5)", fixed = TRUE)
  expect_match(text, "move m = nuts([mu, tau, eta]", fixed = TRUE)
  recipe = suppressWarnings(stan_compile(program, d))
  P = list(mu = 1.5, tau = 3.2, eta = seq(-1, 1, length = 8))
  theta = P$mu + P$tau*P$eta
  truth = dnorm(P$mu, 0, 5, log = TRUE) + log(2) + dcauchy(P$tau, 0, 5, log = TRUE) + sum(dnorm(P$eta, log = TRUE)) +
    sum(dnorm(d$y, theta, d$sigma, log = TRUE))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
})

test_that("a regression with a matrix-vector product, logistic and negative binomial models", {
  set.seed(1)
  X = matrix(round(rnorm(12), 2), 6); y = round(rnorm(6), 2)
  program = '
data { int N; int K; matrix[N, K] X; vector[N] y; }
parameters { real alpha; vector[K] beta; real<lower=0> sigma; }
model {
  alpha ~ normal(0, 10); beta ~ normal(0, 2); sigma ~ exponential(1);
  y ~ normal(alpha + X * beta, sigma);
}'
  recipe = stan_compile(program, list(N = 6, K = 2, X = X, y = y))
  P = list(alpha = 0.3, beta = c(0.5, -0.7), sigma = 1.3)
  truth = dnorm(0.3, 0, 10, log = TRUE) + sum(dnorm(P$beta, 0, 2, log = TRUE)) + dexp(1.3, 1, log = TRUE) +
    sum(dnorm(y, P$alpha + X %*% P$beta, P$sigma, log = TRUE))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
  # bernoulli_logit with target += ..._lpmf, and neg_binomial_2
  program = '
data { int N; array[N] int z; vector[N] x; array[N] int counts; }
parameters { real a; real b; real<lower=0> phi; }
model {
  a ~ normal(0, 2);
  b ~ student_t(3, 0, 1);
  phi ~ gamma(2, 0.5);
  target += bernoulli_logit_lpmf(z | a + b * x);
  for (n in 1:N) counts[n] ~ neg_binomial_2(exp(a), phi);
}'
  d = list(N = 5, z = c(1, 0, 1, 1, 0), x = c(-1, 0.5, 2, 0.1, -0.3), counts = c(3, 0, 7, 2, 1))
  recipe = stan_compile(program, d)
  P = list(a = 0.4, b = -0.8, phi = 2.5)
  truth = dnorm(0.4, 0, 2, log = TRUE) + dt(-0.8, 3, log = TRUE) + dgamma(2.5, 2, 0.5, log = TRUE) +
    sum(dbinom(d$z, 1, plogis(P$a + P$b*d$x), log = TRUE)) + sum(dnbinom(d$counts, size = P$phi, mu = exp(P$a), log = TRUE))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
})

test_that("truncation, bounded parameters without priors, target += expressions and a simplex", {
  program = '
data { int N; array[N] int c; real w; }
parameters { real<lower=-1, upper=2> u; real<lower=0> s; simplex[3] p; }
model {
  s ~ normal(1, 2);
  w ~ normal(u, s) T[0, ];
  target += -0.5 * square(u);
  c ~ categorical(p);
}'
  d = list(N = 4, c = c(1, 3, 3, 2), w = 0.7)
  text = ilike::convert_stan_model(program, d)
  expect_match(text, "u ~ unif(-1, 2)", fixed = TRUE)
  expect_match(text, "p ~ dirichlet([1, 1, 1])", fixed = TRUE)
  expect_match(text, "s ~ tnorm(1, 2, 0, Inf)", fixed = TRUE)
  expect_match(text, "move m = nuts([u, s, p]", fixed = TRUE)
  recipe = stan_compile(program, d)
  P = list(u = 0.4, s = 1.1, p = c(0.2, 0.5, 0.3))
  truth = dunif(0.4, -1, 2, log = TRUE) + log(gamma(3)) + dnorm(1.1, 1, 2, log = TRUE) - pnorm(0, 1, 2, lower.tail = FALSE, log.p = TRUE) +
    dnorm(0.7, 0.4, 1.1, log = TRUE) - pnorm(0, 0.4, 1.1, lower.tail = FALSE, log.p = TRUE) - 0.5*0.4^2 + sum(log(P$p[d$c]))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
})

test_that("what cannot be translated is refused with a message", {
  d = list(N = 2, y = c(1, 2))
  refuse = function(program, message) expect_error(ilike::convert_stan_model(program, d), message, fixed = TRUE)
  refuse('data { int N; vector[N] y; } parameters { real mu; } model { y ~ normal(mu, 1); }',
         "the parameter mu has no prior")
  refuse('data { int N; vector[N] y; } parameters { real<lower=0> s; } model { s ~ cauchy(1, 2); y ~ normal(0, s); }',
         "the prior of s is cauchy on a constrained range")
  refuse('data { int N; vector[N] y; } parameters { real mu; } model { mu ~ normal(0, 1); if (mu > 0) y ~ normal(mu, 1); }',
         "\"if\" statements cannot be translated")
  refuse('functions { real f(real x) { return x; } } data { int N; } parameters { real mu; } model { mu ~ normal(0, 1); }',
         "functions blocks")
  refuse('data { int N; vector[N] y; } parameters { real mu; } model { mu ~ normal(0, 1); target += sum(y); }',
         "the function sum cannot be translated")
  refuse('data { int N; vector[N] y; vector[N] z; } parameters { real mu; } model { mu ~ normal(0, 1); y ~ normal(mu, 1); }',
         "the data z is not in the data given")
})
