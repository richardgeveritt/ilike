# Distributions added to model statements (dev/plan.md C214): invgamma, t,
# cauchy, halfnorm, halfcauchy, halft. Each is checked against R: its log
# density, its gradients (finite differences), its sampler (quantiles), and
# the transform moves choose for it.

distributions_model = function(with_data)
{
  c("ilike 1.0",
    if (with_data) "data { y = [0.7, -1.3, 2.4]; u = [0.2, 1.9]; e = [0.5, 3.1] }",
    "model {",
    "  nu ~ gamma(4, 1)",
    "  m ~ cauchy(0.5, 2)",
    "  s ~ halfnorm(1.5)",
    "  c ~ halfcauchy(0.7)",
    "  h ~ halft(3, 1.2)",
    "  v ~ invgamma(3, 2)",
    "  w ~ t(5, -1, 0.8)",
    if (with_data) c("  for (i in 1:3) y[i] ~ t(nu, m, s)",
                     "  for (j in 1:2) u[j] ~ cauchy(w, c)",
                     "  for (j in 1:2) e[j] ~ halft(nu, h)",
                     "  for (j in 1:2) u[j] ~ invgamma(v, s)"),
    "}")
}

distributions_point = list(nu = 4.3, m = 0.4, s = 1.1, c = 0.6, h = 0.9, v = 1.7, w = -0.6)

test_that("the new distributions' log densities and gradients agree with R and with finite differences", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(distributions_model(TRUE), f)
  recipe = ilike::compile(f)
  P = distributions_point
  dt_ls = function(x, nu, m, s) dt((x - m)/s, nu, log = TRUE) - log(s)
  dinvgamma = function(x, a, b) a*log(b) - lgamma(a) - (a + 1)*log(x) - b/x
  y = c(0.7, -1.3, 2.4); u = c(0.2, 1.9); e = c(0.5, 3.1)
  truth = dgamma(P$nu, 4, 1, log = TRUE) + dcauchy(P$m, 0.5, 2, log = TRUE) + log(2) + dnorm(P$s, 0, 1.5, log = TRUE) +
    log(2) + dcauchy(P$c, 0, 0.7, log = TRUE) + log(2) + dt_ls(P$h, 3, 0, 1.2) + dinvgamma(P$v, 3, 2) + dt_ls(P$w, 5, -1, 0.8) +
    sum(dt_ls(y, P$nu, P$m, P$s)) + sum(dcauchy(u, P$w, P$c, log = TRUE)) + sum(log(2) + dt_ls(e, P$nu, 0, P$h)) +
    sum(dinvgamma(u, P$v, P$s))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
  # outside the support: density 0
  expect_equal(ilike::evaluate_log(recipe, parameters = replace(P, "c", -0.1)), -Inf)
  expect_equal(ilike::evaluate_log(recipe, parameters = replace(P, "v", -0.1)), -Inf)
  # gradients with respect to every variable, including the distributions' arguments
  for (name in names(P))
  {
    g = ilike:::evaluate_gradient_log(recipe, P, name)
    h = 1e-5
    fd = (ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] + h)) -
          ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] - h)))/(2*h)
    expect_lt(abs(g[1] - fd)/(1 + abs(fd)), 1e-6, label = name)
  }
})

test_that("the new distributions are simulated from their priors", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "distributions_prior"), recursive = TRUE), add = TRUE)
  writeLines(distributions_model(FALSE), f)
  unlink("distributions_prior", recursive = TRUE)
  ilike::IS(ilike::compile(f), results_name = "distributions_prior", number_of_importance_points = 20000, seed = 1)
  o = ilike::load_smc_output("distributions_prior")
  p = c(0.25, 0.5, 0.75)
  quantiles = list(m = qcauchy(p, 0.5, 2), s = 1.5*qnorm((1 + p)/2), c = 0.7*qcauchy((1 + p)/2), h = 1.2*qt((1 + p)/2, 3),
                   v = 1/qgamma(1 - p, 3, 2), w = -1 + 0.8*qt(p, 5))
  for (name in names(quantiles))
  {
    x = o$Value[o$ParameterName == name]
    expect_equal(length(x), 20000)
    target = quantiles[[name]]
    # (within 3% of the interquartile range)
    expect_lt(max(abs(quantile(x, p) - target))/(target[3] - target[1]), 0.03, label = name)
  }
  expect_true(all(o$Value[o$ParameterName %in% c("s", "c", "h", "v")] > 0))
})

test_that("a variable with a half or inverse-gamma distribution moves on the log scale", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "distributions_nuts"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  s ~ halfcauchy(1)", "  v ~ invgamma(3, 2)", "}", "algorithm mcmc {",
               "  move m = nuts([s, v], step = 0.5)", "  termination = iterations(4000)", "}"), f)
  unlink("distributions_nuts", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "distributions_nuts", seed = 1, initial_values = list(list(s = 1, v = 1)))
  o = ilike::load_mcmc_output("distributions_nuts")
  s = o$Value[o$ParameterName == "s" & o$Iteration > 400]
  v = o$Value[o$ParameterName == "v" & o$Iteration > 400]
  expect_true(all(s > 0) && all(v > 0))
  # median of halfcauchy(1) is 1; mean of invgamma(3, 2) is 1
  expect_equal(median(s), 1, tolerance = 0.1)
  expect_equal(mean(v), 1, tolerance = 0.1)
})

# tnorm in the registry (prior and independent proposal) and the truncated
# random walk tnorm_rw (C217)
test_that("the registry's tnorm prior has the right density, also far in a tail, and is simulated correctly", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "tnorm_prior"), recursive = TRUE), add = TRUE)
  # a one-sided limit, and limits far in the upper tail (the erf difference would cancel)
  writeLines(c("ilike 1.0", "model {", "  a ~ tnorm(1, 2, 0, Inf)", "  b ~ tnorm(0, 1, 8, 9)", "}"), f)
  recipe = ilike::compile(f)
  log_tnorm = function(x, m, s, l, u) dnorm(x, m, s, log = TRUE) - log(pnorm((u - m)/s) - pnorm((l - m)/s))
  expect_equal(ilike::evaluate_log(recipe, parameters = list(a = 0.7, b = 8.2)),
               log_tnorm(0.7, 1, 2, 0, Inf) + dnorm(8.2, log = TRUE) - pnorm(8, lower.tail = FALSE, log.p = TRUE) -
                 log1p(-exp(pnorm(9, lower.tail = FALSE, log.p = TRUE) - pnorm(8, lower.tail = FALSE, log.p = TRUE))),
               tolerance = 1e-10)
  expect_equal(ilike::evaluate_log(recipe, parameters = list(a = -0.1, b = 8.2)), -Inf)
  unlink("tnorm_prior", recursive = TRUE)
  ilike::IS(recipe, results_name = "tnorm_prior", number_of_importance_points = 20000, seed = 1)
  o = ilike::load_smc_output("tnorm_prior")
  p = c(0.25, 0.5, 0.75)
  qa = 1 + 2*qnorm(pnorm(-0.5) + p*(1 - pnorm(-0.5)))
  a = o$Value[o$ParameterName == "a"]
  expect_lt(max(abs(quantile(a, p) - qa))/(qa[3] - qa[1]), 0.03)
  b = o$Value[o$ParameterName == "b"]
  expect_true(all(b >= 8 & b <= 9))
})

test_that("tnorm_rw and independent_mh(tnorm(...)) target a bounded posterior", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "tnorm_moves"), recursive = TRUE), add = TRUE)
  # p ~ unif(0, 1), 7 successes in 10: p | z ~ beta(8, 4)
  for (move in c("metropolis_hastings(tnorm_rw(p, 0.3, 0, 1))", "independent_mh(tnorm(p, 0.65, 0.2, 0, 1))"))
  {
    writeLines(c("ilike 1.0", "data { z = 7 }", "model {", "  p ~ unif(0, 1)", "  z ~ binom(10, p)", "}", "algorithm mcmc {",
                 paste0("  move m = ", move), "  termination = iterations(20000)", "}"), f)
    unlink("tnorm_moves", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "tnorm_moves", seed = 1, initial_values = list(list(p = 0.5)))
    x = ilike::load_mcmc_output("tnorm_moves")$Value[-(1:1000)]
    expect_true(all(x > 0 & x < 1), label = move)
    expect_equal(mean(x), 8/12, tolerance = 0.01, label = move)
    expect_equal(sd(x), sqrt(8*4/(12^2*13)), tolerance = 0.05, label = move)
  }
})

test_that("tnorm_rw is refused in a metropolis move (it is not symmetric)", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { z = 7 }", "model {", "  p ~ unif(0, 1)", "  z ~ binom(10, p)", "}", "algorithm mcmc {",
               "  move m = metropolis(tnorm_rw(p, 0.3, 0, 1))", "  termination = iterations(10)", "}"), f)
  expect_error(ilike::compile(f), "tnorm_rw")
})

test_that("a tnorm statement's density and gradient stay right with the mean far below the limit", {
  # (C218) y ~ tnorm(mu, 1, 0, Inf) with mu 9 to 40 sds below 0: the
  # normalising constant was an erf difference, which cancelled to 0 (log
  # density +Inf, which Metropolis-Hastings would accept and stick at)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.3, 1.2] }", "model {", "  mu ~ norm(0, 10)", "  for (i in 1:2) y[i] ~ tnorm(mu, 1, 0, Inf)", "}"), f)
  recipe = ilike::compile(f)
  y = c(0.3, 1.2)
  for (mu in c(-9, -20, -40))
  {
    truth = dnorm(mu, 0, 10, log = TRUE) + sum(dnorm(y, mu, 1, log = TRUE) - pnorm(0, mu, 1, lower.tail = FALSE, log.p = TRUE))
    expect_equal(ilike::evaluate_log(recipe, parameters = list(mu = mu)), truth, tolerance = 1e-10, label = mu)
    h = 1e-5
    fd = (ilike::evaluate_log(recipe, parameters = list(mu = mu + h)) - ilike::evaluate_log(recipe, parameters = list(mu = mu - h)))/(2*h)
    expect_lt(abs(ilike:::evaluate_gradient_log(recipe, list(mu = mu), "mu")[1] - fd), 1e-5, label = mu)
  }
})
