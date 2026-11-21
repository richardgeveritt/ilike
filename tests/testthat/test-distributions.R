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

# (C227) negbinom, geom, categorical, vonmises; multinom, dirichlet, mvt,
# wishart, invwishart, lkj
test_that("the second group of distributions have the right log densities", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  S = matrix(c(2, 0.5, 0.5, 1), 2)
  writeLines(c("ilike 1.0",
               "data { n1 = 3; g1 = 2; c1 = 2; a1 = 0.4; m1 = [1, 3, 2]; d1 = [0.2, 0.5, 0.3]; t1 = [0.5, -1]; W = [[1.5, 0.3], [0.3, 0.8]]; R = [[1, 0.3], [0.3, 1]] }",
               "model {", "  q ~ unif(0, 1)", "  k ~ gamma(2, 1)", "  mu ~ norm(0, 1)",
               "  n1 ~ negbinom(k, q)", "  g1 ~ geom(q)", "  c1 ~ categorical([0.2, 0.5, 0.3])", "  a1 ~ vonmises(mu, k)",
               "  m1 ~ multinom(6, [0.2, 0.5, 0.3])", "  d1 ~ dirichlet([2, 3, 1.5])", "  t1 ~ mvt(4, [0.2, 0], [[2, 0.5], [0.5, 1]])",
               "  W ~ wishart(5, [[2, 0.5], [0.5, 1]])", "  W ~ invwishart(6, [[2, 0.5], [0.5, 1]])", "  R ~ lkj(2.5)", "}"), f)
  recipe = ilike::compile(f)
  P = list(q = 0.35, k = 1.7, mu = 0.2)
  ldirichlet = function(x, a) lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1)*log(x))
  lmvt = function(x, nu, m, S) { p = length(x); q = sum((x - m)*solve(S, x - m))
    lgamma((nu + p)/2) - lgamma(nu/2) - p/2*log(nu*pi) - 0.5*determinant(S)$modulus[1] - (nu + p)/2*log1p(q/nu) }
  lmvgamma = function(a, p) p*(p - 1)/4*log(pi) + sum(lgamma(a + (1 - seq_len(p))/2))
  lwishart = function(X, nu, S) { p = nrow(X); (nu - p - 1)/2*determinant(X)$modulus[1] - 0.5*sum(diag(solve(S, X))) -
    nu*p/2*log(2) - nu/2*determinant(S)$modulus[1] - lmvgamma(nu/2, p) }
  linvwishart = function(X, nu, Psi) { p = nrow(X); nu/2*determinant(Psi)$modulus[1] - nu*p/2*log(2) - lmvgamma(nu/2, p) -
    (nu + p + 1)/2*determinant(X)$modulus[1] - 0.5*sum(diag(Psi %*% solve(X))) }
  W = matrix(c(1.5, 0.3, 0.3, 0.8), 2)
  # lkj in two dimensions: (1 - r^2)^(eta - 1)/(2^(2 eta - 1) B(eta, eta))
  truth = 0 + dgamma(P$k, 2, 1, log = TRUE) + dnorm(P$mu, log = TRUE) +
    dnbinom(3, P$k, P$q, log = TRUE) + dgeom(2, P$q, log = TRUE) + log(0.5) +
    (P$k*cos(0.4 - P$mu) - log(2*pi*besselI(P$k, 0))) + dmultinom(c(1, 3, 2), 6, c(0.2, 0.5, 0.3), log = TRUE) +
    ldirichlet(c(0.2, 0.5, 0.3), c(2, 3, 1.5)) + lmvt(c(0.5, -1), 4, c(0.2, 0), S) +
    lwishart(W, 5, S) + linvwishart(W, 6, S) + (1.5*log(1 - 0.09) - ((2*2.5 - 1)*log(2) + lbeta(2.5, 2.5)))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
  # gradients with respect to the parameters in the arguments, against finite differences
  for (name in names(P))
  {
    h = 1e-6
    fd = (ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] + h)) -
          ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] - h)))/(2*h)
    if (name == "k" || name == "q" || name == "mu")
      expect_equal(ilike:::evaluate_gradient_log(recipe, P, name)[1], fd, tolerance = 1e-5, label = name)
  }
})

test_that("the gradients of the vector distributions", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { m1 = [1, 3, 2]; c1 = 2; t1 = [0.5, -1] }", "model {",
               "  p ~ dirichlet([2, 3, 1.5])", "  a ~ mvnorm([1, 1, 1], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])",
               "  x ~ dirichlet(a)", "  m ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
               "  m1 ~ multinom(6, p)", "  c1 ~ categorical(p)", "  t1 ~ mvt(4, m, [[2, 0.5], [0.5, 1]])", "}"), f)
  recipe = ilike::compile(f)
  P = list(p = c(0.2, 0.5, 0.3), a = c(1.5, 2, 2.5), x = c(0.3, 0.3, 0.4), m = c(0.2, -0.4))
  # p and x: on the simplex, so the formula rather than finite differences
  expect_equal(as.vector(ilike:::evaluate_gradient_log(recipe, P, "p")),
               (c(2, 3, 1.5) - 1)/P$p + c(1, 3, 2)/P$p + c(0, 1/0.5, 0))
  expect_equal(as.vector(ilike:::evaluate_gradient_log(recipe, P, "x")), (P$a - 1)/P$x)
  for (name in c("a", "m"))
  {
    g = as.vector(ilike:::evaluate_gradient_log(recipe, P, name))
    for (j in seq_along(P[[name]]))
    {
      h = 1e-6; Pp = P; Pm = P
      Pp[[name]][j] = Pp[[name]][j] + h; Pm[[name]][j] = Pm[[name]][j] - h
      fd = (ilike::evaluate_log(recipe, parameters = Pp) - ilike::evaluate_log(recipe, parameters = Pm))/(2*h)
      expect_equal(g[j], fd, tolerance = 1e-5, label = paste(name, j))
    }
  }
})

test_that("the second group of distributions are simulated from their priors", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "distributions_prior2"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  nb ~ negbinom(3, 0.4)", "  ge ~ geom(0.3)", "  ca ~ categorical([0.2, 0.5, 0.3])",
               "  vm ~ vonmises(1, 2)", "  mn ~ multinom(10, [0.2, 0.5, 0.3])", "  di ~ dirichlet([2, 3, 1.5])",
               "  mt ~ mvt(6, [1, -1], [[2, 0.5], [0.5, 1]])", "  wi ~ wishart(5, [[2, 0.5], [0.5, 1]])",
               "  iw ~ invwishart(7, [[2, 0.5], [0.5, 1]])", "  lk ~ lkj(2, 3)", "}"), f)
  unlink("distributions_prior2", recursive = TRUE)
  ilike::IS(ilike::compile(f), results_name = "distributions_prior2", number_of_importance_points = 20000, seed = 1)
  o = ilike::load_smc_output("distributions_prior2")
  v = function(name, k = 1) o$Value[o$ParameterName == name & o$Dimension == k]
  near = function(x, target, sd, label) expect_lt(abs(x - target)/sd, 5, label = label) # (5 Monte Carlo sds)
  n = 20000
  near(mean(v("nb")), 3*0.6/0.4, sqrt(3*0.6/0.4^2/n), "negbinom mean")
  near(mean(v("ge")), 0.7/0.3, sqrt(0.7/0.3^2/n), "geom mean")
  expect_equal(as.vector(table(v("ca")))/n, c(0.2, 0.5, 0.3), tolerance = 0.03)
  # the mean resultant: E cos(x - mu) = I1(kappa)/I0(kappa)
  near(mean(cos(v("vm") - 1)), besselI(2, 1)/besselI(2, 0), 0.5/sqrt(n), "vonmises")
  expect_equal(c(mean(v("mn", 1)), mean(v("mn", 2)), mean(v("mn", 3))), 10*c(0.2, 0.5, 0.3), tolerance = 0.02)
  expect_true(all(abs(v("mn", 1) + v("mn", 2) + v("mn", 3) - 10) < 1e-12))
  expect_equal(c(mean(v("di", 1)), mean(v("di", 2)), mean(v("di", 3))), c(2, 3, 1.5)/6.5, tolerance = 0.02)
  expect_equal(c(mean(v("mt", 1)), mean(v("mt", 2))), c(1, -1), tolerance = 0.05)
  expect_equal(var(v("mt", 1)), 2*6/4, tolerance = 0.1)
  # wishart: mean nu S; inverse wishart: Psi/(nu - p - 1) (elements column by column)
  expect_equal(c(mean(v("wi", 1)), mean(v("wi", 2)), mean(v("wi", 4))), 5*c(2, 0.5, 1), tolerance = 0.03)
  expect_equal(c(mean(v("iw", 1)), mean(v("iw", 2)), mean(v("iw", 4))), c(2, 0.5, 1)/4, tolerance = 0.05)
  # lkj(2) in 3 dimensions: unit diagonal; each correlation 2 Beta(b, b) - 1, b = 2.5: variance 1/6
  expect_true(all(abs(v("lk", 1) - 1) < 1e-12) && all(abs(v("lk", 5) - 1) < 1e-12))
  expect_equal(c(var(v("lk", 2)), var(v("lk", 6))), c(1/6, 1/6), tolerance = 0.05)
})

# (C229) chisq, weibull, laplace, logis
test_that("chisq, weibull, laplace and logis: log densities, gradients and simulation", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "distributions_prior3"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "data { a = 2.3; b = 1.1; c = -0.4; e = 0.9 }", "model {", "  k ~ gamma(3, 1)", "  m ~ norm(0, 1)", "  s ~ gamma(2, 2)",
               "  a ~ chisq(k)", "  b ~ weibull(k, s)", "  c ~ laplace(m, s)", "  e ~ logis(m, s)", "}"), f)
  recipe = ilike::compile(f)
  P = list(k = 2.6, m = 0.3, s = 0.8)
  truth = dgamma(P$k, 3, 1, log = TRUE) + dnorm(P$m, log = TRUE) + dgamma(P$s, 2, 2, log = TRUE) +
    dchisq(2.3, P$k, log = TRUE) + dweibull(1.1, P$k, P$s, log = TRUE) + (-log(2*P$s) - abs(-0.4 - P$m)/P$s) +
    dlogis(0.9, P$m, P$s, log = TRUE)
  expect_equal(ilike::evaluate_log(recipe, parameters = P), truth, tolerance = 1e-10)
  for (name in names(P))
  {
    h = 1e-6
    fd = (ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] + h)) -
          ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] - h)))/(2*h)
    expect_equal(ilike:::evaluate_gradient_log(recipe, P, name)[1], fd, tolerance = 1e-5, label = name)
  }
  # their gradients with respect to the target, with the target as a parameter
  writeLines(c("ilike 1.0", "model {", "  a ~ chisq(3)", "  b ~ weibull(1.5, 2)", "  c ~ laplace(0.5, 1.2)", "  e ~ logis(-1, 0.7)", "}"), f)
  recipe = ilike::compile(f)
  Q = list(a = 1.7, b = 0.9, c = -0.2, e = 0.4)
  for (name in names(Q))
  {
    h = 1e-6
    fd = (ilike::evaluate_log(recipe, parameters = replace(Q, name, Q[[name]] + h)) -
          ilike::evaluate_log(recipe, parameters = replace(Q, name, Q[[name]] - h)))/(2*h)
    expect_equal(ilike:::evaluate_gradient_log(recipe, Q, name)[1], fd, tolerance = 1e-5, label = name)
  }
  # simulation: quartiles against R's quantile functions (laplace's by its inverse cdf)
  unlink("distributions_prior3", recursive = TRUE)
  ilike::IS(recipe, results_name = "distributions_prior3", number_of_importance_points = 20000, seed = 1)
  o = ilike::load_smc_output("distributions_prior3")
  p = c(0.25, 0.5, 0.75)
  qlaplace = function(p, m, b) ifelse(p < 0.5, m + b*log(2*p), m - b*log(2*(1 - p)))
  quantiles = list(a = qchisq(p, 3), b = qweibull(p, 1.5, 2), c = qlaplace(p, 0.5, 1.2), e = qlogis(p, -1, 0.7))
  for (name in names(quantiles))
  {
    x = o$Value[o$ParameterName == name]
    target = quantiles[[name]]
    expect_lt(max(abs(quantile(x, p) - target))/(target[3] - target[1]), 0.03, label = name)
  }
})

test_that("a division of whole numbers in a statement is not integer division", {
  # (C231) 1/2 in a statement was C++'s integer division (0)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = 0.3 }", "model {", "  mu ~ norm(0, 1)", "  y ~ norm(mu, 1/2 + 3/4)", "}"), f)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = list(mu = 0.1)),
               dnorm(0.1, log = TRUE) + dnorm(0.3, 0.1, 1.25, log = TRUE), tolerance = 1e-12)
})

test_that("a seed written in scientific notation by as.character is accepted", {
  # (C247) a clock-like seed ending in zeros: as.character(1234567000000) is
  # "1.234567e+12", whose first nine characters set.seed used to reject
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  theta ~ norm(0, 1)", "}"), f)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = 0.5), seed = 1234567000000),
               dnorm(0.5, log = TRUE))
})
