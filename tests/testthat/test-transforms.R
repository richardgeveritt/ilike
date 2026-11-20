# Built-in constraining transforms (dev/plan.md 3.18 H2, C209): moves on a
# transformed scale (transform = auto, a list, or none; auto by default for
# hmc), the SMC adaptor on that scale, and IENKI with transform = … .

transform_scale_model = function(move, algorithm = "mcmc", extra = "  termination = iterations(20000)")
{
  y = c(0.3, -1.1, 0.8, 2.2, -0.4)
  c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  sigma ~ gamma(2, 1)",
    "  for (i in 1:5) y[i] ~ norm(0, sigma)", "}", paste0("algorithm ", algorithm, " {"), paste0("  move m = ", move), extra, "}")
}

transform_scale_posterior = function()
{
  y = c(0.3, -1.1, 0.8, 2.2, -0.4)
  density = function(s) sapply(s, function(v) dgamma(v, 2, 1)*prod(dnorm(y, 0, v)))
  grid = seq(1e-4, 20, length = 2e5)
  w = density(grid)
  m = sum(grid*w)/sum(w)
  list(mean = m, sd = sqrt(sum(grid^2*w)/sum(w) - m^2), log_z = log(integrate(density, 0, 50)$value))
}

transform_run = function(lines, name, start = list(sigma = 1))
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  writeLines(lines, f)
  unlink(name, recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(start))
  ilike::load_mcmc_output(name)$Value
}

test_that("moves on the log scale of a positive parameter target its posterior", {
  skip_on_cran()
  exact = transform_scale_posterior()
  for (move in c("metropolis(norm_rw(sigma, 0.5), transform = auto)",
                 "metropolis_hastings(langevin(sigma, [[1]], 0.3), transform = [lower(sigma, 0)])",
                 "metropolis_hastings(hmc(sigma, step = 0.3, steps = 5))"))
  {
    v = transform_run(transform_scale_model(move), "transform_scale")[-(1:2000)]
    expect_equal(mean(v), exact$mean, tolerance = 0.03, label = move)
    expect_equal(sd(v), exact$sd, tolerance = 0.06, label = move)
    expect_true(all(v > 0), label = move)
  }
  # transform = [log(sigma)] and transform = auto give the same chain
  expect_identical(transform_run(transform_scale_model("metropolis(norm_rw(sigma, 0.5), transform = auto)", extra = "  termination = iterations(300)"), "transform_auto"),
                   transform_run(transform_scale_model("metropolis(norm_rw(sigma, 0.5), transform = [log(sigma)])", extra = "  termination = iterations(300)"), "transform_list"))
})

test_that("bounded parameters move on the logit scale", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "transform_bounded"), recursive = TRUE), add = TRUE)
  # beta-binomial (exact Beta(17, 17)); a uniform(0, 10) prior on a scale
  writeLines(c("ilike 1.0", "data { k = [3, 7, 5]; n = 10; z = [1.5, -0.7, 2.4, 0.2] }", "model {", "  p ~ beta(2, 2)",
               "  for (i in 1:3) k[i] ~ binom(n, p)", "  tau ~ unif(0, 10)", "  for (j in 1:4) z[j] ~ norm(0, tau)", "}",
               "algorithm mcmc {", "  move a = metropolis_hastings(hmc(p, step = 0.5, steps = 4))",
               "  move b = metropolis(norm_rw(tau, 0.6), transform = auto)", "  termination = iterations(15000)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "transform_bounded", seed = 1, initial_values = list(list(p = 0.5, tau = 1)))
  o = ilike::load_mcmc_output("transform_bounded")
  p = o$Value[o$ParameterName == "p" & o$Iteration > 1500]
  tau = o$Value[o$ParameterName == "tau" & o$Iteration > 1500]
  expect_equal(mean(p), 0.5, tolerance = 0.05)
  expect_equal(sd(p), sqrt(17*17/(34^2*35)), tolerance = 0.1)
  density = function(t) sapply(t, function(v) prod(dnorm(c(1.5, -0.7, 2.4, 0.2), 0, v)))
  grid = seq(1e-4, 10, length = 1e5)
  w = density(grid)
  expect_equal(mean(tau), sum(grid*w)/sum(w), tolerance = 0.05)
  expect_true(all(tau > 0 & tau < 10))
})

test_that("SMC with moves and the adaptor on the transformed scale gives the evidence", {
  skip_on_cran()
  exact = transform_scale_posterior()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "transform_smc"), recursive = TRUE), add = TRUE)
  writeLines(transform_scale_model("metropolis(norm_rw(sigma, 1), transform = auto, smc_adapt = gaussian(2.38))", "smc",
                                   c("  termination = iterations(5)", "  sequence = annealing(power, [0, 0.01, 0.1, 0.4, 0.7, 1])")), f)
  z = sapply(1:3, function(s) { unlink("transform_smc", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "transform_smc", number_of_particles = 2000, seed = s) })
  expect_lt(abs(mean(z) - exact$log_z), 0.05)
})

test_that("IENKI moves the ensemble on the transformed scale", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "transform_enki"), recursive = TRUE), add = TRUE)
  run = function(transform)
  {
    writeLines(c("ilike 1.0", "data { y = 0.2 }", "model {", "  sigma ~ lnorm(0, 0.8)",
                 "  factor lik = linear_gaussian_data_model(y, sigma, 1, 0.09)", "}", "algorithm enki {", transform,
                 "  sequence = annealing(power, [0, 0.1, 0.3, 0.6, 1])", "}"), f)
    unlink("transform_enki", recursive = TRUE)
    ilike::IENKI(ilike::compile(f), results_name = "transform_enki", number_of_ensemble_members = 2000, seed = 1)
    o = ilike::load_enk_output("transform_enki")
    o$Value[o$Target == max(o$Target)]
  }
  # without a transform members go below 0; on the log scale they stay
  # positive, and the mean is closer to the exact one (0.43; EnKI is
  # approximate for this nonlinear problem)
  expect_lt(min(run(character(0))), 0)
  logged = run("  transform = auto")
  expect_true(all(logged > 0))
  expect_lt(abs(mean(logged) - 0.43), 0.08)
})

test_that("transforms are refused where they cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "transform_refused"), recursive = TRUE), add = TRUE)
  check = function(lines, message)
  {
    writeLines(lines, f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check(transform_scale_model("metropolis(norm_rw(sigma, 0.5), transform = [log(tau)])"), "transform: tau is not moved by this move")
  check(transform_scale_model("metropolis(norm_rw(sigma, 0.5), transform = [logit(sigma, 2, 1)])"), "logit(x, a, b): a must be below b")
  check(transform_scale_model("metropolis(norm_rw(sigma, 0.5), transform = sqrt)"), "transform is auto, none, or a list of")
  check(transform_scale_model("independent_mh(gamma(sigma, 2, 1), transform = auto)"), "an independent proposal does not move on a transformed scale")
  writeLines(c(transform_scale_model("metropolis(norm_rw(sigma, 0.5))", extra = c("  transform = auto", "  termination = iterations(10)"))), f)
  expect_error(ilike::MCMC(ilike::compile(f), results_name = "transform_refused", seed = 1, initial_values = list(list(sigma = 1))),
               "transform = ... in an algorithm section is for IENKI", fixed = TRUE)
})

# (C232) the simplex, cov and corr transforms, and the gradients of the
# matrix distributions
test_that("the gradients of wishart, invwishart and lkj targets and of mvnorm and mvt matrices", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.4, -0.7]; t = [1.1, 0.2] }", "model {", "  W ~ wishart(5, [[2, 0.5], [0.5, 1]])",
               "  V ~ invwishart(6, [[2, 0.5], [0.5, 1]])", "  R ~ lkj(2.5, 2)", "  y ~ mvnorm([0, 0], V)", "  t ~ mvt(4, [0, 0], W)", "}"), f)
  recipe = ilike::compile(f)
  P = list(W = matrix(c(1.5, 0.3, 0.3, 0.8), 2), V = matrix(c(1.2, -0.2, -0.2, 0.9), 2), R = matrix(c(1, 0.4, 0.4, 1), 2))
  # diagonal elements by finite differences (an off-diagonal change alone leaves the matrices non-symmetric)
  for (name in c("W", "V"))
  {
    g = ilike:::evaluate_gradient_log(recipe, P, name)
    for (k in c(1, 4))
    {
      h = 1e-6; Pp = P; Pm = P
      Pp[[name]][k] = Pp[[name]][k] + h; Pm[[name]][k] = Pm[[name]][k] - h
      fd = (ilike::evaluate_log(recipe, parameters = Pp) - ilike::evaluate_log(recipe, parameters = Pm))/(2*h)
      expect_equal(g[k], fd, tolerance = 1e-5, label = paste(name, k))
    }
  }
  # and the formula: lkj's (eta - 1) R^-1
  expect_equal(ilike:::evaluate_gradient_log(recipe, P, "R"), 1.5*solve(P$R), tolerance = 1e-10)
})

test_that("nuts moves a probability vector on the simplex, and a correlation matrix", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "transforms_whole"), recursive = TRUE), add = TRUE)
  run = function(lines, move, initial, iterations)
  {
    writeLines(c("ilike 1.0", lines, "algorithm mcmc {", paste0("  move m = ", move), paste0("  termination = iterations(", iterations, ")"), "}"), f)
    unlink("transforms_whole", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "transforms_whole", seed = 1, initial_values = list(initial))
    ilike::load_mcmc_output("transforms_whole")
  }
  # the Dirichlet posterior: dirichlet([2, 3, 1.5]) with counts (2, 3, 1)
  o = run(c("data { c = [1, 2, 2, 3, 2, 1] }", "model {", "  p ~ dirichlet([2, 3, 1.5])", "  for (i in 1:6) c[i] ~ categorical(p)", "}"),
          "nuts(p, step = 0.3, adapt = dual_averaging(until = iterations(500)))", list(p = c(1, 1, 1)/3), 3000)
  a = c(4, 6, 2.5)
  x = sapply(1:3, function(k) o$Value[o$ParameterName == "p" & o$Dimension == k & o$Iteration > 500])
  expect_true(all(abs(rowSums(x) - 1) < 1e-12) && all(x > 0))
  sds = sqrt(a*(sum(a) - a)/(sum(a)^2*(sum(a) + 1)))
  expect_lt(max(abs(colMeans(x) - a/sum(a))/sds), 0.15)
  expect_equal(apply(x, 2, sd), sds, tolerance = 0.1)
  # an explicit transform list
  expect_silent(ilike::compile(f))
  skip_on_cran()
  # the LKJ prior: each correlation 2 Beta(b, b) - 1, b = eta + (d - 2)/2, variance 1/6 for eta = 2, d = 3
  o = run(c("model {", "  R ~ lkj(2, 3)", "}"), "nuts(R, step = 0.3, transform = [corr(R)], adapt = dual_averaging(until = iterations(1000)))",
          list(R = diag(3)), 8000)
  r = sapply(c(2, 3, 6), function(k) o$Value[o$ParameterName == "R" & o$Dimension == k & o$Iteration > 1000])
  expect_true(all(abs(o$Value[o$ParameterName == "R" & o$Dimension == 1] - 1) < 1e-10))
  expect_equal(apply(r, 2, var), rep(1/6, 3), tolerance = 0.08)
})

test_that("nuts moves a covariance matrix: the inverse Wishart posterior of normal data", {
  skip_on_cran()
  set.seed(2); Y = matrix(round(rnorm(16), 2), 8); S0 = matrix(c(2, 0.5, 0.5, 1), 2)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "transforms_cov"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", paste0("y", 1:8, " = [", Y[, 1], ", ", Y[, 2], "]", collapse = "; "), " }"), "model {",
               "  Sigma ~ invwishart(6, [[2, 0.5], [0.5, 1]])", paste0("  y", 1:8, " ~ mvnorm([0, 0], Sigma)"), "}",
               "algorithm mcmc {", "  move m = nuts(Sigma, step = 0.3, adapt = dual_averaging(until = iterations(1000)))",
               "  termination = iterations(11000)", "}"), f)
  unlink("transforms_cov", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "transforms_cov", seed = 1, initial_values = list(list(Sigma = diag(2))))
  o = ilike::load_mcmc_output("transforms_cov")
  # IW(14, S0 + Y'Y): mean Psi/(nu - p - 1), var of a diagonal element 2 Psi_kk^2/((nu - p - 1)^2 (nu - p - 3))
  Psi = S0 + crossprod(Y); nu = 14
  for (k in c(1, 4))
  {
    x = o$Value[o$ParameterName == "Sigma" & o$Dimension == k & o$Iteration > 1000]
    sd_k = sqrt(2*Psi[k]^2/((nu - 3)^2*(nu - 5)))
    expect_lt(abs(mean(x) - Psi[k]/(nu - 3))/sd_k, 0.1, label = k)
  }
})
