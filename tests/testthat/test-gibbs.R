# Gibbs updates from normal full conditionals (dev/plan.md C206): gibbs(mu)
# on a scalar and gibbs(theta[j]) on elements, when the model statements make
# the full conditional normal; and the evidence of SMC with gibbs moves
# (discrete and normal), which was wrong before C206.

gibbs_file = function(lines)
{
  f = tempfile(fileext = ".ilike")
  writeLines(lines, f)
  f
}

test_that("gibbs on a normal mean draws from its full conditional", {
  skip_on_cran()
  y = c(1.2, 0.4, 2.1, 1.7, 0.9)
  f = gibbs_file(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  mu ~ norm(0, 10)",
                   "  for (i in 1:5) y[i] ~ norm(mu, 1)", "}", "algorithm mcmc {", "  move m = gibbs(mu)",
                   "  termination = iterations(10000)", "}"))
  on.exit(unlink(c(f, "gibbs_mean"), recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "gibbs_mean", seed = 1, initial_values = list(list(mu = 0)))
  mu = ilike::load_mcmc_output("gibbs_mean")$Value
  P = 1/100 + 5
  expect_equal(mean(mu[-1]), sum(y)/P, tolerance = 0.02)
  expect_equal(sd(mu[-1]), sqrt(1/P), tolerance = 0.03)
  # an exact draw from the conditional: every proposal accepted
  expect_equal(mean(diff(mu) != 0), 1)
})

test_that("gibbs on elements gives the exact posterior of a hierarchical normal model", {
  skip_on_cran()
  # eight schools with tau fixed: (mu, theta) are jointly normal a posteriori
  ys = c(28, 8, -3, 7, -1, 1, 18, 12); s = c(15, 10, 16, 11, 9, 11, 10, 18); tau = 6
  f = gibbs_file(c("ilike 1.0", paste0("data { y = [", paste(ys, collapse = ", "), "]; s = [", paste(s, collapse = ", "), "] }"),
                   "model {", "  mu ~ norm(0, 20)", paste0("  for (j in 1:8) theta[j] ~ norm(mu, ", tau, ")"),
                   "  for (j in 1:8) y[j] ~ norm(theta[j], s[j])", "}", "algorithm mcmc {",
                   "  for (j in 1:8) move g = gibbs(theta[j])", "  move m = gibbs(mu)", "  termination = iterations(20000)", "}"))
  on.exit(unlink(c(f, "gibbs_schools"), recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "gibbs_schools", seed = 1, initial_values = list(list(mu = 0, theta = rep(0, 8))))
  o = ilike::load_mcmc_output("gibbs_schools")
  keep = o$Iteration > 1000
  theta = sapply(1:8, function(j) o$Value[keep & o$ParameterName == "theta" & o$Dimension == j])
  mu = o$Value[keep & o$ParameterName == "mu"]
  # the exact posterior of (mu, theta): precision Q, mean Q^-1 b
  Q = matrix(0, 9, 9); b = rep(0, 9)
  Q[1, 1] = 1/20^2 + 8/tau^2
  Q[1, 2:9] = Q[2:9, 1] = -1/tau^2
  diag(Q)[2:9] = 1/tau^2 + 1/s^2
  b[2:9] = ys/s^2
  S = solve(Q); m = c(S %*% b)
  expect_lt(max(abs(c(mean(mu), colMeans(theta)) - m)/sqrt(diag(S))), 0.1)
  expect_equal(c(sd(mu), apply(theta, 2, sd)), sqrt(diag(S)), tolerance = 0.05)
})

test_that("SMC with gibbs moves gives the evidence", {
  skip_on_cran()
  # normal mean: exact evidence
  y = c(1.2, 0.4, 2.1, 1.7, 0.9)
  f = gibbs_file(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  mu ~ norm(0, 10)",
                   "  for (i in 1:5) y[i] ~ norm(mu, 1)", "}", "algorithm smc {", "  move m = gibbs(mu)",
                   "  termination = iterations(5)", "  sequence = annealing(power, [0, 0.001, 0.01, 0.1, 0.4, 1])", "}"))
  on.exit(unlink(c(f, "gibbs_smc"), recursive = TRUE), add = TRUE)
  run = function(s) { unlink("gibbs_smc", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "gibbs_smc", number_of_particles = 2000, seed = s) }
  S = diag(5) + 100
  exact = c(-0.5*determinant(2*pi*S)$modulus - 0.5*t(y) %*% solve(S, y))
  expect_lt(abs(mean(sapply(1:3, run)) - exact), 0.1)
  # a discrete gibbs move (enumeration over the states)
  y = c(0.8, 1.3)
  writeLines(c("ilike 1.0", "data { y = [0.8, 1.3] }", "model {", "  s ~ bern(0.3)", "  for (i in 1:2) y[i] ~ norm(s, 1)", "}",
               "algorithm smc {", "  move m = gibbs(s)", "  termination = iterations(4)", "  sequence = annealing(power, [0, 0.1, 0.4, 0.7, 1])", "}"), f)
  exact = log(0.7*prod(dnorm(y, 0, 1)) + 0.3*prod(dnorm(y, 1, 1)))
  expect_lt(abs(mean(sapply(1:3, run)) - exact), 0.05)
})

test_that("gibbs is refused where the full conditional is not normal", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(lines, move, message)
  {
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", lines, "}", "algorithm mcmc {", paste0("  ", move),
                 "  termination = iterations(10)", "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check(c("  mu ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ norm(exp(mu), 1)"), "move m = gibbs(mu)",
        "gibbs(mu): mu is not discrete (bern, binom, or values = [...]), and its full conditional is not normal: mu enters the mean nonlinearly: mu is inside exp(…)")
  check(c("  mu ~ norm(0, 1)", "  s ~ lnorm(0, 1)", "  for (i in 1:2) y[i] ~ norm(mu, s)"), "move m = gibbs(s)",
        "s's distribution is lnorm")
  check(c("  mu ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ norm(mu, abs(mu) + 1)"), "move m = gibbs(mu)", "mu is in the standard deviation")
  check(c("  mu ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ pois(mu * mu)"), "move m = gibbs(mu)", "this term's distribution is pois")
  check(c("  for (k in 1:2) mu[k] ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ norm(mu[i], 1)"), "move m = gibbs(mu)",
        "gibbs(mu) moves a scalar")
  # affine through a node and with a coefficient that is another parameter: accepted
  writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", "  mu ~ norm(0, 1)", "  b ~ norm(0, 1)", "  m2 <- 2 * mu + 1",
               "  for (i in 1:2) y[i] ~ norm(b * m2 / 3 - mu, 1)", "}", "algorithm mcmc {", "  move m = gibbs(mu)",
               "  termination = iterations(10)", "}"), f)
  expect_silent(ilike::compile(f))
})
