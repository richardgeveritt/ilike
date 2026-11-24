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
        "gibbs(mu): mu is not discrete (bern, binom, or values = [...]), and its full conditional is not normal: mu enters the mean nonlinearly: mu is inside exp(") # (ASCII only: the message continues with an ellipsis)
  check(c("  mu ~ norm(0, 1)", "  s ~ lnorm(0, 1)", "  for (i in 1:2) y[i] ~ norm(mu, s)"), "move m = gibbs(s)",
        "s's distribution is lnorm")
  check(c("  mu ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ norm(mu, abs(mu) + 1)"), "move m = gibbs(mu)", "mu is in the standard deviation")
  check(c("  mu ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ pois(mu * mu)"), "move m = gibbs(mu)", "this term's distribution is pois")
  # (gibbs(mu) on a vector is the block update, C216: accepted)
  writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", "  for (k in 1:2) mu[k] ~ norm(0, 1)",
               "  for (i in 1:2) y[i] ~ norm(mu[i], 1)", "}", "algorithm mcmc {", "  move m = gibbs(mu)",
               "  termination = iterations(10)", "}"), f)
  expect_silent(ilike::compile(f))
  # affine through a node and with a coefficient that is another parameter: accepted
  writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", "  mu ~ norm(0, 1)", "  b ~ norm(0, 1)", "  m2 <- 2 * mu + 1",
               "  for (i in 1:2) y[i] ~ norm(b * m2 / 3 - mu, 1)", "}", "algorithm mcmc {", "  move m = gibbs(mu)",
               "  termination = iterations(10)", "}"), f)
  expect_silent(ilike::compile(f))
})

# Conjugate gamma, inverse gamma and beta full conditionals (C215)
gibbs_conjugate_model = c("ilike 1.0",
  "data { y = [0.3, 1.9, -0.4, 1.2, 0.8]; z = [3, 5, 4]; cnt = [2, 7, 1, 4]; g = [1, 2, 1, 2]; e = [0.5, 2, 1.5, 1] }",
  "model {",
  "  tau ~ gamma(2, 1)",
  "  v ~ invgamma(3, 2)",
  "  p ~ beta(2, 3)",
  "  for (j in 1:2) lambda[j] ~ gamma(2, 1)",
  "  sd_y <- 1 / sqrt(tau)",
  "  for (i in 1:5) y[i] ~ norm(0.5, sd_y)",
  "  for (i in 1:5) y[i] ~ norm(0.5, 2 * sqrt(v))",
  "  for (i in 1:3) z[i] ~ binom(10, p)",
  "  for (i in 1:4) cnt[i] ~ pois(e[i] * lambda[g[i]])",
  "}")

test_that("gibbs draws gamma, inverse gamma and beta full conditionals exactly", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gibbs_conjugate"), recursive = TRUE), add = TRUE)
  writeLines(c(gibbs_conjugate_model, "algorithm mcmc {", "  move a = gibbs(tau)", "  move b = gibbs(v)", "  move c = gibbs(p)",
               "  for (j in 1:2) move d = gibbs(lambda[j])", "  termination = iterations(4000)", "}"), f)
  unlink("gibbs_conjugate", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "gibbs_conjugate", seed = 1,
              initial_values = list(list(tau = 1, v = 1, p = 0.5, lambda = c(1, 1))))
  o = ilike::load_mcmc_output("gibbs_conjugate")
  draws = function(name, k = 1) o$Value[o$ParameterName == name & o$Dimension == k]
  y = c(0.3, 1.9, -0.4, 1.2, 0.8); ss = sum((y - 0.5)^2)
  cnt = c(2, 7, 1, 4); g = c(1, 2, 1, 2); e = c(0.5, 2, 1.5, 1)
  # exact conditionals (each variable's terms involve no other parameter)
  exact = list(tau = c(2 + 5/2, 1 + ss/2),                  # gamma(shape, rate)
               v = c(3 + 5/2, 2 + ss/8),                    # invgamma(shape, scale): sd 2 sqrt(v)
               p = c(2 + 12, 3 + 30 - 12),                  # beta
               lambda1 = c(2 + sum(cnt[g == 1]), 1 + sum(e[g == 1])),
               lambda2 = c(2 + sum(cnt[g == 2]), 1 + sum(e[g == 2])))
  check = function(x, mean, sd, label)
  {
    expect_lt(abs(mean(x) - mean)/sd, 0.1, label = label)
    expect_equal(sd(x), sd, tolerance = 0.06, label = label)
    # an exact draw: every proposal accepted
    expect_true(all(diff(x) != 0), label = label)
  }
  a = exact$tau; check(draws("tau"), a[1]/a[2], sqrt(a[1])/a[2], "tau")
  a = exact$v; check(draws("v"), a[2]/(a[1] - 1), a[2]/((a[1] - 1)*sqrt(a[1] - 2)), "v")
  a = exact$p; check(draws("p"), a[1]/sum(a), sqrt(a[1]*a[2]/(sum(a)^2*(sum(a) + 1))), "p")
  a = exact$lambda1; check(draws("lambda", 1), a[1]/a[2], sqrt(a[1])/a[2], "lambda[1]")
  a = exact$lambda2; check(draws("lambda", 2), a[1]/a[2], sqrt(a[1])/a[2], "lambda[2]")
})

test_that("gibbs is refused where the full conditional is not conjugate", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(lines, move, message)
  {
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", lines, "}", "algorithm mcmc {", paste0("  ", move),
                 "  termination = iterations(10)", "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check(c("  tau ~ gamma(2, 1)", "  for (i in 1:2) y[i] ~ norm(0, tau)"), "move m = gibbs(tau)",
        "its full conditional is not gamma: the standard deviation must be c/sqrt(tau) (tau a precision)")
  check(c("  p ~ beta(1, 1)", "  for (i in 1:2) y[i] ~ binom(3, inv_logit(p))"), "move m = gibbs(p)",
        "the probability must be p itself")
  check(c("  l ~ gamma(1, 1)", "  for (i in 1:2) y[i] ~ pois(l + 1)"), "move m = gibbs(l)",
        "the rate must be c*l: l is in a sum or difference")
  check(c("  v ~ invgamma(1, 1)", "  for (i in 1:2) y[i] ~ t(3, 0, sqrt(v))"), "move m = gibbs(v)",
        "this term's distribution is t")
  check(c("  v ~ invgamma(1, 1)", "  for (i in 1:2) y[i] ~ norm(v, 1)"), "move m = gibbs(v)", "v is in the mean")
})

test_that("SMC with a conjugate gibbs move gives the evidence", {
  skip_on_cran()
  # gamma-Poisson: the marginal of the counts is negative binomial
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gibbs_conjugate_smc"), recursive = TRUE), add = TRUE)
  cnt = c(2, 7, 1, 4); e = c(0.5, 2, 1.5, 1)
  writeLines(c("ilike 1.0", "data { cnt = [2, 7, 1, 4]; e = [0.5, 2, 1.5, 1] }", "model {", "  l ~ gamma(2, 1)",
               "  for (i in 1:4) cnt[i] ~ pois(e[i] * l)", "}", "algorithm smc {", "  move m = gibbs(l)",
               "  termination = iterations(5)", "  sequence = annealing(power, [0, 0.05, 0.2, 0.5, 0.8, 1])", "}"), f)
  unlink("gibbs_conjugate_smc", recursive = TRUE)
  z = ilike::SMC_with_MCMC(ilike::compile(f), results_name = "gibbs_conjugate_smc", number_of_particles = 2000, seed = 1)
  exact = log(integrate(function(l) sapply(l, function(v) dgamma(v, 2, 1)*prod(dpois(cnt, e*v))), 0, 60)$value)
  expect_lt(abs(z - exact), 0.05)
})

# gibbs on a whole vector: the block's normal full conditional (C216)
gibbs_regression = function()
{
  set.seed(2)
  n = 8
  X = cbind(1, round(rnorm(n), 2), round(runif(n), 2))
  y = round(drop(X %*% c(0.5, -1, 2) + rnorm(n, 0, 0.7)), 2)
  rows = paste0("[", apply(X, 1, function(r) paste(r, collapse = ", ")), "]", collapse = ", ")
  list(X = X, y = y, sigma = 0.7,
       data = paste0("data { y = [", paste(y, collapse = ", "), "]; X = [", rows, "] }"),
       likelihood = "  for (i in 1:8) y[i] ~ norm(X[i, 1] * beta[1] + X[i, 2] * beta[2] + X[i, 3] * beta[3], 0.7)")
}

test_that("gibbs on a vector draws regression coefficients from their normal full conditional", {
  r = gibbs_regression()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gibbs_block"), recursive = TRUE), add = TRUE)
  run = function(prior, iterations = 3000)
  {
    writeLines(c("ilike 1.0", r$data, "model {", prior, r$likelihood, "}", "algorithm mcmc {", "  move m = gibbs(beta)",
                 paste0("  termination = iterations(", iterations, ")"), "}"), f)
    unlink("gibbs_block", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "gibbs_block", seed = 1, initial_values = list(list(beta = c(0, 0, 0))))
    o = ilike::load_mcmc_output("gibbs_block")
    sapply(1:3, function(k) o$Value[o$Dimension == k])
  }
  # independent normal elements (sd 3), and a multivariate normal prior
  for (case in list(list(prior = "  for (k in 1:3) beta[k] ~ norm(0, 3)", S0 = diag(9, 3)),
                    list(prior = "  beta ~ mvnorm([0, 0, 0], [[4, 1, 0], [1, 4, 0], [0, 0, 4]])",
                         S0 = matrix(c(4, 1, 0, 1, 4, 0, 0, 0, 4), 3))))
  {
    beta = run(case$prior)
    Q = solve(case$S0) + crossprod(r$X)/r$sigma^2
    V = solve(Q); m = drop(V %*% crossprod(r$X, r$y))/r$sigma^2
    expect_lt(max(abs(colMeans(beta) - m)/sqrt(diag(V))), 0.1, label = case$prior)
    expect_equal(apply(beta, 2, sd), sqrt(diag(V)), tolerance = 0.06, label = case$prior)
    expect_lt(max(abs(cor(beta) - cov2cor(V))), 0.05, label = case$prior)
    # exact draws: every proposal accepted, consecutive draws independent
    expect_true(all(diff(beta[, 1]) != 0), label = case$prior)
    expect_lt(abs(cor(beta[-1, 1], beta[-nrow(beta), 1])), 0.06, label = case$prior)
  }
})

test_that("gibbs on a vector is refused where the block's conditional is not normal", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(lines, message)
  {
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", lines, "}", "algorithm mcmc {", "  move m = gibbs(beta)",
                 "  termination = iterations(10)", "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check(c("  for (k in 1:2) beta[k] ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ norm(beta[1] * beta[2], 1)"),
        "the full conditional of the vector beta is not normal: beta enters the mean nonlinearly: two factors of a product depend on beta")
  check(c("  for (k in 1:2) beta[k] ~ bern(0.5)", "  for (i in 1:2) y[i] ~ norm(beta[i], 1)"),
        "for the elements of a discrete beta, write for (j in 1:J) move m = gibbs(beta[j])")
  # its proposal is not one to write
  writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", "  for (k in 1:2) beta[k] ~ norm(0, 1)",
               "  for (i in 1:2) y[i] ~ norm(beta[i], 1)", "}", "algorithm mcmc {",
               "  move m = metropolis_hastings(normal_block(beta))", "  termination = iterations(10)", "}"), f)
  expect_error(ilike::compile(f), "normal_block is not a proposal to write: write move m = gibbs(", fixed = TRUE)
})

test_that("SMC with gibbs on a vector gives the evidence of a linear regression", {
  skip_on_cran()
  r = gibbs_regression()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gibbs_block_smc"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", r$data, "model {", "  for (k in 1:3) beta[k] ~ norm(0, 3)", r$likelihood, "}", "algorithm smc {",
               "  move m = gibbs(beta)", "  termination = iterations(5)",
               "  sequence = annealing(power, [0, 0.02, 0.1, 0.3, 0.6, 1])", "}"), f)
  unlink("gibbs_block_smc", recursive = TRUE)
  z = ilike::SMC_with_MCMC(ilike::compile(f), results_name = "gibbs_block_smc", number_of_particles = 1000, seed = 1)
  S = r$sigma^2*diag(8) + 9*tcrossprod(r$X)
  exact = -0.5*8*log(2*pi) - 0.5*determinant(S)$modulus[1] - 0.5*sum(r$y*solve(S, r$y))
  expect_lt(abs(z - exact), 0.1)
})

# (C228) Dirichlet-categorical/multinomial and normal-inverse-gamma updates
test_that("gibbs draws a Dirichlet full conditional, and categorical elements with Dirichlet probabilities", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gibbs_dirichlet"), recursive = TRUE), add = TRUE)
  run = function(lines, moves, iterations = 4000, initial = list(list(p = c(1, 1, 1)/3)))
  {
    writeLines(c("ilike 1.0", "data { c = [1, 2, 2, 3, 2, 1]; m = [2, 5, 3] }", "model {", lines, "}", "algorithm mcmc {",
                 paste0("  ", moves), paste0("  termination = iterations(", iterations, ")"), "}"), f)
    unlink("gibbs_dirichlet", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "gibbs_dirichlet", seed = 1, initial_values = initial)
    ilike::load_mcmc_output("gibbs_dirichlet")
  }
  o = run(c("  p ~ dirichlet([1, 2, 1.5])", "  for (i in 1:6) c[i] ~ categorical(p)", "  m ~ multinom(10, p)"), "move g = gibbs(p)")
  a = c(1, 2, 1.5) + c(2, 3, 1) + c(2, 5, 3) # the posterior Dirichlet(a)
  for (k in 1:3)
  {
    x = o$Value[o$ParameterName == "p" & o$Dimension == k]
    expect_lt(abs(mean(x) - a[k]/sum(a))/sqrt(a[k]*(sum(a) - a[k])/(sum(a)^2*(sum(a) + 1))), 0.1, label = k)
    expect_equal(sd(x), sqrt(a[k]*(sum(a) - a[k])/(sum(a)^2*(sum(a) + 1))), tolerance = 0.06, label = k)
    expect_true(all(diff(x) != 0))
  }
  # a categorical parameter with Dirichlet probabilities, both by gibbs: z's
  # marginal is E[p] = alpha/sum(alpha)
  o = run(c("  p ~ dirichlet([1, 2, 1.5])", "  z ~ categorical(p)"), c("move g = gibbs(p)", "  move h = gibbs(z)"),
          iterations = 8000, initial = list(list(p = c(1, 1, 1)/3, z = 1)))
  z = o$Value[o$ParameterName == "z"]
  expect_equal(as.vector(table(factor(z, levels = 1:3)))/length(z), c(1, 2, 1.5)/4.5, tolerance = 0.06)
})

test_that("gibbs draws a mean and a variance from their normal-inverse-gamma full conditional", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gibbs_nig"), recursive = TRUE), add = TRUE)
  y = c(1.3, 0.2, 2.1, 1.7, 0.9, 1.4); n = length(y)
  writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  v ~ invgamma(3, 2)",
               "  mu ~ norm(1, 2 * sqrt(v))", "  for (i in 1:6) y[i] ~ norm(mu, sqrt(v))", "}", "algorithm mcmc {",
               "  move g = gibbs([mu, v])", "  termination = iterations(4000)", "}"), f)
  unlink("gibbs_nig", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "gibbs_nig", seed = 1, initial_values = list(list(mu = 0, v = 1)))
  o = ilike::load_mcmc_output("gibbs_nig")
  mu = o$Value[o$ParameterName == "mu"]; v = o$Value[o$ParameterName == "v"]
  # the exact posterior: kappa0 = 1/4 (the prior variance of mu is 4 v)
  k0 = 0.25; kn = k0 + n; mn = (k0*1 + n*mean(y))/kn; an = 3 + n/2
  bn = 2 + 0.5*sum((y - mean(y))^2) + k0*n*(mean(y) - 1)^2/(2*kn)
  expect_lt(abs(mean(v) - bn/(an - 1))/(bn/((an - 1)*sqrt(an - 2))), 0.1)
  expect_equal(sd(v), bn/((an - 1)*sqrt(an - 2)), tolerance = 0.08)
  expect_lt(abs(mean(mu) - mn)/sqrt(bn/((an - 1)*kn)), 0.1)
  expect_equal(sd(mu), sqrt(bn/((an - 1)*kn)), tolerance = 0.06)
  expect_true(all(diff(v) != 0))
})

test_that("gibbs is refused where the Dirichlet or normal-inverse-gamma conditional does not hold", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(lines, move, message)
  {
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", lines, "}", "algorithm mcmc {", paste0("  ", move),
                 "  termination = iterations(10)", "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check(c("  p ~ dirichlet([1, 1])", "  for (i in 1:2) y[i] ~ pois(p[1])"), "move m = gibbs(p)",
        "the full conditional of p is not Dirichlet: this term's distribution is pois")
  check(c("  v ~ invgamma(2, 1)", "  mu ~ norm(0, 1)", "  for (i in 1:2) y[i] ~ norm(mu, sqrt(v))"), "move m = gibbs([mu, v])",
        "the standard deviation must be c*sqrt(v) (v a variance)")
  check(c("  v ~ gamma(2, 1)", "  mu ~ norm(0, sqrt(v))", "  for (i in 1:2) y[i] ~ norm(mu, sqrt(v))"), "move m = gibbs([mu, v])",
        "one of the two needs an invgamma distribution")
})
