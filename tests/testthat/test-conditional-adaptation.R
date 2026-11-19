# Adaptation of moves on part of the parameters (dev/plan.md, "adaptation of
# element moves, and conditional scaling of partial moves"): the conditional
# covariance of the moved block given the rest, shrunk towards its diagonal
# when singular; element moves with a step per element.

test_that("the covariance numerics agree with direct computations", {
  set.seed(1)
  n = 200; D = 5
  A = matrix(rnorm(D*D), D); X = matrix(rnorm(n*D), n) %*% chol(crossprod(A) + diag(D))
  w = rexp(n); wn = w/sum(w)
  r = ilike:::covariance_tools_check(X, w, 2)
  expect_equal(r$covariance, cov.wt(X, wt = wn, method = "unbiased")$cov, tolerance = 1e-12)
  # Schafer & Strimmer's intensity (target "D") with weights, computed naively
  m = colSums(X*wn); Z = sweep(X, 2, m); cc = 1 - sum(wn^2)
  Sb = t(Z) %*% (Z*wn); S = Sb/cc; num = 0
  for (i in 1:D) for (j in 1:D) if (i != j) num = num + sum(wn^2*(Z[, i]*Z[, j] - Sb[i, j])^2)/cc^2
  expect_equal(r$lambda, min(1, max(0, num/(sum(S^2) - sum(diag(S)^2)))), tolerance = 1e-10)
  expect_equal(r$block_covariance, S[1:2, 1:2] - S[1:2, 3:D] %*% solve(S[3:D, 3:D]) %*% S[3:D, 1:2], tolerance = 1e-10)
  expect_equal(as.numeric(r$variances), 1/diag(solve(S)), tolerance = 1e-10)
  # fewer points than dimensions: shrunk, and a conditional variance below the marginal one
  r2 = ilike:::covariance_tools_check(X[1:4, ], rep(1, 4), 1)
  expect_true(r2$ok && (r2$lambda > 0) && (r2$block_covariance[1, 1] < r2$covariance[1, 1]))
  # a column that does not vary is left out
  r3 = ilike:::covariance_tools_check(cbind(X, 3), w, 2)
  expect_equal(r3$block_covariance, r$block_covariance, tolerance = 1e-10)
  expect_true(is.na(r3$variances[D + 1]))
})

test_that("element moves adapt a step per element to its acceptance rate", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.5, -1.0, 3.0]; s = [0.1, 1, 10] }", "model {",
               "  for (j in 1:3) theta[j] ~ norm(0, 10)", "  for (j in 1:3) y[j] ~ norm(theta[j], s[j])", "}",
               "algorithm mcmc {", "  for (j in 1:3) move m = metropolis(norm_rw(theta[j], 1), adapt = scale(averaged))",
               "  termination = iterations(20000)", "}"), f)
  run = function(name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(theta = c(0, 0, 0))))
    list(output = ilike::load_mcmc_output(name), adaptation = ilike::load_mcmc_adaptation(name))
  }
  adapted = run("element_adapt")
  s = c(0.1, 1, 10); y = c(0.5, -1, 3); v = 1/(1/100 + 1/s^2)
  for (j in 1:3)
  {
    x = adapted$output$Value[adapted$output$Dimension == j & adapted$output$Iteration > 5000]
    expect_equal(mean(diff(x) != 0), 0.44, tolerance = 0.05/0.44)
    expect_equal(sd(x), sqrt(v[j]), tolerance = 0.05)
  }
  # steps of about 2.4 posterior sds, whatever the element's spread
  final = adapted$adaptation$final
  steps = final$Value[match(paste0("theta[", 1:3, "].scale"), final$Parameter)]
  expect_equal(steps/sqrt(v), rep(2.4, 3), tolerance = 0.15)
  # the same chain through the copy-based moves
  old = options(ilike.in_place_moves = FALSE)
  on.exit(options(old), add = TRUE)
  expect_identical(run("element_adapt_copies")$output$Value, adapted$output$Value)
})

test_that("SMC sets element steps and named variables' covariances from conditional moments", {
  skip_on_cran()
  # smoothing model (Gaussian, exact evidence) with element moves set from the particles
  y = c(1.2, 0.4, -0.3, 0.8, 2.1, 1.7); K = 6
  L = matrix(0, K, K)
  for (k in 1:(K - 1)) L[c(k, k + 1), c(k, k + 1)] = L[c(k, k + 1), c(k, k + 1)] + matrix(c(1, -1, -1, 1), 2)
  Q = diag(K)/100 + 2*L + diag(K)
  log_z = -K/2*log(2*pi) - K/2*log(100) - 0.5*as.numeric(determinant(Q)$modulus) - 0.5*(sum(y^2) - sum(y*solve(Q, y)))
  lines = readLines(system.file("extdata", "lang_model_smooth.ilike", package = "ilike"))
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(lines[seq_len(grep("^algorithm", lines) - 1)], "algorithm smc {",
               "  for (k in 1:K) move m = metropolis(norm_rw(theta[k], 5), smc_adapt = gaussian(2.38))",
               "  termination = iterations(5)", "  sequence = annealing(power, [0, 0.01, 0.03, 0.1, 0.3, 1])", "}"), f)
  z = sapply(1:3, function(s) {
    name = paste0("smc_element_adapt_", s)
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = name, number_of_particles = 2000, seed = s) })
  expect_equal(mean(z), log_z, tolerance = 0.02)

  # two strongly correlated named parameters, each moved alone
  writeLines(c("ilike 1.0", "data { y = 1 }", "model {", "  a ~ norm(0, 1)", "  b ~ norm(0, 1)", "  y ~ norm(a + b, 0.1)", "}",
               "algorithm smc {", "  move m1 = metropolis(norm_rw(a, 1), smc_adapt = gaussian(2.38))",
               "  move m2 = metropolis(norm_rw(b, 1), smc_adapt = gaussian(2.38))", "  termination = iterations(5)",
               "  sequence = annealing(power, [0, 0.01, 0.1, 0.3, 1])", "}"), f)
  log_z = dnorm(1, 0, sqrt(2 + 0.01), log = TRUE)
  z = sapply(1:3, function(s) {
    name = paste0("smc_named_adapt_", s)
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = name, number_of_particles = 2000, seed = s) })
  expect_equal(mean(z), log_z, tolerance = 0.02)
})

test_that("element adaptation is refused where there is no step", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  write_model = function(move)
    writeLines(c("ilike 1.0", "model { for (k in 1:3) s[k] ~ bern(0.3) }", "algorithm mcmc {", move, "}"), f)
  write_model("  for (k in 1:3) move m = metropolis(flip(s[k]), adapt = scale())")
  expect_error(ilike::compile(f), "flip moves have no step to adapt")
  write_model("  for (k in 1:3) move m = gibbs(s[k], adapt = scale())")
  expect_error(ilike::compile(f), "gibbs moves have no step to adapt")
  write_model("  for (k in 1:3) move m = metropolis(norm_rw(s[k], 1), adapt = scaled_am())")
  expect_error(ilike::compile(f), "adapt each element's step with adapt = scale")
})

test_that("adaptive Metropolis on one of several parameters uses its conditional covariance", {
  skip_on_cran()
  # a, b ~ N(0, 1), y = 1 ~ N(a + b, 0.1): the conditional variance of a given b
  # is 1/101, its marginal variance about 0.5
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  run = function(moves, name)
  {
    writeLines(c("ilike 1.0", "data { y = 1 }", "model {", "  a ~ norm(0, 1)", "  b ~ norm(0, 1)", "  y ~ norm(a + b, 0.1)", "}",
                 "algorithm mcmc {", moves, "  termination = iterations(20000)", "}"), f)
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(a = 0.5, b = 0.5)))
    o = ilike::load_mcmc_output(name)
    list(a = o$Value[o$ParameterName == "a" & o$Iteration > 5000], b = o$Value[o$ParameterName == "b" & o$Iteration > 5000],
         final = ilike::load_mcmc_adaptation(name)$final)
  }
  # am keeps the scale 2.38^2/d: with the conditional covariance the acceptance
  # rate is that of a well scaled one-dimensional walk (about 0.44)
  am = run(c("  move m1 = metropolis(norm_rw(a, 1), adapt = am())", "  move m2 = metropolis(norm_rw(b, 1), adapt = am())"),
           "conditional_am")
  conditional = am$final$Value[am$final$Parameter == "move1.a.conditional_covariance[1,1]"]
  expect_equal(conditional, 1/101, tolerance = 0.2)
  expect_equal(mean(diff(am$a) != 0), 0.44, tolerance = 0.1/0.44)
  # the mirror proposal, centred on the conditional mean of a given b
  mirror = run(c("  move m1 = metropolis(mirror(a, 0, 1, 1), adapt = scaled_am())",
                 "  move m2 = metropolis(mirror(b, 0, 1, 1), adapt = scaled_am())"), "conditional_mirror")
  S = solve(diag(2) + 100*matrix(1, 2, 2))
  expect_equal(mean(mirror$a), (S %*% c(100, 100))[1], tolerance = 0.1)
  expect_equal(sd(mirror$a), sqrt(S[1, 1]), tolerance = 0.1)
  expect_equal(cor(mirror$a, mirror$b), S[1, 2]/S[1, 1], tolerance = 0.01)
})
