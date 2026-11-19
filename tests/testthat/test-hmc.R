# Hamiltonian Monte Carlo as a proposal (dev/plan.md 3.18 stage H1, C208):
# metropolis_hastings(hmc(…)) and unadjusted(hmc(…)), with the gradients of
# model statements.

hmc_smooth = function()
{
  y = c(1.2, 0.4, -0.3, 0.8, 2.1, 1.7); K = 6
  L = matrix(0, K, K)
  for (k in 1:(K - 1)) L[c(k, k + 1), c(k, k + 1)] = L[c(k, k + 1), c(k, k + 1)] + matrix(c(1, -1, -1, 1), 2)
  Q = diag(K)/100 + 2*L + diag(K)
  S = solve(Q)
  list(K = K, S = S, m = solve(Q, y),
       covariance = paste0("[", paste(apply(S, 1, function(r) paste0("[", paste(signif(r, 8), collapse = ", "), "]")), collapse = ", "), "]"),
       lines = readLines(system.file("extdata", "lang_model_smooth.ilike", package = "ilike")))
}

hmc_run = function(model, move, name, iterations = 4000, algorithm = "mcmc", extra = character(0))
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  writeLines(c(model$lines[seq_len(grep("^algorithm", model$lines) - 1)], paste0("algorithm ", algorithm, " {"),
               paste0("  move m = ", move), paste0("  termination = iterations(", iterations, ")"), extra, "}"), f)
  unlink(name, recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(theta = rep(0, model$K))))
  o = ilike::load_mcmc_output(name)
  list(theta = sapply(1:model$K, function(k) o$Value[o$Dimension == k]))
}

expect_posterior = function(model, theta, burn = 1000, mean_tolerance = 0.1, sd_tolerance = 0.05)
{
  theta = theta[-(1:burn), ]
  expect_lt(max(abs(colMeans(theta) - model$m)/sqrt(diag(model$S))), mean_tolerance)
  expect_equal(apply(theta, 2, sd), sqrt(diag(model$S)), tolerance = sd_tolerance)
}

test_that("HMC targets the posterior, with and without a mass matrix", {
  skip_on_cran()
  model = hmc_smooth()
  identity = hmc_run(model, "metropolis_hastings(hmc(theta, step = 0.2, steps = 10))", "hmc_identity")
  expect_posterior(model, identity$theta)
  expect_gt(mean(diff(identity$theta[, 1]) != 0), 0.9)
  mass = hmc_run(model, paste0("metropolis_hastings(hmc(theta, step = 0.5, length = 2.5, mass = ", model$covariance, "))"),
                 "hmc_mass")
  expect_posterior(model, mass$theta)
})

test_that("HMC with one leapfrog step is the Langevin proposal", {
  skip_on_cran()
  # (eps C L^-T z = eps L z with C = L L': the same random numbers)
  model = hmc_smooth()
  one = hmc_run(model, paste0("metropolis_hastings(hmc(theta, step = 0.9, steps = 1, mass = ", model$covariance, "))"),
                "hmc_one", 300)
  langevin = hmc_run(model, paste0("metropolis_hastings(langevin(theta, ", model$covariance, ", 0.81))"), "hmc_langevin", 300)
  expect_equal(one$theta, langevin$theta, tolerance = 1e-12)
})

test_that("unadjusted HMC, and adaptation of the step size to the target 0.8", {
  skip_on_cran()
  model = hmc_smooth()
  unadjusted = hmc_run(model, "unadjusted(hmc(theta, step = 0.05, steps = 10))", "hmc_unadjusted")
  expect_posterior(model, unadjusted$theta)
  adapted = hmc_run(model, "metropolis_hastings(hmc(theta, step = 0.1, steps = 10), adapt = scaled_am())", "hmc_adapted", 6000)
  expect_equal(mean(diff(adapted$theta[-(1:2000), 1]) != 0), 0.8, tolerance = 0.05)
  expect_posterior(model, adapted$theta, burn = 2000)
})

test_that("HMC moves a list of variables as their stacked vector", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "hmc_list"), recursive = TRUE), add = TRUE)
  run = function(named)
  {
    x = if (named) c("a", "v[1]", "v[2]") else paste0("theta[", 1:3, "]")
    writeLines(c("ilike 1.0", "data { y = [0.5, -1.2, 2.0] }", "model {", paste0("  ", x, " ~ norm(0, 1)"),
                 paste0("  y[1] ~ norm(", x[1], " + ", x[2], ", 1)"), paste0("  y[2] ~ norm(", x[3], ", 0.5)"),
                 paste0("  y[3] ~ norm(", x[1], " - ", x[3], ", 1)"), "}", "algorithm mcmc {",
                 paste0("  move m = metropolis_hastings(hmc(", if (named) "[a, v]" else "theta", ", step = 0.3, steps = 5))"),
                 "  termination = iterations(1000)", "}"), f)
    unlink("hmc_list", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "hmc_list", seed = 2,
                initial_values = list(if (named) list(a = 0, v = c(0, 0)) else list(theta = rep(0, 3))))
    matrix(ilike::load_mcmc_output("hmc_list")$Value, ncol = 3, byrow = TRUE)
  }
  expect_identical(run(TRUE), run(FALSE))
})

test_that("SMC with HMC moves gives the evidence", {
  skip_on_cran()
  y = c(1.2, 0.4, 2.1, 1.7, 0.9)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "hmc_smc"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  mu ~ norm(0, 10)",
               "  for (i in 1:5) y[i] ~ norm(mu, 1)", "}", "algorithm smc {", "  move m = metropolis_hastings(hmc(mu, step = 0.3, steps = 5))",
               "  termination = iterations(5)", "  sequence = annealing(power, [0, 0.001, 0.01, 0.1, 0.4, 1])", "}"), f)
  z = sapply(1:3, function(s) { unlink("hmc_smc", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "hmc_smc", number_of_particles = 2000, seed = s) })
  S = diag(5) + 100
  expect_lt(abs(mean(z) - c(-0.5*determinant(2*pi*S)$modulus - 0.5*t(y) %*% solve(S, y))), 0.1)
})

test_that("HMC is refused where it cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "hmc_refused"), recursive = TRUE), add = TRUE)
  run = function(move)
  {
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", "  for (k in 1:2) theta[k] ~ norm(0, 1)",
                 "  for (k in 1:2) y[k] ~ norm(theta[k], 1)", "}", "algorithm mcmc {", paste0("  move m = ", move),
                 "  termination = iterations(10)", "}"), f)
    unlink("hmc_refused", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "hmc_refused", seed = 1, initial_values = list(list(theta = c(0, 0))))
  }
  expect_error(run("metropolis_hastings(hmc(theta, step = 0.1))"), "give the number of leapfrog steps (steps = L) or the trajectory length", fixed = TRUE)
  expect_error(run("metropolis_hastings(hmc(theta, step = 0.1, steps = 2, length = 1))"), "one of them", fixed = TRUE)
  expect_error(run("metropolis_hastings(hmc(theta, step = 0.1, steps = 2, mass = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))"),
               "The proposal's mass (covariance) for theta is 3 x 3, but theta has 2 element(s).", fixed = TRUE)
  expect_error(run("metropolis(hmc(theta, step = 0.1, steps = 2))"), "cannot be used as a Metropolis proposal", fixed = TRUE)
})
