# The No-U-Turn sampler (dev/plan.md 3.18 stage H3, C210): move m = nuts(…).

nuts_smooth = function()
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

nuts_run = function(model, move, name, iterations = 3000)
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  writeLines(c(model$lines[seq_len(grep("^algorithm", model$lines) - 1)], "algorithm mcmc {", paste0("  move m = ", move),
               paste0("  termination = iterations(", iterations, ")"), "}"), f)
  unlink(name, recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(theta = rep(0, model$K))))
  o = ilike::load_mcmc_output(name)
  sapply(1:model$K, function(k) o$Value[o$Dimension == k])
}

test_that("NUTS targets the posterior, with and without a mass matrix", {
  skip_on_cran()
  model = nuts_smooth()
  for (move in c("nuts(theta, step = 0.2)", paste0("nuts(theta, step = 0.8, mass = ", model$covariance, ")"),
                 "nuts(theta, step = 0.05, adapt = scale())"))
  {
    theta = nuts_run(model, move, "nuts_smooth")[-(1:300), ]
    expect_lt(max(abs(colMeans(theta) - model$m)/sqrt(diag(model$S))), 0.1, label = move)
    expect_equal(apply(theta, 2, sd), sqrt(diag(model$S)), tolerance = 0.06, label = move)
  }
})

test_that("NUTS moves a positive parameter on the log scale, and in SMC gives the evidence", {
  skip_on_cran()
  y = c(0.3, -1.1, 0.8, 2.2, -0.4)
  density = function(s) sapply(s, function(v) dgamma(v, 2, 1)*prod(dnorm(y, 0, v)))
  grid = seq(1e-4, 20, length = 2e5)
  w = density(grid)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_scale"), recursive = TRUE), add = TRUE)
  model = function(algorithm, settings)
    c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  sigma ~ gamma(2, 1)",
      "  for (i in 1:5) y[i] ~ norm(0, sigma)", "}", paste0("algorithm ", algorithm, " {"), "  move m = nuts(sigma, step = 0.5)",
      settings, "}")
  writeLines(model("mcmc", "  termination = iterations(6000)"), f)
  ilike::MCMC(ilike::compile(f), results_name = "nuts_scale", seed = 1, initial_values = list(list(sigma = 1)))
  v = ilike::load_mcmc_output("nuts_scale")$Value[-(1:600)]
  expect_equal(mean(v), sum(grid*w)/sum(w), tolerance = 0.03)
  expect_true(all(v > 0))
  writeLines(model("smc", c("  termination = iterations(5)", "  sequence = annealing(power, [0, 0.01, 0.1, 0.4, 0.7, 1])")), f)
  z = sapply(1:2, function(s) { unlink("nuts_scale", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "nuts_scale", number_of_particles = 1000, seed = s) })
  expect_lt(abs(mean(z) - log(integrate(density, 0, 50)$value)), 0.1)
})

test_that("NUTS explores a 50-dimensional normal, and moves a list as its stacked vector", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_normal"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  for (k in 1:50) x[k] ~ norm(0, 1)", "}", "algorithm mcmc {", "  move m = nuts(x, step = 0.3)",
               "  termination = iterations(600)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "nuts_normal", seed = 1, initial_values = list(list(x = rep(0, 50))))
  X = matrix(ilike::load_mcmc_output("nuts_normal")$Value, ncol = 50, byrow = TRUE)[-(1:100), ]
  expect_lt(abs(mean(colMeans(X))), 0.05)
  expect_equal(mean(apply(X, 2, var)), 1, tolerance = 0.05)
  # [a, b] and a vector theta: the same chain
  run = function(named)
  {
    x = if (named) c("a", "b") else c("theta[1]", "theta[2]")
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", paste0("  ", x, " ~ norm(0, 1)"),
                 paste0("  y[", 1:2, "] ~ norm(", x, ", 1)"), "}", "algorithm mcmc {",
                 paste0("  move m = nuts(", if (named) "[a, b]" else "theta", ", step = 0.4)"), "  termination = iterations(200)", "}"), f)
    unlink("nuts_normal", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "nuts_normal", seed = 2,
                initial_values = list(if (named) list(a = 0, b = 0) else list(theta = c(0, 0))))
    matrix(ilike::load_mcmc_output("nuts_normal")$Value, ncol = 2, byrow = TRUE)
  }
  expect_identical(run(TRUE), run(FALSE))
})

test_that("nuts is refused where it cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(move, message)
  {
    writeLines(c("ilike 1.0", "data { y = 1 }", "model {", "  theta ~ norm(0, 1)", "  y ~ norm(theta, 1)", "}", "algorithm mcmc {",
                 paste0("  move m = ", move), "  termination = iterations(10)", "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check("metropolis_hastings(nuts(theta, step = 0.1))", "nuts is a move of its own, not a proposal: move m = nuts(x, step = …).")
  check("metropolis(nuts(theta, step = 0.1))", "nuts is a move of its own, not a proposal")
})
