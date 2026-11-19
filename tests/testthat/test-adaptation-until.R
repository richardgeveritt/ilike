# Stopping adaptation (dev/plan.md 3.18 item 10, C207): adapt = …(until =
# criterion), with iterations(n), stable(…), acceptance(…), all(…), any(…).
# After the criterion holds the proposal's parameters stay fixed.

until_model = c("ilike 1.0", "data { y = 1 }", "model {", "  a ~ norm(0, 1)", "  b ~ norm(0, 1)", "  y ~ norm(a + b, 0.3)", "}")

until_run = function(move, name, iterations = 3000)
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  writeLines(c(until_model, "algorithm mcmc {", paste0("  ", move), paste0("  termination = iterations(", iterations, ")"), "}"), f)
  unlink(name, recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(a = 0, b = 0)))
  list(output = ilike::load_mcmc_output(name), adaptation = ilike::load_mcmc_adaptation(name))
}

stopped_at = function(r) r$adaptation$final$Value[grepl("adaptation_stopped_at$", r$adaptation$final$Parameter)]

scale_trace = function(r, parameter = "a;b.scale")
{
  tr = r$adaptation$trace
  tr$Value[tr$Parameter == parameter]
}

test_that("adaptation stops after a number of iterations, and the chain stays right", {
  skip_on_cran()
  r = until_run("move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), adapt = scaled_am(until = iterations(500)))",
                "until_iterations", 6000)
  expect_equal(stopped_at(r), 500)
  s = scale_trace(r)
  expect_true(length(unique(s[1:500])) > 100)
  expect_equal(length(unique(s[-(1:500)])), 1)
  # the posterior (exact: mean 1/2.09 for each of a and b)
  o = r$output
  for (p in c("a", "b"))
    expect_equal(mean(o$Value[o$ParameterName == p & o$Iteration > 1000]), 1/2.09, tolerance = 0.1)
})

test_that("stable, acceptance, all and any stop adaptation", {
  skip_on_cran()
  # stable: stops before the end, and the parameters stay fixed after
  r = until_run(paste0("move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), ",
                       "adapt = scaled_am(until = stable(0.2, 50)))"), "until_stable")
  k = stopped_at(r)
  expect_true(k > 50 && k < 3000)
  expect_equal(length(unique(scale_trace(r)[-(1:k)])), 1)
  # any: the first criterion that holds; all: the last
  r = until_run(paste0("move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), ",
                       "adapt = scaled_am(until = any(iterations(300), stable(0, 50))))"), "until_any", 700)
  expect_equal(stopped_at(r), 300)
  r = until_run(paste0("move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), ",
                       "adapt = scaled_am(until = all(iterations(300), iterations(600))))"), "until_all", 700)
  expect_equal(stopped_at(r), 600)
  # watching chosen quantities, with absolute changes
  r = until_run(paste0("move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), ",
                       "adapt = scaled_am(until = stable(0.5, 50, watch = [scale, mean], change = absolute)))"), "until_watch", 1000)
  expect_true(stopped_at(r) > 50)
  # acceptance: the average acceptance near the target (scale's default, 0.234
  # in dimension 2 is not used: scale's target is given here)
  r = until_run(paste0("move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), ",
                       "adapt = scale(target = 0.3, until = acceptance(0.05, 100)))"), "until_acceptance", 2000)
  k = stopped_at(r)
  expect_true(k >= 100 && k < 2000)
  expect_equal(length(unique(scale_trace(r)[-(1:k)])), 1)
})

test_that("element moves stop adapting each element on its own", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "until_elements"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.5, 1.5, -1] }", "model {", "  for (k in 1:3) theta[k] ~ norm(0, 1)",
               "  for (k in 1:3) y[k] ~ norm(theta[k], 0.5)", "}", "algorithm mcmc {",
               "  for (k in 1:3) move m = metropolis(norm_rw(theta[k], 1), adapt = scale(until = iterations(200)))",
               "  termination = iterations(1000)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "until_elements", seed = 1, initial_values = list(list(theta = c(0, 0, 0))))
  final = ilike::load_mcmc_adaptation("until_elements")$final
  expect_equal(final$Value[grepl("adaptation_stopped_at$", final$Parameter)], c(200, 200, 200))
  expect_true(all(c("theta[1].adaptation_stopped_at", "theta[3].adaptation_stopped_at") %in% final$Parameter))
})

test_that("criteria that cannot apply are refused clearly", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "until_refused"), recursive = TRUE), add = TRUE)
  compile_move = function(move)
  {
    writeLines(c(until_model, "algorithm mcmc {", paste0("  ", move), "  termination = iterations(50)", "}"), f)
    ilike::compile(f)
  }
  rw = "move m = metropolis(mvnorm_rw([a, b], [[1, 0], [0, 1]]), adapt = "
  expect_error(compile_move(paste0(rw, "scaled_am(until = forever()))")), "until is a criterion: iterations(n)", fixed = TRUE)
  expect_error(compile_move(paste0(rw, "scaled_am(until = stable(0.1, 0)))")), "the window is a whole number from 1", fixed = TRUE)
  expect_error(compile_move(paste0(rw, "scaled_am(until = stable(0.1, 10, watch = speed)))")), "stable: watch is default", fixed = TRUE)
  expect_error(compile_move(paste0(rw, "scaled_am(until = iterations(0)))")), "n is a whole number from 1", fixed = TRUE)
  run = function(move)
  {
    unlink("until_refused", recursive = TRUE)
    ilike::MCMC(compile_move(move), results_name = "until_refused", seed = 1, initial_values = list(list(a = 0, b = 0)))
  }
  expect_error(run(paste0(rw, "am(until = acceptance(0.05, 10)))")), "has no target acceptance rate", fixed = TRUE)
  expect_error(run(paste0(rw, "scaled_am(until = stable(0.1, 10, watch = component_scales)))")),
               "cannot watch component_scales (it can watch: scale, covariance, covariance_elements, mean)", fixed = TRUE)
})
