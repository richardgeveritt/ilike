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
      "  for (i in 1:5) y[i] ~ norm(0, sigma)", "}", paste0("algorithm ", algorithm, " {"), "  move m = nuts(sigma, step = 0.3)",
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
  check("metropolis_hastings(nuts(theta, step = 0.1))", "nuts is a move of its own, not a proposal: move m = nuts(x, step = ") # (ASCII only)
  check("metropolis(nuts(theta, step = 0.1))", "nuts is a move of its own, not a proposal")
})

test_that("nuts records its diagnostics, thinned with the points, and warns of divergences after warm-up", {
  # (3.18 item 9, C212)
  model = nuts_smooth()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_diagnostics"), recursive = TRUE), add = TRUE)
  run = function(move, iterations = 300, thin = 3)
  {
    writeLines(c(model$lines[seq_len(grep("^algorithm", model$lines) - 1)], "algorithm mcmc {", paste0("  move m = ", move),
                 paste0("  termination = iterations(", iterations, ")"), "}"), f)
    unlink("nuts_diagnostics", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "nuts_diagnostics", seed = 1, thin = thin,
                initial_values = list(list(theta = rep(0, model$K)), list(theta = rep(1, model$K))))
    ilike::load_mcmc_diagnostics("nuts_diagnostics")
  }
  expect_no_warning(d <- run("nuts(theta, step = 0.2)"))
  expect_equal(names(d), c("Chain", "Iteration", "tree_depth", "leapfrog_steps", "divergent", "accept_stat"))
  # one row per stored point: 300/3 per chain
  expect_equal(as.vector(table(d$Chain)), c(100, 100))
  # a tree of depth j has at most 2^j - 1 leapfrog steps (fewer if a subtree stopped early)
  expect_true(all(d$leapfrog_steps >= d$tree_depth & d$leapfrog_steps <= 2^d$tree_depth - 1))
  expect_true(all(d$accept_stat >= 0 & d$accept_stat <= 1) && !any(d$divergent))
  # a step far too large diverges, counted over every iteration (not only the stored ones)
  expect_warning(run("nuts(theta, step = 3)"), "of 600 transitions after warm-up were divergent")
  # divergences while the step adapts (until it stops) are not warned of
  d = expect_no_warning(run("nuts(theta, step = 3, adapt = scale(until = iterations(200)))", iterations = 400, thin = 1))
  expect_gt(sum(d$divergent[d$Iteration <= 200]), 0)
  expect_equal(sum(d$divergent[d$Iteration > 200]), 0)
})

test_that("dual averaging finds the step of nuts and hmc from a poor start, then fixes it", {
  # (3.18 stage H4, C219)
  model = nuts_smooth()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_dual"), recursive = TRUE), add = TRUE)
  run = function(move)
  {
    writeLines(c(model$lines[seq_len(grep("^algorithm", model$lines) - 1)], "algorithm mcmc {", paste0("  move m = ", move),
                 "  termination = iterations(3000)", "}"), f)
    unlink("nuts_dual", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "nuts_dual", seed = 1, initial_values = list(list(theta = rep(0, model$K))))
    list(o = ilike::load_mcmc_output("nuts_dual"), a = ilike::load_mcmc_adaptation("nuts_dual"),
         d = tryCatch(ilike::load_mcmc_diagnostics("nuts_dual"), error = function(e) NULL))
  }
  r = run("nuts(theta, step = 0.001, adapt = dual_averaging(until = iterations(500)))")
  final = r$a$final
  # the step after stopping is the average of the iterates, and stays there
  expect_equal(final$Value[final$Parameter == "theta.step"], final$Value[final$Parameter == "theta.step_average"])
  expect_equal(final$Value[final$Parameter == "theta.adaptation_stopped_at"], 500)
  trace = r$a$trace
  expect_true(all(trace$Value[trace$Parameter == "theta.step" & trace$Iteration > 500] == final$Value[final$Parameter == "theta.step"]))
  # far from the 0.001 given, and the acceptance statistic near its target 0.8
  expect_gt(final$Value[final$Parameter == "theta.step"], 0.1)
  expect_lt(abs(mean(r$d$accept_stat[r$d$Iteration > 500]) - 0.8), 0.1)
  theta = sapply(1:model$K, function(k) r$o$Value[r$o$Dimension == k & r$o$Iteration > 500])
  expect_lt(max(abs(colMeans(theta) - model$m)/sqrt(diag(model$S))), 0.15)
  # hmc, from a step far too large
  r = run("metropolis_hastings(hmc(theta, step = 2, steps = 5), adapt = dual_averaging(until = iterations(500)))")
  expect_lt(r$a$final$Value[r$a$final$Parameter == "theta.step"], 1)
})

test_that("dual averaging is refused without until, and for other proposals", {
  model = nuts_smooth()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(move, message)
  {
    writeLines(c(model$lines[seq_len(grep("^algorithm", model$lines) - 1)], "algorithm mcmc {", paste0("  move m = ", move),
                 "  termination = iterations(10)", "}"), f)
    expect_error(ilike::MCMC(ilike::compile(f), initial_values = list(list(theta = rep(0, model$K)))), message, fixed = TRUE)
  }
  check("nuts(theta, step = 0.1, adapt = dual_averaging())", "dual_averaging needs until = criterion")
  check("metropolis(norm_rw(theta, 0.1), adapt = dual_averaging(until = iterations(5)))",
        "applies to the proposals ilike::hmc and ilike::nuts")
})

test_that("dual averaging with mass windows estimates the inverse mass matrix of a badly scaled target", {
  # (3.18 stage H5, C220) a 10-d normal with sds from 0.01 to 10 and two correlated pairs
  d = 10; sds = exp(seq(log(0.01), log(10), length = d))
  R = diag(d); R[1, 2] = R[2, 1] = 0.9; R[5, 6] = R[6, 5] = -0.8
  S = diag(sds) %*% R %*% diag(sds)
  text = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(signif(r, 10), collapse = ", "), "]")), collapse = ", "), "]")
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_mass"), recursive = TRUE), add = TRUE)
  run = function(adapt)
  {
    writeLines(c("ilike 1.0", "model {", paste0("  theta ~ mvnorm([", paste(rep(0, d), collapse = ", "), "], ", text(S), ")"), "}",
                 "algorithm mcmc {", paste0("  move m = nuts(theta, step = 0.1, adapt = ", adapt, ")"),
                 "  termination = iterations(3000)", "}"), f)
    unlink("nuts_mass", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "nuts_mass", seed = 1, initial_values = list(list(theta = rep(0, d))))
    o = ilike::load_mcmc_output("nuts_mass")
    list(x = sapply(1:d, function(k) o$Value[o$Dimension == k & o$Iteration > 1000]),
         steps = mean(ilike::load_mcmc_diagnostics("nuts_mass")$leapfrog_steps[-(1:1000)]))
  }
  diagonal = run("dual_averaging(mass = diagonal, until = iterations(1000))")
  dense = run("dual_averaging(mass = dense, until = iterations(1000))")
  for (r in list(diagonal, dense))
  {
    expect_equal(apply(r$x, 2, sd), sds, tolerance = 0.1)
    expect_lt(abs(cor(r$x[, 1], r$x[, 2]) - 0.9), 0.03)
  }
  # (with the step alone, about 800 leapfrog steps per iteration: the
  # scales differ by 1000 times)
  expect_lt(diagonal$steps, 30)
  expect_lt(dense$steps, diagonal$steps)
  # the windows need the warm-up's length
  expect_error(run("dual_averaging(mass = dense, until = stable())"), "until = iterations(n)", fixed = TRUE)
})

test_that("the paper's variant of NUTS targets the posterior", {
  # (C221: Stan's variant, multinomial sampling and the generalised U-turn
  # criterion, is the default, used by the other tests; the paper's
  # Algorithm 3, slice sampling, is variant = paper)
  skip_on_cran()
  model = nuts_smooth()
  theta = nuts_run(model, "nuts(theta, step = 0.2, variant = paper)", "nuts_paper", iterations = 6000)[-(1:600), ]
  expect_lt(max(abs(colMeans(theta) - model$m)/sqrt(diag(model$S))), 0.1)
  expect_equal(apply(theta, 2, sd), sqrt(diag(model$S)), tolerance = 0.05)
})

test_that("Stan's variant records its diagnostics, and an unknown variant is refused", {
  model = nuts_smooth()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_stan_diagnostics"), recursive = TRUE), add = TRUE)
  run = function(move)
  {
    writeLines(c(model$lines[seq_len(grep("^algorithm", model$lines) - 1)], "algorithm mcmc {", paste0("  move m = ", move),
                 "  termination = iterations(200)", "}"), f)
    unlink("nuts_stan_diagnostics", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "nuts_stan_diagnostics", seed = 1,
                initial_values = list(list(theta = rep(0, model$K))))
    ilike::load_mcmc_diagnostics("nuts_stan_diagnostics")
  }
  d = run("nuts(theta, step = 0.2, variant = stan)")
  expect_true(all(d$leapfrog_steps >= 1 & d$leapfrog_steps <= 2^d$tree_depth - 1))
  expect_true(all(d$accept_stat >= 0 & d$accept_stat <= 1))
  expect_error(run("nuts(theta, step = 0.2, variant = other)"), "variant is paper or stan", fixed = TRUE)
})

test_that("hmc and moves inside scans record their diagnostics, and the mass matrix is in the adaptation output", {
  # (C226)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "nuts_c226"), recursive = TRUE), add = TRUE)
  run = function(moves, extra = character(0), iterations = 200)
  {
    writeLines(c("ilike 1.0", "data { y = [0.3, 1.2] }", "model {", "  mu ~ norm(0, 1)", "  s ~ gamma(2, 2)",
                 "  for (i in 1:2) y[i] ~ norm(mu, s)", "}", "algorithm mcmc {", paste0("  move ", moves), extra,
                 paste0("  termination = iterations(", iterations, ")"), "}"), f)
    unlink("nuts_c226", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "nuts_c226", seed = 1, initial_values = list(list(mu = 0, s = 1)))
    ilike::load_mcmc_diagnostics("nuts_c226")
  }
  d = expect_no_warning(run("m = metropolis_hastings(hmc([mu, s], step = 0.2, steps = 5))"))
  expect_equal(names(d), c("Chain", "Iteration", "leapfrog_steps", "accept_stat", "energy_error", "divergent"))
  expect_true(all(d$leapfrog_steps == 5) && !any(d$divergent))
  # the acceptance probability is min(1, exp(-energy error))
  expect_equal(d$accept_stat, pmin(1, exp(-d$energy_error)))
  expect_warning(run("m = metropolis_hastings(hmc([mu, s], step = 3, steps = 5))"), "of 200 transitions after warm-up were divergent")
  # a scan: each move's columns, prefixed
  d = run(c("a = nuts(mu, step = 0.3)", "b = metropolis(norm_rw(s, 0.3))"))
  expect_equal(names(d)[-(1:2)], c("move1.tree_depth", "move1.leapfrog_steps", "move1.divergent", "move1.accept_stat"))
  expect_true(all(d$move1.tree_depth >= 1))
  # a stochastic scan: the columns of a move not chosen are NA
  d = run(c("a = nuts(mu, step = 0.3)", "b = metropolis(norm_rw(s, 0.3))"), "  weights = [0.5, 0.5]")
  expect_true(any(is.na(d$move1.tree_depth)) && any(!is.na(d$move1.tree_depth)))
  # the inverse mass matrix estimated in windows, at the end of the chain
  run("m = nuts([mu, s], step = 0.3, adapt = dual_averaging(mass = dense, until = iterations(150)))", iterations = 300)
  final = ilike::load_mcmc_adaptation("nuts_c226")$final
  expect_true(all(c("mu;s.mass[1,1]", "mu;s.mass[1,2]", "mu;s.mass[2,2]") %in% final$Parameter) ||
              all(c("[mu, s].mass[1,1]", "[mu, s].mass[1,2]", "[mu, s].mass[2,2]") %in% final$Parameter))
})
