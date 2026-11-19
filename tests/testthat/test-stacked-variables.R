# Proposals over several named variables stacked into one block, e.g.
# mvnorm_rw([a, v, c], cov) (dev/plan.md C203). The stack is the column vector
# of the variables' elements in the order written, so a move on [a, v, c] is
# the same move as on a vector theta = (a, v[1], v[2], c): with the same seed
# the chains are identical.

stacked_model = function(named, move, iterations = 2000)
{
  x = if (named) c("a", "v[1]", "v[2]", "c") else paste0("theta[", 1:4, "]")
  c("ilike 1.0", "data { y = [0.5, -1.2, 2.0] }", "model {", paste0("  ", x, " ~ norm(0, 1)"),
    paste0("  y[1] ~ norm(", x[1], " + ", x[2], ", 1)"), paste0("  y[2] ~ norm(", x[3], " - ", x[4], ", 0.5)"),
    paste0("  y[3] ~ norm(", x[1], " + ", x[4], ", 1)"), "}", "algorithm mcmc {",
    paste0("  move m = ", gsub("VARS", if (named) "[a, v, c]" else "theta", move)),
    paste0("  termination = iterations(", iterations, ")"), "}")
}

run_stacked = function(lines, initial, name)
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  writeLines(lines, f)
  ilike::MCMC(ilike::compile(f), results_name = name, seed = 3, initial_values = list(initial))
  o = ilike::load_mcmc_output(name)
  matrix(o$Value, ncol = length(unlist(initial)), byrow = TRUE)
}

I4 = "[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]"
named_start = list(a = 0, v = c(0, 0), c = 0)
vector_start = list(theta = rep(0, 4))

test_that("a stacked random walk is the random walk on the vector of the variables", {
  move = paste0("metropolis(mvnorm_rw(VARS, ", I4, ", scale = 0.5))")
  named = run_stacked(stacked_model(TRUE, move), named_start, "stacked_rw_named")
  vector = run_stacked(stacked_model(FALSE, move), vector_start, "stacked_rw_vector")
  expect_identical(named, vector)
})

test_that("stacked gradient-based and other proposals equal their vector versions", {
  skip_on_cran()
  moves = c(langevin = paste0("metropolis_hastings(langevin(VARS, ", I4, ", 0.5))"),
            barker = paste0("metropolis_hastings(barker(VARS, ", I4, ", 0.5))"),
            mirror = paste0("metropolis_hastings(mirror(VARS, [0.5, 0, -0.3, 0.9], ", I4, ", 0.3))"),
            componentwise = paste0("metropolis(componentwise_rw(VARS, ", I4, ", 1))"),
            pc = paste0("metropolis(pc_rw(VARS, ", I4, ", 1))"),
            scaled_am = paste0("metropolis(mvnorm_rw(VARS, ", I4, ", 0.5), adapt = scaled_am())"))
  for (k in names(moves))
  {
    named = run_stacked(stacked_model(TRUE, moves[k]), named_start, paste0("stacked_", k, "_named"))
    vector = run_stacked(stacked_model(FALSE, moves[k]), vector_start, paste0("stacked_", k, "_vector"))
    expect_identical(named, vector, info = k)
  }
  # and the Langevin chain targets the exact Gaussian posterior
  X = rbind(c(1, 1, 0, 0), c(0, 0, 1, -1), c(1, 0, 0, 1))
  W = diag(c(1, 4, 1))
  P = diag(4) + t(X) %*% W %*% X
  mu = c(solve(P, t(X) %*% W %*% c(0.5, -1.2, 2.0)))
  named = run_stacked(stacked_model(TRUE, moves["langevin"], 10000), named_start, "stacked_langevin_long")[-(1:1000), ]
  expect_lt(max(abs(colMeans(named) - mu)/sqrt(diag(solve(P)))), 0.1)
})

test_that("a stacked matrix variable keeps its shape", {
  skip_on_cran()
  I5 = "[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]"
  matrix_model = c("ilike 1.0", "data { y = [0.3, -0.4, 1.1, 0.8] }", "model {",
                   "  for (i in 1:2) for (j in 1:2) M[i, j] ~ norm(0, 1)", "  a ~ norm(0, 1)",
                   "  for (j in 1:2) for (i in 1:2) y[i + 2*(j - 1)] ~ norm(M[i, j] + a, 1)", "}", "algorithm mcmc {",
                   paste0("  move m = metropolis(mvnorm_rw([M, a], ", I5, ", 0.3))"), "  termination = iterations(500)", "}")
  vector_model = c("ilike 1.0", "data { y = [0.3, -0.4, 1.1, 0.8] }", "model {", "  for (k in 1:5) t[k] ~ norm(0, 1)",
                   "  for (k in 1:4) y[k] ~ norm(t[k] + t[5], 1)", "}", "algorithm mcmc {",
                   paste0("  move m = metropolis(mvnorm_rw(t, ", I5, ", 0.3))"), "  termination = iterations(500)", "}")
  stacked = run_stacked(matrix_model, list(M = matrix(1:4, 2, 2), a = 5), "stacked_matrix")
  vector = run_stacked(vector_model, list(t = 1:5), "stacked_matrix_vector")
  # the stack and the output both take M column by column (M starts at
  # matrix(1:4, 2, 2) and t at 1:5, so the chains are identical only then)
  expect_identical(stacked, vector)
})

test_that("adaptation on a stacked block conditions on the variables outside it", {
  skip_on_cran()
  I3 = "[[1, 0, 0], [0, 1, 0], [0, 0, 1]]"
  model = function(named, algorithm, adapt)
  {
    x = if (named) c("a", "v[1]", "v[2]") else paste0("u[", 1:3, "]")
    c("ilike 1.0", "data { y = [0.5, -1.2, 2.0] }", "model {", paste0("  ", x, " ~ norm(0, 1)"), "  c ~ norm(0, 1)",
      paste0("  y[1] ~ norm(", x[1], " + ", x[2], ", 0.3)"), paste0("  y[2] ~ norm(", x[3], " - c, 0.5)"),
      paste0("  y[3] ~ norm(", x[1], " + c, 0.2)"), "}", paste0("algorithm ", algorithm, " {"),
      paste0("  move m1 = metropolis(mvnorm_rw(", if (named) "[a, v]" else "u", ", ", I3, ", 0.5), ", adapt, ")"),
      paste0("  move m2 = metropolis(norm_rw(c, 1), ", adapt, ")"),
      if (algorithm == "mcmc") "  termination = iterations(3000)"
      else c("  termination = iterations(5)", "  sequence = annealing(power, [0, 0.01, 0.1, 0.3, 1])"), "}")
  }
  # conditional adaptive Metropolis: the block [a, v] given c is the vector u given c
  named = run_stacked(model(TRUE, "mcmc", "adapt = am()"), list(a = 0, v = c(0, 0), c = 0), "stacked_am_named")
  vector = run_stacked(model(FALSE, "mcmc", "adapt = am()"), list(u = rep(0, 3), c = 0), "stacked_am_vector")
  expect_identical(named, vector)
  # SMC with the adaptor on the block: the same evidence, which is the exact one
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "stacked_smc"), recursive = TRUE), add = TRUE)
  log_z = sapply(c(TRUE, FALSE), function(named) {
    writeLines(model(named, "smc", "smc_adapt = gaussian(2.38)"), f)
    unlink("stacked_smc", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "stacked_smc", number_of_particles = 2000, seed = 3) })
  expect_identical(log_z[1], log_z[2])
  X = rbind(c(1, 1, 0, 0), c(0, 0, 1, -1), c(1, 0, 0, 1))
  S = diag(c(0.3, 0.5, 0.2)^2) + X %*% t(X)
  y = c(0.5, -1.2, 2.0)
  exact = c(-0.5*determinant(2*pi*S)$modulus - 0.5*t(y) %*% solve(S, y))
  expect_lt(abs(log_z[1] - exact), 0.1)
})

test_that("every proposal takes a list, and a scalar proposal acts on each element", {
  skip_on_cran()
  # (C205) on a block of several elements, a one-dimensional proposal applies
  # to each element independently: a list [a, b] and a vector theta give the
  # same chain; the posterior of a, b is N(y/2, 1/2) (or a gamma prior's)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "stacked_scalar"), recursive = TRUE), add = TRUE)
  run = function(move, named, prior = "norm(0, 1)")
  {
    x = if (named) c("a", "b") else c("theta[1]", "theta[2]")
    writeLines(c("ilike 1.0", "data { y = [1, 2] }", "model {", paste0("  ", x, " ~ ", prior), paste0("  y[", 1:2, "] ~ norm(", x, ", 1)"), "}",
                 "algorithm mcmc {", paste0("  move m = ", gsub("VARS", if (named) "[a, b]" else "theta", move)),
                 "  termination = iterations(4000)", "}"), f)
    unlink("stacked_scalar", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "stacked_scalar", seed = 1,
                initial_values = list(if (named) list(a = 0.5, b = 0.5) else list(theta = c(0.5, 0.5))))
    matrix(ilike::load_mcmc_output("stacked_scalar")$Value, ncol = 2, byrow = TRUE)
  }
  moves = c("metropolis(norm_rw(VARS, 0.8))", "metropolis(unif_rw(VARS, 1))", "metropolis(mvunif_rw(VARS, [1, 0.8]))",
            "metropolis(mvnorm_rw(VARS, [[0.5]]))", "metropolis(norm_rw(VARS, 0.8), adapt = am())",
            "independent_mh(norm(VARS, 0.8, 1))", "independent_mh(mvnorm(VARS, [0.5, 1], [[1, 0], [0, 1]]))")
  for (move in moves)
  {
    named = run(move, TRUE)
    expect_identical(named, run(move, FALSE), info = move)
    expect_lt(max(abs(colMeans(named[-(1:500), ]) - c(0.5, 1))), 0.1, label = move)
  }
  posterior_mean = function(y) { x = seq(1e-4, 15, length = 2e5); w = dgamma(x, 2, 1)*dnorm(y, x, 1); sum(x*w)/sum(w) }
  for (move in c("independent_mh(gamma(VARS, 2, 1))", "independent_mh(lnorm(VARS, 0, 0.7))"))
  {
    named = run(move, TRUE, "gamma(2, 1)")
    expect_identical(named, run(move, FALSE, "gamma(2, 1)"), info = move)
    expect_lt(max(abs(colMeans(named[-(1:500), ]) - c(posterior_mean(1), posterior_mean(2)))), 0.1, label = move)
  }
})

test_that("a matrix variable keeps its shape when it is moved", {
  skip_on_cran()
  # the model reads M[1, 2] after each move: before C205 a Gaussian random
  # walk returned M as a column, and this failed
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "stacked_shape"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.3, -0.4, 1.1, 0.8] }", "model {", "  for (i in 1:2) for (j in 1:2) M[i, j] ~ norm(0, 1)",
               "  for (j in 1:2) for (i in 1:2) y[i + 2*(j - 1)] ~ norm(M[i, j], 1)", "}", "algorithm mcmc {",
               "  move m = metropolis(mvnorm_rw(M, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], 0.5))",
               "  termination = iterations(4000)", "}"), f)
  unlink("stacked_shape", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "stacked_shape", seed = 1, initial_values = list(list(M = matrix(0, 2, 2))))
  m = matrix(ilike::load_mcmc_output("stacked_shape")$Value, ncol = 4, byrow = TRUE)[-(1:500), ]
  expect_lt(max(abs(colMeans(m) - c(0.3, -0.4, 1.1, 0.8)/2)), 0.1)
})

test_that("a proposal of the wrong size, and lists where they cannot be used, are refused clearly", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "stacked_size"), recursive = TRUE), add = TRUE)
  model = c("ilike 1.0", "data { y = 1 }", "model {", "  a ~ norm(0, 1)", "  b ~ norm(0, 1)", "  y ~ norm(a + b, 1)", "}")
  check = function(move, message)
  {
    writeLines(c(model, "algorithm mcmc {", paste0("  move m = ", move), "  termination = iterations(10)", "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check("metropolis(mvnorm_rw([a, a], [[1, 0], [0, 1]]))", "a appears twice in the list of variables")
  check("metropolis(mvnorm_rw([a, 2], [[1, 0], [0, 1]]))", "each element of the list is the name of a variable")
  check("metropolis(mvnorm_rw([], [[1, 0], [0, 1]]))", "the list of variables is empty")
  writeLines(c(model, "algorithm is {", "  proposal = mvnorm([a, b], [0, 0], [[1, 0], [0, 1]])", "}"), f)
  expect_error(ilike::compile(f), "a list of variables can be moved together by the proposal of a move", fixed = TRUE)
  run = function(move)
  {
    writeLines(c(model, "algorithm mcmc {", paste0("  move m = ", move), "  termination = iterations(10)", "}"), f)
    unlink("stacked_size", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "stacked_size", seed = 1, initial_values = list(list(a = 0, b = 0)))
  }
  expect_error(run("metropolis(mvnorm_rw([a, b], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))"),
               "The proposal's covariance for [a, b] is 3 x 3, but [a, b] has 2 element(s).", fixed = TRUE)
  expect_error(run("metropolis_hastings(langevin([a, b], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 0.5))"),
               "The proposal's covariance for [a, b] is 3 x 3", fixed = TRUE)
  expect_error(run("metropolis(mvunif_rw([a, b], [1, 1, 1]))"),
               "The proposal has 3 half-width(s) for [a, b], but [a, b] has 2 element(s).", fixed = TRUE)
  expect_error(run("independent_mh(mvnorm([a, b], [0, 0, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))"),
               "The proposal has 3 mean value(s) for [a, b]", fixed = TRUE)
})
