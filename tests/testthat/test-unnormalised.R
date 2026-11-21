# Unnormalised distributions: x ~ ising(beta, nb, nn), x ~ potts(beta, q, nb, nn)
# (dev/plan.md 3.20, S3; C242), against enumeration on a 3 x 3 lattice

unnormalised_edges = function(g)
{
  e = do.call(rbind, lapply(seq_len(g$n), function(i) cbind(i, g$nb_list[[i]])))
  e[e[, 1] < e[, 2], ]
}

test_that("the simulators of ising and potts fields draw from them", {
  g = lattice(3); e = unnormalised_edges(g); beta = 0.6
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  # ising: spins s = 2x - 1, density proportional to exp(beta * sum over edges of s_a s_b)
  states = as.matrix(expand.grid(rep(list(0:1), 9)))
  lw = apply(states, 1, function(x) { s = 2*x - 1; beta * sum(s[e[, 1]] * s[e[, 2]]) })
  w = exp(lw - max(lw)); w = w / sum(w)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", "  x ~ ising(0.6, nb, nn)", "}"), f)
  m = ilike::compile(f)
  for (method in c("exact", "sweeps"))
  {
    d = ilike:::simulate_unnormalised(m, "x", draws = 5000, method = method, sweeps = 30, seed = 2)
    expect_equal(mean(d[, 1] == d[, 2]), sum(w * (states[, 1] == states[, 2])), tolerance = 0.03)
    expect_equal(mean(d[, 1] == d[, 9]), sum(w * (states[, 1] == states[, 9])), tolerance = 0.03)
  }
  # potts, q = 3: density proportional to exp(beta * the number of edges whose ends agree)
  states = as.matrix(expand.grid(rep(list(1:3), 9)))
  lw = apply(states, 1, function(x) beta * sum(x[e[, 1]] == x[e[, 2]]))
  w = exp(lw - max(lw)); w = w / sum(w)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", "  x ~ potts(0.6, 3, nb, nn)", "}"), f)
  d = ilike:::simulate_unnormalised(ilike::compile(f), "x", draws = 5000, method = "sweeps", sweeps = 30, seed = 2)
  expect_equal(mean(d[, 1] == d[, 2]), sum(w * (states[, 1] == states[, 2])), tolerance = 0.05)
})

# a noisily observed Potts field: the states are inferred, with beta fixed
noisy_potts = function()
{
  g = lattice(3); e = unnormalised_edges(g); beta = 0.6
  y = c(1.2, 0.8, 2.1, 1.1, 2.9, 2.2, 3.1, 0.7, 2.6)
  states = as.matrix(expand.grid(rep(list(1:3), 9)))
  lw = 9 * log(1/3) + apply(states, 1, function(x) beta * sum(x[e[, 1]] == x[e[, 2]])) +
       apply(states, 1, function(x) sum(dnorm(y, x, 1, log = TRUE)))
  w = exp(lw - max(lw)); w = w / sum(w)
  list(head = c("ilike 1.0", paste0("data { ", g$data, "; y = [", paste(y, collapse = ", "), "] }"), "model {",
                "  x ~ potts(0.6, 3, nb, nn)", "  for (i in 1:9) y[i] ~ norm(x[i], 1)", "}"),
       p12 = sum(w * (states[, 1] == states[, 2])), p5_3 = sum(w * (states[, 5] == 3)),
       log_z = max(lw) + log(sum(exp(lw - max(lw)))))
}

test_that("the states of a noisily observed potts field are inferred by gibbs moves and by SMC from the base measure", {
  skip_on_cran()
  p = noisy_potts()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_mcmc", "unnormalised_smc"), recursive = TRUE), add = TRUE)
  writeLines(c(p$head, "algorithm mcmc {", "  for (j in 1:9) move m = gibbs(x[j])", "  termination = iterations(20000)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "unnormalised_mcmc", seed = 1)
  o = ilike::load_mcmc_output("unnormalised_mcmc")
  X = matrix(o$Value[o$ParameterName == "x"], ncol = 9, byrow = TRUE)
  expect_equal(mean(X[, 1] == X[, 2]), p$p12, tolerance = 0.03)
  expect_equal(mean(X[, 5] == 3), p$p5_3, tolerance = 0.05)
  # SMC annealing the potential and the likelihood from the base measure (the
  # prior's simulation); its evidence is with respect to the unnormalised
  # prior (base measure times exp(potential): Z(beta) is not included)
  writeLines(c(p$head, "algorithm smc {", "  for (j in 1:9) move m = gibbs(x[j])", "  termination = iterations(5)",
               "  sequence = annealing(power, [0, 1])", "  target = cess(0.9)", "  particles = 5000", "}"), f)
  z = ilike::SMC_with_MCMC(ilike::compile(f), results_name = "unnormalised_smc", seed = 1)
  expect_lt(abs(z - p$log_z), 0.1)
  o = ilike::load_smc_output("unnormalised_smc")
  last = o[(o$Target == max(o$Target)) & (o$ParameterName == "x"), ]
  last = last[last$Iteration == max(last$Iteration), ]
  X = matrix(last$Value, ncol = 9, byrow = TRUE)
  lw = matrix(last$LogWeight, ncol = 9, byrow = TRUE)[, 1]
  w = exp(lw - max(lw)); w = w / sum(w)
  expect_equal(sum(w * (X[, 1] == X[, 2])), p$p12, tolerance = 0.03)
})

test_that("moves and algorithms that would need an unknown normalising constant are refused", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  head = c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
           "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn)", "}")
  writeLines(c(head, "algorithm mcmc {", "  move m = metropolis(norm_rw(beta, 0.1))", "  termination = iterations(10)", "}"), f)
  expect_error(ilike::compile(f), "a parameter of the normalising constant of x's distribution", fixed = TRUE)
  writeLines(c(head, "algorithm is {", "  importance_points = 50", "}"), f)
  expect_error(ilike::IS(ilike::compile(f), seed = 1), "vary between particles", fixed = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", "  x ~ potts(0.5, 2.5, nb, nn)", "}"), f)
  expect_error(ilike::compile(f), "q, the number of colours, is a whole number", fixed = TRUE)
})

# (C243, dev/plan.md 3.20 S4) the single and multiple auxiliary variable
# methods: x ~ ising(beta, ..., normaliser = sav(points, at)), so beta can be inferred
test_that("beta of an observed ising field is inferred with the single and multiple auxiliary variable methods", {
  skip_on_cran()
  g = lattice(4)
  e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(b) { a = b * S_all; b * stat(x) - (max(a) + log(sum(exp(a - max(a))))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  exact_mean = sum(p * grid); exact_sd = sqrt(sum(p * grid^2) - exact_mean^2)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_sav"), recursive = TRUE), add = TRUE)
  for (points in c(1, 5))
  {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", paste0("  x ~ ising(beta, nb, nn, normaliser = sav(points = ", points, ", at = 0.3))"), "}",
                 "algorithm mcmc {", "  move m = metropolis(norm_rw(beta, 1), transform = auto)", "  termination = iterations(20000)", "}"), f)
    unlink("unnormalised_sav", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "unnormalised_sav", seed = 1, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_sav")
    b = o$Value[o$ParameterName == "beta"][-(1:1000)]
    expect_equal(mean(b), exact_mean, tolerance = 0.05)
    expect_equal(sd(b), exact_sd, tolerance = 0.1)
  }
})

test_that("beta of a potts field is inferred with draws made by sweeps (approximate, with a note)", {
  skip_on_cran()
  g = lattice(3)
  e = unnormalised_edges(g)
  x = c(1, 1, 2, 1, 1, 2, 3, 3, 2)
  S_all = apply(as.matrix(expand.grid(rep(list(1:3), 9))), 1, function(z) sum(z[e[, 1]] == z[e[, 2]]))
  grid = seq(0.001, 1.999, by = 0.002)
  lp = sapply(grid, function(b) { a = b * S_all; b * sum(x[e[, 1]] == x[e[, 2]]) - (max(a) + log(sum(exp(a - max(a))))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_potts_sav"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
               "  beta ~ unif(0, 2)", "  x ~ potts(beta, 3, nb, nn, normaliser = sav(points = 5, at = 0.8, simulate = sweeps(30)))", "}",
               "algorithm mcmc {", "  move m = metropolis(norm_rw(beta, 1), transform = auto)", "  termination = iterations(20000)", "}"), f)
  expect_message(m <- ilike::compile(f), "heat-bath sweeps")
  ilike::MCMC(m, results_name = "unnormalised_potts_sav", seed = 1, initial_values = list(list(beta = 0.8)))
  o = ilike::load_mcmc_output("unnormalised_potts_sav")
  b = o$Value[o$ParameterName == "beta"][-(1:1000)]
  expect_equal(mean(b), sum(p * grid), tolerance = 0.05)
})

test_that("normaliser = sav(...) is refused where it does not apply", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  field = function(statement) c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
                                "  beta ~ unif(0, 1)", statement, "}")
  writeLines(field("  x ~ ising(0.4, nb, nn, normaliser = sav(at = 0.3))"), f)
  expect_error(ilike::compile(f), "depends on no parameter", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(points = 2))"), f)
  expect_error(ilike::compile(f), "at = value", fixed = TRUE)
  writeLines(field("  x ~ potts(beta, 2, nb, nn, normaliser = sav(at = 0.3, simulate = exact))"), f)
  expect_error(ilike::compile(f), "exact simulation is available for ising", fixed = TRUE)
})

test_that("proposals outside the prior's range are rejected before the field is drawn (early rejection)", {
  # (C244, dev/plan.md 3.21) beta ~ unif(0, 1) and a wide random walk: most
  # proposals leave [0, 1], where exact draws could not be made (beta < 0) or
  # would take very long (beta >> 1); the prior rejects them first
  g = lattice(4)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0] }"), "model {",
               "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3))", "}",
               "algorithm mcmc {", "  move m = metropolis(norm_rw(beta, 3))", "  termination = iterations(300)", "}"), f)
  expect_silent(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(beta = 0.3))))
})
