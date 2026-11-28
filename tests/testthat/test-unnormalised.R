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
  # (C299, dev/plan.md 3.22 U4 step 3) exact potts draws, by coupling from the
  # past on the random-cluster representation, below and above the critical
  # value (log(1 + sqrt(3)) for q = 3)
  for (b in c(0.6, 1.4))
  {
    lw = apply(states, 1, function(x) b * sum(x[e[, 1]] == x[e[, 2]]))
    w = exp(lw - max(lw)); w = w / sum(w)
    writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", paste0("  x ~ potts(", b, ", 3, nb, nn)"), "}"), f)
    d = ilike:::simulate_unnormalised(ilike::compile(f), "x", draws = 5000, method = "exact", seed = 2)
    expect_equal(mean(d[, 1] == d[, 2]), sum(w * (states[, 1] == states[, 2])), tolerance = 0.03)
    expect_equal(mean(d[, 1] == d[, 9]), sum(w * (states[, 1] == states[, 9])), tolerance = 0.05)
    expect_equal(mean(apply(d, 1, function(r) all(r == r[1]))), sum(w * apply(states, 1, function(r) all(r == r[1]))), tolerance = 0.1)
  }
})

test_that("beta of an observed potts field by exchange with exact draws, against enumeration", {
  skip_on_cran()
  g = lattice(3); e = unnormalised_edges(g)
  states = as.matrix(expand.grid(rep(list(1:3), 9))); A = apply(states, 1, function(x) sum(x[e[, 1]] == x[e[, 2]]))
  y = c(1, 1, 2, 1, 1, 2, 3, 3, 2); ay = sum(y[e[, 1]] == y[e[, 2]])
  grid = seq(0.001, 2.999, by = 0.002)
  lp = sapply(grid, function(b) { a = b * A; b * ay - (max(a) + log(sum(exp(a - max(a))))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_potts"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(y, collapse = ", "), "] }"), "model {", "  beta ~ unif(0, 3)",
               "  x ~ potts(beta, 3, nb, nn)", "}", "algorithm mcmc {", "  move b = exchange(norm_rw(beta, 1), simulate = exact, transform = auto)",
               "  termination = iterations(20000)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "unnormalised_potts", seed = 1, initial_values = list(list(beta = 0.5)))
  o = ilike::load_mcmc_output("unnormalised_potts")
  b = o$Value[o$ParameterName == "beta"][-(1:1000)]
  expect_equal(mean(b), sum(p * grid), tolerance = 0.05)
  expect_equal(sd(b), sqrt(sum(p * grid^2) - sum(p * grid)^2), tolerance = 0.1)
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
  # (C283, dev/plan.md 3.22 U0) and the annealed estimate, with 5 bridges
  for (method in c("points = 1", "points = 5", "points = 1, bridges = 5"))
  {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", paste0("  x ~ ising(beta, nb, nn, normaliser = sav(", method, ", at = 0.3))"), "}",
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
  # (C299, dev/plan.md 3.22 U4 step 3: potts fields have exact draws too)
  writeLines(field("  x ~ potts(beta, 2, nb, nn, normaliser = sav(at = 0.3, simulate = exact))"), f)
  expect_no_error(suppressMessages(ilike::compile(f)))
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, bridge_sweeps = 2))"), f)
  expect_error(ilike::compile(f), "goes with bridges = K", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, bridges = 1.5))"), f)
  expect_error(ilike::compile(f), "bridges is a whole number of at least 0", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(points = 5, at = 0.3, resampling = ess(0.5)))"), f)
  expect_error(ilike::compile(f), "needs bridges = K (K ≥ 1) and points = M (M ≥ 2)", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(points = 5, at = 0.3, bridges = 2, resampling = ess(2)))"), f)
  expect_error(ilike::compile(f), "resampling = ess(p), with 0 < p ≤ 1", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, normalise = is()))"), f)
  expect_error(ilike::compile(f), "normalise = smc(particles = P", fixed = TRUE)
})

# (C283, dev/plan.md 3.22 U0) the annealed estimate (Murray, Ghahramani &
# MacKay 2006): exp(estimate) is unbiased for Z(at)/Z(beta) for any number of
# bridges, with a variance that falls as they are added. Checked against
# enumeration on a 4 x 4 lattice, through evaluate_log (whose seed now seeds
# the draws: before C283 every seed gave the same value)
test_that("the annealed estimate of a normalising constant ratio is unbiased, and less variable with more bridges", {
  g = lattice(4)
  e = unnormalised_edges(g)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  lZ = function(b) { a = b * S_all; max(a) + log(sum(exp(a - max(a)))) }
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  variances = c()
  # (C293, dev/plan.md 3.22 U3c) and with the annealed runs resampled (resampling = ess(p), C294)
  for (bridges in c("0", "10", "10, points = 5, resampling = ess(0.5)"))
  {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {", "  beta ~ unif(0, 1)",
                 paste0("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, bridges = ", bridges, "))"), "}"), f)
    r = ilike::compile(f)
    # the estimate: the log target less the prior (0), the potential and the base measure
    v = sapply(1:1000, function(s) ilike::evaluate_log(r, parameters = list(beta = 0.6), seed = s)) - 0.6 * stat(x) + 16 * log(2)
    w = exp(v)
    expect_lt(abs(mean(w) - exp(lZ(0.3) - lZ(0.6))), 4 * sd(w) / sqrt(length(w)))
    variances = c(variances, var(v))
  }
  expect_lt(variances[2], variances[1] / 2)
  expect_lt(variances[3], variances[1] / 2)
})

# (C283, dev/plan.md 3.22 U0) the single and multiple auxiliary variable
# methods with a latent field: y_i ~ N(2 x_i - 1, 1) observing a 3 x 3 ising
# field x, with beta ~ U(0, 1). The estimate reads only beta, so the gibbs
# moves on x leave it alone; the exact posterior is by enumeration (2^9 fields
# on a grid of beta). Also SMC from the prior, whose evidence includes log Z(at)
test_that("beta and a latent ising field are inferred with the auxiliary variable methods, by MCMC and SMC", {
  skip_on_cran()
  g = lattice(3)
  e = unnormalised_edges(g)
  y = c(0.9, 1.3, -0.2, 0.7, 1.1, -0.8, 0.4, -1.2, -0.6)
  states = as.matrix(expand.grid(rep(list(0:1), 9)))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) })
  L = apply(states, 1, function(z) sum(dnorm(y, 2*z - 1, 1, log = TRUE)))
  grid = seq(0.0005, 0.9995, by = 0.001)
  lj = sapply(grid, function(b) { a = b * S; a + L - (max(a) + log(sum(exp(a - max(a))))) })
  top = max(lj); w = exp(lj - top); w = w / sum(w)
  pb = colSums(w); exact_mean = sum(pb * grid); exact_sd = sqrt(sum(pb * grid^2) - exact_mean^2)
  p12 = sum(rowSums(w) * (states[, 1] == states[, 2]))
  a = 0.3 * S
  log_evidence = top + log(sum(exp(lj - top))) + log(0.001) - 9 * log(2) + max(a) + log(sum(exp(a - max(a))))
  head = function(points) c("ilike 1.0", paste0("data { ", g$data, "; y = [", paste(y, collapse = ", "), "] }"), "model {",
                            "  beta ~ unif(0, 1)", paste0("  x ~ ising(beta, nb, nn, normaliser = sav(points = ", points, ", at = 0.3))"),
                            "  for (i in 1:9) y[i] ~ norm(2 * x[i] - 1, 1)", "}")
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_latent", "unnormalised_latent_smc"), recursive = TRUE), add = TRUE)
  writeLines(c(head(1), "algorithm mcmc {", "  for (j in 1:9) move g = gibbs(x[j])",
               "  move m = metropolis(norm_rw(beta, 1), transform = auto)", "  termination = iterations(40000)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "unnormalised_latent", seed = 1, initial_values = list(list(beta = 0.3, x = rep(1, 9))))
  o = ilike::load_mcmc_output("unnormalised_latent")
  keep = o$Iteration > 2000
  b = o$Value[o$ParameterName == "beta" & keep]
  X = matrix(o$Value[o$ParameterName == "x" & keep], ncol = 9, byrow = TRUE)
  expect_equal(mean(b), exact_mean, tolerance = 0.03)
  expect_equal(sd(b), exact_sd, tolerance = 0.05)
  expect_equal(mean(X[, 1] == X[, 2]), p12, tolerance = 0.02)
  # (C284, dev/plan.md 3.22 U1) the exchange algorithm on the same model, with
  # no estimate of the normalising constant
  writeLines(c(sub(", normaliser = sav(points = 1, at = 0.3)", "", head(1), fixed = TRUE), "algorithm mcmc {",
               "  for (j in 1:9) move g = gibbs(x[j])", "  move b = exchange(norm_rw(beta, 1), transform = auto)",
               "  termination = iterations(40000)", "}"), f)
  unlink("unnormalised_latent", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "unnormalised_latent", seed = 1, initial_values = list(list(beta = 0.3, x = rep(1, 9))))
  o = ilike::load_mcmc_output("unnormalised_latent")
  keep = o$Iteration > 2000
  b = o$Value[o$ParameterName == "beta" & keep]
  expect_equal(mean(b), exact_mean, tolerance = 0.03)
  expect_equal(sd(b), exact_sd, tolerance = 0.05)
  writeLines(c(head(5), "algorithm smc {", "  for (j in 1:9) move g = gibbs(x[j])", "  move m = metropolis(norm_rw(beta, 0.3), transform = auto)",
               "  sequence = annealing(power, [0, 1])", "  target = cess(0.9)", "  particles = 2000", "  termination = iterations(5)", "}"), f)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), results_name = "unnormalised_latent_smc", seed = 1))
  expect_lt(abs(z - log_evidence), 0.1)
  o = ilike::load_smc_output("unnormalised_latent_smc")
  last = o[o$Target == max(o$Target), ]
  last = last[(last$Iteration == max(last$Iteration)) & (last$ParameterName == "beta"), ]
  lw = exp(last$LogWeight - max(last$LogWeight)); lw = lw / sum(lw)
  expect_equal(sum(lw * last$Value), exact_mean, tolerance = 0.03)
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

# (C284, dev/plan.md 3.22 U1) the exchange algorithm (Murray, Ghahramani &
# MacKay 2006): beta of an observed ising field, against the exact posterior
# by enumeration (the C243 model); with exact draws it is exact
test_that("beta of an observed ising field is inferred with the exchange algorithm", {
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
  on.exit(unlink(c(f, "unnormalised_exchange"), recursive = TRUE), add = TRUE)
  # (C285, dev/plan.md 3.22 U2) and with the averaged acceptance ratio, points = 5;
  # (C292, U3) and annealed, bridges = 5
  for (points in c("1", "5", "1, bridges = 5"))
  {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn)", "}",
                 "algorithm mcmc {", paste0("  move b = exchange(norm_rw(beta, 1), points = ", points, ", transform = auto)"),
                 "  termination = iterations(20000)", "}"), f)
    unlink("unnormalised_exchange", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "unnormalised_exchange", seed = 1, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_exchange")
    b = o$Value[o$ParameterName == "beta"][-(1:1000)]
    expect_equal(mean(b), exact_mean, tolerance = 0.05)
    expect_equal(sd(b), exact_sd, tolerance = 0.1)
  }
})

# (C285, dev/plan.md 3.22 U2) the averaged acceptance ratio (MHAAR; Andrieu,
# Yildirim, Doucet & Chopin, arXiv 2101.01253, Algorithm 1) as ratio_estimating_mcmc.cpp
# computes it, checked exactly: the transition matrix on beta in {0.2, 0.5, 0.8}
# (a symmetric proposal to one of the other two values, an observed 2 x 2 ising
# field), enumerating the N fields drawn, leaves the posterior invariant and is
# reversible; the plain average (always the first direction: Alquier et al.'s
# noisy exchange) is not
test_that("the averaged acceptance ratio of the exchange move leaves the posterior invariant", {
  states = as.matrix(expand.grid(rep(list(0:1), 4)))
  edges = rbind(c(1, 2), c(3, 4), c(1, 3), c(2, 4))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[edges[, 1]] * s[edges[, 2]]) })
  thetas = c(0.2, 0.5, 0.8); prior = c(0.5, 0.3, 0.2); y = 1
  lg = function(b, i) b * S[i]
  f = function(b) { w = exp(b * S); w / sum(w) }
  post = prior * sapply(thetas, function(b) exp(lg(b, y)) / sum(exp(b * S))); post = post / sum(post)
  kernel = function(N, noisy = FALSE) {
    P = matrix(0, 3, 3)
    tuples = as.matrix(expand.grid(rep(list(1:16), N)))
    for (a in 1:3) for (b in setdiff(1:3, a)) {
      th = thetas[a]; vt = thetas[b]
      lbase = log(prior[b]) - log(prior[a]) + lg(vt, y) - lg(th, y)
      # c = 1: all N fields at vartheta
      a1 = sum(apply(tuples, 1, function(t) prod(f(vt)[t]) * min(1, exp(lbase) * mean(exp(lg(th, t) - lg(vt, t))))))
      if (noisy) { P[a, b] = 0.5 * a1; next }
      # c = 2: the first at vartheta, the others at theta
      a2 = sum(apply(tuples, 1, function(t) f(vt)[t[1]] * prod(f(th)[t[-1]]) * min(1, exp(lbase) / mean(exp(lg(vt, t) - lg(th, t))))))
      P[a, b] = 0.5 * (0.5 * a1 + 0.5 * a2)
    }
    diag(P) = 1 - rowSums(P)
    P
  }
  acceptance = c()
  for (N in 1:3)
  {
    P = kernel(N)
    expect_lt(max(abs(post %*% P - post)), 1e-12)
    expect_lt(max(abs(post * P - t(post * P))), 1e-12)
    acceptance = c(acceptance, sum(post * (1 - diag(P))))
  }
  expect_true(all(diff(acceptance) > 0))
  expect_gt(max(abs(post %*% kernel(3, noisy = TRUE) - post)), 0.01)
})

test_that("exchange(...) is refused where it does not apply, and notes approximate draws", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  model = function(statement, move) c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
                                      "  beta ~ unif(0, 1)", "  mu ~ norm(0, 1)", statement, "}", "algorithm mcmc {", move,
                                      "  termination = iterations(10)", "}")
  # (C290: exchange with normaliser = sav(…) is allowed, the estimate left out of its ratio)
  writeLines(model("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3))", "  move b = exchange(norm_rw(beta, 0.2))"), f)
  expect_no_error(ilike::compile(f))
  writeLines(model("  x ~ ising(beta, nb, nn)", "  move b = exchange(norm_rw(mu, 0.2))"), f)
  expect_error(ilike::compile(f), "changes no parameter of the normalising constant", fixed = TRUE)
  writeLines(model("  x ~ potts(beta, 2, nb, nn)", "  move b = exchange(norm_rw(beta, 0.2), simulate = exact)"), f)
  expect_no_error(suppressMessages(ilike::compile(f)))
  writeLines(model("  x ~ potts(beta, 2, nb, nn)", "  move b = exchange(norm_rw(beta, 0.2), simulate = sweeps(20))"), f)
  expect_message(ilike::compile(f), "draws x by 20 heat-bath sweeps", fixed = TRUE)
  writeLines(model("  x ~ ising(beta, nb, nn)", "  move b = exchange(norm_rw(beta, 0.2), bridge_sweeps = 2)"), f)
  expect_error(ilike::compile(f), "goes with bridges = K", fixed = TRUE)
  writeLines(model("  x ~ ising(beta, nb, nn)", "  move b = exchange(norm_rw(beta, 0.2), points = 0)"), f)
  expect_error(ilike::compile(f), "points is a whole number of at least 1", fixed = TRUE)
})

# (C286, dev/plan.md 3.22 U2) the N draws of the averaged acceptance ratio
# use their own random-number streams (a seed from the chain, jumped once per
# draw), so making them in parallel gives exactly the sequential chain
test_that("exchange(..., points = N, parallel = true) gives the same chain as without parallel", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_parallel"), recursive = TRUE), add = TRUE)
  chain = function(parallel) {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
                 "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn)", "}", "algorithm mcmc {",
                 paste0("  move b = exchange(norm_rw(beta, 0.5), points = 5", if (parallel) ", parallel = true", ", transform = auto)"),
                 "  termination = iterations(300)", "}"), f)
    unlink("unnormalised_parallel", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "unnormalised_parallel", seed = 3, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_parallel")
    o$Value[o$ParameterName == "beta"]
  }
  expect_identical(chain(TRUE), chain(FALSE))
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
               "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn)", "}", "algorithm mcmc {",
               "  move b = exchange(norm_rw(beta, 0.5), parallel = true)", "  termination = iterations(10)", "}"), f)
  expect_error(ilike::compile(f), "nothing to share out", fixed = TRUE)
})

# (C287, dev/plan.md 3.22 U1b) the estimated gradient of log 1/Z(beta),
# gradient = estimate(points = N): its mean is the exact gradient of the log
# posterior, S(x) - E_beta[S(u)] (uniform prior), by enumeration (4 x 4 lattice)
test_that("the estimated gradient of an ising field's log density in beta is unbiased", {
  g = lattice(4)
  e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
               "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn, gradient = estimate(points = 5))", "}"), f)
  r = ilike::compile(f)
  w = exp(0.6 * S_all - max(0.6 * S_all))
  v = sapply(1:500, function(s) ilike:::evaluate_gradient_log(r, parameters = list(beta = 0.6), variable = "beta", seed = s))
  expect_lt(abs(mean(v) - (stat(x) - sum(w * S_all) / sum(w))), 4 * sd(v) / sqrt(length(v)))
})

# (C287, dev/plan.md 3.22 U1b) the noisy methods of Alquier, Friel, Everitt &
# Boland (2016) and the exact MALA-exchange, on the C243 model: MALA-exchange
# against the exact posterior; the noisy exchange, noisy Langevin and noisy
# MALA-exchange (approximate) close to it
test_that("MALA-exchange and the noisy exchange, Langevin and MALA-exchange moves infer beta", {
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
  on.exit(unlink(c(f, "unnormalised_noisy"), recursive = TRUE), add = TRUE)
  moves = c(exact = "move b = exchange(langevin(beta, [[1]], 0.8), transform = auto)",
            noisy = "move b = exchange(norm_rw(beta, 1), points = 5, noisy = true, transform = auto)",
            noisy = "move b = unadjusted(langevin(beta, [[1]], 0.05), transform = auto)",
            noisy = "move b = exchange(langevin(beta, [[1]], 0.8), points = 5, noisy = true, transform = auto)")
  for (k in seq_along(moves))
  {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn, gradient = estimate(points = 5))", "}",
                 "algorithm mcmc {", paste0("  ", moves[[k]]), "  termination = iterations(20000)", "}"), f)
    unlink("unnormalised_noisy", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_noisy", seed = 1, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_noisy")
    b = o$Value[o$ParameterName == "beta"][-(1:1000)]
    expect_equal(mean(b), exact_mean, tolerance = 0.05)
    expect_equal(sd(b), exact_sd, tolerance = 0.1)
  }
})

# (C287) why MALA-exchange draws its exchange field afresh: on a discrete
# analogue (the state is beta and one 2 x 2 field drawn at it, kept with the
# state as the fields behind the gradient are; the proposal for beta depends
# on that field), enumerating the kernel, an independent exchange draw leaves
# the extended target invariant, and reusing the kept proposed field as the
# exchange variable (as in Algorithm 9 of Alquier et al. 2016) does not
test_that("an exchange move whose proposal depends on kept fields needs an independent exchange draw", {
  states = as.matrix(expand.grid(rep(list(0:1), 4)))
  edges = rbind(c(1, 2), c(3, 4), c(1, 3), c(2, 4))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[edges[, 1]] * s[edges[, 2]]) })
  th = c(0.2, 0.5, 0.8); prior = c(0.5, 0.3, 0.2); obs = 1; lambda = 0.7
  lg = function(b, i) b * S[i]
  f = function(b) { w = exp(b * S); w / sum(w) }
  Z = sapply(th, function(b) sum(exp(b * S)))
  pi_ext = matrix(0, 3, 16)
  for (a in 1:3) pi_ext[a, ] = prior[a] * exp(lg(th[a], obs)) / Z[a] * f(th[a])
  p = as.vector(t(pi_ext / sum(pi_ext)))
  h = function(b, a, y) { others = setdiff(1:3, a); w = exp(lambda * (th[others] - th[a]) * S[y]); w[others == b] / sum(w) }
  kernel = function(reuse) {
    P = matrix(0, 48, 48)
    for (a in 1:3) for (y in 1:16) for (b in setdiff(1:3, a)) for (yp in 1:16) {
      base = prior[b] * exp(lg(th[b], obs)) / (prior[a] * exp(lg(th[a], obs))) * h(a, b, yp) / h(b, a, y)
      acc = if (reuse) min(1, base * exp(lg(th[a], yp) - lg(th[b], yp)))
            else sum(f(th[b]) * pmin(1, base * exp(lg(th[a], 1:16) - lg(th[b], 1:16))))
      P[(a - 1)*16 + y, (b - 1)*16 + yp] = P[(a - 1)*16 + y, (b - 1)*16 + yp] + h(b, a, y) * f(th[b])[yp] * acc
    }
    diag(P) = diag(P) + 1 - rowSums(P)
    P
  }
  expect_lt(max(abs(p %*% kernel(FALSE) - p)), 1e-12)
  expect_gt(max(abs(p %*% kernel(TRUE) - p)), 0.01)
})

test_that("gradient = estimate(...), noisy exchange and noisy Langevin are refused where they do not apply", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  model = function(statement, move) c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
                                      "  beta ~ unif(0, 1)", "  lb ~ norm(0, 1)", statement, "}", "algorithm mcmc {", move,
                                      "  termination = iterations(10)", "}")
  writeLines(model("  x ~ ising(beta, nb, nn)", "  move b = exchange(langevin(beta, [[1]], 0.5))"), f)
  expect_error(ilike::compile(f), "gradient = estimate(points = N)", fixed = TRUE)
  writeLines(model("  x ~ ising(beta, nb, nn)", "  move b = unadjusted(langevin(beta, [[1]], 0.05))"), f)
  expect_error(ilike::compile(f), "a parameter of the normalising constant", fixed = TRUE)
  writeLines(model("  x ~ ising(beta, nb, nn)", "  move b = exchange(norm_rw(beta, 0.2), noisy = true)"), f)
  expect_error(ilike::compile(f), "averages the draws of points = N", fixed = TRUE)
  writeLines(model("  x ~ ising(exp(lb), nb, nn, gradient = estimate(points = 5))", "  move b = exchange(norm_rw(lb, 0.2))"), f)
  expect_error(ilike::compile(f), "beta is a parameter itself", fixed = TRUE)
  writeLines(model("  x ~ ising(beta, nb, nn, gradient = estimate(points = 5))", "  move b = unadjusted(langevin(beta, [[1]], 0.05))"), f)
  expect_message(ilike::compile(f), "noisy Langevin algorithm", fixed = TRUE)
})

# (C289, dev/plan.md 3.22 U3b) sav(..., normalise = smc(...)): log Z(at) is
# estimated before the run (SMC on the field from its base measure) and
# removed, so that importance sampling (SAVIS/MAVIS, Everitt, Johansen, Rowing
# & Evdemon-Hogan 2017) and SMC report the evidence itself; against
# enumeration and quadrature on a 4 x 4 lattice
test_that("normalise = smc(...) makes importance sampling with the auxiliary variable methods estimate the evidence", {
  g = lattice(4)
  e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  lZ = function(b) { a = b * S_all; max(a) + log(sum(exp(a - max(a)))) - 16 * log(2) }
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  model = c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
            "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn, normaliser = sav(points = 5, at = 0.3, normalise = smc(particles = 1000)))", "}")
  writeLines(model, f)
  r = ilike::compile(f)
  # at beta = at the estimate of Z(at)/Z(beta) is exactly 1, so the log target is the rest less log Zhat(at)
  v = ilike::evaluate_log(r, parameters = list(beta = 0.3), seed = 1)
  expect_lt(abs(-(v - (-16 * log(2) + 0.3 * stat(x))) - lZ(0.3)), 0.05)
  skip_on_cran()
  grid = seq(0.00005, 0.99995, by = 0.0001)
  ll = sapply(grid, function(b) -16 * log(2) + b * stat(x) - lZ(b))
  log_evidence = max(ll) + log(mean(exp(ll - max(ll))))
  writeLines(c(model, "algorithm is {", "  importance_points = 5000", "}"), f)
  expect_lt(abs(suppressMessages(ilike::IS(ilike::compile(f), seed = 1)) - log_evidence), 0.1)
})

# (C290, dev/plan.md 3.22 U3b) replicated fields, x ~ ising(..., replicates = R)
replicated_fields = function()
{
  g = lattice(3); e = unnormalised_edges(g)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  Y = rbind(c(1, 1, 0, 1, 0, 0, 1, 0, 0), c(0, 0, 1, 0, 1, 1, 0, 1, 1), c(1, 1, 1, 1, 0, 0, 1, 1, 1))
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 9))), 1, stat)
  lZ = function(b) { a = b * S_all; max(a) + log(sum(exp(a - max(a)))) - 9 * log(2) }
  grid = seq(0.00005, 0.99995, by = 0.0001)
  ll = sapply(grid, function(b) 3 * (-9 * log(2)) + b * sum(apply(Y, 1, stat)) - 3 * lZ(b))
  p = exp(ll - max(ll)); p = p / sum(p)
  list(head = function(extra) c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste0("[", apply(Y, 1, paste, collapse = ", "), "]", collapse = ", "), "] }"),
                                "model {", "  beta ~ unif(0, 1)", paste0("  x ~ ising(beta, nb, nn, replicates = 3", extra, ")"), "}"),
       stat_total = sum(apply(Y, 1, stat)), mean = sum(p * grid), sd = sqrt(sum(p * grid^2) - sum(p * grid)^2),
       log_evidence = max(ll) + log(mean(exp(ll - max(ll)))))
}

test_that("replicated ising fields: one potential term per replicate (the base measure and the potentials summed)", {
  z = replicated_fields()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(z$head(""), f)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = list(beta = 0.4)), 27 * log(0.5) + 0.4 * z$stat_total, tolerance = 1e-12)
})

test_that("replicated ising fields: beta by exchange MCMC, and the evidence by SMC with exchange moves and random weights", {
  skip_on_cran()
  z = replicated_fields()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_replicates"), recursive = TRUE), add = TRUE)
  writeLines(c(z$head(""), "algorithm mcmc {", "  move b = exchange(norm_rw(beta, 1), transform = auto)", "  termination = iterations(40000)", "}"), f)
  ilike::MCMC(ilike::compile(f), results_name = "unnormalised_replicates", seed = 1, initial_values = list(list(beta = 0.3)))
  o = ilike::load_mcmc_output("unnormalised_replicates")
  b = o$Value[o$ParameterName == "beta"][-(1:1000)]
  expect_equal(mean(b), z$mean, tolerance = 0.03)
  expect_equal(sd(b), z$sd, tolerance = 0.05)
  # Everitt, Johansen, Rowing & Evdemon-Hogan (2017), Algorithm 1: data tempering
  # adds one replicate per step, weighted by its potential times a fresh SAV
  # estimate of 1/Z (normalised by log Z(at)), with exchange moves
  writeLines(c(z$head(", normaliser = sav(points = 5, at = 0.5, normalise = smc(particles = 1000))"), "algorithm smc {",
               "  move b = exchange(norm_rw(beta, 0.5), transform = auto)", "  termination = iterations(5)",
               "  sequence = data(fraction, [0, 0.3333333333, 0.6666666667, 1], order = given, anneal = false)", "  particles = 2000", "}"), f)
  unlink("unnormalised_replicates", recursive = TRUE)
  log_z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), results_name = "unnormalised_replicates", seed = 13))
  expect_lt(abs(log_z - z$log_evidence), 0.15)
  o = ilike::load_smc_output("unnormalised_replicates")
  last = o[(o$Target == max(o$Target)) & (o$ParameterName == "beta"), ]
  last = last[last$Iteration == max(last$Iteration), ]
  w = exp(last$LogWeight - max(last$LogWeight)); w = w / sum(w)
  expect_equal(sum(w * last$Value), z$mean, tolerance = 0.06)
})

test_that("replicated fields: what does not apply is refused", {
  z = replicated_fields()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(sub("replicates = 3", "replicates = 2.5", z$head(""), fixed = TRUE), f)
  expect_error(ilike::compile(f), "replicates is a whole number", fixed = TRUE)
  writeLines(c(z$head(", normaliser = sav(at = 0.5)"), "algorithm mcmc {", "  move b = exchange(norm_rw(beta, 1))",
               "  move c = metropolis(norm_rw(beta, 1))", "  termination = iterations(10)", "}"), f)
  expect_error(ilike::compile(f), "use exchange moves only", fixed = TRUE)
  writeLines(c(z$head(", normaliser = sav(at = 0.5)"), "algorithm smc {", "  move b = exchange(norm_rw(beta, 0.5), transform = auto)",
               "  termination = iterations(2)", "  sequence = data(fraction, [0, 0.5, 1], order = given, anneal = true)", "  particles = 50", "}"), f)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1)), "the replicates are added whole and in the given order", fixed = TRUE)
})

# (C291, dev/plan.md 3.22 U3b) likelihood = abc(...) | sl(...): an observed
# field's likelihood by ABC or synthetic likelihood with its own simulator
# and sufficient statistic. With a uniform kernel and a tolerance below the
# statistic's spacing, ABC matches the statistic exactly, so (the statistic
# being sufficient) it targets the exact posterior
test_that("an observed ising field's beta by ABC (exact match on the sufficient statistic) and by synthetic likelihood", {
  skip_on_cran()
  g = lattice(3)
  e = unnormalised_edges(g)
  x = c(1, 1, 0, 1, 1, 0, 0, 0, 1)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 9))), 1, stat)
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(b) { a = b * S_all; b * stat(x) - (max(a) + log(sum(exp(a - max(a))))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  exact_mean = sum(p * grid); exact_sd = sqrt(sum(p * grid^2) - exact_mean^2)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_abc"), recursive = TRUE), add = TRUE)
  for (likelihood in c("abc(points = 20, tolerance = 0.5)", "sl(points = 20)"))
  {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", paste0("  x ~ ising(beta, nb, nn, likelihood = ", likelihood, ")"), "}",
                 "algorithm mcmc {", "  move b = metropolis(norm_rw(beta, 1), transform = auto)", "  termination = iterations(20000)", "}"), f)
    unlink("unnormalised_abc", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "unnormalised_abc", seed = 1, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_abc")
    b = o$Value[o$ParameterName == "beta"][-(1:1000)]
    # (synthetic likelihood assumes the statistic Gaussian: approximate)
    tolerance = if (startsWith(likelihood, "abc")) 0.05 else 0.1
    expect_equal(mean(b), exact_mean, tolerance = tolerance)
  }
})

test_that("likelihood = abc/sl is refused where it does not apply", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  field = function(statement) c("ilike 1.0", paste0("data { ", g$data, "; y = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
                                "  beta ~ unif(0, 1)", statement, "}")
  writeLines(field("  x ~ ising(beta, nb, nn, likelihood = sl(points = 10))"), f)
  expect_error(ilike::compile(f), "is an observed field (a data variable)", fixed = TRUE)
  writeLines(field("  y ~ ising(beta, nb, nn, likelihood = sl(points = 10), normaliser = sav(at = 0.3))"), f)
  expect_error(ilike::compile(f), "replaces the field's likelihood by simulation", fixed = TRUE)
  writeLines(field("  y ~ ising(beta, nb, nn, likelihood = abc(points = 10))"), f)
  expect_error(ilike::compile(f), "the tolerance is missing", fixed = TRUE)
})

# (C292, dev/plan.md 3.22 U3) annealed exchange (Murray, Ghahramani & MacKay
# 2006, section 5) is exact only with reversible updates at the bridges: on an
# asymmetric 4-site graph, enumerating the kernel with one bridge, a sweep
# then a reverse sweep leaves the posterior invariant, and a sweep in one
# order alone does not (on a lattice symmetric under reversing the order the
# two would agree, which is why the graph is asymmetric)
test_that("annealed exchange is exact with a sweep then a reverse sweep at the bridge, and not with one sweep", {
  states = as.matrix(expand.grid(rep(list(0:1), 4)))
  edges = rbind(c(1, 2), c(2, 3), c(3, 4), c(1, 3))
  nbs = lapply(1:4, function(i) c(edges[edges[, 1] == i, 2], edges[edges[, 2] == i, 1]))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[edges[, 1]] * s[edges[, 2]]) })
  idx = function(z) sum(z * 2^(0:3)) + 1
  th = c(0.2, 0.5, 0.8); prior = c(0.5, 0.3, 0.2); obs = 1
  f = function(b) { w = exp(b * S); w / sum(w) }
  post = prior * sapply(th, function(b) exp(b * S[obs]) / sum(exp(b * S))); post = post / sum(post)
  site_kernel = function(b, i) {
    K = matrix(0, 16, 16)
    for (a in 1:16) {
      z = states[a, ]; p1 = 1 / (1 + exp(-2 * b * sum(2 * z[nbs[[i]]] - 1)))
      z1 = z; z1[i] = 1; z0 = z; z0[i] = 0
      K[a, idx(z1)] = K[a, idx(z1)] + p1; K[a, idx(z0)] = K[a, idx(z0)] + 1 - p1
    }
    K
  }
  forward = function(b) Reduce(`%*%`, lapply(1:4, function(i) site_kernel(b, i)))
  reversible = function(b) forward(b) %*% Reduce(`%*%`, lapply(4:1, function(i) site_kernel(b, i)))
  kernel = function(sweep) {
    P = matrix(0, 3, 3)
    for (a in 1:3) for (b in setdiff(1:3, a)) {
      t = th[a]; v = th[b]; m = (t + v) / 2; T = sweep(m)
      lbase = log(prior[b] / prior[a]) + (v - t) * S[obs]
      lw = outer((m - v) * S, (t - m) * S, `+`)
      P[a, b] = 0.5 * sum(outer(f(v), rep(1, 16)) * T * pmin(1, exp(lbase + lw)))
    }
    diag(P) = 1 - rowSums(P)
    P
  }
  expect_lt(max(abs(post %*% kernel(reversible) - post)), 1e-12)
  expect_gt(max(abs(post %*% kernel(forward) - post)), 1e-4)
})

# (C295, dev/plan.md 3.22 U4) unnormalised distributions written in C++:
# x ~ unnormalised(log = cpp::potential, given = [...], simulate = exact(cpp::draw) | sweeps(cpp::sweep, m), ...)
user_ising_cpp = function(n) c("cpp {",
  paste0("double pot(const arma::mat &field, const Parameters &parameters, const Data &data) { const arma::mat &nb = data[\"nb\"]; const arma::mat &nn = data[\"nn\"]; const double beta = (parameters[\"beta\"])(0); double total = 0.0; ",
         "for (size_t i = 0; i < ", n, "; ++i) for (size_t k = 0; k < (size_t)nn(i); ++k) total += (2.0 * field(i) - 1.0) * (2.0 * field((size_t)nb(i, k) - 1) - 1.0); return 0.5 * beta * total; }"),
  paste0("void sw(RandomNumberGenerator &rng, Parameters &parameters, const Data &data) { const arma::mat &nb = data[\"nb\"]; const arma::mat &nn = data[\"nn\"]; const double beta = (parameters[\"beta\"])(0); const size_t n = ", n, "; ",
         "arma::mat &x = parameters[\"x\"]; if (x.n_elem != n) { x = arma::mat(n, 1); for (size_t i = 0; i < n; ++i) x(i) = (runif(rng) < 0.5) ? 1.0 : 0.0; } ",
         "for (size_t i = 0; i < n; ++i) { double s = 0.0; for (size_t k = 0; k < (size_t)nn(i); ++k) s += 2.0 * x((size_t)nb(i, k) - 1) - 1.0; x(i) = (runif(rng) < 1.0 / (1.0 + std::exp(-2.0 * beta * s))) ? 1.0 : 0.0; } }"),
  "}")

test_that("an ising field written in C++ gives the built-in's exchange chain (the same sweep and potential)", {
  g = lattice(4)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_user"), recursive = TRUE), add = TRUE)
  chain = function(statement, move, extra = NULL) {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", statement, "}", extra, "algorithm mcmc {", move, "  termination = iterations(300)", "}"), f)
    unlink("unnormalised_user", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_user", seed = 2, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_user")
    o$Value[o$ParameterName == "beta"]
  }
  expect_identical(chain("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 30))",
                         "  move b = exchange(norm_rw(beta, 1), transform = auto)", user_ising_cpp(16)),
                   chain("  x ~ ising(beta, nb, nn)", "  move b = exchange(norm_rw(beta, 1), simulate = sweeps(30), transform = auto)"))
})

# the exponential toy of the original exchange plan: y_i >= 0, gamma(y | theta)
# = exp(-theta sum y), exact draws y_i ~ Exp(theta) (Z(theta) = theta^-n is not
# given), theta ~ gamma(2, 1), so the posterior is gamma(2 + n, 1 + sum y)
test_that("a continuous distribution written in C++: exchange and SAV against the conjugate posterior", {
  skip_on_cran()
  set.seed(3); y = round(rexp(20, 1.5), 3)
  a = 2 + 20; b = 1 + sum(y)
  cpp = c("cpp {",
          "double potential(const arma::mat &field, const Parameters &parameters, const Data &data) { return -parameters[\"theta\"][0] * arma::accu(field); }",
          "Parameters draw(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data) { double t = parameters[\"theta\"][0]; arma::mat u(20, 1); if (!(t > 0.0)) { u.fill(arma::datum::nan); return Parameters(\"y\", u); } for (int i = 0; i < 20; ++i) u(i) = rexp(rng, t); return Parameters(\"y\", u); }",
          "}")
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_toy"), recursive = TRUE), add = TRUE)
  for (setup in list(c("", "move m = exchange(norm_rw(theta, 0.5), transform = auto)"),
                     c(", normaliser = sav(points = 5, at = 1.4)", "move m = metropolis(norm_rw(theta, 0.5), transform = auto)")))
  {
    writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  theta ~ gamma(2, 1)",
                 paste0("  y ~ unnormalised(log = cpp::potential, given = [theta], simulate = exact(cpp::draw)", setup[1], ")"), "}", cpp,
                 "algorithm mcmc {", paste0("  ", setup[2]), "  termination = iterations(30000)", "}"), f)
    unlink("unnormalised_toy", recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = "unnormalised_toy", seed = 1, initial_values = list(list(theta = 1)))
    o = ilike::load_mcmc_output("unnormalised_toy")
    th = o$Value[o$ParameterName == "theta"][-(1:1000)]
    expect_equal(mean(th), a / b, tolerance = 0.03)
    expect_equal(sd(th), sqrt(a) / b, tolerance = 0.06)
  }
})

test_that("unnormalised(...) written in C++ refuses what does not apply", {
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  model = function(statement, move = NULL) c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {",
                                             "  beta ~ unif(0, 1)", statement, "}", user_ising_cpp(9),
                                             if (!is.null(move)) c("algorithm mcmc {", move, "  termination = iterations(10)", "}"))
  writeLines(model("  z ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 10))"), f)
  expect_error(ilike::compile(f), "needs its base measure", fixed = TRUE)
  writeLines(model("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 10), gradient = estimate(points = 2))"), f)
  expect_error(ilike::compile(f), "give statistic = cpp::f", fixed = TRUE)
  writeLines(model("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 10))", "  move b = exchange(norm_rw(beta, 0.2), bridges = 2)"), f)
  expect_error(ilike::compile(f), "reverse = cpp::reverse_sweep", fixed = TRUE)
  writeLines(model("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 10))", "  move b = exchange(norm_rw(beta, 0.2), simulate = exact)"), f)
  expect_error(ilike::compile(f), "has no exact draw", fixed = TRUE)
  writeLines(model("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 10), likelihood = sl(points = 5))"), f)
  expect_error(ilike::compile(f), "not written yet", fixed = TRUE)
  writeLines(model("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 10))", "  move b = metropolis(norm_rw(beta, 0.2))"), f)
  expect_error(ilike::compile(f), "a parameter of the normalising constant", fixed = TRUE)
})

# (C296, dev/plan.md 3.22 U4 step 2) x ~ unnormalised(…) { … }: a distribution
# written as statements, its base measure and the terms of its log density,
# run by the general engine (the package's SMC on its inner model)
block_ising = function(extra = "", sweeps = 30)
{
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  lZ = function(b) { a = b * S_all; max(a) + log(sum(exp(a - max(a)))) - 16 * log(2) }
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(b) b * stat(x) - lZ(b))
  p = exp(lp - max(lp)); p = p / sum(p)
  list(model = c("ilike 1.0", paste0("data { x = [", paste(x, collapse = ", "), "]; from = [", paste(e[, 1], collapse = ", "),
                                     "]; to = [", paste(e[, 2], collapse = ", "), "]; E = ", nrow(e), " }"),
                 "model {", "  beta ~ unif(0, 1)",
                 paste0("  x ~ unnormalised(simulate = sweeps(", sweeps, ")", extra, ") {"),
                 "    for (i in 1:16) x[i] ~ bern(0.5)",
                 "    for (k in 1:E) log += beta * (2 * x[from[k]] - 1) * (2 * x[to[k]] - 1)",
                 "  }", "}"),
       stat = stat(x), lZ = lZ, mean = sum(p * grid), sd = sqrt(sum(p * grid^2) - sum(p * grid)^2))
}

test_that("an ising field written as a block: normalise through the general engine estimates log Z", {
  z = block_ising(", normaliser = sav(points = 5, at = 0.3, normalise = smc(particles = 1000))")
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(z$model, f)
  r = suppressMessages(ilike::compile(f))
  # at beta = at the path's estimate of Z(at)/Z(beta) is exactly 1
  v = ilike::evaluate_log(r, parameters = list(beta = 0.3), seed = 1)
  expect_lt(abs(-(v - (-16 * log(2) + 0.3 * z$stat)) - z$lZ(0.3)), 0.05)
})

test_that("an ising field written as a block: beta by exchange and by SAV through the general engine, against enumeration", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_block"), recursive = TRUE), add = TRUE)
  # exchange (one particle: no bridges); SAV with annealed runs resampled, as
  # waste-free SMC (sav(…, waste_free = 2), the general engine's)
  for (setup in list(c("", "move b = exchange(norm_rw(beta, 1), transform = auto)"),
                     c(", normaliser = sav(points = 8, at = 0.35, bridges = 3, resampling = ess(0.5), waste_free = 2)",
                       "move b = metropolis(norm_rw(beta, 0.5))")))
  {
    z = block_ising(setup[1])
    writeLines(c(z$model, "algorithm mcmc {", paste0("  ", setup[2]), "  termination = iterations(10000)", "}"), f)
    unlink("unnormalised_block", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_block", seed = 1, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_block")
    b = o$Value[o$ParameterName == "beta"][-(1:500)]
    expect_equal(mean(b), z$mean, tolerance = 0.05)
    expect_equal(sd(b), z$sd, tolerance = 0.1)
  }
})

test_that("a latent field written as a block: its states by gibbs moves and beta by exchange, against enumeration", {
  skip_on_cran()
  g = lattice(3); e = unnormalised_edges(g)
  y = c(1.2, 0.8, -0.3, 1.1, 0.9, -0.2, 1.3, -0.7, 0.1)
  states = as.matrix(expand.grid(rep(list(0:1), 9)))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) })
  ly = apply(states, 1, function(z) sum(dnorm(y, z, 1, log = TRUE)))
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(b) { a = b * S; lZ = max(a) + log(sum(exp(a - max(a)))); l = a - lZ + ly; max(l) + log(sum(exp(l - max(l)))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_block"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "]; from = [", paste(e[, 1], collapse = ", "),
                                   "]; to = [", paste(e[, 2], collapse = ", "), "]; E = ", nrow(e), " }"),
               "model {", "  beta ~ unif(0, 1)",
               "  x ~ unnormalised(simulate = sweeps(30)) {",
               "    for (i in 1:9) x[i] ~ bern(0.5)",
               "    for (k in 1:E) log += beta * (2 * x[from[k]] - 1) * (2 * x[to[k]] - 1)",
               "  }", "  for (i in 1:9) y[i] ~ norm(x[i], 1)", "}",
               "algorithm mcmc {", "  for (j in 1:9) move g = gibbs(x[j])", "  move b = exchange(norm_rw(beta, 0.5), transform = auto)",
               "  termination = iterations(60000)", "}"), f)
  # (beta and the states mix slowly: 60000 iterations give a Monte Carlo
  # error of about 0.009 in the mean; 200000 gave 0.4533 +/- 0.005)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_block", seed = 2,
              initial_values = list(list(beta = 0.3, x = matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 0), 9, 1))))
  o = ilike::load_mcmc_output("unnormalised_block")
  b = o$Value[o$ParameterName == "beta"][-(1:1000)]
  expect_equal(mean(b), sum(p * grid), tolerance = 0.06)
  expect_equal(sd(b), sqrt(sum(p * grid^2) - sum(p * grid)^2), tolerance = 0.1)
})

test_that("the base measure of a distribution written in C++ can be given as a block", {
  g = lattice(3)
  y = c(1.2, 0.8, -0.3, 1.1, 0.9, -0.2, 1.3, -0.7, 0.1)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_user"), recursive = TRUE), add = TRUE)
  chain = function(statement) {
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; y = [", paste(y, collapse = ", "), "] }"), "model {",
                 "  beta ~ unif(0, 1)", statement, "  for (i in 1:9) y[i] ~ norm(x[i], 1)", "}", user_ising_cpp(9),
                 "algorithm mcmc {", "  for (j in 1:9) move g = gibbs(x[j])", "  move b = exchange(norm_rw(beta, 0.5), transform = auto)",
                 "  termination = iterations(300)", "}"), f)
    unlink("unnormalised_user", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_user", seed = 2,
                initial_values = list(list(beta = 0.3, x = matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 0), 9, 1))))
    o = ilike::load_mcmc_output("unnormalised_user")
    o$Value[o$ParameterName == "beta"]
  }
  expect_identical(chain("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 30), base = bern(0.5), n = 9)"),
                   chain(c("  x ~ unnormalised(log = cpp::pot, given = [beta], simulate = sweeps(cpp::sw, 30)) {",
                           "    for (i in 1:9) x[i] ~ bern(0.5)", "  }")))
})

test_that("x ~ unnormalised(…) { … } refuses what does not apply", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  head = c("ilike 1.0", "data { x = [1, 0, 1]; from = [1, 2]; to = [2, 3]; E = 2 }", "model {", "  beta ~ unif(0, 1)", "  mu ~ norm(0, 1)")
  block = function(..., extra = "", move = NULL)
    c(head, paste0("  x ~ unnormalised(simulate = sweeps(10)", extra, ") {"), paste0("    ", c(...)), "  }", "}",
      if (!is.null(move)) c("algorithm mcmc {", move, "  termination = iterations(10)", "}"))
  term = "for (k in 1:E) log += beta * x[from[k]] * x[to[k]]"
  writeLines(block(term), f)
  expect_error(ilike::compile(f), "needs x's base measure", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ bern(0.5)"), f)
  expect_error(ilike::compile(f), "has no terms", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ bern(0.5)", "for (k in 1:E) log += x[from[k]] * x[to[k]]"), f)
  expect_error(ilike::compile(f), "read no parameters", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ bern(beta)", term), f)
  expect_error(ilike::compile(f), "the base measure is fixed", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ norm(0, 1)", term), f)
  expect_error(ilike::compile(f), "a discrete distribution per element", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ bern(0.5)", "mu ~ norm(0, 1)"), f)
  expect_error(ilike::compile(f), "holds x's base measure", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ bern(0.5)", term, move = "  move b = exchange(norm_rw(beta, 0.2), bridges = 2)"), f)
  expect_error(ilike::compile(f), "needs reversible moves at each bridge", fixed = TRUE)
  writeLines(block("for (i in 1:3) x[i] ~ bern(0.5)", term, move = "  move b = metropolis(norm_rw(beta, 0.2))"), f)
  expect_error(ilike::compile(f), "a parameter of the normalising constant", fixed = TRUE)
  g = lattice(3)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {", "  beta ~ unif(0, 1)",
               "  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, waste_free = 2))", "}"), f)
  expect_error(ilike::compile(f), "the fast engine has no waste-free SMC", fixed = TRUE)
})

# (C298, dev/plan.md 3.22 U4 step 2b) domain declarations, x[i] in {0, 1}:
# the terms that read one element make the base measure the independence
# model (normalised exactly over the domain), the others are the potential
declared_field = function(field_line, extra = "", h_data = NULL)
{
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
  c("ilike 1.0", paste0("data { x = [", paste(x, collapse = ", "), "]; from = [", paste(e[, 1], collapse = ", "),
                        "]; to = [", paste(e[, 2], collapse = ", "), "]; E = ", nrow(e), if (!is.null(h_data)) paste0("; h = ", h_data), " }"),
    "model {", "  beta ~ unif(0, 1)", if (is.null(h_data)) "  h ~ norm(0, 1)",
    paste0("  x ~ unnormalised(simulate = sweeps(30)", extra, ") {"), "    for (i in 1:16) x[i] in {0, 1}",
    "    for (k in 1:E) log += beta * (2 * x[from[k]] - 1) * (2 * x[to[k]] - 1)", field_line, "  }", "}")
}

test_that("with domain declarations and only terms that read one element, the distribution is the independence model, exactly", {
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { x = [", paste(x, collapse = ", "), "] }"), "model {", "  h ~ norm(0, 1)",
               "  x ~ unnormalised() {", "    for (i in 1:16) x[i] in {0, 1}", "    for (i in 1:16) log += h * x[i]", "  }", "}"), f)
  v = ilike::evaluate_log(suppressMessages(ilike::compile(f)), parameters = list(h = 0.7))
  expect_equal(v, dnorm(0.7, log = TRUE) + sum(0.7 * x - log(1 + exp(0.7))), tolerance = 1e-10)
})

test_that("domain declarations: normalise from the independence model estimates log Z with less variance than from the uniform base", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
  states = as.matrix(expand.grid(rep(list(0:1), 16)))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) })
  a = 0.3 * S + 3 * rowSums(states)
  exact = 0.3 * sum((2*x[e[, 1]] - 1) * (2*x[e[, 2]] - 1)) + 3 * sum(x) - (max(a) + log(sum(exp(a - max(a)))))
  spread = function(model) {
    writeLines(model, f)
    m = suppressMessages(ilike::compile(f))
    sapply(1:8, function(s) ilike::evaluate_log(m, parameters = list(beta = 0.3), seed = s)) - exact
  }
  normalise = ", normaliser = sav(points = 5, at = 0.3, normalise = smc(particles = 100, steps = 10))"
  declared = spread(declared_field("    for (i in 1:16) log += h * x[i]", normalise, h_data = 3))
  explicit_model = sub("    for (i in 1:16) x[i] in {0, 1}", "    for (i in 1:16) x[i] ~ bern(0.5)",
                       declared_field("    for (i in 1:16) log += h * x[i]", normalise, h_data = 3), fixed = TRUE)
  explicit = spread(explicit_model)
  expect_lt(abs(mean(declared)), 0.05)
  expect_lt(sd(declared), sd(explicit) / 2)
})

test_that("domain declarations: beta and the field h by exchange through the general engine, against enumeration", {
  skip_on_cran()
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
  states = as.matrix(expand.grid(rep(list(0:1), 16)))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }); M = rowSums(states)
  sx = sum((2*x[e[, 1]] - 1) * (2*x[e[, 2]] - 1)); mx = sum(x)
  bg = seq(0.005, 0.995, by = 0.01); hg = seq(-3, 3, by = 0.02)
  lp = outer(bg, hg, Vectorize(function(b, h) { a = b * S + h * M; b * sx + h * mx - (max(a) + log(sum(exp(a - max(a))))) + dnorm(h, log = TRUE) }))
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_declared"), recursive = TRUE), add = TRUE)
  writeLines(c(declared_field("    for (i in 1:16) log += h * x[i]"), "algorithm mcmc {",
               "  move b = exchange(norm_rw(beta, 0.5), transform = auto)", "  move m = exchange(norm_rw(h, 0.5))",
               "  termination = iterations(20000)", "}"), f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_declared", seed = 1, initial_values = list(list(beta = 0.3, h = 0)))
  o = ilike::load_mcmc_output("unnormalised_declared")
  b = o$Value[o$ParameterName == "beta"][-(1:1000)]; h = o$Value[o$ParameterName == "h"][-(1:1000)]
  expect_equal(mean(b), sum(p * bg), tolerance = 0.05)
  expect_lt(abs(mean(h) - sum(t(p) * hg)), 0.06)
})

test_that("domain declarations refuse what does not apply", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  model = function(...) c("ilike 1.0", "data { x = [1, 0, 1]; y = [0.2, 0.1, 0.3]; from = [1, 2]; to = [2, 3]; E = 2 }", "model {",
                          "  beta ~ unif(0, 1)", ..., "}")
  pair = "    for (k in 1:E) log += beta * x[from[k]] * x[to[k]]"
  writeLines(model("  z ~ unnormalised() {", "    for (i in 1:3) z[i] in {0, 1}", "    for (k in 1:E) log += beta * z[from[k]] * z[to[k]]", "  }",
                   "  for (i in 1:3) y[i] ~ norm(z[i], 1)"), f)
  expect_error(ilike::compile(f), "domain declarations are for data for now", fixed = TRUE)
  writeLines(model("  x ~ unnormalised() {", "    for (i in 1:3) x[i] in {0, 1}", "    x[1] ~ bern(0.5)", pair, "  }"), f)
  expect_error(ilike::compile(f), "not both, for now", fixed = TRUE)
  writeLines(model("  x ~ unnormalised() {", "    for (i in 1:3) x[i] in {0, 1}", "    log += beta", "  }"), f)
  expect_error(ilike::compile(f), "reads no element of x", fixed = TRUE)
  writeLines(model("  x ~ unnormalised() {", "    for (i in 1:2) x[i] in {0, 1}", "    x[3] in {0, 1, 2}", pair, "  }"), f)
  expect_error(ilike::compile(f), "have one domain for now", fixed = TRUE)
  writeLines(model("  x ~ unnormalised() {", "    for (i in 1:3) x[i] in {1}", pair, "  }"), f)
  expect_error(ilike::compile(f), "at least two values", fixed = TRUE)
  writeLines(model("  x ~ unnormalised(normaliser = sav(at = 0.3)) {", "    for (i in 1:3) x[i] in {0, 1}", "    for (i in 1:3) log += beta * x[i]", "  }"), f)
  expect_error(ilike::compile(f), "normaliser = … is not needed", fixed = TRUE)
})

# (C300, dev/plan.md 3.22 U9 R1) normaliser = estimate(smc(…)): log Z(beta)
# estimated by SMC from the base measure at each value, and plugged in
test_that("normaliser = estimate(smc(...)) plugs in an SMC estimate of log Z through either engine", {
  z = block_ising()
  g = lattice(4)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  est = "normaliser = estimate(smc(particles = 1000, steps = 20))"
  # the fast engine (the built-in field) and the general engine (a block):
  # the value at beta = 0.3 is log gamma(x | 0.3) - log Zhat(0.3)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
               "  beta ~ unif(0, 1)", paste0("  x ~ ising(beta, nb, nn, ", est, ")"), "}"), f)
  v = ilike::evaluate_log(suppressMessages(ilike::compile(f)), parameters = list(beta = 0.3), seed = 1)
  expect_lt(abs(-(v - (-16 * log(2) + 0.3 * z$stat)) - z$lZ(0.3)), 0.05)
  writeLines(block_ising(paste0(", ", est))$model, f)
  v = ilike::evaluate_log(suppressMessages(ilike::compile(f)), parameters = list(beta = 0.3), seed = 1)
  expect_lt(abs(-(v - (-16 * log(2) + 0.3 * z$stat)) - z$lZ(0.3)), 0.05)
  # refusals
  writeLines(block_ising(", normaliser = estimate(smc(particles = 10), debias = roulette)")$model, f)
  expect_error(ilike::compile(f), "debias = none", fixed = TRUE)
  writeLines(block_ising(", normaliser = estimate(is(points = 10))")$model, f)
  expect_error(ilike::compile(f), "normaliser = estimate(smc(particles = P", fixed = TRUE)
})

test_that("normaliser = estimate(smc(...)): the posterior of beta approaches the exact one as the SMC grows", {
  skip_on_cran()
  z = block_ising()
  g = lattice(4)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_plugin"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {",
               "  beta ~ unif(0, 1)", "  x ~ ising(beta, nb, nn, normaliser = estimate(smc(particles = 200, steps = 20)))", "}",
               "algorithm mcmc {", "  move b = metropolis(norm_rw(beta, 0.5))", "  termination = iterations(5000)", "}"), f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_plugin", seed = 1, initial_values = list(list(beta = 0.3)))
  o = ilike::load_mcmc_output("unnormalised_plugin")
  b = o$Value[o$ParameterName == "beta"][-(1:500)]
  expect_equal(mean(b), z$mean, tolerance = 0.05)
  expect_equal(sd(b), z$sd, tolerance = 0.1)
})

# (C301) sav(…) for normalising constants of several parameters, at = [beta = …, h = …]
test_that("sav(...) estimates Z(at)/Z(theta) for a normalising constant of several parameters", {
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
  states = as.matrix(expand.grid(rep(list(0:1), 16)))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }); M = rowSums(states)
  lZ = function(b, h) { a = b * S + h * M; max(a) + log(sum(exp(a - max(a)))) }
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  model = function(at) c("ilike 1.0", paste0("data { x = [", paste(x, collapse = ", "), "]; from = [", paste(e[, 1], collapse = ", "),
                                             "]; to = [", paste(e[, 2], collapse = ", "), "]; E = ", nrow(e), " }"),
                         "model {", "  beta ~ unif(0, 1)", "  h ~ norm(0, 1)",
                         paste0("  x ~ unnormalised(simulate = sweeps(100), normaliser = sav(at = ", at, ")) {"), "    for (i in 1:16) x[i] ~ bern(0.5)",
                         "    for (k in 1:E) log += beta * (2 * x[from[k]] - 1) * (2 * x[to[k]] - 1)", "    for (i in 1:16) log += h * x[i]", "  }", "}")
  writeLines(model("[h = 0.16, beta = 0.43]"), f)
  m = suppressMessages(ilike::compile(f))
  sx = sum((2*x[e[, 1]] - 1) * (2*x[e[, 2]] - 1))
  w = sapply(1:400, function(s) ilike::evaluate_log(m, parameters = list(beta = 0.6, h = 0.16), seed = s)) -
      (dnorm(0.16, log = TRUE) - 16 * log(2) + 0.6 * sx + 0.16 * sum(x))
  expect_lt(abs(log(mean(exp(w))) - (lZ(0.43, 0.16) - lZ(0.6, 0.16))), 0.2)
  writeLines(model("0.43"), f)
  expect_error(ilike::compile(f), "by name, at = [beta = …, h = …]", fixed = TRUE)
  writeLines(model("[beta = 0.43]"), f)
  expect_error(ilike::compile(f), "the value of h is missing", fixed = TRUE)
  writeLines(model("[beta = 0.43, h = 0.1, mu = 1]"), f)
  expect_error(ilike::compile(f), "mu is not a parameter of the normalising constant", fixed = TRUE)
})

# (C302) ising(…, field = h): an external field (spin convention), a number,
# an expression or one value per site, with exact draws by coupling from the past
test_that("an ising field with an external field: exact draws against enumeration", {
  g = lattice(3); e = unnormalised_edges(g)
  states = as.matrix(expand.grid(rep(list(0:1), 9)))
  hv = c(0.5, -0.2, 0.1, -0.6, 0.3, 0, 0.2, -0.4, 0.7)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  for (case in list(list(0.6, "-0.4", rep(-0.4, 9)), list(0.9, "hv", hv)))
  {
    lw = apply(states, 1, function(x) { s = 2*x - 1; case[[1]] * sum(s[e[, 1]] * s[e[, 2]]) + sum(case[[3]] * s) })
    w = exp(lw - max(lw)); w = w / sum(w)
    writeLines(c("ilike 1.0", paste0("data { ", g$data, "; hv = [", paste(hv, collapse = ", "), "] }"), "model {",
                 paste0("  x ~ ising(", case[[1]], ", nb, nn, field = ", case[[2]], ")"), "}"), f)
    d = ilike:::simulate_unnormalised(suppressMessages(ilike::compile(f)), "x", draws = 5000, method = "exact", seed = 2)
    # (absolute tolerances: some of these probabilities are small, and the
    # Monte Carlo standard error of 5000 draws is at most 0.007)
    expect_lt(abs(mean(d[, 1]) - sum(w * states[, 1])), 0.02)
    expect_lt(abs(mean(d[, 4]) - sum(w * states[, 4])), 0.02)
    expect_lt(abs(mean(d[, 1] == d[, 2]) - sum(w * (states[, 1] == states[, 2]))), 0.02)
  }
  # a per-site field of the wrong length, and what is not written yet
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; hv = [1, 2] }"), "model {", "  x ~ ising(0.5, nb, nn, field = hv)", "}"), f)
  expect_error(ilike:::simulate_unnormalised(suppressMessages(ilike::compile(f)), "x", draws = 1, method = "exact"), "hv has 2 values for 9 sites", fixed = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 1, 0, 1, 0, 0, 1, 0, 0] }"), "model {", "  beta ~ unif(0, 1)", "  h ~ norm(0, 1)",
               "  x ~ ising(beta, nb, nn, field = h, gradient = estimate(points = 2))", "}"), f)
  expect_error(ilike::compile(f), "with an external field it is not written yet", fixed = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [1, 2, 1, 3, 1, 2, 1, 1, 2] }"), "model {", "  x ~ potts(0.5, 3, nb, nn, field = 0.2)", "}"), f)
  expect_error(ilike::compile(f), "field", fixed = TRUE)
})

test_that("beta and an external field h by SAV with exact draws, against enumeration", {
  skip_on_cran()
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
  states = as.matrix(expand.grid(rep(list(0:1), 16)))
  S = apply(states, 1, function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }); M = rowSums(2*states - 1)
  sx = sum((2*x[e[, 1]] - 1) * (2*x[e[, 2]] - 1)); mx = sum(2*x - 1)
  bg = seq(0.005, 0.995, by = 0.01); hg = seq(-2, 2, by = 0.01)
  lp = outer(bg, hg, Vectorize(function(b, h) { a = b * S + h * M; b * sx + h * mx - (max(a) + log(sum(exp(a - max(a))))) + dnorm(h, log = TRUE) }))
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_field"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {", "  beta ~ unif(0, 1)", "  h ~ norm(0, 1)",
               "  x ~ ising(beta, nb, nn, field = h, normaliser = sav(at = [beta = 0.43, h = 0.08]))", "}",
               "algorithm mcmc {", "  move b = metropolis(norm_rw(beta, 0.3))", "  move m = metropolis(norm_rw(h, 0.3))", "  termination = iterations(20000)", "}"), f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_field", seed = 1, initial_values = list(list(beta = 0.3, h = 0)))
  o = ilike::load_mcmc_output("unnormalised_field")
  b = o$Value[o$ParameterName == "beta"][-(1:1000)]; h = o$Value[o$ParameterName == "h"][-(1:1000)]
  expect_equal(mean(b), sum(p * bg), tolerance = 0.06)
  expect_lt(abs(mean(h) - sum(t(p) * hg)), 0.04)
})

# (C303, dev/plan.md 3.22 U9 R2) debias = roulette(…): Lyne et al.'s unbiased,
# possibly negative, estimate of 1/Z; signed pseudo-marginal MCMC records each
# iteration's sign (the Sign column), and signed_mean weights by it
roulette_model = function(setting, prior = "unif(0, 1)", iterations = 300)
{
  g = lattice(4)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  c("ilike 1.0", paste0("data { ", g$data, "; x = [", paste(x, collapse = ", "), "] }"), "model {", paste0("  beta ~ ", prior),
    paste0("  x ~ ising(beta, nb, nn, normaliser = estimate(smc(", setting, "))"), "}",
    "algorithm mcmc {", "  move b = metropolis(norm_rw(beta, 0.5))", paste0("  termination = iterations(", iterations, ")"), "}")
}

test_that("signed pseudo-marginal MCMC: the Sign column, signed_mean, and the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_signed"), recursive = TRUE), add = TRUE)
  writeLines(roulette_model("particles = 20, steps = 5), debias = roulette(continue = 0.8, scale = 0.6)"), f)
  m = suppressMessages(ilike::compile(f))
  ilike::MCMC(m, results_name = "unnormalised_signed", seed = 1, initial_values = list(list(beta = 0.3)))
  o = ilike::load_mcmc_output("unnormalised_signed")
  expect_true("Sign" %in% names(o))
  expect_true(all(o$Sign %in% c(-1, 1)))
  s = ilike::signed_mean(o, "beta")
  b = o[o$ParameterName == "beta", ]
  expect_equal(s$estimate, sum(b$Sign * b$Value) / sum(b$Sign))
  expect_equal(s$mean_sign, mean(b$Sign))
  # a run without estimates that may be negative has no Sign column
  writeLines(roulette_model("particles = 20, steps = 5)"), f)
  unlink("unnormalised_signed", recursive = TRUE)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_signed", seed = 1, initial_values = list(list(beta = 0.3)))
  expect_false("Sign" %in% names(ilike::load_mcmc_output("unnormalised_signed")))
  # importance sampling and SMC refuse it
  writeLines(c(roulette_model("particles = 20, steps = 5), debias = roulette()")[1:6], "algorithm is {", "  importance_points = 10", "}"), f)
  # (C309) importance sampling weights by the absolute values and returns the signed evidence
  expect_false(is.null(attr(suppressMessages(ilike::IS(ilike::compile(f), seed = 1)), "mean_sign")))
  writeLines(roulette_model("particles = 20, steps = 5), debias = roulette(continue = 1.5)"), f)
  expect_error(ilike::compile(f), "between 0 and 1", fixed = TRUE)
})

test_that("Russian roulette estimates of 1/Z: the posterior of beta against enumeration, with and without negative estimates", {
  skip_on_cran()
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  exact_mean = function(top) {
    grid = seq(0.0005, top - 0.0005, by = 0.001)
    lp = sapply(grid, function(b) { a = b * S_all; b * stat(x) - (max(a) + log(sum(exp(a - max(a))))) })
    p = exp(lp - max(lp)); p = p / sum(p)
    sum(p * grid)
  }
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_roulette"), recursive = TRUE), add = TRUE)
  run = function(model) {
    writeLines(model, f)
    unlink("unnormalised_roulette", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_roulette", seed = 1, initial_values = list(list(beta = 0.3)))
    ilike::load_mcmc_output("unnormalised_roulette")
  }
  # k = 2: the terms 1 - c Zhat_i mostly in (0, 1), no negative estimates
  o = run(roulette_model("particles = 50, steps = 10), debias = roulette(continue = 0.75, scale = 2)", iterations = 5000))
  expect_equal(ilike::signed_mean(o, "beta", burn_in = 500)$estimate, exact_mean(1), tolerance = 0.05)
  # k = 0.65 (beta below 0.6, where the SMC's estimates are not too variable):
  # some estimates negative, and the signs correct the absolute value's bias
  o = run(roulette_model("particles = 100, steps = 10), debias = roulette(continue = 0.8, scale = 0.65)", prior = "unif(0, 0.6)", iterations = 10000))
  s = ilike::signed_mean(o, "beta", burn_in = 500)
  expect_gt(mean(o$Sign < 0), 0.01)
  expect_equal(s$estimate, exact_mean(0.6), tolerance = 0.05)
})

# (C304, dev/plan.md 3.22 U9 R3a) debias = block_poisson(…): the block-Poisson
# estimate of exp(-nu Z) with the auxiliary nu (2210.02734, Algorithm 1, made
# by any Metropolis move); correlated refresh keeps all blocks but one of the
# current estimate's random numbers
test_that("block-Poisson estimates: signs, correlated refresh, and the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_bp"), recursive = TRUE), add = TRUE)
  acceptance = function(setting) {
    writeLines(sub("0.5)", "0.05)", roulette_model(setting, iterations = 2000), fixed = TRUE), f)
    unlink("unnormalised_bp", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_bp", seed = 1, initial_values = list(list(beta = 0.3)))
    o = ilike::load_mcmc_output("unnormalised_bp")
    expect_true(all(o$Sign %in% c(-1, 1)))
    mean(diff(o$Value[o$ParameterName == "beta"]) != 0)
  }
  # a noisy estimate of Z and small steps: with correlated refresh successive
  # estimates are closer, so more proposals are accepted
  correlated = acceptance("particles = 5, steps = 3), debias = block_poisson()")
  fresh = acceptance("particles = 5, steps = 3), debias = block_poisson(correlated = false)")
  expect_gt(correlated, fresh + 0.04)
  writeLines(roulette_model("particles = 5, steps = 3), debias = block_poisson(blocks = 0)"), f)
  expect_error(ilike::compile(f), "whole number of at least 1", fixed = TRUE)
  writeLines(roulette_model("particles = 5, steps = 3), debias = block_poisson(correlated = 1)"), f)
  expect_error(ilike::compile(f), "correlated = true | false", fixed = TRUE)
  writeLines(c(roulette_model("particles = 5, steps = 3), debias = block_poisson()")[1:6], "algorithm is {", "  importance_points = 10", "}"), f)
  # (C309) importance sampling weights by the absolute values and returns the signed evidence
  expect_false(is.null(attr(suppressMessages(ilike::IS(ilike::compile(f), seed = 1)), "mean_sign")))
})

test_that("block-Poisson estimates with correlated refresh: the posterior of beta against enumeration", {
  skip_on_cran()
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(b) { a = b * S_all; b * stat(x) - (max(a) + log(sum(exp(a - max(a))))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_bp"), recursive = TRUE), add = TRUE)
  writeLines(sub("0.5)", "0.2)", roulette_model("particles = 20, steps = 5), debias = block_poisson()", iterations = 20000), fixed = TRUE), f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_bp", seed = 1, initial_values = list(list(beta = 0.3)))
  o = ilike::load_mcmc_output("unnormalised_bp")
  expect_equal(ilike::signed_mean(o, "beta", burn_in = 2000)$estimate, sum(p * grid), tolerance = 0.03)
})

# (C308, dev/plan.md 3.22 U9 R4) debias = taylor(…): the randomly truncated
# Taylor-series estimate of Z^-n (Chopin et al., arXiv 2403.20313; jointly for
# replicates, arXiv 2510.03587), tuned from pilot estimates at each value
test_that("Taylor-series estimates: signs, and the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_taylor"), recursive = TRUE), add = TRUE)
  writeLines(roulette_model("particles = 10, steps = 3), debias = taylor()", iterations = 100), f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_taylor", seed = 1, initial_values = list(list(beta = 0.3)))
  expect_true(all(ilike::load_mcmc_output("unnormalised_taylor")$Sign %in% c(-1, 1)))
  writeLines(roulette_model("particles = 10, steps = 3), debias = taylor(joint = true)"), f)
  expect_error(ilike::compile(f), "has no replicates", fixed = TRUE)
  writeLines(roulette_model("particles = 10, steps = 3), debias = taylor(tuning = scale(alpha = 3))"), f)
  expect_error(ilike::compile(f), "between 0 and 2", fixed = TRUE)
  writeLines(c(roulette_model("particles = 10, steps = 3), debias = taylor()")[1:6], "algorithm is {", "  importance_points = 10", "}"), f)
  # (C309) importance sampling weights by the absolute values and returns the signed evidence
  expect_false(is.null(attr(suppressMessages(ilike::IS(ilike::compile(f), seed = 1)), "mean_sign")))
})

test_that("Taylor-series estimates: the posterior of beta against enumeration, one field and three replicates", {
  skip_on_cran()
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(b) { a = b * S_all; b * stat(x) - (max(a) + log(sum(exp(a - max(a))))) })
  p = exp(lp - max(lp)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_taylor"), recursive = TRUE), add = TRUE)
  run = function(lines) {
    writeLines(lines, f)
    unlink("unnormalised_taylor", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "unnormalised_taylor", seed = 1, initial_values = list(list(beta = 0.3)))
    ilike::load_mcmc_output("unnormalised_taylor")
  }
  # the cycling estimate with the bootstrap rule (the default)
  o = run(sub("0.5)", "0.3)", roulette_model("particles = 20, steps = 5), debias = taylor()", iterations = 5000), fixed = TRUE))
  expect_equal(ilike::signed_mean(o, "beta", burn_in = 500)$estimate, sum(p * grid), tolerance = 0.03)
  # three replicates: one series for Z^-3, and (C308's fix) the plug-in of all three
  z = replicated_fields()
  mcmc = c("algorithm mcmc {", "  move b = metropolis(norm_rw(beta, 0.3))", "  termination = iterations(5000)", "}")
  o = run(c(z$head(", normaliser = estimate(smc(particles = 20, steps = 5), debias = taylor(joint = true))"), mcmc))
  expect_equal(ilike::signed_mean(o, "beta", burn_in = 500)$estimate, z$mean, tolerance = 0.03)
  o = run(c(z$head(", normaliser = estimate(smc(particles = 20, steps = 5))"), mcmc))
  expect_equal(mean(o$Value[o$ParameterName == "beta"][-(1:500)]), z$mean, tolerance = 0.03)
})

# (C309, dev/plan.md 3.22 U9 R5) importance sampling and SMC with factor
# estimates that may be negative: weights by the absolute values, the signs
# written with the points, and the evidence of the model with its sign
test_that("signed importance sampling and SMC: the posterior of beta and the evidence against enumeration", {
  skip_on_cran()
  g = lattice(4); e = unnormalised_edges(g)
  x = c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0)
  stat = function(z) { s = 2*z - 1; sum(s[e[, 1]] * s[e[, 2]]) }
  S_all = apply(as.matrix(expand.grid(rep(list(0:1), 16))), 1, stat)
  exact = function(top) {
    grid = seq(top/2000, top - top/2000, by = top/1000)
    lp = sapply(grid, function(b) { a = b * S_all; b * stat(x) - (max(a) + log(sum(exp(a - max(a))))) })
    p = exp(lp - max(lp)); p = p / sum(p)
    list(mean = sum(p * grid), log_evidence = max(lp) + log(mean(exp(lp - max(lp)))))
  }
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "unnormalised_signed_is"), recursive = TRUE), add = TRUE)
  run = function(debias, algorithm, prior = "unif(0, 1)", smc = FALSE) {
    writeLines(c(roulette_model(paste0("particles = 20, steps = 5), debias = ", debias), prior = prior)[1:6], algorithm), f)
    unlink("unnormalised_signed_is", recursive = TRUE)
    m = suppressMessages(ilike::compile(f))
    z = if (smc) suppressMessages(ilike::SMC_with_MCMC(m, results_name = "unnormalised_signed_is", seed = 1))
        else suppressMessages(ilike::IS(m, results_name = "unnormalised_signed_is", seed = 1))
    o = ilike::load_smc_output("unnormalised_signed_is")
    expect_true(all(o$Sign %in% c(-1, 1)))
    list(z = z, mean = ilike::signed_mean(o, "beta")$estimate)
  }
  is = c("algorithm is {", "  importance_points = 3000", "}")
  smc = c("algorithm smc {", "  move m = metropolis(norm_rw(beta, 0.2))", "  termination = iterations(5)",
          "  sequence = annealing(power, [0, 0.1, 0.3, 0.6, 1])", "  particles = 1000", "}")
  # block-Poisson over the whole prior (exact: mean 0.365, log evidence -10.008)
  ex = exact(1)
  r = run("block_poisson()", is)
  expect_equal(r$mean, ex$mean, tolerance = 0.03)
  expect_lt(abs(r$z - ex$log_evidence), 0.15)
  expect_equal(attr(r$z, "sign"), 1)
  r = run("block_poisson()", smc, smc = TRUE)
  expect_equal(r$mean, ex$mean, tolerance = 0.03)
  expect_lt(abs(r$z - ex$log_evidence), 0.15)
  # the Taylor series where the SMC's estimates of Z are precise enough (beta below 0.6)
  ex = exact(0.6)
  r = run("taylor()", is, prior = "unif(0, 0.6)")
  expect_equal(r$mean, ex$mean, tolerance = 0.03)
  expect_lt(abs(r$z - ex$log_evidence), 0.15)
})
