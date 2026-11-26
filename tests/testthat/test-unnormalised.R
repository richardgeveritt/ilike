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
  writeLines(field("  x ~ potts(beta, 2, nb, nn, normaliser = sav(at = 0.3, simulate = exact))"), f)
  expect_error(ilike::compile(f), "exact simulation is available for ising", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, bridge_sweeps = 2))"), f)
  expect_error(ilike::compile(f), "goes with bridges = K", fixed = TRUE)
  writeLines(field("  x ~ ising(beta, nb, nn, normaliser = sav(at = 0.3, bridges = 1.5))"), f)
  expect_error(ilike::compile(f), "bridges is a whole number of at least 0", fixed = TRUE)
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
  for (bridges in c(0, 10))
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
  # (C285, dev/plan.md 3.22 U2) and with the averaged acceptance ratio, points = 5
  for (points in c(1, 5))
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
  expect_error(ilike::compile(f), "exact simulation is available for ising", fixed = TRUE)
  writeLines(model("  x ~ potts(beta, 2, nb, nn)", "  move b = exchange(norm_rw(beta, 0.2), simulate = sweeps(20))"), f)
  expect_message(ilike::compile(f), "draws x by 20 heat-bath sweeps", fixed = TRUE)
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
