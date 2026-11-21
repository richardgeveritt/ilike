# Simulators written as x = H(u, theta) (dev/plan.md 3.21 (2), C248): C++
# deterministic nodes, x <- cpp::f(…), and ABC kernels as distributions of
# observed data, y ~ abc_uniform(s, epsilon) / abc_gaussian(s, epsilon)

reparameterised_model = function(kernel, extra = character(0))
{
  c("ilike 1.0", "data { ybar = [0.9] }", "model {", "  theta ~ norm(0, 3)",
    "  u ~ mvnorm([0, 0, 0, 0, 0], [[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])",
    "  s <- cpp::summary_of(u, theta)", paste0("  ybar ~ ", kernel), "}",
    "cpp {", "arma::mat summary_of(const arma::mat &u, const arma::mat &theta)",
    "{ arma::mat s(1, 1); s(0, 0) = theta(0, 0) + arma::mean(arma::vectorise(u)); return s; }", "}", extra)
}

test_that("ABC kernel distributions and C++ nodes give the right log density", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  u = c(0.3, -0.2, 0.5, 0.1, -0.4); theta = 0.7; s = theta + mean(u)
  prior = dnorm(theta, 0, 3, log = TRUE) + sum(dnorm(u, log = TRUE))
  P = list(theta = theta, u = u)
  # uniform, p = 2 (n = 1: 1/(2 epsilon) inside)
  writeLines(reparameterised_model("abc_uniform(s, 0.5)"), f)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = P), prior - log(2*0.5), tolerance = 1e-10)
  writeLines(reparameterised_model("abc_uniform(s, 0.1)"), f) # (|0.9 - 0.76| > 0.1: outside)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = P), -Inf)
  # Gaussian with a scale (a standard deviation)
  writeLines(reparameterised_model("abc_gaussian(s, 0.5, scale = 2)"), f)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = P), prior + dnorm(0.9, s, 0.5*2, log = TRUE), tolerance = 1e-10)
})

test_that("the uniform kernel's distance and normalising constant are right for p = 1", {
  # (the older ABC kernels' distance has no absolute value: wrong for odd p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [1.0, 2.0] }", "model {", "  m ~ norm(0, 3)", "  s <- cpp::both(m)",
               "  y ~ abc_uniform(s, 1.5, p = 1)", "}", "cpp {",
               "arma::mat both(const arma::mat &m) { arma::mat s(2, 1); s(0) = m(0) - 0.5; s(1) = m(0) + 0.5; return s; }", "}"), f)
  r = ilike::compile(f)
  # m = 1: s = (0.5, 1.5), |y - s| = (0.5, 0.5), L1 distance 1 <= 1.5: 1 / (area 2 epsilon^2)
  expect_equal(ilike::evaluate_log(r, parameters = list(m = 1)), dnorm(1, 0, 3, log = TRUE) - log(2*1.5^2), tolerance = 1e-10)
  # m = 2.2: |y - s| = (0.7, 0.3) (signed: -0.7 and 0.3, which would cancel): distance 1.0: inside
  expect_equal(ilike::evaluate_log(r, parameters = list(m = 2.2)), dnorm(2.2, 0, 3, log = TRUE) - log(2*1.5^2), tolerance = 1e-10)
  # m = 3: |y - s| = (1.5, 0.5): distance 2 > 1.5: outside
  expect_equal(ilike::evaluate_log(r, parameters = list(m = 3)), -Inf)
})

test_that("joint MCMC and joint SMC on (theta, u) target the ABC posterior", {
  skip_on_cran()
  eps = 0.5; sd = 1/sqrt(5); ybar = 0.9
  grid = seq(-5, 6, length.out = 20001)
  p = dnorm(grid, 0, 3) * (pnorm((ybar + eps - grid)/sd) - pnorm((ybar - eps - grid)/sd)); p = p / sum(p)
  exact_mean = sum(p * grid); exact_sd = sqrt(sum(p * grid^2) - exact_mean^2)
  exact_log_z = log(pnorm((ybar + eps)/sqrt(9.2)) - pnorm((ybar - eps)/sqrt(9.2))) - log(2*eps)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "reparameterised_mcmc", "reparameterised_smc"), recursive = TRUE), add = TRUE)
  moves = c("  move a = metropolis(norm_rw(theta, 0.5))",
            "  move b = metropolis(mvnorm_rw(u, [[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], 0.1))")
  writeLines(reparameterised_model("abc_uniform(s, 0.5)", c("algorithm mcmc {", moves, "  termination = iterations(100000)", "}")), f)
  ilike::MCMC(ilike::compile(f), results_name = "reparameterised_mcmc", seed = 1, initial_values = list(list(theta = 0.9, u = rep(0, 5))))
  o = ilike::load_mcmc_output("reparameterised_mcmc")
  b = o$Value[o$ParameterName == "theta"][-(1:5000)]
  expect_equal(mean(b), exact_mean, tolerance = 0.05)
  expect_equal(sd(b), exact_sd, tolerance = 0.1)
  writeLines(reparameterised_model("abc_uniform(s, epsilon)", c("algorithm smc {", moves, "  termination = iterations(10)",
                                   "  sequence = parameter(epsilon, [Inf, 4, 2, 1, 0.5])", "  particles = 2000", "}")), f)
  z = ilike::SMC_with_MCMC(ilike::compile(f), results_name = "reparameterised_smc", seed = 1)
  expect_lt(abs(z - exact_log_z), 0.1)
  o = ilike::load_smc_output("reparameterised_smc")
  m = smc_weighted_moments_by_target(o[o$ParameterName == "theta", ])
  expect_equal(unname(m[nrow(m), "mean"]), exact_mean, tolerance = 0.05)
})

test_that("C++ nodes and ABC kernels refuse what they cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  m = reparameterised_model("abc_uniform(s)")
  writeLines(m, f)
  expect_error(ilike::compile(f), "the tolerance epsilon are needed", fixed = TRUE)
  m = reparameterised_model("abc_uniform(s, 0.5)")
  m[m == "  s <- cpp::summary_of(u, theta)"] = "  for (i in 1:1) s[i] <- cpp::summary_of(u, theta)"
  writeLines(m, f)
  expect_error(ilike::compile(f), "is a whole variable: no index, no loop", fixed = TRUE)
  m = reparameterised_model("abc_uniform(s, 0.5)", c("algorithm mcmc {", "  move a = metropolis_hastings(langevin(theta, 0.1, 1))",
                                                     "  termination = iterations(10)", "}"))
  writeLines(m, f)
  expect_error(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(theta = 0.9, u = rep(0, 5)))), "Gradients are not available")
})

# integrate u with is(points = N) (dev/plan.md 3.21, step (b), C249): u's prior
# and the data that depend on u become an inner model, whose likelihood an
# importance sampler over u (proposal: u's prior, given theta) estimates

# theta ~ N(0, 3^2), u | theta ~ N(theta, 1), ybar ~ abc_uniform(u, eps): the
# integrated likelihood is P(|ybar - u| <= eps | theta) / (2 eps)
integrated_model = function(points, extra = character(0), method = "is")
{
  c("ilike 1.0", "data { ybar = [0.9] }", "model {", "  theta ~ norm(0, 3)",
    paste0("  integrate u with ", method, "(points = ", points, ")"), "  u ~ norm(theta, 1)",
    "  ybar ~ abc_uniform(u, 0.5)", "}", extra)
}

test_that("integrate u with is estimates the integrated likelihood", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  ybar = 0.9; eps = 0.5; theta = 0.7
  # (u's prior reads theta: it is simulated given the outer parameters)
  writeLines(integrated_model(100000), f)
  exact = dnorm(theta, 0, 3, log = TRUE) + log((pnorm(ybar + eps - theta) - pnorm(ybar - eps - theta))/(2*eps))
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = theta)), exact, tolerance = 0.01) # (Monte Carlo sd about 0.005)
  # u a vector read by a C++ node, with abc(…), the shorthand for is(…)
  m = reparameterised_model("abc_uniform(s, 0.5)")
  m = append(m, "  integrate u with abc(points = 100000)", after = 3)
  writeLines(m, f)
  r = ilike::compile(f)
  expect_equal(length(r$factor), 2) # (the integrated factor and theta's prior)
  sd = 1/sqrt(5)
  exact = dnorm(theta, 0, 3, log = TRUE) + log((pnorm((ybar + eps - theta)/sd) - pnorm((ybar - eps - theta)/sd))/(2*eps))
  expect_equal(ilike::evaluate_log(r, parameters = list(theta = theta)), exact, tolerance = 0.01) # (Monte Carlo sd about 0.005)
})

test_that("MCMC, annealed SMC and ABC-SMC on theta with u integrated out target the ABC posterior", {
  skip_on_cran()
  ybar = 0.9; eps = 0.5
  grid = seq(-8, 9, length.out = 20001)
  p = dnorm(grid, 0, 3) * (pnorm(ybar + eps - grid) - pnorm(ybar - eps - grid)); p = p / sum(p)
  exact_mean = sum(p * grid); exact_sd = sqrt(sum(p * grid^2) - exact_mean^2)
  exact_log_z = log(pnorm((ybar + eps)/sqrt(10)) - pnorm((ybar - eps)/sqrt(10))) - log(2*eps)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "integrated_mcmc"), recursive = TRUE), add = TRUE)
  writeLines(integrated_model(50, c("algorithm mcmc {", "  move a = metropolis(norm_rw(theta, 1.5))",
                                    "  termination = iterations(20000)", "}",
                                    "algorithm smc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(10)",
                                    "  sequence = annealing(power, [0, 0.01, 0.03, 0.1, 0.3, 1])", "  particles = 2000", "}")), f)
  r = ilike::compile(f)
  ilike::MCMC(r, results_name = "integrated_mcmc", seed = 1, initial_values = list(list(theta = 0.9)), algorithm = "mcmc")
  o = ilike::load_mcmc_output("integrated_mcmc")
  expect_setequal(unique(o$ParameterName), "theta")
  b = o$Value[o$ParameterName == "theta"][-(1:2000)]
  expect_equal(mean(b), exact_mean, tolerance = 0.05)
  expect_equal(sd(b), exact_sd, tolerance = 0.1)
  z = suppressMessages(ilike::SMC_with_MCMC(r, seed = 1, algorithm = "smc"))
  expect_lt(abs(z - exact_log_z), 0.1)
  # a decreasing tolerance (ABC-SMC, arXiv:2211.02172 alg. 3): each particle's
  # points are kept and weighted again at each tolerance (redrawing them, or
  # keeping the first estimate, would give a wrong evidence)
  m = integrated_model(50, c("algorithm smc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(10)",
                             "  sequence = parameter(epsilon, [Inf, 4, 2, 1, 0.5])", "  particles = 2000", "}"))
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(u, epsilon)"
  writeLines(m, f)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1))
  expect_lt(abs(z - exact_log_z), 0.1)
})

test_that("integrate refuses what it cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  m = integrated_model(10, c("algorithm mcmc {", "  move a = metropolis(norm_rw(u, 1))", "  termination = iterations(10)", "}"))
  writeLines(m, f)
  expect_error(ilike::compile(f), "changes u, which is integrated out", fixed = TRUE)
  m = integrated_model(10)
  m[m == "  theta ~ norm(0, 3)"] = "  theta ~ norm(u, 3)"
  writeLines(m, f)
  expect_error(ilike::compile(f), "the prior of theta depends on u, which is integrated out", fixed = TRUE)
  m = integrated_model(10)
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(theta, 0.5)"
  writeLines(m, f)
  expect_error(ilike::compile(f), "no data depend on u", fixed = TRUE)
  m = integrated_model(10)
  m[m == "  integrate u with is(points = 10)"] = "  integrate u with is(10)"
  writeLines(m, f)
  expect_error(ilike::compile(f), "its settings are named", fixed = TRUE)
  m = integrated_model(10)
  m[m == "  integrate u with is(points = 10)"] = "  integrate u with pf(particles = 10)"
  writeLines(m, f)
  expect_error(ilike::compile(f), "integrate u with is(points = N)", fixed = TRUE)
})

# integrate u with smc(…) / rare_event_abc(…) (dev/plan.md 3.21, (c1), C250):
# an SMC sampler over u for each theta, through decreasing tolerances, with
# MCMC moves on u (Prangle, Everitt & Kypraios's rare-event ABC)
I5 = "[[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]]"
rare_event_setting = function(particles, move, sequence = "[Inf, 2, 1, 0.5]", method = "smc", extra = "")
  paste0("  integrate u with ", method, "(particles = ", particles, ", move = ", move, ", sequence = parameter(epsilon, ",
         sequence, "), termination = iterations(3)", extra, ")")

test_that("integrate u with smc estimates the integrated likelihood", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  ybar = 0.9; eps = 0.5; theta = 0.7
  # u | theta ~ N(theta, 1): a random walk on u
  m = integrated_model(10)
  m[m == "  integrate u with is(points = 10)"] = rare_event_setting(2000, "metropolis(norm_rw(u, 0.5))")
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(u, epsilon)"
  writeLines(m, f)
  exact = dnorm(theta, 0, 3, log = TRUE) + log((pnorm(ybar + eps - theta) - pnorm(ybar - eps - theta))/(2*eps))
  expect_equal(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = theta), seed = 1)),
               exact, tolerance = 0.01)
  # u a 5-vector with a N(0, I) prior, moved by pCN (which leaves the prior
  # invariant); and an adaptive random walk with an adaptive tolerance schedule
  sd = 1/sqrt(5)
  exact = dnorm(theta, 0, 3, log = TRUE) + log((pnorm((ybar + eps - theta)/sd) - pnorm((ybar - eps - theta)/sd))/(2*eps))
  m = reparameterised_model("abc_uniform(s, epsilon)")
  writeLines(append(m, rare_event_setting(1000, paste0("metropolis_hastings(pcn(u, [0, 0, 0, 0, 0], ", I5, ", 0.8))"),
                                          method = "rare_event_abc"), after = 3), f)
  expect_equal(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = theta), seed = 1)),
               exact, tolerance = 0.02)
  writeLines(append(m, rare_event_setting(1000, paste0("metropolis(mvnorm_rw(u, ", I5, ", 0.1), smc_adapt = gaussian(2.38))"),
                                          sequence = "[Inf, 0.5]", extra = ", target = cess(0.9)"), after = 3), f)
  expect_equal(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = theta), seed = 1)),
               exact, tolerance = 0.02)
})

test_that("rare-event ABC-MCMC targets the ABC posterior", {
  skip_on_cran()
  ybar = 0.9; eps = 0.5
  grid = seq(-8, 9, length.out = 20001)
  p = dnorm(grid, 0, 3) * (pnorm(ybar + eps - grid) - pnorm(ybar - eps - grid)); p = p / sum(p)
  exact_mean = sum(p * grid); exact_sd = sqrt(sum(p * grid^2) - exact_mean^2)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rare_event_mcmc"), recursive = TRUE), add = TRUE)
  m = integrated_model(10, c("algorithm mcmc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(10000)", "}"))
  m[m == "  integrate u with is(points = 10)"] = rare_event_setting(50, "metropolis(norm_rw(u, 0.5))", method = "rare_event_abc")
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(u, epsilon)"
  writeLines(m, f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "rare_event_mcmc", seed = 1, initial_values = list(list(theta = 0.9))))
  o = ilike::load_mcmc_output("rare_event_mcmc")
  b = o$Value[o$ParameterName == "theta"][-(1:1000)]
  expect_equal(mean(b), exact_mean, tolerance = 0.06)
  expect_equal(sd(b), exact_sd, tolerance = 0.1)
})

# rare-event ABC-SMC^2 (dev/plan.md 3.21 (c2), C251; Kerama, Thorne & Everitt,
# arXiv:2211.02172): SMC on theta through decreasing tolerances, each particle
# with its own inner SMC on u that follows the outer tolerances, one step per
# outer step (reweight, resample, move u at the outer move stage)
test_that("rare-event ABC-SMC^2 targets the ABC posterior and evidence", {
  skip_on_cran()
  ybar = 0.9; eps = 0.5
  grid = seq(-8, 9, length.out = 20001)
  p = dnorm(grid, 0, 3) * (pnorm(ybar + eps - grid) - pnorm(ybar - eps - grid)); p = p / sum(p)
  exact_mean = sum(p * grid)
  exact_log_z = log(pnorm((ybar + eps)/sqrt(10)) - pnorm((ybar - eps)/sqrt(10))) - log(2*eps)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rare_event_smc2"), recursive = TRUE), add = TRUE)
  m = integrated_model(10, c("algorithm smc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(5)",
                             "  sequence = parameter(epsilon, [Inf, 4, 2, 1, 0.5])", "  particles = 1000", "}"))
  m[m == "  integrate u with is(points = 10)"] =
    "  integrate u with rare_event_abc(particles = 50, move = metropolis(norm_rw(u, 0.5)), termination = iterations(3))"
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(u, epsilon)"
  writeLines(m, f)
  r = ilike::compile(f)
  z = suppressMessages(ilike::SMC_with_MCMC(r, results_name = "rare_event_smc2", seed = 1))
  expect_lt(abs(z - exact_log_z), 0.1)
  o = ilike::load_smc_output("rare_event_smc2")
  mm = smc_weighted_moments_by_target(o[o$ParameterName == "theta", ])
  expect_equal(unname(mm[nrow(mm), "mean"]), exact_mean, tolerance = 0.06)
  # an adaptive outer schedule: the inner states are tried at each candidate tolerance
  m[m == "  sequence = parameter(epsilon, [Inf, 4, 2, 1, 0.5])"] = "  sequence = parameter(epsilon, [Inf, 0.5])"
  writeLines(append(m, "  target = cess(0.9)", after = which(m == "  particles = 1000")), f)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1))
  expect_lt(abs(z - exact_log_z), 0.1)
})

test_that("integrate with smc refuses what it cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  m = integrated_model(10)
  m[m == "  integrate u with is(points = 10)"] = "  integrate u with smc(particles = 10, move = metropolis(norm_rw(u, 0.5)))"
  writeLines(m, f)
  expect_error(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = 0.7), seed = 1)),
               "the SMC's sequence is needed", fixed = TRUE)
  # annealing outside: the factor is fixed along it, so the inner SMC runs its own sequence
  m = integrated_model(10, c("algorithm smc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(1)",
                             "  sequence = annealing(power, [0, 1])", "  particles = 20", "}"))
  m[m == "  integrate u with is(points = 10)"] = rare_event_setting(10, "metropolis(norm_rw(u, 0.5))", sequence = "[Inf, 1]")
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(u, epsilon)"
  writeLines(m, f)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1)), NA)
  # an outer sequence that changes the factor: the inner SMC follows it, so it has no sequence of its own
  m[m == "  sequence = annealing(power, [0, 1])"] = "  sequence = parameter(delta, [1, 0.5])"
  writeLines(m, f)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1)), "remove sequence (and target) from smc(...)", fixed = TRUE)
  # pCN's settings
  m = reparameterised_model("abc_uniform(s, epsilon)")
  writeLines(append(m, rare_event_setting(10, "metropolis_hastings(pcn(u, [0, 0, 0], [[1,0,0],[0,1,0],[0,0,1]], 0.8))"), after = 3), f)
  expect_error(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = 0.7), seed = 1)),
               "the variable has 5 element(s), the mean 3", fixed = TRUE)
  writeLines(append(m, rare_event_setting(10, paste0("metropolis_hastings(pcn(u, [0, 0, 0, 0, 0], ", I5, ", 1))")), after = 3), f)
  expect_error(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = 0.7), seed = 1)),
               "rho must be in [0, 1)", fixed = TRUE)
})
