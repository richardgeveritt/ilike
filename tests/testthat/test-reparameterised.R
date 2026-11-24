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

# termination = adaptive_iterations(c, max) (dev/plan.md 3.21 (c3), C252):
# each SMC step's number of MCMC iterations from the acceptance rate p of a
# pilot iteration, ceil(log c / log(1 - p)) (South, Pettitt & Drovandi 2019)
test_that("adaptive_iterations keeps the SMC exact as the acceptance rate falls", {
  skip_on_cran()
  ybar = 0.9; eps = 0.05
  exact_log_z = log(pnorm((ybar + eps)/sqrt(10)) - pnorm((ybar - eps)/sqrt(10))) - log(2*eps)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  # (a poor random walk on (theta, u): its acceptance rate falls from about 0.3 to 0.01 as epsilon falls)
  writeLines(c("ilike 1.0", "data { ybar = [0.9] }", "model {", "  theta ~ norm(0, 3)", "  u ~ norm(theta, 1)",
               "  ybar ~ abc_uniform(u, epsilon)", "}", "algorithm smc {",
               "  move a = metropolis(mvnorm_rw([theta, u], [[1, 0], [0, 1]], 4))", "  termination = adaptive_iterations(0.2)",
               "  sequence = parameter(epsilon, [Inf, 0.05])", "  target = cess(0.9)", "  particles = 2000", "}"), f)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1))
  expect_lt(abs(z - exact_log_z), 0.1)
})

test_that("adaptive_iterations refuses what it cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  theta ~ norm(0, 1)", "}", "algorithm mcmc {",
               "  move a = metropolis(norm_rw(theta, 1))", "  termination = adaptive_iterations(0.2)", "}"), f)
  expect_error(suppressMessages(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(theta = 0)))),
               "it is for SMC moves", fixed = TRUE)
  writeLines(c("ilike 1.0", "model {", "  theta ~ norm(0, 1)", "}", "algorithm smc {",
               "  move a = metropolis(norm_rw(theta, 1))", "  termination = adaptive_iterations(1.5)",
               "  sequence = annealing(power, [0, 1])", "  particles = 20", "}"), f)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(f), seed = 1)), "c must be in (0, 1)", fixed = TRUE)
})

# move m = slice(x, width, lower, upper, max_steps) (dev/plan.md 3.21 (c3),
# C252): slice sampling of each element in turn (Neal 2003); with u on
# [0, 1]^m the interval is the whole of [0, 1] (Prangle et al.'s rare-event ABC)
uniform_inputs_model = function(extra = character(0), kernel = "abc_uniform(s, 0.5)")
{
  c("ilike 1.0", "data { ybar = [0.9] }", "model {", "  theta ~ norm(0, 3)", "  for (i in 1:5) u[i] ~ unif(0, 1)",
    "  s <- cpp::summary_of(u, theta)", paste0("  ybar ~ ", kernel), "}", "cpp {",
    "arma::mat summary_of(const arma::mat &u, const arma::mat &theta)",
    "{ arma::mat s(1, 1); double m = 0.0; for (arma::uword k = 0; k < u.n_elem; ++k) m += R::qnorm(u(k), 0.0, 1.0, 1, 0); s(0, 0) = theta(0, 0) + m/u.n_elem; return s; }",
    "}", extra)
}

test_that("slice sampling targets the posterior", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "slice_mcmc"), recursive = TRUE), add = TRUE)
  # a conjugate model, stepping out (posterior N(0.6, 0.5))
  writeLines(c("ilike 1.0", "data { y = [1.2] }", "model {", "  theta ~ norm(0, 1)", "  y ~ norm(theta, 1)", "}",
               "algorithm mcmc {", "  move a = slice(theta, width = 1)", "  termination = iterations(20000)", "}"), f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "slice_mcmc", seed = 1, initial_values = list(list(theta = 0))))
  b = ilike::load_mcmc_output("slice_mcmc")$Value[-(1:1000)]
  expect_equal(mean(b), 0.6, tolerance = 0.05)
  expect_equal(var(b), 0.5, tolerance = 0.05)
  # u on [0, 1]^5 (u = Phi(z) with z of C248's model, so the same ABC posterior), whole interval
  unlink("slice_mcmc", recursive = TRUE)
  writeLines(uniform_inputs_model(c("algorithm mcmc {", "  move a = metropolis(norm_rw(theta, 0.5))",
                                    "  move b = slice(u, lower = 0, upper = 1)", "  termination = iterations(20000)", "}")), f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "slice_mcmc", seed = 1, initial_values = list(list(theta = 0.9, u = rep(0.5, 5)))))
  o = ilike::load_mcmc_output("slice_mcmc")
  b = o$Value[o$ParameterName == "theta"][-(1:2000)]
  expect_equal(mean(b), 0.8725, tolerance = 0.05)
  expect_equal(sd(b), 0.5245, tolerance = 0.1)
})

test_that("rare-event ABC with slice moves on uniform inputs estimates the likelihood", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  m = uniform_inputs_model(kernel = "abc_uniform(s, epsilon)")
  writeLines(append(m, rare_event_setting(1000, "slice(u, lower = 0, upper = 1)", method = "rare_event_abc"), after = 4), f)
  theta = 0.7; eps = 0.5; ybar = 0.9; sd = 1/sqrt(5)
  exact = dnorm(theta, 0, 3, log = TRUE) + log((pnorm((ybar + eps - theta)/sd) - pnorm((ybar - eps - theta)/sd))/(2*eps))
  expect_equal(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = theta), seed = 1)),
               exact, tolerance = 0.02)
})

test_that("slice refuses what it cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  # an estimated target (u integrated out)
  writeLines(integrated_model(10, c("algorithm mcmc {", "  move a = slice(theta, width = 1)", "  termination = iterations(5)", "}")), f)
  expect_error(suppressMessages(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(theta = 0.9)))),
               "slice sampling needs the exact target", fixed = TRUE)
  # as a proposal
  writeLines(c("ilike 1.0", "model {", "  theta ~ norm(0, 1)", "}", "algorithm mcmc {",
               "  move a = metropolis_hastings(slice(theta, width = 1))", "  termination = iterations(5)", "}"), f)
  expect_error(ilike::compile(f), "slice is a move of its own", fixed = TRUE)
  # an unbounded interval needs a width
  writeLines(c("ilike 1.0", "model {", "  theta ~ norm(0, 1)", "}", "algorithm mcmc {",
               "  move a = slice(theta)", "  termination = iterations(5)", "}"), f)
  expect_error(suppressMessages(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(theta = 0)))),
               "a positive width is needed", fixed = TRUE)
})

# pcn_mala(u, mean, cov, delta) (dev/plan.md 3.21 (c3), C252): pCN with a step
# along the gradient of the other factors (Cotter, Roberts, Stuart & White
# 2013); abc_gaussian's gradient through an expression node
Z5 = "[0, 0, 0, 0, 0]"
gaussian_kernel_model = function(extra = character(0), integrate = character(0), kernel = "abc_gaussian(s, 0.5)")
  c("ilike 1.0", "data { ybar = [0.9] }", "model {", "  theta ~ norm(0, 3)", integrate, paste0("  u ~ mvnorm(", Z5, ", ", I5, ")"),
    "  s <- theta + sum(k in 1:5, u[k])/5", paste0("  ybar ~ ", kernel), "}", extra)

test_that("abc_gaussian's gradient and pcn_mala in rare-event ABC", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(gaussian_kernel_model(), f)
  g = ilike::check_gradient(ilike::compile(f), parameters = list(theta = 0.7, u = c(0.3, -0.2, 0.5, 0.1, -0.4)))
  expect_lt(max(abs(g$error)), 1e-6)
  # ybar | theta ~ N(theta, 1/5 + 0.25)
  writeLines(gaussian_kernel_model(integrate = rare_event_setting(1000, paste0("metropolis_hastings(pcn_mala(u, ", Z5, ", ", I5, ", 0.5))"),
                                                                 method = "rare_event_abc"), kernel = "abc_gaussian(s, epsilon)"), f)
  expect_equal(suppressMessages(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = 0.7), seed = 1)),
               dnorm(0.7, 0, 3, log = TRUE) + dnorm(0.9, 0.7, sqrt(0.45), log = TRUE), tolerance = 0.02)
})

test_that("joint MCMC with pcn_mala targets the posterior", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "pcn_mala_mcmc"), recursive = TRUE), add = TRUE)
  writeLines(gaussian_kernel_model(c("algorithm mcmc {", "  move a = metropolis(norm_rw(theta, 0.8))",
                                     paste0("  move b = metropolis_hastings(pcn_mala(u, ", Z5, ", ", I5, ", 0.5))"),
                                     "  termination = iterations(30000)", "}")), f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "pcn_mala_mcmc", seed = 1, initial_values = list(list(theta = 0.9, u = rep(0, 5)))))
  o = ilike::load_mcmc_output("pcn_mala_mcmc")
  b = o$Value[o$ParameterName == "theta"][-(1:2000)]
  # posterior N(0.9*9/9.45, 9*0.45/9.45)
  expect_equal(mean(b), 0.9*9/9.45, tolerance = 0.05)
  expect_equal(sd(b), sqrt(9*0.45/9.45), tolerance = 0.1)
})

test_that("pcn_mala refuses what it cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  # a C++ node has no gradient
  m = reparameterised_model("abc_gaussian(s, 0.5)", c("algorithm mcmc {",
        paste0("  move b = metropolis_hastings(pcn_mala(u, ", Z5, ", ", I5, ", 0.5))"), "  termination = iterations(5)", "}"))
  writeLines(m, f)
  expect_error(suppressMessages(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(theta = 0.9, u = rep(0, 5))))),
               "Gradients are not available")
  writeLines(gaussian_kernel_model(c("algorithm mcmc {", paste0("  move b = metropolis_hastings(pcn_mala(u, ", Z5, ", ", I5, ", 0))"),
                                     "  termination = iterations(5)", "}")), f)
  expect_error(suppressMessages(ilike::MCMC(ilike::compile(f), seed = 1, initial_values = list(list(theta = 0.9, u = rep(0, 5))))), "delta")
})

# The older ABC kernels of sbi(…) (dev/plan.md C253): the L^p distance with
# absolute values (and p = Inf), the unit L^p ball's volume for the kernel's
# p, the Gaussian kernel's scale as a standard deviation, and a scale given
# without a variable's name (it was ignored). With a simulator that has no
# randomness and one simulation, the estimate is the kernel itself.
test_that("the older ABC kernels give the right values", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  sim = c("cpp {", "Data both(RandomNumberGenerator &rng, const Parameters &parameters)",
          "{ arma::colvec s(2); s(0) = parameters[\"m\"][0] - 0.5; s(1) = parameters[\"m\"][0] + 0.5; return Data(\"y\", s); }", "}")
  value = function(method, m)
  {
    writeLines(c("ilike 1.0", "data { y = [1.0, 2.0] }", "model {", "  m ~ norm(0, 3)",
                 paste0("  factor lik = sbi(", method, ", simulate = cpp::both, data = y)"), "}", sim), f)
    ilike::evaluate_log(ilike::compile(f), parameters = list(m = m), seed = 1) - dnorm(m, 0, 3, log = TRUE)
  }
  y = c(1, 2)
  s = function(m) c(m - 0.5, m + 0.5)
  for (m in c(1, 2.2))
  {
    # p = 1: |y - s| = (0.5, 0.5) at m = 1, (0.7, 0.3) at m = 2.2 (signed: -0.7 and 0.3)
    expect_equal(value("lp_uniform_abc(1, epsilon, 1.5, 1)", m), -log(2*1.5^2), tolerance = 1e-10)
    expect_equal(value("lp_uniform_abc(1, epsilon, 1.5, Inf)", m), -log(3^2), tolerance = 1e-10)
    expect_equal(value("lp_uniform_abc(1, epsilon, 1.5, 1, scale = [2, 3])", m), -log(2*1.5^2*6), tolerance = 1e-10)
    expect_equal(value("gaussian_abc(1, epsilon, 0.5, scale = [2, 3])", m), sum(dnorm(y, s(m), 0.5*c(2, 3), log = TRUE)), tolerance = 1e-10)
  }
  expect_equal(value("lp_uniform_abc(1, epsilon, 1.5, 1)", 3), -Inf) # (|y - s| = (1.5, 0.5): outside)
})

# Early stopping in rare-event ABC-MCMC (dev/plan.md 3.21 (c4), C254): the
# generated factors' maxima (each term's density at its mode); a proposal's
# inner SMC stops when its estimate so far plus log max K(eps_T) - log max
# K(eps_s) is below the acceptance threshold (ABC kernels, tolerances not
# increasing); other inner factors run to the end
test_that("rare-event ABC-MCMC with early stopping targets the ABC posterior", {
  skip_on_cran()
  ybar = 0.9
  posterior = function(eps)
  {
    grid = seq(-8, 9, length.out = 20001)
    p = dnorm(grid, 0, 3) * (pnorm(ybar + eps - grid) - pnorm(ybar - eps - grid)); p = p / sum(p)
    c(sum(p * grid), sqrt(sum(p * grid^2) - sum(p * grid)^2))
  }
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "early_stopping_mcmc"), recursive = TRUE), add = TRUE)
  m = integrated_model(10, c("algorithm mcmc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(10000)", "}"))
  m[m == "  integrate u with is(points = 10)"] = rare_event_setting(50, "metropolis(norm_rw(u, 0.3))", sequence = "[Inf, 2, 1, 0.5, 0.25, 0.1]",
                                                                    method = "rare_event_abc")
  m[m == "  ybar ~ abc_uniform(u, 0.5)"] = "  ybar ~ abc_uniform(u, epsilon)"
  writeLines(m, f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "early_stopping_mcmc", seed = 1, initial_values = list(list(theta = 0.9))))
  b = ilike::load_mcmc_output("early_stopping_mcmc")$Value[-(1:1000)]
  exact = posterior(0.1)
  expect_equal(mean(b), exact[1], tolerance = 0.06)
  expect_equal(sd(b), exact[2], tolerance = 0.1)
  # an inner factor that is not an ABC kernel: no early stopping (ybar | theta ~ N(theta, 1 + 0.25))
  unlink("early_stopping_mcmc", recursive = TRUE)
  m[m == "  ybar ~ abc_uniform(u, epsilon)"] = "  ybar ~ norm(u, 0.5*epsilon)"
  m[grepl("integrate u with", m)] = rare_event_setting(50, "metropolis(norm_rw(u, 0.3))", sequence = "[4, 2, 1]", method = "rare_event_abc")
  writeLines(m, f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "early_stopping_mcmc", seed = 1, initial_values = list(list(theta = 0.9))))
  b = ilike::load_mcmc_output("early_stopping_mcmc")$Value[-(1:1000)]
  expect_equal(mean(b), 0.9*9/10.25, tolerance = 0.06)
  expect_equal(sd(b), sqrt(9*1.25/10.25), tolerance = 0.1)
})

# A term written in C++ with maximum = cpp::g and kernel = true (C254): its
# factor's maximum, and early stopping in rare-event ABC as for abc_uniform
user_kernel_model = function(extra, mcmc = character(0))
  c("ilike 1.0", "model {", "  theta ~ norm(0, 3)",
    rare_event_setting(50, "metropolis(norm_rw(u, 0.3))", sequence = "[Inf, 2, 1, 0.5, 0.25, 0.1]", method = "rare_event_abc"),
    "  u ~ norm(theta, 1)", paste0("  factor lik += cpp::kern(reads = [u, epsilon]", extra, ")"), "}", "cpp {",
    "double kern(const Parameters &parameters, const Data &data) { double e = parameters[\"epsilon\"][0]; if (std::isinf(e)) return 0.0;",
    "  return (std::abs(0.9 - parameters[\"u\"][0]) <= e) ? -std::log(2.0*e) : -std::numeric_limits<double>::infinity(); }",
    "double kmax(const Parameters &parameters, const Data &data) { double e = parameters[\"epsilon\"][0]; return std::isinf(e) ? 0.0 : -std::log(2.0*e); }",
    "}", mcmc)

test_that("a kernel written in C++ with its maximum stops early and targets the ABC posterior", {
  skip_on_cran()
  ybar = 0.9; eps = 0.1
  grid = seq(-8, 9, length.out = 20001)
  p = dnorm(grid, 0, 3) * (pnorm(ybar + eps - grid) - pnorm(ybar - eps - grid)); p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "user_kernel_mcmc"), recursive = TRUE), add = TRUE)
  writeLines(user_kernel_model(", maximum = cpp::kmax, kernel = true",
                               c("algorithm mcmc {", "  move a = metropolis(norm_rw(theta, 1.5))", "  termination = iterations(10000)", "}")), f)
  suppressMessages(ilike::MCMC(ilike::compile(f), results_name = "user_kernel_mcmc", seed = 1, initial_values = list(list(theta = 0.9))))
  b = ilike::load_mcmc_output("user_kernel_mcmc")$Value[-(1:1000)]
  expect_equal(mean(b), sum(p * grid), tolerance = 0.06)
  expect_equal(sd(b), sqrt(sum(p * grid^2) - sum(p * grid)^2), tolerance = 0.1)
})

test_that("a C++ term's maximum and kernel settings refuse what they cannot do", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(user_kernel_model(", kernel = true"), f)
  expect_error(ilike::compile(f), "kernel = true needs the term's maximum", fixed = TRUE)
  writeLines(user_kernel_model(", maximum = cpp::missing_function"), f)
  expect_error(ilike::compile(f), "cpp::missing_function is not defined", fixed = TRUE)
  writeLines(user_kernel_model(", maximum = cpp::kmax, kernel = 2.5"), f)
  expect_error(ilike::compile(f), NA) # (a number other than 0: true)
})
