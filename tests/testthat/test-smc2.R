# SMC^2 (dev/plan.md 3.21 (c2'), C256; Chopin, Jacob & Papaspiliopoulos): an
# outer SMC taking in one observation at a time (data tempering, whole
# observations) with a particle filter as a factor; each outer particle's
# filter is kept between the steps (weighted by the new observation at the
# reweighting, resampled and propagated at the move stage), and a PMMH
# proposal runs a new filter through the observations so far

smc2_files = function(n, outer_particles, inner_particles, anneal = "false", extra = character(0))
{
  set.seed(4); x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  d = tempfile(); dir.create(d)
  inner = function(meas, alg) c("ilike 1.0", "model {", "  x ~ norm(0, 1)",
                                "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)", meas,
                                paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
                                "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
                                "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }",
                                "double meas(const Parameters &parameters, const Data &data) { return dnorm(data[\"y\"][0], parameters[\"x\"][0], 1.0); }",
                                "}", alg)
  writeLines(inner("  factor obs = density(log = cpp::meas)", c("algorithm filter {", paste0("  particles = ", inner_particles), "}")),
             file.path(d, "inner.ilike"))
  writeLines(inner("  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)", c("algorithm kalman {", "}")),
             file.path(d, "inner_kf.ilike"))
  outer = function(import, factor, alg) c("ilike 1.0", paste0("import \"", import, "\" as inner"),
                                          paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  a ~ unif(-1, 1)",
                                          paste0("  factor filter = ", factor), "}", alg)
  schedule = paste0("[", paste(format((0:n)/n, digits = 15), collapse = ", "), "]")
  writeLines(outer("inner.ilike", "pf(inner)", c("algorithm smc {", "  move m = metropolis(norm_rw(a, 0.3))", "  termination = iterations(2)",
                                                 paste0("  sequence = data(fraction, ", schedule, ", factors = [filter], anneal = ", anneal, ")"),
                                                 paste0("  particles = ", outer_particles), extra, "}")), file.path(d, "smc2.ilike"))
  writeLines(outer("inner_kf.ilike", "kf(inner)", character(0)), file.path(d, "kf.ilike"))
  d
}

test_that("SMC^2 gives the evidence and the posterior of a state-space model", {
  skip_on_cran()
  d = smc2_files(20, 400, 100)
  on.exit(unlink(c(d, "smc2_output"), recursive = TRUE), add = TRUE)
  # exact: the Kalman filter's likelihood integrated over a ~ unif(-1, 1)
  rk = ilike::compile(file.path(d, "kf.ilike"))
  grid = seq(-0.999, 0.999, length.out = 400)
  ll = sapply(grid, function(a) suppressMessages(ilike::evaluate_log(rk, parameters = list(a = a))))
  m = max(ll)
  exact_log_z = m + log(sum(exp(ll - m))*(grid[2] - grid[1]))
  w = exp(ll - m); exact_mean = sum(w*grid)/sum(w)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), results_name = "smc2_output", seed = 1))
  expect_lt(abs(z - exact_log_z), 0.2)
  o = ilike::load_smc_output("smc2_output")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight))
  expect_equal(sum(w*o$Value)/sum(w), exact_mean, tolerance = 0.1)
  # moving only at the steps with resampling, as in the SMC^2 paper (C257)
  d2 = smc2_files(20, 400, 100, extra = c("  resampling = ess(0.5)", "  move_only_after_resampling = true"))
  on.exit(unlink(d2, recursive = TRUE), add = TRUE)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), seed = 1))
  expect_lt(abs(z - exact_log_z), 0.2)
})

test_that("SMC^2 refuses what it cannot do yet", {
  d = smc2_files(5, 20, 10, anneal = "true")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), seed = 1)),
               "takes in whole observations for now", fixed = TRUE)
})
