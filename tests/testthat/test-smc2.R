# SMC^2 (dev/plan.md 3.21 (c2'), C256; Chopin, Jacob & Papaspiliopoulos): an
# outer SMC taking in one observation at a time (data tempering, whole
# observations) with a particle filter as a factor; each outer particle's
# filter is kept between the steps (weighted by the new observation at the
# reweighting, resampled and propagated at the move stage), and a PMMH
# proposal runs a new filter through the observations so far. With
# data(..., anneal = true) (3.21 (c2''), F1) a step may take in a share r of
# the next observation: the target is the geometric bridge between the two
# extended targets (the filter's estimate of the observation's ratio to the
# power r). With pf(inner, increase = exchange(...)) (F2) the filters get more
# particles when the moves' acceptance rate falls (the exchange step)

smc2_files = function(n, outer_particles, inner_particles, anneal = "false", extra = character(0), schedule = NULL,
                      pf = "pf(inner)", algorithm = NULL, inner_file = "inner.ilike")
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
  if (is.null(schedule))
    schedule = paste0("[", paste(format((0:n)/n, digits = 15), collapse = ", "), "]")
  if (is.null(algorithm))
    algorithm = c("algorithm smc {", "  move m = metropolis(norm_rw(a, 0.3))", "  termination = iterations(2)",
                  paste0("  sequence = data(fraction, ", schedule, ", factors = [filter], anneal = ", anneal, ")"),
                  paste0("  particles = ", outer_particles), extra, "}")
  writeLines(outer(inner_file, pf, algorithm), file.path(d, "smc2.ilike"))
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
  # each observation annealed in, with an adaptive target from [0, 1]: the
  # whole observations become points of the schedule (F1)
  d3 = smc2_files(20, 400, 100, anneal = "true", extra = "  target = cess(0.9)", schedule = "[0, 1]")
  on.exit(unlink(c(d3, "smc2_anneal_output"), recursive = TRUE), add = TRUE)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d3, "smc2.ilike")), results_name = "smc2_anneal_output", seed = 1))
  expect_lt(abs(z - exact_log_z), 0.2)
  o = ilike::load_smc_output("smc2_anneal_output")
  expect_gt(max(o$Target), 20) # (steps within observations)
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight))
  expect_equal(sum(w*o$Value)/sum(w), exact_mean, tolerance = 0.1)
  # the exchange step (F2): filters of 10 particles, doubled while the moves'
  # acceptance rate is below 0.6, up to 400
  d4 = smc2_files(20, 400, 10, pf = "pf(inner, increase = exchange(below = 0.6, max = 400))")
  on.exit(unlink(c(d4, "smc2_exchange_output"), recursive = TRUE), add = TRUE)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d4, "smc2.ilike")), results_name = "smc2_exchange_output", seed = 1))
  expect_lt(abs(z - exact_log_z), 0.2)
  x = ilike::load_smc_inner_particles("smc2_exchange_output")
  expect_equal(x$particles_before[1], 10)
  expect_equal(tail(x$particles_after, 1), 400)
  expect_true(all(x$acceptance_rate[x$particles_after > x$particles_before] < 0.6))
  o = ilike::load_smc_output("smc2_exchange_output")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight))
  expect_equal(sum(w*o$Value)/sum(w), exact_mean, tolerance = 0.1)
})

test_that("SMC^2's exchange step, never triggered, is the run without it", {
  d1 = smc2_files(6, 30, 10)
  d2 = smc2_files(6, 30, 10, pf = "pf(inner, increase = exchange(below = 0))")
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  z1 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d1, "smc2.ilike")), seed = 1))
  z2 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), seed = 1))
  expect_identical(z1, z2)
})

test_that("SMC^2's exchange step refuses what it cannot do", {
  d = smc2_files(5, 20, 10, pf = "pf(inner, increase = exchange(below = 0.2))",
                 algorithm = c("algorithm mcmc {", "  move m = metropolis(norm_rw(a, 0.3))", "  termination = iterations(5)", "}"))
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_error(suppressMessages(ilike::MCMC(ilike::compile(file.path(d, "smc2.ilike")), seed = 1)),
               "is for SMC^2", fixed = TRUE)
  d2 = smc2_files(5, 20, 10, pf = "pf(inner, increase = exchange(below = 0.2, multiply = 1))")
  on.exit(unlink(d2, recursive = TRUE), add = TRUE)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), seed = 1)),
               "must be more than 1", fixed = TRUE)
  d3 = smc2_files(5, 20, 10, pf = "pf(inner, increase = exchange(multiply = 3))")
  on.exit(unlink(d3, recursive = TRUE), add = TRUE)
  expect_error(ilike::compile(file.path(d3, "smc2.ilike")), "the acceptance rate below which the filters get more particles", fixed = TRUE)
})

test_that("SMC^2 with anneal = true on whole observations is the run without it", {
  d1 = smc2_files(6, 30, 10)
  d2 = smc2_files(6, 30, 10, anneal = "true")
  on.exit(unlink(c(d1, d2, "smc2_id1", "smc2_id2"), recursive = TRUE), add = TRUE)
  z1 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d1, "smc2.ilike")), results_name = "smc2_id1", seed = 1))
  z2 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), results_name = "smc2_id2", seed = 1))
  expect_identical(z1, z2)
  o1 = ilike::load_smc_output("smc2_id1"); o2 = ilike::load_smc_output("smc2_id2")
  expect_identical(o1$Value, o2$Value)
  expect_identical(o1$LogWeight, o2$LogWeight)
})

test_that("SMC^2 refuses a fixed schedule that passes a whole observation", {
  d = smc2_files(5, 20, 10, anneal = "true", schedule = "[0, 0.1, 1]")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_error(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), seed = 1)),
               "a fixed schedule must contain every whole observation", fixed = TRUE)
})

# (dev/plan.md 3.21 (c2'''), C265) SMC^2 with a Kalman filter: the increments
# are the exact predictive densities, so this is IBIS with an exact
# likelihood; with anneal = true, exact geometric bridges between posteriors
smc2_kf_exact = function(d)
{
  rk = ilike::compile(file.path(d, "kf.ilike"))
  grid = seq(-0.999, 0.999, length.out = 400)
  ll = sapply(grid, function(a) suppressMessages(ilike::evaluate_log(rk, parameters = list(a = a))))
  m = max(ll); w = exp(ll - m)
  c(log_z = m + log(sum(w)*(grid[2] - grid[1])), mean = sum(w*grid)/sum(w))
}

test_that("SMC^2 with a Kalman filter gives the exact evidence and posterior", {
  d = smc2_files(20, 400, 100, pf = "kf(inner)", inner_file = "inner_kf.ilike")
  on.exit(unlink(c(d, "smc2_kf_output"), recursive = TRUE), add = TRUE)
  exact = smc2_kf_exact(d)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), results_name = "smc2_kf_output", seed = 1))
  expect_lt(abs(z - exact[["log_z"]]), 0.15)
  o = ilike::load_smc_output("smc2_kf_output")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight))
  expect_equal(sum(w*o$Value)/sum(w), exact[["mean"]], tolerance = 0.05)
  # anneal = true on whole observations is the same run
  d2 = smc2_files(20, 400, 100, anneal = "true", pf = "kf(inner)", inner_file = "inner_kf.ilike")
  on.exit(unlink(d2, recursive = TRUE), add = TRUE)
  expect_identical(suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), seed = 1)), z)
})

test_that("SMC^2 with a Kalman filter anneals each observation in, with an adaptive target", {
  skip_on_cran()
  d = smc2_files(20, 1000, 100, anneal = "true", extra = "  target = cess(0.9)", schedule = "[0, 1]",
                 pf = "kf(inner)", inner_file = "inner_kf.ilike")
  on.exit(unlink(c(d, "smc2_kf_anneal_output"), recursive = TRUE), add = TRUE)
  exact = smc2_kf_exact(d)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), results_name = "smc2_kf_anneal_output", seed = 1))
  expect_lt(abs(z - exact[["log_z"]]), 0.15)
  o = ilike::load_smc_output("smc2_kf_anneal_output")
  expect_gt(max(o$Target), 20)
})

# (dev/plan.md 3.21 (c2'''), C266) SMC^2 with an ensemble Kalman filter:
# each outer particle keeps its ensemble between steps (moved on at the move
# stage, as a particle filter's particles are); the same approximate target
# as MCMC or annealed SMC with the same enkf(...) factor, close to the exact
# one on this linear-Gaussian model with 100 members
smc2_enkf_factor = "enkf(inner, ensemble_members = 100, shifter = \"stochastic\")"

test_that("SMC^2 with an ensemble Kalman filter gives the evidence and posterior", {
  skip_on_cran()
  d = smc2_files(20, 200, 100, pf = smc2_enkf_factor, inner_file = "inner_kf.ilike")
  on.exit(unlink(c(d, "smc2_enkf_output"), recursive = TRUE), add = TRUE)
  exact = smc2_kf_exact(d)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), results_name = "smc2_enkf_output", seed = 1))
  expect_lt(abs(z - exact[["log_z"]]), 0.25)
  o = ilike::load_smc_output("smc2_enkf_output")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight))
  expect_equal(sum(w*o$Value)/sum(w), exact[["mean"]], tolerance = 0.05)
})

test_that("SMC^2 with an ensemble Kalman filter and anneal = true on whole observations is the run without it", {
  factor = "enkf(inner, ensemble_members = 20, shifter = \"stochastic\")"
  d1 = smc2_files(6, 30, 10, pf = factor, inner_file = "inner_kf.ilike")
  d2 = smc2_files(6, 30, 10, anneal = "true", pf = factor, inner_file = "inner_kf.ilike")
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  z1 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d1, "smc2.ilike")), seed = 1))
  z2 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), seed = 1))
  expect_identical(z1, z2)
})

# (dev/plan.md 3.21 (c2'''), C267) the exchange step for the ensemble filter's
# size: ensembles of 10 members, doubled while the moves' acceptance rate is
# below 0.6, up to 200. The ensemble filter's estimate is biased, so each
# exchange is an importance-sampling step from the target with the old size
# to the one with the new size; the run ends on the 200-member target.
test_that("SMC^2 with an ensemble Kalman filter grows the ensembles by the exchange step", {
  skip_on_cran()
  factor = "enkf(inner, ensemble_members = 10, shifter = \"stochastic\", increase = exchange(below = 0.6, max = 200))"
  d = smc2_files(20, 400, 100, pf = factor, inner_file = "inner_kf.ilike")
  on.exit(unlink(c(d, "smc2_enkf_exchange_output"), recursive = TRUE), add = TRUE)
  exact = smc2_kf_exact(d)
  z = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d, "smc2.ilike")), results_name = "smc2_enkf_exchange_output", seed = 1))
  expect_lt(abs(z - exact[["log_z"]]), 0.3)
  x = ilike::load_smc_inner_particles("smc2_enkf_exchange_output")
  expect_equal(x$particles_before[1], 10)
  expect_equal(tail(x$particles_after, 1), 200)
  expect_true(all(x$acceptance_rate[x$particles_after > x$particles_before] < 0.6))
  o = ilike::load_smc_output("smc2_enkf_exchange_output")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight))
  expect_equal(sum(w*o$Value)/sum(w), exact[["mean"]], tolerance = 0.05)
})

test_that("SMC^2's exchange step for an ensemble Kalman filter, never triggered, is the run without it", {
  base = "enkf(inner, ensemble_members = 20, shifter = \"stochastic\""
  d1 = smc2_files(6, 30, 10, pf = paste0(base, ")"), inner_file = "inner_kf.ilike")
  d2 = smc2_files(6, 30, 10, pf = paste0(base, ", increase = exchange(below = 0))"), inner_file = "inner_kf.ilike")
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  z1 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d1, "smc2.ilike")), seed = 1))
  z2 = suppressMessages(ilike::SMC_with_MCMC(ilike::compile(file.path(d2, "smc2.ilike")), seed = 1))
  expect_identical(z1, z2)
})
