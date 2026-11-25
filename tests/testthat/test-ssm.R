test_that("KF gives the exact log-likelihood of a linear-Gaussian SSM", {

  on.exit(unlink("kf_linear_gaussian", recursive=TRUE), add=TRUE)

  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_kf_fixture("test_kf_linear_gaussian.ilike", "kf_linear_gaussian"), truth, tolerance = 1e-6)
})

test_that("KF explains that the prior must be a built-in Gaussian", {

  on.exit(unlink("kf_user_prior", recursive=TRUE), add=TRUE)

  # Used to stop with "Invalid factor."
  expect_error(run_kf_fixture("test_kf_user_prior.ilike", "kf_user_prior"), "built-in Gaussian")
})

test_that("PF with a user-defined transition model (simulate and evaluate) gives the log-likelihood", {

  on.exit(unlink("pf_user_transition", recursive=TRUE), add=TRUE)

  # evaluate_log_transition_model used to be generated with one argument
  # instead of two, so the model did not compile (B15).
  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_pf_fixture("test_pf_user_transition.ilike", "pf_user_transition"), truth, tolerance = 0.01)
})

test_that("PF with a user-defined transition proposal taking the data gives the log-likelihood", {

  on.exit(unlink("pf_transition_proposal", recursive=TRUE), add=TRUE)

  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_pf_fixture("test_pf_transition_proposal.ilike", "pf_transition_proposal"), truth, tolerance = 0.01)
})

test_that("the one-line linear-Gaussian transition and data models work in KF and PF", {

  on.exit({
    unlink("kf_one_line", recursive=TRUE)
    unlink("pf_one_line", recursive=TRUE)
  }, add=TRUE)

  # These built-ins were parsed but nothing used them (B17).
  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_kf_fixture("test_ssm_one_line.ilike", "kf_one_line"), truth, tolerance = 1e-6)
  expect_equal(run_pf_fixture("test_ssm_one_line.ilike", "pf_one_line"), truth, tolerance = 0.01)
})

test_that("PF and EnKF on a linear-Gaussian SSM give the log-likelihood and the final filtering distribution", {

  on.exit({
    unlink("pf_linear_gaussian", recursive=TRUE)
    unlink("enkf_linear_gaussian", recursive=TRUE)
  }, add=TRUE)

  truth = kf_log_likelihood(ssm_y)
  filter = kf_final_filter(ssm_y)
  recipe = ilike::compile(system.file("extdata", "test_kf_linear_gaussian.ilike", package = "ilike"))

  z = ilike::PF(recipe, results_name = "pf_linear_gaussian", number_of_particles = 5000, seed = 1,
                smc_iterations_to_store = 10)
  expect_equal(z, truth, tolerance = 0.01)
  o = ilike::load_smc_output("pf_linear_gaussian")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight)); w = w/sum(w)
  m = sum(w*o$Value)
  expect_equal(m, unname(filter["mean"]), tolerance = 0.1)
  expect_equal(sum(w*(o$Value - m)^2), unname(filter["var"]), tolerance = 0.1)

  z = ilike::EnKF(recipe, results_name = "enkf_linear_gaussian", number_of_ensemble_members = 5000, seed = 1)
  expect_equal(z, truth, tolerance = 0.01)
})

test_that("PF and EnKF work with nonlinear-Gaussian transition and data models given in C++", {

  on.exit({
    unlink("pf_nonlinear", recursive=TRUE)
    unlink("enkf_nonlinear", recursive=TRUE)
  }, add=TRUE)

  truth = kf_log_likelihood(ssm_y)
  recipe = ilike::compile(system.file("extdata", "test_ssm_nonlinear.ilike", package = "ilike"))
  expect_equal(ilike::PF(recipe, results_name = "pf_nonlinear", number_of_particles = 5000, seed = 1), truth, tolerance = 0.01)
  expect_equal(ilike::EnKF(recipe, results_name = "enkf_nonlinear", number_of_ensemble_members = 5000, seed = 1), truth, tolerance = 0.01)
})

test_that("a user function with the wrong arguments gives a clear message", {

  # Used to fail inside the compiler ("different number of parameters (3 vs 2)").
  filename = system.file("extdata", "test_pf_transition_proposal_wrong_arguments.ilike", package = "ilike")
  expect_error(ilike::compile(filename), "simulate_transition_proposal function must have arguments")
})

test_that("KF and the stochastic and square-root EnKF handle a 2-D measurement with correlated noise", {
  # Until C82-C84 the KF stopped on vector measurements (a 1 x 2 data row was
  # assigned to a column vector), and the square-root EnKF used Cholesky
  # factors where the update needs symmetric square roots (about 0.45 too low
  # on this model). Both were correct only for scalar measurements.
  skip_on_cran()
  pl = list(c(0, 0), diag(2))
  model = system.file("extdata", "test_enkf_bivariate.ilike", package = "ilike")
  truth = kf_bivariate_log_likelihood()
  expect_equal(ilike::KF(ilike::compile(model, model_parameter_list = pl), model_parameter_list = pl),
               truth, tolerance = 1e-6)
  enkf = ilike::compile(model, model_parameter_list = pl)
  sqrt_enkf = ilike::compile(c(model, system.file("extdata", "test_enkf_sqrt.ilike", package = "ilike")),
                             model_parameter_list = pl)
  expect_lt(abs(ilike::EnKF(enkf, model_parameter_list = pl, number_of_ensemble_members = 20000, seed = 1) - truth), 0.15)
  expect_lt(abs(ilike::EnKF(sqrt_enkf, model_parameter_list = pl, number_of_ensemble_members = 20000, seed = 1) - truth), 0.15)
  # the ensemble adjustment update (before C183 it used the full SVD with V
  # in place of V^T, and a J x J inverse, too slow at this ensemble size)
  adjustment = tempfile(fileext = ".ilike")
  writeLines("/***enk_shifter,ilike::adjustment()***/", adjustment)
  on.exit(unlink(adjustment), add = TRUE)
  adjustment_enkf = ilike::compile(c(model, adjustment), model_parameter_list = pl)
  expect_lt(abs(ilike::EnKF(adjustment_enkf, model_parameter_list = pl, number_of_ensemble_members = 20000, seed = 1) - truth), 0.15)
})

# (C260) a measurement written as a model statement in a state-space model:
# with its data declared (data { vector y }) the filter uses it; undeclared,
# it would be a prior on a parameter y depending on the state (a filter that
# ignored the measurement, silently), which is refused
test_that("a state-space model's measurement can be a model statement on declared data", {
  set.seed(4); n = 20; x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  ssm = function(meas, first, last) c(first, "model {", "  x ~ norm(0, 1)",
    "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)", meas,
    paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
    "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
    "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }", "}", last)
  outer = function(factor, inner) c("ilike 1.0", paste0("import \"", inner, "\" as inner"), paste0("data { y = [", paste(y, collapse = ", "), "] }"),
    "model {", "  a ~ unif(-1, 1)", paste0("  factor filter = ", factor), "}")
  d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines(ssm("  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)", "ilike 1.0", c("algorithm kalman {", "}")),
             file.path(d, "inner_kf.ilike"))
  writeLines(outer("kf(inner)", "inner_kf.ilike"), file.path(d, "kf.ilike"))
  writeLines(ssm("  y ~ norm(x, 1)", c("ilike 1.0", "data { vector y }"), c("algorithm filter {", "  particles = 2000", "}")),
             file.path(d, "inner.ilike"))
  writeLines(outer("pf(inner)", "inner.ilike"), file.path(d, "pf.ilike"))
  exact = suppressMessages(ilike::evaluate_log(ilike::compile(file.path(d, "kf.ilike")), parameters = list(a = 0.7)))
  estimate = suppressMessages(ilike::evaluate_log(ilike::compile(file.path(d, "pf.ilike")), parameters = list(a = 0.7), seed = 1))
  expect_lt(abs(estimate - exact), 0.3)
  # undeclared: refused
  writeLines(ssm("  y ~ norm(x, 1)", "ilike 1.0", c("algorithm filter {", "  particles = 20", "}")), file.path(d, "inner.ilike"))
  expect_error(ilike::compile(file.path(d, "pf.ilike")), "y is not data, so this would be a prior on a parameter y that depends on the state x",
               fixed = TRUE)
})

# (C264) a Kalman or ensemble Kalman filter as the likelihood in MCMC stores
# each iteration compactly (nothing of the filter is stored by default): the
# chain and the filters' output files are as before
test_that("KF and EnKF likelihoods in MCMC write their output and the chain loads", {
  n = 10
  set.seed(4); x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  x ~ norm(0, 1)", "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)",
               "  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }", "}",
               "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
  for (factor in c("kf(inner)", "enkf(inner, ensemble_members = 20, shifter = \"stochastic\")")) {
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {",
                 "  a ~ unif(-1, 1)", paste0("  factor filter = ", factor), "}",
                 "algorithm mcmc {", "  move m = metropolis(norm_rw(a, 0.2))", "  termination = iterations(50)", "}"), file.path(d, "outer.ilike"))
    out = "kf_enkf_mcmc_output"; unlink(out, recursive = TRUE)
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    suppressMessages(ilike::MCMC(ilike::compile(file.path(d, "outer.ilike")), results_name = out, seed = 1, initial_values = list(list(a = 0.5))))
    o = ilike::load_mcmc_output(out)
    expect_equal(length(unique(o$Iteration)), 50)
    expect_true(any(grepl(if (grepl("^kf", factor)) "factor2_kf/output.h5" else "factor2_enk/output.h5",
                          list.files(out, recursive = TRUE), fixed = TRUE)))
  }
})
