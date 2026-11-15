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
})
