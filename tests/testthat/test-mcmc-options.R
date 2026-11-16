# MCMC options from dev/plan.md section 1.4 (T1, T3, T5, T7, T8, T11, T12).
# Posterior of test_gaussian_posterior.ilike: N(1, 1/2).

test_that("MCMC works with a built-in prior", {
  on.exit(unlink("mcmc_builtin_prior", recursive = TRUE), add = TRUE)
  o = run_mcmc_fixture("test_is_builtin_norm_prior.ilike", "test_m_proposal_norm_rw.ilike",
                       "mcmc_builtin_prior", list(theta = 1))
  expect_equal(mean(o$Value), 1, tolerance = 0.1)
  expect_equal(var(o$Value), 0.5, tolerance = 0.1)
})

test_that("a C++ function can give the mcmc_weights of a stochastic sweep", {
  on.exit(unlink("mcmc_weights_cpp", recursive = TRUE), add = TRUE)
  o = run_mcmc_fixture("test_gaussian_posterior.ilike", "test_mcmc_weights_cpp.ilike",
                       "mcmc_weights_cpp", list(theta = 1))
  expect_equal(mean(o$Value), 1, tolerance = 0.1)
  expect_equal(var(o$Value), 0.5, tolerance = 0.1)
})

test_that("independent MH with a built-in multivariate log-normal proposal", {
  on.exit(unlink("imh_mvlnorm_2d", recursive = TRUE), add = TRUE)
  truth = positive_posterior_moments(function(t) dlnorm(t, 0, 1))
  o = run_mcmc_fixture("test_imh_mvlnorm_posterior_2d.ilike", character(0),
                       "imh_mvlnorm_2d", list(theta = c(1, 1)),
                       model_parameter_list = list(c(0, 0), diag(2)))
  for (d in 1:2)
  {
    v = o$Value[o$Dimension == d]
    expect_equal(mean(v), unname(truth["mean"]), tolerance = 0.1)
    expect_equal(var(v), unname(truth["var"]), tolerance = 0.1)
  }
})

test_that("proposals can use proposal parameters (algorithm_parameter_list) and data", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally
  cases = c("test_m_proposal_proposal_parameters.ilike",
            "test_m_proposal_data.ilike",
            "test_mh_proposal_parameters_and_data.ilike",
            "test_imh_proposal_parameters.ilike",
            "test_imh_proposal_data.ilike")
  on.exit(for (k in seq_along(cases)) unlink(paste("mcmc_proposal_args_", k, sep = ""), recursive = TRUE),
          add = TRUE)
  for (k in seq_along(cases))
  {
    o = run_mcmc_fixture("test_gaussian_posterior.ilike", cases[k],
                         paste("mcmc_proposal_args_", k, sep = ""), list(theta = 1),
                         algorithm_parameter_list = list(sd = 1, mean = 1))
    expect_equal(mean(o$Value), 1, tolerance = 0.1, label = cases[k])
    expect_equal(var(o$Value), 0.5, tolerance = 0.1, label = cases[k])
  }
})

test_that("fixed_parameter_list conditions on a variable the moves do not change", {
  on.exit(unlink("mcmc_fixed_tau", recursive = TRUE), add = TRUE)
  # Likelihood raised to the power tau = 1/2: posterior N(1, 2/3)
  o = run_mcmc_fixture("test_gaussian_posterior_tau.ilike", "test_m_proposal_norm_rw.ilike",
                       "mcmc_fixed_tau", list(theta = 1),
                       fixed_parameter_list = list(tau = 0.5))
  expect_equal(mean(o$Value), 1, tolerance = 0.1)
  expect_equal(var(o$Value), 2/3, tolerance = 0.1)
})

test_that("unfinished MCMC options stop with a clear message", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally
  on.exit(unlink(c("mcmc_hmc", "mcmc_se"), recursive = TRUE), add = TRUE)
  expect_error(run_mcmc_fixture("test_gaussian_posterior.ilike", "test_mh_proposal_hamiltonian.ilike",
                                "mcmc_hmc", list(theta = 1)),
               "hamiltonian \\(HMC\\) is not written yet")
  expect_error(run_mcmc_fixture("test_gaussian_posterior.ilike", "test_mcmc_termination_se.ilike",
                                "mcmc_se", list(theta = 1)),
               "ilike::se is not written yet")
  # Parallel runs exist for IS, SMC_with_MCMC, MCMC and PF (dev/plan.md 2.6a-b)
  expect_error(ilike::EnKF(system.file("extdata", "test_ssm_one_line.ilike", package = "ilike"),
                           number_of_ensemble_members = 10, parallel = TRUE),
               "parallel = TRUE is not available yet")
})
