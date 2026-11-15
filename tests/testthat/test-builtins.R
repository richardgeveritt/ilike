test_that("evaluate_log gives the built-in prior densities", {

  m = c(1, 2)
  S = matrix(c(2, 0.5, 0.5, 1), 2, 2)
  pl = list(m, S)
  recipe = ilike::compile(system.file("extdata", "test_builtin_priors.ilike", package = "ilike"),
                          model_parameter_list = pl)
  point = list(a = 0.3, b = 1.3, c = 1.3, d = c(0.5, 2.5), e = c(1.3, 0.7))
  truth = c(dnorm(0.3, 1, 2, log = TRUE),
            dlnorm(1.3, 0, 1, log = TRUE),
            dgamma(1.3, 3, 2, log = TRUE),
            log_dmvnorm(c(0.5, 2.5), m, S),
            log_dmvnorm(log(c(1.3, 0.7)), m, S) - sum(log(c(1.3, 0.7))))
  for (i in 1:5)
    expect_equal(ilike::evaluate_log(recipe, parameters = point, index = i, model_parameter_list = pl),
                 truth[i], tolerance = 1e-8)
  expect_equal(ilike::evaluate_log(recipe, parameters = point, model_parameter_list = pl), sum(truth), tolerance = 1e-8)
})

test_that("IS works with a built-in prior and with a built-in importance proposal", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("is_builtin_prior", recursive=TRUE)
    unlink("is_builtin_proposal", recursive=TRUE)
  }, add=TRUE)

  truth = dnorm(1, 1, sqrt(2), log = TRUE)
  for (case in list(list("test_is_builtin_norm_prior.ilike", "is_builtin_prior"),
                    list(c("test_gaussian_posterior.ilike", "test_is_proposal_norm.ilike"), "is_builtin_proposal")))
  {
    expect_equal(run_is_fixture(case[[1]], case[[2]]), truth, tolerance = 0.02)
    o = ilike::load_smc_output(case[[2]])
    w = exp(o$LogWeight - max(o$LogWeight)); w = w/sum(w)
    m = sum(w*o$Value)
    expect_equal(m, 1, tolerance = 0.1)
    expect_equal(sum(w*(o$Value - m)^2), 0.5, tolerance = 0.1)
  }
})

test_that("independent MH with built-in log-normal, gamma and multivariate normal proposals", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    for (d in c("imh_lnorm", "imh_gamma", "imh_mvnorm_2d")) unlink(d, recursive=TRUE)
  }, add=TRUE)

  # Prior and proposal are the same distribution, so the chain is the
  # posterior of a positive parameter (moments by quadrature).
  for (case in list(list("lnorm", function(t) dlnorm(t, 0, 1)),
                    list("gamma", function(t) dgamma(t, 3, 2))))
  {
    truth = positive_posterior_moments(case[[2]])
    o = run_mcmc_fixture(paste("test_imh_", case[[1]], "_posterior.ilike", sep = ""), character(0),
                         paste("imh_", case[[1]], sep = ""), list(theta = 1))
    expect_equal(mean(o$Value), unname(truth["mean"]), tolerance = 0.1)
    expect_equal(var(o$Value), unname(truth["var"]), tolerance = 0.1)
  }

  o = run_mcmc_fixture("test_gaussian_posterior_vector.ilike", "test_imh_proposal_mvnorm_2d.ilike",
                       "imh_mvnorm_2d", list(theta = c(1, 1)),
                       model_parameter_list = list(c(1, 1), diag(2)))
  for (d in 1:2)
  {
    v = o$Value[o$Dimension == d]
    expect_equal(mean(v), 1, tolerance = 0.1)
    expect_equal(var(v), 0.5, tolerance = 0.1)
  }
})

test_that("uniform random walk and mirror proposals work on a 2-D parameter", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  cases = c("test_m_proposal_mvunif_rw_2d.ilike", "test_mh_proposal_mvunif_rw_2d.ilike",
            "test_m_proposal_mirror_2d.ilike")
  on.exit(for (k in seq_along(cases)) unlink(paste("builtin_2d_", k, sep = ""), recursive=TRUE), add=TRUE)

  pl = list(c(1, 1), 0.5*diag(2), c(1, 1), 1)
  for (k in seq_along(cases))
  {
    o = run_mcmc_fixture("test_gaussian_posterior_vector.ilike", cases[k],
                         paste("builtin_2d_", k, sep = ""), list(theta = c(1, 1)),
                         model_parameter_list = pl)
    for (d in 1:2)
    {
      v = o$Value[o$Dimension == d]
      expect_equal(mean(v), 1, tolerance = 0.1)
      expect_equal(var(v), 0.5, tolerance = 0.1)
    }
  }
})

test_that("MCMC runs several chains and load_mcmc_output gives the ggmcmc format", {

  on.exit(unlink("mcmc_two_chains", recursive=TRUE), add=TRUE)

  recipe = ilike::compile(system.file("extdata", c("test_gaussian_posterior.ilike", "test_m_proposal_norm_rw.ilike"),
                                      package = "ilike"))
  ilike::MCMC(recipe, results_name = "mcmc_two_chains", seed = 1,
              initial_values = list(list(theta = 1), list(theta = -3)))
  o = ilike::load_mcmc_output("mcmc_two_chains")
  expect_equal(sort(unique(o$Chain)), c(1, 2))
  for (chain in 1:2)
  {
    v = o$Value[o$Chain == chain]
    v = v[-(1:1000)]
    expect_equal(mean(v), 1, tolerance = 0.1)
  }

  g = ilike::load_mcmc_output("mcmc_two_chains", ggmcmc = TRUE)
  expect_true(all(c("Iteration", "Chain", "Parameter", "value") %in% names(g)))
  expect_equal(length(unique(g$Chain)), 2)
})
