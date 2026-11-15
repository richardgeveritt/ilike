test_that("IS works when using the prior as proposal, and when using a custom proposal", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("is_prior_as_proposal", recursive=TRUE)
    unlink("is_custom", recursive=TRUE)
  }, add=TRUE)

  # Test IS where prior is used as proposal
  expect_equal(test_is_prior_as_proposal_z(), dnorm(1,1,sqrt(1+1),log=TRUE), tolerance = 0.05)
  sampler_output = ilike::load_smc_output("is_prior_as_proposal")
  m = sum(sampler_output$Value*exp(sampler_output$LogWeight))
  expect_equal(m, 1, tolerance = 0.1)
  expect_equal(sum((sampler_output$Value-m)^2*exp(sampler_output$LogWeight)), 1/(1+1), tolerance = 0.1)

  # Test IS where custom proposal is used
  expect_equal(test_is_custom_proposal_z(), dnorm(1,1,sqrt(1+1),log=TRUE), tolerance = 0.05)
  sampler_output = ilike::load_smc_output("is_custom")
  m = sum(sampler_output$Value*exp(sampler_output$LogWeight))
  expect_equal(m, 1, tolerance = 0.1)
  expect_equal(sum((sampler_output$Value-m)^2*exp(sampler_output$LogWeight)), 1/(1+1), tolerance = 0.1)
})

test_that("IS with built-in gamma and log-normal priors and proposals gives the right evidence", {

  on.exit({
    unlink("is_gamma_prior", recursive=TRUE)
    unlink("is_lnorm_prior", recursive=TRUE)
  }, add=TRUE)

  # Prior and proposal differ, so the gamma log density does not cancel in the
  # weights (it did cancel in the earlier tests, which hid B8).
  truth = log(integrate(function(t) dgamma(t,3,2)*dnorm(1,t,1), 0, Inf)$value)
  expect_equal(test_is_builtin_z("test_is_gamma_prior_lnorm_proposal.ilike", "is_gamma_prior"), truth, tolerance = 0.02)

  truth = log(integrate(function(t) dlnorm(t,0,1)*dnorm(1,t,1), 0, Inf)$value)
  expect_equal(test_is_builtin_z("test_is_lnorm_prior_gamma_proposal.ilike", "is_lnorm_prior"), truth, tolerance = 0.02)
})
