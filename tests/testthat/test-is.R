test_that("IS works when using the prior as proposal, and when using a custom proposal", {

  # Test IS where prior is used as proposal
  expect_equal(signif(test_is_prior_as_proposal_z(), digits = 2), signif(dnorm(1,1,sqrt(1+1),log=TRUE), digits = 2))
  sampler_output = ilike::load_smc_output("is_prior_as_proposal")
  m = sum(sampler_output$Value*exp(sampler_output$LogWeight))
  expect_equal(signif(m, digits = 2), 1)
  expect_equal(signif(sum((sampler_output$Value-m)^2*exp(sampler_output$LogWeight)), digits = 2), signif(1/(1+1), digits = 2))
  unlink("is_prior_as_proposal", recursive=TRUE)

  # Test IS where custom proposal is used
  expect_equal(signif(test_is_custom_proposal_z(), digits = 2), signif(dnorm(1,1,sqrt(1+1),log=TRUE), digits = 2))
  sampler_output = ilike::load_smc_output("is_custom")
  m = sum(sampler_output$Value*exp(sampler_output$LogWeight))
  expect_equal(signif(m, digits = 2), 1)
  expect_equal(signif(sum((sampler_output$Value-m)^2*exp(sampler_output$LogWeight)), digits = 2), signif(1/(1+1), digits = 2))
  unlink("is_custom", recursive=TRUE)
})
