# Adaptation of MCMC proposals from the SMC particles (dev/plan.md 3.2, B2):
# /***smc_adaptor,ilike::gaussian(theta[, c])***/ after a proposal sets its
# covariance to the particles' weighted covariance before each move (and a
# mirror's or an independent proposal's mean to their weighted mean), and
# with c the random walk's scale to c^2/d.

smc_adaptor_run <- function(proposal_lines, seed, ...)
{
  m = c(0.5, -0.5); S = 4*matrix(c(1, 0.95, 0.95, 1), 2)
  algorithm = tempfile(fileext = ".ilike")
  writeLines(c(proposal_lines, "",
               "/***mcmc_termination,ilike::iterations(5)***/", "",
               "/***smc_sequence,ilike::annealing(power,c(0,1))***/", "",
               "/***adaptive_target,ilike::cess(0.9)***/"), algorithm)
  pl = list(m, S, 100*diag(2), 1, c(3, -3))
  recipe = suppressMessages(ilike::compile(c(system.file("extdata", "test_mvnorm_correlated_prior.ilike", package = "ilike"),
                                             algorithm), model_parameter_list = pl))
  suppressMessages(ilike::SMC_with_MCMC(recipe, model_parameter_list = pl, number_of_particles = 1000,
                                        seed = seed, ...))
}

# exact evidence: x = (1, 1) ~ N(m, S + I)
smc_adaptor_log_z <- function()
{
  m = c(0.5, -0.5); C = 4*matrix(c(1, 0.95, 0.95, 1), 2) + diag(2); d = c(1, 1) - m
  -log(2*pi) - 0.5*log(det(C)) - 0.5*sum(d*solve(C, d))
}

test_that("adapting a random walk and an independent proposal from the particles reduces the variance of the evidence", {
  skip_on_cran()
  for (proposal in list(c("/***m_proposal,ilike::mvnorm_rw(theta,p3,p4)***/", "/***smc_adaptor,ilike::gaussian(theta,2.38)***/"),
                        c("/***imh_proposal,ilike::mvnorm(theta,p5,p3)***/", "/***smc_adaptor,ilike::gaussian(theta)***/")))
  {
    plain = sapply(1:8, function(seed) smc_adaptor_run(proposal[1], seed))
    adapted = sapply(1:8, function(seed) smc_adaptor_run(proposal, seed))
    expect_lt(abs(mean(adapted) - smc_adaptor_log_z()), 0.05)
    expect_lt(sd(adapted), sd(plain))
  }
})

test_that("a run with an SMC adaptor is the same in parallel, and works with waste-free SMC", {
  skip_on_cran()
  proposal = c("/***m_proposal,ilike::mvnorm_rw(theta,p3,p4)***/", "/***smc_adaptor,ilike::gaussian(theta,2.38)***/")
  expect_identical(smc_adaptor_run(proposal, 3, parallel = TRUE), smc_adaptor_run(proposal, 3))
  z = sapply(1:4, function(seed) smc_adaptor_run(proposal, seed, waste_free_chains = 50))
  expect_lt(abs(mean(z) - smc_adaptor_log_z()), 0.15)
})

test_that("an SMC adaptor must follow a proposal with a Gaussian covariance", {
  skip_on_cran()
  expect_error(smc_adaptor_run(c("/***m_proposal,ilike::unif_rw(theta,1)***/", "/***smc_adaptor,ilike::gaussian(theta)***/"), 1),
               "has no Gaussian covariance")
  expect_error(smc_adaptor_run(c("/***mcmc_termination,ilike::iterations(5)***/", "/***smc_adaptor,ilike::gaussian(theta)***/"), 1),
               "must follow the MCMC proposal it adapts")
})
