test_that("MCMC works with all available proposals", {

  on.exit({
    unlink("mcmc_m_proposal", recursive=TRUE)
    unlink("mcmc_m_proposal_prior", recursive=TRUE)
    unlink("mcmc_m_proposal_norm_rw", recursive=TRUE)
    unlink("mcmc_m_proposal_unif_rw", recursive=TRUE)
    unlink("mcmc_m_proposal_mirror", recursive=TRUE)
    unlink("mcmc_imh_proposal", recursive=TRUE)
    unlink("mcmc_imh_proposal_norm", recursive=TRUE)
    unlink("mcmc_mh_proposal_custom", recursive=TRUE)
    unlink("mcmc_mh_proposal_norm_rw", recursive=TRUE)
    unlink("mcmc_mh_proposal_unif_rw", recursive=TRUE)
    unlink("mcmc_mh_proposal_langevin", recursive=TRUE)
    unlink("mcmc_mh_proposal_barker", recursive=TRUE)
    unlink("mcmc_mh_proposal_mirror", recursive=TRUE)
    unlink("mcmc_unadjusted_proposal", recursive=TRUE)
    unlink("mcmc_unadjusted_proposal_langevin", recursive=TRUE)
  }, add=TRUE)

  # Test MCMC with custom Metropolis proposal
  run_mcmc_m_proposal()
  sampler_output = ilike::load_mcmc_output("mcmc_m_proposal")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test custom Metropolis proposal in MCMC where only the first factor (the prior) is used
  run_mcmc_m_proposal_prior()
  sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_prior")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis proposal in MCMC
  # ilike::norm_rw
  run_mcmc_m_proposal_norm_rw()
  sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_norm_rw")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis proposal in MCMC
  # ilike::unif_rw
  run_mcmc_m_proposal_unif_rw()
  sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_unif_rw")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis proposal in MCMC
  # ilike::mirror
  run_mcmc_m_proposal_mirror()
  sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_mirror")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test custom independent Metropolis-Hastings proposal in MCMC
  run_mcmc_imh_proposal()
  sampler_output = ilike::load_mcmc_output("mcmc_imh_proposal")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in independent Metropolis-Hastings proposal in MCMC
  # ilike::norm
  run_mcmc_imh_proposal_norm()
  sampler_output = ilike::load_mcmc_output("mcmc_imh_proposal_norm")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test custom Metropolis-Hastings proposal in MCMC
  run_mcmc_mh_proposal_custom()
  sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_custom")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis-Hastings proposal in MCMC
  # ilike::norm_rw
  run_mcmc_mh_proposal_norm_rw()
  sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_norm_rw")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis-Hastings proposal in MCMC
  # ilike::unif_rw
  run_mcmc_mh_proposal_unif_rw()
  sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_unif_rw")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis-Hastings proposal in MCMC
  # ilike::langevin
  run_mcmc_mh_proposal_langevin()
  sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_langevin")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis-Hastings proposal in MCMC
  # ilike::barker
  run_mcmc_mh_proposal_barker()
  sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_barker")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # Test built-in Metropolis-Hastings proposal in MCMC
  # ilike::mirror
  run_mcmc_mh_proposal_mirror()
  sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_mirror")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  # The unadjusted proposals below are (unadjusted) Langevin moves with step
  # size h = 0.5. There is no accept-reject correction, so the chain does not
  # target the posterior N(1, 1/2) exactly. For a Gaussian target with variance
  # s2, the chain is the AR(1) process
  #   theta' = theta - (h/(2*s2)) (theta - 1) + N(0, h),
  # whose stationary distribution is N(1, h / (1 - (1 - h/(2*s2))^2)).
  # With s2 = 1/2 this variance is 1/(2-h) = 2/3 when h = 0.5.

  # Test custom unadjusted proposal in MCMC
  run_mcmc_unadjusted_proposal()
  sampler_output = ilike::load_mcmc_output("mcmc_unadjusted_proposal")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(2-0.5), tolerance = 0.1)

  # Test built-in unadjusted proposal in MCMC
  # ilike::langevin
  run_mcmc_unadjusted_proposal_langevin()
  sampler_output = ilike::load_mcmc_output("mcmc_unadjusted_proposal_langevin")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(2-0.5), tolerance = 0.1)

})

test_that("MCMC works with deterministic and stochastic sweeps", {

  on.exit({
    unlink("mcmc_sweep_deterministic", recursive=TRUE)
    unlink("mcmc_sweep_stochastic", recursive=TRUE)
  }, add=TRUE)

  # The second move is gradient-based, so this also checks that every move's
  # proposal is registered with the particles (bug B1 in dev/plan.md), and the
  # run returning at all checks that the moves no longer share one Index (B2).
  run_mcmc_sweep_deterministic()
  sampler_output = ilike::load_mcmc_output("mcmc_sweep_deterministic")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

  run_mcmc_sweep_stochastic()
  sampler_output = ilike::load_mcmc_output("mcmc_sweep_stochastic")
  expect_equal(mean(sampler_output$Value), 1, tolerance = 0.1)
  expect_equal(var(sampler_output$Value), 1/(1+1), tolerance = 0.1)

})
