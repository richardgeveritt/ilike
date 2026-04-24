test_that("MCMC works with all available proposals", {

  # Test MCMC with custom Metropolis proposal
  run_mcmc_m_proposal()
  sampler_output = ilike::load_mcmc_output("mcmc_m_proposal")
  expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  unlink("mcmc_m_proposal", recursive=TRUE)

  # # Test custom Metropolis proposal in MCMC where only the first factor (the prior) is used
  #
  # run_mcmc_m_proposal_prior()
  # sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_prior")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_m_proposal_prior", recursive=TRUE)
  #
  # # Test built-in Metropolis proposal in MCMC
  # # ilike::norm_rw
  #
  # run_mcmc_m_proposal_norm_rw()
  # sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_norm_rw")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_m_proposal_norm_rw", recursive=TRUE)
  #
  # # Test built-in Metropolis proposal in MCMC
  # # ilike::unif_rw
  #
  # run_mcmc_m_proposal_unif_rw()
  # sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_unif_rw")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_m_proposal_unif_rw", recursive=TRUE)
  #
  # # Test built-in Metropolis proposal in MCMC
  # # ilike::mirror
  #
  # run_mcmc_m_proposal_mirror()
  # sampler_output = ilike::load_mcmc_output("mcmc_m_proposal_mirror")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_m_proposal_mirror", recursive=TRUE)
  #
  # # Test custom independent Metropolis-Hastings proposal in MCMC
  #
  # run_mcmc_m_proposal_independent()
  # sampler_output = ilike::load_mcmc_output("mcmc_imh_proposal")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_imh_proposal", recursive=TRUE)
  #
  # # Test built-in independent Metropolis-Hastings proposal in MCMC
  # # ilike::norm
  #
  # run_mcmc_imh_proposal_norm()
  # sampler_output = ilike::load_mcmc_output("mcmc_imh_proposal_norm")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_imh_proposal_norm", recursive=TRUE)
  #
  # # Test custom Metropolis-Hastings proposal in MCMC
  #
  # run_mcmc_mh_proposal_custom()
  # sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_custom")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_mh_proposal_custom", recursive=TRUE)
  #
  # # Test built-in Metropolis-Hastings proposal in MCMC
  # # ilike::norm_rw
  #
  # run_mcmc_mh_proposal_norm_rw()
  # sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_mh_norm_rw")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_mh_proposal_mh_norm_rw", recursive=TRUE)
  #
  # # Test built-in Metropolis-Hastings proposal in MCMC
  # # ilike::unif_rw
  #
  # run_mcmc_mh_proposal_unif_rw()
  # sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_mh_unif_rw")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_mh_proposal_mh_unif_rw", recursive=TRUE)
  #
  # # Test built-in Metropolis-Hastings proposal in MCMC
  # # ilike::langevin
  #
  # run_mcmc_mh_proposal_langevin()
  # sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_langevin")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_mh_proposal_langevin", recursive=TRUE)
  #
  # # Test built-in Metropolis-Hastings proposal in MCMC
  # # ilike::barker
  #
  # run_mcmc_mh_proposal_barker()
  # sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_barker")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_mh_proposal_barker", recursive=TRUE)
  #
  # # Test built-in Metropolis-Hastings proposal in MCMC
  # # ilike::mirror
  #
  # run_mcmc_mh_proposal_mirror()
  # sampler_output = ilike::load_mcmc_output("mcmc_mh_proposal_mirror")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_mh_proposal_mirror", recursive=TRUE)
  #
  # # Test custom unadjusted proposal in MCMC
  #
  # run_mcmc_unadjusted()
  # sampler_output = ilike::load_mcmc_output("mcmc_unadjusted_proposal")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_unadjusted_proposal", recursive=TRUE)
  #
  # # Test built-in unadjusted proposal in MCMC
  # # ilike::langevin
  #
  # run_mcmc_unadjusted_proposal_langevin()
  # sampler_output = ilike::load_mcmc_output("mcmc_unadjusted_proposal_langevin")
  # expect_equal(signif(mean(sampler_output$Value), digits = 2), 1)
  # expect_equal(signif(var(sampler_output$Value), digits = 2), signif(1/(1+1), digits = 2))
  # unlink("mcmc_unadjusted_proposal_langevin", recursive=TRUE)

})
