test_that("SMC with MCMC moves and fixed annealing gives the right evidence and weighted particles", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("smc_annealing", recursive=TRUE), add=TRUE)

  expect_equal(test_smc_annealing_z("smc_annealing"), annealed_log_evidence(1), tolerance = 0.05)

  # Each target after the first is reached by reweighting the particles from
  # the previous one and storing several MCMC iterations per particle. The
  # weights must be attached to every row of their particle (B9).
  output = ilike::load_smc_output("smc_annealing")
  moments = smc_weighted_moments_by_target(output)
  # The power-0 target (the prior) is not written out: target 1 is the prior
  # sample reweighted to power 0.25.
  powers = c(0.25,0.5,0.75,1)
  target_parameters = unique(output[,c("Target","TargetParameters")])
  expect_equal(target_parameters$TargetParameters, paste("power=",powers,";",sep=""))
  expect_equal(nrow(moments), length(powers))
  expect_equal(moments[,"mean"], rep(1,length(powers)), tolerance = 0.1)
  expect_equal(moments[,"var"], 1/(1+powers), tolerance = 0.1)
})

test_that("SMC with MCMC moves and no smc_sequence defaults to importance sampling", {

  on.exit(unlink("smc_default_is", recursive=TRUE), add=TRUE)

  # Used to segfault when writing the output (B5).
  expect_equal(test_smc_z("test_smc_default_is.ilike", "smc_default_is"), annealed_log_evidence(1), tolerance = 0.05)
  output = ilike::load_smc_output("smc_default_is")
  moments = smc_weighted_moments_by_target(output)
  expect_equal(unname(moments[nrow(moments),"mean"]), 1, tolerance = 0.1)
  expect_equal(unname(moments[nrow(moments),"var"]), 0.5, tolerance = 0.1)
})

test_that("ilike::stable stops SMC once the target parameter changes little for several iterations", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("smc_stable", recursive=TRUE), add=TRUE)

  # Used to stop after the first step (B12): stops at power 0.54.
  expect_equal(test_smc_z("test_smc_stable.ilike", "smc_stable"), annealed_log_evidence(0.54), tolerance = 0.05)
  output = ilike::load_smc_output("smc_stable")
  expect_equal(unique(output$TargetParameters), paste("power=",c(0.1,0.5,0.52,0.54),";",sep=""))
})

test_that("ilike::annealing with a factor argument anneals only those factors", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("smc_annealing_factor", recursive=TRUE)
    unlink("smc_annealing_factor_is", recursive=TRUE)
    unlink("smc_annealing_prior", recursive=TRUE)
  }, add=TRUE)

  # Was read as the schedule, annealing to power 2 (B11). Annealing only the
  # likelihood gives the targets prior x likelihood^power.
  powers = c(0.25,0.5,0.75,1)
  for (case in list(c("test_smc_annealing_factor.ilike", "smc_annealing_factor"),
                    c("test_smc_annealing_factor_is_proposal.ilike", "smc_annealing_factor_is")))
  {
    expect_equal(test_smc_z(case[1], case[2]), annealed_log_evidence(1), tolerance = 0.05)
    output = ilike::load_smc_output(case[2])
    moments = smc_weighted_moments_by_target(output)
    expect_equal(unique(output$TargetParameters), paste("power=",powers,";",sep=""))
    expect_equal(moments[,"mean"], rep(1,length(powers)), tolerance = 0.1)
    expect_equal(moments[,"var"], 1/(1+powers), tolerance = 0.1)
  }

  # The prior is the proposal here, so it cancels and cannot be annealed.
  expect_error(test_smc_z("test_smc_annealing_prior_as_proposal.ilike", "smc_annealing_prior"), "cannot be annealed")
})

test_that("adaptive_target ilike::ess chooses each power to give the requested effective sample size", {

  on.exit(unlink("smc_adaptive_ess", recursive=TRUE), add=TRUE)

  # Used to run indefinitely, the bisection reading a stale ESS (B10).
  expect_equal(test_smc_z("test_smc_adaptive_ess.ilike", "smc_adaptive_ess"), annealed_log_evidence(1), tolerance = 0.05)
  output = ilike::load_smc_output("smc_adaptive_ess")
  targets = unique(output[,c("Target","TargetParameters","ISESS")])
  expect_true(nrow(targets) > 1)
  expect_equal(targets$TargetParameters[nrow(targets)], "power=1;")
  # Every intermediate target has ESS 0.95 x 2000 particles; the last may be higher.
  expect_equal(targets$ISESS[-nrow(targets)], rep(1900, nrow(targets)-1), tolerance = 0.01)
  expect_true(targets$ISESS[nrow(targets)] >= 1880)
})

test_that("ilike::parameter sequences work with and without the factor argument", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("smc_parameter", recursive=TRUE)
    unlink("smc_parameter_factors", recursive=TRUE)
  }, add=TRUE)

  # The factor argument used to be required; without it every factor is
  # re-evaluated when the sequence variable changes.
  for (n in c("smc_parameter", "smc_parameter_factors"))
  {
    z = test_smc_z(paste("test_", n, ".ilike", sep=""), n, model_file = "test_gaussian_posterior_tau.ilike")
    expect_equal(z, annealed_log_evidence(1), tolerance = 0.05)
    moments = smc_weighted_moments_by_target(ilike::load_smc_output(n))
    expect_equal(unname(moments[nrow(moments),"mean"]), 1, tolerance = 0.1)
    expect_equal(unname(moments[nrow(moments),"var"]), 0.5, tolerance = 0.1)
  }
})

test_that("tempering, an adaptive CESS target and ESS-based resampling give the evidence and posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  cases = c("test_smc_tempering.ilike", "test_smc_adaptive_cess.ilike", "test_smc_resample_ess.ilike")
  on.exit(for (k in seq_along(cases)) unlink(paste("smc_option_", k, sep=""), recursive=TRUE), add=TRUE)

  for (k in seq_along(cases))
  {
    name = paste("smc_option_", k, sep = "")
    expect_equal(test_smc_z(cases[k], name), annealed_log_evidence(1), tolerance = 0.05)
    moments = smc_weighted_moments_by_target(ilike::load_smc_output(name))
    expect_equal(unname(moments[nrow(moments), "mean"]), 1, tolerance = 0.1)
    expect_equal(unname(moments[nrow(moments), "var"]), 0.5, tolerance = 0.1)
  }
})

test_that("SMC with MCMC moves runs without writing results", {

  recipe = ilike::compile(system.file("extdata", c("test_gaussian_posterior.ilike", "test_smc_annealing.ilike"),
                                      package = "ilike"))
  z = ilike::SMC_with_MCMC(recipe, number_of_particles = 2000, seed = 1)
  expect_equal(z, annealed_log_evidence(1), tolerance = 0.05)
})
