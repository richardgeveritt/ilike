test_that("Langevin and Barker MH proposals use the scale s^2 correctly", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("mcmc_langevin_scaled", recursive=TRUE)
    unlink("mcmc_barker_scaled", recursive=TRUE)
  }, add=TRUE)

  # Before the fix the proposal densities used precision Sigma^-1/s instead of
  # Sigma^-1/s^2 (Langevin) or a z without the 1/s (Barker); with s^2 = 0.25
  # the posterior variance came out as 0.34 and 0.39 instead of 0.5.
  for (k in c("langevin", "barker"))
  {
    o = run_mcmc_fixture("test_gaussian_posterior.ilike",
                         paste("test_mh_proposal_", k, "_scaled.ilike", sep=""),
                         paste("mcmc_", k, "_scaled", sep=""), list(theta = 1))
    expect_equal(mean(o$Value), 1, tolerance = 0.1)
    expect_equal(var(o$Value), 0.5, tolerance = 0.1)
  }
})

test_that("the mirror proposal uses its scale argument", {

  on.exit(unlink("mcmc_mirror_scaled", recursive=TRUE), add=TRUE)

  # With s^2 = 0.0001 about the posterior mean, each proposal is almost the
  # reflection 2 - theta (always accepted here), so theta two iterations later
  # is almost unchanged. The scale used to be ignored (s^2 = 1), which gives a
  # median squared lag-2 change of about 1. (The median, because the small
  # noise lets theta drift, and an occasional rejection then gives one large
  # lag-2 change.)
  o = run_mcmc_fixture("test_gaussian_posterior.ilike", "test_m_proposal_mirror_scaled.ilike",
                       "mcmc_mirror_scaled", list(theta = 1.5))
  v = o$Value
  n = length(v)
  expect_lt(median((v[3:n] - v[1:(n-2)])^2), 0.01)
})

test_that("mvnorm_rw, Langevin and Barker MH proposals work on a 2-D parameter", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    for (k in c("mvnorm_rw", "langevin", "barker"))
      unlink(paste("mcmc_", k, "_2d", sep=""), recursive=TRUE)
  }, add=TRUE)

  # A correlated preconditioning matrix, so that the Cholesky factor in the
  # Barker proposal (y = x + C^T z) matters; s^2 = 0.5. Used to stop with
  # "Mat::init(): requested size is not compatible with column vector layout".
  pl = list(matrix(c(1,0.3,0.3,1),2,2), 0.5)
  for (k in c("mvnorm_rw", "langevin", "barker"))
  {
    o = run_mcmc_fixture("test_gaussian_posterior_vector.ilike",
                         paste("test_mh_proposal_", k, "_2d.ilike", sep=""),
                         paste("mcmc_", k, "_2d", sep=""), list(theta = c(1,1)),
                         model_parameter_list = pl)
    for (d in 1:2)
    {
      v = o$Value[o$Dimension == d]
      expect_equal(mean(v), 1, tolerance = 0.1)
      expect_equal(var(v), 0.5, tolerance = 0.1)
    }
  }
})

test_that("mvnorm_rw accepts (theta, Sigma) and (theta, Sigma, s^2) as M and MH proposals", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  cases = c("test_m_proposal_mvnorm_rw_2d.ilike", "test_m_proposal_mvnorm_rw_2d_short.ilike",
            "test_mh_proposal_mvnorm_rw_2d_short.ilike")
  on.exit(for (k in seq_along(cases)) unlink(paste("mcmc_mvnorm_rw_form", k, sep=""), recursive=TRUE), add=TRUE)

  # m_proposal accepted only the short form and mh_proposal only the long one.
  pl = list(matrix(c(1,0.3,0.3,1),2,2), 0.5)
  for (k in seq_along(cases))
  {
    o = run_mcmc_fixture("test_gaussian_posterior_vector.ilike", cases[k],
                         paste("mcmc_mvnorm_rw_form", k, sep=""), list(theta = c(1,1)),
                         model_parameter_list = pl)
    for (d in 1:2)
    {
      v = o$Value[o$Dimension == d]
      expect_equal(mean(v), 1, tolerance = 0.1)
      expect_equal(var(v), 0.5, tolerance = 0.1)
    }
  }
})
