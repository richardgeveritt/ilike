test_that("IS with synthetic likelihood and ABC likelihoods gives the evidence and posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  kinds = c("sl", "gaussian_abc", "euclidean_uniform_abc", "lp_uniform_abc")
  on.exit(for (k in kinds) unlink(paste("sbi_is_", k, sep=""), recursive=TRUE), add=TRUE)

  for (k in kinds)
  {
    truth = sbi_truth(k)
    name = paste("sbi_is_", k, sep = "")
    z = run_is_fixture(paste("test_sbi_", k, ".ilike", sep = ""), name)
    expect_equal(z, unname(truth["log_evidence"]), tolerance = 0.05)
    o = ilike::load_smc_output(name)
    w = exp(o$LogWeight - max(o$LogWeight)); w = w/sum(w)
    m = sum(w*o$Value)
    expect_equal(m, unname(truth["mean"]), tolerance = 0.1)
    expect_equal(sum(w*(o$Value - m)^2), unname(truth["var"]), tolerance = 0.1)
  }
})

test_that("pseudo-marginal MCMC with simulation-based likelihoods targets the posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  kinds = c("sl", "gaussian_abc")
  on.exit(for (k in kinds) unlink(paste("sbi_mcmc_", k, sep=""), recursive=TRUE), add=TRUE)

  for (k in kinds)
  {
    truth = sbi_truth(k)
    o = run_mcmc_fixture(paste("test_sbi_", k, ".ilike", sep = ""), "test_m_proposal_norm_rw.ilike",
                         paste("sbi_mcmc_", k, sep = ""), list(theta = 1))
    expect_equal(mean(o$Value), unname(truth["mean"]), tolerance = 0.1)
    expect_equal(var(o$Value), unname(truth["var"]), tolerance = 0.15)
  }
})
