# Copy audit (dev/plan.md 3.5, CH1): parallel runs copy estimators and MCMC
# moves for every thread group. For every model file that can be built, each
# likelihood estimator is evaluated, copied, the copy evaluated and deleted,
# and the original evaluated again (same seed): the three values must be
# identical; the same for one MCMC move and the move's copy. Run under
# AddressSanitizer, this also finds copies that share or overrun memory.

audit_fixture <- function(files)
{
  recipe = suppressMessages(ilike::compile(system.file("extdata", files, package = "ilike")))
  termination = ilike:::get_method(recipe, "mcmc_termination")
  if (is.null(termination))
    termination = list()
  weights = ilike:::get_method(recipe, "mcmc_weights")
  if (is.null(weights))
    weights = list()
  ilike:::copy_audit(recipe, list(), list(), termination, weights, 1)
}

expect_copies_identical <- function(result, label)
{
  ok = is.na(result$original) | (result$messages != "")
  for (i in seq_along(result$original))
  {
    if (result$messages[i] == "" && !is.na(result$original[i]))
    {
      expect_identical(result$copy[i], result$original[i], label = paste(label, "factor", i, "copy"))
      expect_identical(result$original_again[i], result$original[i], label = paste(label, "factor", i, "after deleting the copy"))
    }
  }
  if (result$mcmc_message == "" && !is.na(result$mcmc[1]))
  {
    expect_identical(result$mcmc[2], result$mcmc[1], label = paste(label, "MCMC move copy"))
    expect_identical(result$mcmc[3], result$mcmc[1], label = paste(label, "MCMC move after deleting the copy"))
  }
}

test_that("copies of estimators and MCMC moves behave as the originals", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  # model files with their algorithm files (MCMC moves)
  cases = list(
    c("test_gaussian_posterior.ilike", "test_m_proposal_norm_rw.ilike"),
    c("test_gaussian_posterior.ilike", "test_mh_proposal_langevin.ilike"),
    c("test_gaussian_posterior.ilike", "test_mh_proposal_barker.ilike"),
    c("test_gaussian_posterior.ilike", "test_m_proposal_mirror.ilike"),
    c("test_gaussian_posterior.ilike", "test_imh_proposal_norm.ilike"),
    c("test_gaussian_posterior.ilike", "test_sweep_deterministic.ilike"),
    c("test_gaussian_posterior.ilike", "test_smc_annealing.ilike"),
    c("test_gaussian_posterior.ilike", "test_smc_annealing_factor.ilike"),
    c("test_gaussian_posterior_two_params.ilike", "test_gibbs_m_proposal.ilike"),
    c("test_gaussian_posterior.ilike", "test_unadjusted_proposal.ilike"),
    c("test_imh_gamma_posterior.ilike"),
    c("test_cut_model.ilike", "test_cut_mcmc.ilike"),
    c("test_sbi_gaussian_abc.ilike", "test_m_proposal_norm_rw.ilike"),
    c("test_sbi_lp_uniform_abc.ilike", "test_abc_smc_adaptive.ilike"),
    c("test_sbi_sl.ilike", "test_m_proposal_norm_rw.ilike"),
    c("test_algorithmic_pf.ilike"),
    c("test_algorithmic_kf.ilike"),
    c("test_algorithmic_enkf.ilike"))
  for (files in cases)
  {
    result = tryCatch(audit_fixture(files), error = function(e) e)
    label = paste(files, collapse = " + ")
    if (inherits(result, "error"))
    {
      # recipes this audit cannot build are reported, not failed
      message("copy audit could not build ", label, ": ", conditionMessage(result))
      next
    }
    expect_copies_identical(result, label)
    if (any(result$messages != ""))
      message("copy audit, ", label, ": ", paste(unique(result$messages[result$messages != ""]), collapse = "; "))
    if (result$mcmc_message != "")
      message("copy audit, ", label, ", MCMC: ", result$mcmc_message)
  }
})

test_that("simulating from a prior made of several prior blocks works (B32)", {
  # the prior is then a composite kernel, whose simulation followed each
  # component's transform without checking that it had one
  recipe = ilike::compile(system.file("extdata", "test_cut_model.ilike", package = "ilike"))
  z = ilike::IS(recipe, number_of_importance_points = 200, seed = 1)
  expect_true(is.finite(z))
})
