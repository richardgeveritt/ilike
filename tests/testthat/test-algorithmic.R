test_that("IS with likelihoods computed by KF, PF and EnKF (algorithmic_likelihood)", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(for (k in c("kf","pf","enkf")) unlink(paste("algorithmic_is_", k, sep=""), recursive=TRUE), add=TRUE)

  truth = algorithmic_truth()
  # The fixtures name the inner model by a path relative to the outer file
  # (resolved relative to the working directory before).
  for (method in c("kf", "pf", "enkf"))
  {
    name = paste("algorithmic_is_", method, sep = "")
    z = ilike::IS(ilike::compile(algorithmic_model_file(method)), results_name = name,
                  number_of_importance_points = 1000, seed = 1)
    expect_equal(z, unname(truth["log_evidence"]), tolerance = 0.02)
    o = ilike::load_smc_output(name)
    w = exp(o$LogWeight - max(o$LogWeight)); w = w/sum(w)
    expect_equal(sum(w*o$Value), unname(truth["mean"]), tolerance = 0.1)
  }
})

test_that("MCMC with a likelihood computed by a PF (pseudo-marginal) targets the posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("algorithmic_mcmc_pf", recursive=TRUE), add=TRUE)

  truth = algorithmic_truth()
  algorithm = tempfile(fileext = ".ilike")
  writeLines(c("/***m_proposal,ilike::norm_rw(a,0.5)***/", "", "/***mcmc_termination,ilike::iterations(3000)***/"), algorithm)
  ilike::MCMC(ilike::compile(c(algorithmic_model_file("pf"), algorithm)), results_name = "algorithmic_mcmc_pf",
              initial_values = list(list(a = 0.5)), seed = 1)
  o = ilike::load_mcmc_output("algorithmic_mcmc_pf")
  expect_equal(mean(o$Value), unname(truth["mean"]), tolerance = 0.1)
})

test_that("a missing inner model file gives a clear message", {

  model = tempfile(fileext = ".ilike")
  outer = readLines(algorithmic_model_file("kf"))
  outer = sub("test_algorithmic_inner_ssm.ilike", "no_such_inner_model.ilike", outer, fixed = TRUE)
  writeLines(outer, model)
  expect_error(ilike::compile(model), "no_such_inner_model.ilike not found")
})
