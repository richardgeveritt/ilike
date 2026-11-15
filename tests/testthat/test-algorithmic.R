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

test_that("nested output: the particle filter's final points and weights inside IS and PMMH", {
  # The particle filter inside an algorithmic likelihood keeps its final
  # generation by default and writes it for every outer point; loaded with
  # nesting_level = 1, the inner points carry outer weight + inner weight, so
  # they are a weighted sample from p(x_T | y) (exact: the Kalman filter
  # averaged over the posterior of a).
  skip_on_cran()
  on.exit(unlink(c("nested_is", "nested_pmmh"), recursive = TRUE), add = TRUE)
  truth = algorithmic_final_state_truth()

  ilike::IS(ilike::compile(algorithmic_model_file("pf")), results_name = "nested_is",
            number_of_importance_points = 1000, seed = 1)
  outer = ilike::load_smc_output("nested_is")
  nested = ilike::load_smc_output("nested_is", nesting_level = 1, factor = 2)
  expect_equal(nrow(nested), 1000*500)
  # the inner weights of each outer point add up to its outer weight
  expect_equal(as.vector(tapply(exp(nested$LogWeight), nested$ExternalIndex, sum)), exp(outer$LogWeight), tolerance = 1e-8)
  expect_equal(nested_final_state_moments(nested), truth, tolerance = 0.1)

  algorithm_file = tempfile(fileext = ".ilike")
  writeLines(c("/***m_proposal,ilike::norm_rw(a,0.5)***/", "", "/***mcmc_termination,ilike::iterations(3000)***/"), algorithm_file)
  on.exit(unlink(algorithm_file), add = TRUE)
  ilike::MCMC(ilike::compile(c(algorithmic_model_file("pf"), algorithm_file)), results_name = "nested_pmmh",
              initial_values = list(list(a = 0.5)), seed = 1)
  chain = ilike::load_mcmc_output("nested_pmmh")
  nested = ilike::load_smc_output("nested_pmmh", as.mcmc = TRUE, nesting_level = 1, factor = 2)
  expect_equal(nrow(nested), nrow(chain)*500)
  expect_equal(nested_final_state_moments(nested), truth, tolerance = 0.15)
})

test_that("nested output: every inner generation and its ancestors when asked for", {
  # third argument of ilike::pf larger than the number of generations: all
  # are kept, with ancestor indices (what particle Gibbs with ancestor
  # sampling needs)
  skip_on_cran()
  d = tempfile(); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  file.copy(system.file("extdata", "test_algorithmic_inner_ssm.ilike", package = "ilike"), d)
  outer_file = file.path(d, "outer.ilike")
  writeLines(sub('ilike::pf("test_algorithmic_inner_ssm.ilike",500)', 'ilike::pf("test_algorithmic_inner_ssm.ilike",50,100)',
                 readLines(algorithmic_model_file("pf")), fixed = TRUE), outer_file)
  ilike::IS(ilike::compile(outer_file), results_name = "o", results_path = d, number_of_importance_points = 20, seed = 1)
  nested = ilike::load_smc_output(file.path(d, "o"), nesting_level = 1, factor = 2)
  expect_equal(sort(unique(nested$Target)), 1:5)
  expect_equal(nrow(nested), 5*20*50)
  for (t in 2:5)
  {
    a = nested$AncestorIndex[nested$Target == t]
    expect_true(all(a >= 1 & a <= 50))
  }
})
