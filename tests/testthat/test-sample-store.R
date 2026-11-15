# SampleStore (dev/plan.md C90): the compact store of points behind MCMC
# output. Checked against the existing Parameters functions.

test_that("SampleStore stores, selects and rebuilds points like Parameters", {
  r = ilike:::sample_store_check()
  expect_equal(r$size, 5)
  # columns of the chosen variables in the chosen order, as get_rowvec gives them
  expect_equal(r$matrix, r$expected)
  # fill_transposed writes the transpose from the given column on
  expect_equal(r$transposed[, 3:7], t(r$expected))
  expect_equal(r$transposed[, c(1, 2, 8)], matrix(0, 7, 3))
  # a point rebuilt with its shapes and its non-matrix variable
  expect_equal(r$theta3, r$expected_theta3)
  expect_equal(dim(r$theta3), c(2, 3))
  expect_equal(r$sigma3, matrix(-3, 1, 1))
  expect_identical(r$label3, "point3")
  expect_equal(r$extras, cbind(0:4, 100 + 0:4))
  expect_true(r$compatible_same)
  expect_false(r$compatible_wrong_size)
  expect_false(r$compatible_missing)
  # copies are independent
  expect_equal(r$copy_cleared_size, 0)
  expect_equal(r$original_size_after_copy_cleared, 5)
})

test_that("MCMC with thin stores every thin-th iteration of the same chain", {
  on.exit(unlink(c("thin_1", "thin_10"), recursive = TRUE), add = TRUE)
  files = system.file("extdata", c("test_gaussian_posterior.ilike", "test_m_proposal_norm_rw.ilike"), package = "ilike")
  recipe = ilike::compile(files)
  ilike::MCMC(recipe, results_name = "thin_1", seed = 3, initial_values = list(list(theta = 1)))
  ilike::MCMC(recipe, results_name = "thin_10", seed = 3, initial_values = list(list(theta = 1)), thin = 10)
  all = ilike::load_mcmc_output("thin_1"); thinned = ilike::load_mcmc_output("thin_10")
  expect_equal(nrow(thinned), nrow(all)/10)
  # thinning does not change the random numbers: iterations 10, 20, ... of the full chain
  expect_equal(thinned$Value, all$Value[seq(10, nrow(all), by = 10)])
  expect_error(ilike::MCMC(recipe, initial_values = list(list(theta = 1)), thin = 0), "positive whole number")
})
