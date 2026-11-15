# ABC-SMC: the ABC tolerance is the SMC sequence variable. These cases used to
# give the answer for tolerance 1 whatever the sequence said (the outer value
# never reached the kernel), or to segfault (the simulations were discarded
# after the first estimate, so they could not be reweighted at the next
# tolerance).

test_that("ABC-SMC with a single tolerance set by the SMC sequence targets that tolerance", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(for (e in c("0.25", "0.5")) unlink(paste("abc_smc_single_", e, sep=""), recursive=TRUE), add=TRUE)

  for (epsilon in c(0.25, 0.5))
  {
    result = run_abc_smc("test_sbi_gaussian_abc.ilike", "test_abc_smc_single.ilike",
                         paste("abc_smc_single_", epsilon, sep = ""), model_parameter_list = list(epsilon))
    expect_abc_smc(result, abc_truth(gaussian_abc_likelihood(epsilon)))
  }
})

test_that("ABC-SMC from an infinite tolerance through a fixed sequence of tolerances", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("abc_smc_fixed", recursive=TRUE), add=TRUE)

  result = run_abc_smc("test_sbi_gaussian_abc.ilike", "test_abc_smc_fixed.ilike", "abc_smc_fixed")
  expect_equal(result$tolerances, c("epsilon=2;", "epsilon=1;", "epsilon=0.5;"))
  expect_abc_smc(result, abc_truth(gaussian_abc_likelihood(0.5)))
})

test_that("ABC-SMC with an adaptive tolerance sequence from infinity (Gaussian and uniform kernels, summary statistics)", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  cases = list(list("test_sbi_gaussian_abc.ilike", gaussian_abc_likelihood(0.5)),
               list("test_sbi_lp_uniform_abc.ilike", function(t) pnorm(1.5 - t) - pnorm(0.5 - t)),
               list("test_sbi_gaussian_abc_summary.ilike", gaussian_abc_likelihood(0.5, v = 0.1)))
  on.exit(for (k in seq_along(cases)) unlink(paste("abc_smc_adaptive_", k, sep=""), recursive=TRUE), add=TRUE)

  for (k in seq_along(cases))
  {
    result = run_abc_smc(cases[[k]][[1]], "test_abc_smc_adaptive.ilike", paste("abc_smc_adaptive_", k, sep = ""))
    expect_true(length(result$tolerances) > 1)
    expect_equal(result$tolerances[length(result$tolerances)], "epsilon=0.5;")
    expect_abc_smc(result, abc_truth(cases[[k]][[2]]))
  }
})
