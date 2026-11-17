# MCMC without initial values (dev/plan.md 3.3): each chain starts at an
# independent draw from the recipe's is_proposal, or from the prior if it has
# none. test_gaussian_posterior.ilike: prior N(1, 1), posterior N(1, 1/2).

initial_run <- function(extra_lines, number_of_chains, iterations, sd, seed = 1)
{
  algorithm = tempfile(fileext = ".ilike")
  writeLines(c(sprintf("/***m_proposal,ilike::norm_rw(theta,%s)***/", sd), "",
               sprintf("/***mcmc_termination,ilike::iterations(%d)***/", iterations), "",
               extra_lines), algorithm)
  recipe = suppressMessages(ilike::compile(c(system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike"),
                                             algorithm)))
  results = basename(tempfile("initial"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)
  suppressMessages(ilike::MCMC(recipe, results_name = results, seed = seed, number_of_chains = number_of_chains))
  ilike::load_mcmc_output(results)
}

test_that("without initial values, chains start at independent draws from the prior or the is_proposal", {
  skip_on_cran()
  # a random walk with a huge step almost never accepts: after one iteration
  # the chains are (almost all) still at their starting points
  o = initial_run("", 500, 1, 1000)
  x = o$Value[o$Iteration == 1]
  expect_equal(length(unique(x)), 500)
  expect_equal(mean(x), 1, tolerance = 0.1)
  expect_equal(var(x), 1, tolerance = 0.15)

  o = initial_run("/***is_proposal,ilike::norm(theta,10,0.1)***/", 500, 1, 1000)
  x = o$Value[o$Iteration == 1]
  near = x[abs(x - 10) < 1]
  expect_gt(length(near), 490)
  expect_equal(mean(near), 10, tolerance = 0.01)
  expect_equal(sd(near), 0.1, tolerance = 0.15)
})

test_that("chains started from the prior sample the posterior", {
  skip_on_cran()
  o = initial_run("", 4, 5000, 1, seed = 2)
  expect_equal(sort(unique(o$Chain)), 1:4)
  v = o$Value[o$Iteration > 1000]
  expect_equal(mean(v), 1, tolerance = 0.1)
  expect_equal(var(v), 0.5, tolerance = 0.1)
})
