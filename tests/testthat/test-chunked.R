# Chunked generations (dev/plan.md 2.6a). An importance sampler with a fixed
# target simulates and weights its points grain_size at a time, with reused
# particle objects, and keeps the points as numbers. One chunk (grain_size at
# least the number of points) must give exactly the results of the unchunked
# path (recipe$force_rich_generations = TRUE); other chunk sizes change the
# random numbers but not what is estimated.

chunk_recipe <- function(files, rich = FALSE)
{
  recipe = ilike::compile(system.file("extdata", files, package = "ilike"))
  recipe$force_rich_generations = rich
  recipe
}

# the loaded output without the wall-clock time
loaded_without_time <- function(results_name)
{
  o = ilike::load_smc_output(results_name)
  o[, setdiff(names(o), "Time")]
}

weighted_moments <- function(output)
{
  w = exp(output$LogWeight - max(output$LogWeight)); w = w/sum(w)
  m = sum(w*output$Value)
  c(mean = m, var = sum(w*(output$Value - m)^2))
}

test_that("chunked IS with one chunk gives exactly the unchunked results and output", {
  on.exit(unlink(c("chunk_one", "chunk_rich"), recursive = TRUE), add = TRUE)

  z_chunked = ilike::IS(chunk_recipe("test_gaussian_posterior.ilike"), results_name = "chunk_one",
                        number_of_importance_points = 500, seed = 3, grain_size = 1000)
  z_rich = ilike::IS(chunk_recipe("test_gaussian_posterior.ilike", rich = TRUE), results_name = "chunk_rich",
                     number_of_importance_points = 500, seed = 3, grain_size = 1000)
  expect_identical(z_chunked, z_rich)
  expect_identical(loaded_without_time("chunk_one"), loaded_without_time("chunk_rich"))
  # the factor outputs (prior and likelihood, both exact)
  for (f in c("factor1_exact.h5", "factor2_exact.h5"))
  {
    f1 = hdf5r::H5File$new(file.path("chunk_one", "ilike_smc", "iteration1", f), "r")
    f2 = hdf5r::H5File$new(file.path("chunk_rich", "ilike_smc", "iteration1", f), "r")
    x1 = f1[["log_likelihood"]]$read(); x2 = f2[["log_likelihood"]]$read()
    f1$close_all(); f2$close_all()
    expect_length(x1, 500)
    expect_identical(x1, x2)
  }
})

test_that("chunked IS gives the right evidence and posterior for chunk sizes that do and do not divide the number of points", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("chunk_sizes", recursive = TRUE), add = TRUE)

  recipe = chunk_recipe("test_gaussian_posterior.ilike")
  for (grain_size in c(1, 7, 255, 1999))
  {
    z = ilike::IS(recipe, results_name = "chunk_sizes", number_of_importance_points = 2000,
                  seed = 2, grain_size = grain_size)
    expect_equal(z, dnorm(1, 1, sqrt(2), log = TRUE), tolerance = 0.05)
    moments = weighted_moments(ilike::load_smc_output("chunk_sizes"))
    expect_equal(unname(moments["mean"]), 1, tolerance = 0.1)
    expect_equal(unname(moments["var"]), 0.5, tolerance = 0.1)
    unlink("chunk_sizes", recursive = TRUE)
  }
})

test_that("chunked IS with a proposal that is evaluated, and without factor output", {
  on.exit(unlink(c("chunk_prop", "chunk_prop_rich"), recursive = TRUE), add = TRUE)

  files = c("test_gaussian_posterior.ilike", "test_is_proposal_norm.ilike")
  z_chunked = ilike::IS(chunk_recipe(files), results_name = "chunk_prop", write_factor_output = FALSE,
                        number_of_importance_points = 300, seed = 1, grain_size = 300)
  z_rich = ilike::IS(chunk_recipe(files, rich = TRUE), results_name = "chunk_prop_rich", write_factor_output = FALSE,
                     number_of_importance_points = 300, seed = 1, grain_size = 300)
  expect_identical(z_chunked, z_rich)
  expect_identical(loaded_without_time("chunk_prop"), loaded_without_time("chunk_prop_rich"))
  expect_false(dir.exists(file.path("chunk_prop", "ilike_smc", "iteration1")))
})

test_that("adaptive SMC choosing its targets from adaptive summaries gives exactly the results of the rich weighting", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink(c("chunk_cess", "chunk_cess_rich"), recursive = TRUE), add = TRUE)

  files = c("test_gaussian_posterior.ilike", "test_smc_adaptive_cess.ilike")
  z = ilike::SMC_with_MCMC(chunk_recipe(files), results_name = "chunk_cess", number_of_particles = 500,
                           smc_iterations_to_store = 100, seed = 1)
  z_rich = ilike::SMC_with_MCMC(chunk_recipe(files, rich = TRUE), results_name = "chunk_cess_rich",
                                number_of_particles = 500, smc_iterations_to_store = 100, seed = 1)
  expect_identical(z, z_rich)
  expect_identical(loaded_without_time("chunk_cess"), loaded_without_time("chunk_cess_rich"))
})

# SMC with MCMC moves keeps its moved generations flat (dev/plan.md 2.6a, G3)
# when every factor can be restored from its adaptive summary: the results
# and output must be exactly those of the per-particle path.
smc_flat_and_rich <- function(files, results_name, number_of_particles, extra_algorithm = NULL)
{
  if (!is.null(extra_algorithm))
  {
    algorithm = tempfile(fileext = ".ilike")
    writeLines(extra_algorithm, algorithm)
    files = c(system.file("extdata", files, package = "ilike"), algorithm)
  }
  else
  {
    files = system.file("extdata", files, package = "ilike")
  }
  out = list()
  for (rich in c(FALSE, TRUE))
  {
    recipe = ilike::compile(files)
    recipe$force_rich_generations = rich
    name = paste0(results_name, if (rich) "_rich" else "")
    z = ilike::SMC_with_MCMC(recipe, results_name = name, number_of_particles = number_of_particles,
                             smc_iterations_to_store = 100, seed = 1)
    out[[length(out) + 1]] = list(z = z, output = loaded_without_time(name))
  }
  out
}

test_that("SMC with MCMC moves on flat generations gives exactly the per-particle results (annealing, chains of 5)", {
  on.exit(unlink(c("chunk_smc", "chunk_smc_rich"), recursive = TRUE), add = TRUE)
  out = smc_flat_and_rich(c("test_gaussian_posterior.ilike", "test_smc_annealing.ilike"), "chunk_smc", 300)
  expect_identical(out[[1]]$z, out[[2]]$z)
  expect_identical(out[[1]]$output, out[[2]]$output)
})

test_that("SMC with Langevin moves on flat generations gives exactly the per-particle results", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink(c("chunk_smc_mala", "chunk_smc_mala_rich"), recursive = TRUE), add = TRUE)
  out = smc_flat_and_rich("test_gaussian_posterior.ilike", "chunk_smc_mala", 300,
                          c("/***mh_proposal,ilike::langevin(theta,1,0.5)***/", "",
                            "/***mcmc_termination,ilike::iterations(3)***/", "",
                            "/***smc_sequence,ilike::annealing(power,c(0,0.3,1))***/"))
  expect_identical(out[[1]]$z, out[[2]]$z)
  expect_identical(out[[1]]$output, out[[2]]$output)
})
