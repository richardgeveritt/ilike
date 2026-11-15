# Julia functions called from block headers (ilike-files.Rmd, "Julia
# functions"). Needs Julia on the PATH and the JuliaCall package.

test_that("Julia functions can give the data, prior and likelihood", {
  skip_on_cran()
  skip_if_not_installed("JuliaCall")
  julia = Sys.which("julia")
  skip_if(julia == "", "Julia is not on the PATH")
  recipe = ilike::compile(system.file("extdata", "test_julia_functions.ilike", package = "ilike"),
                          julia_bin_dir = dirname(julia))
  truth = dlnorm(1.3, 0, 1, log = TRUE) + 3*dnorm(1, 1.3, 1, log = TRUE)
  expect_equal(ilike::evaluate_log(recipe, parameters = list(tau = 1.3)), truth, tolerance = 1e-8)
})
