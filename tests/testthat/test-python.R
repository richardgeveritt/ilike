# Python functions called from block headers (ilike-files.Rmd, "Python
# functions"). Needs a Python with NumPy for reticulate (set RETICULATE_PYTHON
# to choose one). Not run on CI, where reticulate might download a Python.

test_that("Python functions can give the data, prior and likelihood", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not(reticulate::py_module_available("numpy"), "no Python with NumPy")
  recipe = ilike::compile(system.file("extdata", "test_python_functions.ilike", package = "ilike"))
  truth = dlnorm(1.3, 0, 1, log = TRUE) + 3*dnorm(1, 1.3, 1, log = TRUE)
  expect_equal(ilike::evaluate_log(recipe, parameters = list(tau = 1.3)), truth, tolerance = 1e-8)
})
