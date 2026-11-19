# Functions in other languages (dev/ilike-language.md, I5): r::, py::, jl::
# calls and code sections; and old files read through the converter (L-8).

foreign_ext = function(f) system.file("extdata", f, package = "ilike")

foreign_mcmc = function(recipe, name)
{
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(recipe, results_name = name, initial_values = list(list("θ" = 1)), seed = 1)
  ilike::load_mcmc_output(name)$Value
}

old_parser = function(files, ...)
{
  old = options(ilike.old_format = "parse")
  on.exit(options(old), add = TRUE)
  suppressMessages(ilike::compile(foreign_ext(files), ...))
}

test_that("R calls and an r section give the same runs as the old format's R expressions", {
  skip_on_cran()
  new = ilike::compile(foreign_ext("lang_rforms.ilike"))
  expect_true(new$uses_r_functions)
  old = old_parser(c("test_rforms_gaussian_posterior.ilike", "test_rforms_mh_proposal.ilike"))
  expect_identical(foreign_mcmc(new, "foreign_new"), foreign_mcmc(old, "foreign_old"))
})

test_that("Python functions imported from a file give the same value as in the old format", {
  skip_on_cran()
  skip_if_not(reticulate::py_available(initialize = TRUE) && reticulate::py_module_available("numpy"))
  truth = dlnorm(1.3, 0, 1, log = TRUE) + 3*dnorm(1, 1.3, 1, log = TRUE)
  expect_equal(ilike::evaluate_log(ilike::compile(foreign_ext("lang_python.ilike")), parameters = list(tau = 1.3)),
               truth, tolerance = 1e-8)
})

test_that("constants used in calls to other languages need values when compiling", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  writeLines(c("ilike 1.0", "constants { real s }", "data { x = 1.0 }",
               "model {", "  prior p = density(log = r::dnorm(theta, 0, s, log = TRUE), simulate = r::list(theta = rnorm(1, 0, s)))",
               "  factor lik = r::dnorm(data$x, theta, 1, log = TRUE)", "}"), f)
  expect_error(ilike::compile(f), "The constant\\(s\\) s are used in a call to another language")
  recipe = ilike::compile(f, model_parameter_list = list(s = 2))
  truth = dnorm(0.5, 0, 2, log = TRUE) + dnorm(1, 0.5, 1, log = TRUE)
  expect_equal(ilike::evaluate_log(recipe, parameters = list(theta = 0.5)), truth, tolerance = 1e-10)
})

test_that("old files read through the converter give the same results as the old parser", {
  skip_on_cran()
  sets = list(c("test_gaussian_posterior.ilike", "test_m_proposal_norm_rw.ilike"),
              c("test_gaussian_posterior_two_params.ilike", "test_gibbs_m_proposal.ilike"),
              c("test_rforms_gaussian_posterior.ilike", "test_rforms_mh_proposal.ilike"))
  for (k in seq_along(sets))
  {
    start = if (k == 2) list(theta1 = 1, theta2 = 1) else if (k == 3) list("θ" = 1) else list(theta = 1)
    run = function(recipe, name)
    {
      on.exit(unlink(name, recursive = TRUE), add = TRUE)
      ilike::MCMC(recipe, results_name = name, initial_values = list(start), seed = 1)
      ilike::load_mcmc_output(name)$Value
    }
    converted = suppressMessages(ilike::compile(foreign_ext(sets[[k]])))
    expect_identical(converted$language, "ilike")
    expect_identical(run(converted, paste0("via_converter_", k)), run(old_parser(sets[[k]]), paste0("via_parser_", k)),
                     label = paste(sets[[k]], collapse = " + "))
  }
})

test_that("calls to R may use names with \".\" and $ after the call; other names may not", {
  skip_on_cran()
  csv = tempfile(fileext = ".csv")
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(csv, f)), add = TRUE)
  utils::write.csv(data.frame(y = c(1.2, 0.4, 2.2)), csv, row.names = FALSE)
  path = normalizePath(csv, winslash = "/")
  write_model = function(data_line, prior = "mu ~ norm(0, 10)")
    writeLines(c("ilike 1.0", paste("data {", data_line, "}"), "model {", paste(" ", prior),
                 "  for (i in 1:3) y[i] ~ norm(mu, 1)", "}",
                 "algorithm mcmc { move m = metropolis(norm_rw(mu, 1)); termination = iterations(500) }"), f)
  chain = function(name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = name, initial_values = list(list(mu = 0)), seed = 1)
    ilike::load_mcmc_output(name)$Value
  }
  write_model("y = [1.2, 0.4, 2.2]")
  literal = chain("dotted_literal")
  write_model(paste("y = r::read.csv(\"", path, "\")$y", sep = ""))
  expect_identical(chain("dotted_csv"), literal)
  write_model(paste("y = r::as.numeric(read.csv(\"", path, "\")$y)", sep = ""))
  expect_identical(chain("dotted_numeric"), literal)

  write_model("y = [1.2, 0.4, 2.2]", "my.mu ~ norm(0, 10)")
  expect_error(ilike::compile(f), "\"my.mu\": names cannot contain \".\"")
  write_model("y = py::np.ones(3L)")
  expect_error(ilike::compile(f), "\"np.ones\": names cannot contain \".\"")
})
