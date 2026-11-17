# convert_ilike_file() (dev/ilike-language.md section 5, I3c): converted
# recipes give the same runs as the old files.

conv_ext = function(f) system.file("extdata", f, package = "ilike")

conv_compile = function(files)
{
  out = file.path(tempfile("converted"), "recipe.ilike")
  dir.create(dirname(out))
  suppressWarnings(ilike::convert_ilike_file(conv_ext(files), out))
  ilike::compile(out)
}

conv_mcmc = function(recipe, name, start, ...)
{
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(recipe, results_name = name, seed = 2, initial_values = list(start), ...)
  ilike::load_mcmc_output(name)$Value
}

test_that("converted recipes give the same MCMC runs (factor names, sweeps, weights, adaptors)", {
  skip_on_cran()
  sets = list(list(c("test_gaussian_posterior.ilike", "test_m_proposal_prior.ilike"), list(theta = 1)),
              list(c("test_gaussian_posterior_two_params.ilike", "test_gibbs_m_proposal.ilike"), list(theta1 = 1, theta2 = 1)),
              list(c("test_gaussian_posterior_two_params.ilike", "test_gibbs_mh_proposal_langevin.ilike"), list(theta1 = 1, theta2 = 1)),
              list(c("test_gaussian_posterior.ilike", "test_sweep_stochastic.ilike"), list(theta = 1)),
              list(c("test_gaussian_posterior.ilike", "test_mh_proposal.ilike"), list(theta = 1)))
  for (k in seq_along(sets))
  {
    files = sets[[k]][[1]]
    expect_identical(conv_mcmc(conv_compile(files), paste0("conv_new_", k), sets[[k]][[2]]),
                     conv_mcmc(suppressMessages(ilike::compile(conv_ext(files))), paste0("conv_old_", k), sets[[k]][[2]]),
                     label = paste(files, collapse = " + "))
  }
})

test_that("converted SMC, state-space and nested recipes give the same results", {
  skip_on_cran()
  smc = function(recipe, name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::SMC_with_MCMC(recipe, results_name = name, number_of_particles = 200, smc_iterations_to_store = 100, seed = 1)
  }
  files = c("test_gaussian_posterior.ilike", "test_smc_annealing_factor.ilike")
  expect_identical(smc(conv_compile(files), "conv_smc_new"), smc(suppressMessages(ilike::compile(conv_ext(files))), "conv_smc_old"))
  expect_identical(ilike::KF(conv_compile("test_kf_linear_gaussian.ilike")),
                   ilike::KF(suppressMessages(ilike::compile(conv_ext("test_kf_linear_gaussian.ilike")))))
  nested = function(recipe, name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::IS(recipe, results_name = name, number_of_importance_points = 50, seed = 1)
  }
  expect_identical(nested(conv_compile("test_algorithmic_pf.ilike"), "conv_nested_new"),
                   nested(suppressMessages(ilike::compile(conv_ext("test_algorithmic_pf.ilike"))), "conv_nested_old"))
  # constants p1, p2 given by position, as in the old format
  pl = list(c(1, 1), diag(2))
  files = c("test_gaussian_posterior_vector.ilike", "test_imh_proposal_mvnorm_2d.ilike")
  expect_identical(conv_mcmc(conv_compile(files), "conv_pl_new", list(theta = c(1, 1)), model_parameter_list = pl),
                   conv_mcmc(suppressMessages(ilike::compile(conv_ext(files), model_parameter_list = pl)), "conv_pl_old",
                             list(theta = c(1, 1)), model_parameter_list = pl))
})

test_that("the converter writes readable files and refuses what the language cannot express", {
  out = file.path(tempfile("converted"), "recipe.ilike")
  dir.create(dirname(out))
  text = suppressMessages(ilike::convert_ilike_file(conv_ext(c("test_gaussian_posterior_two_params.ilike", "test_gibbs_m_proposal.ilike")), out))
  expect_true(any(grepl("^  move m1 = metropolis\\(norm_rw\\(theta1, 1\\), factors = \\[f1, f2\\]\\)$", text)))
  expect_true(any(grepl("^cpp \\{$", text)))
  expect_error(ilike::convert_ilike_file(conv_ext("test_gaussian_posterior.ilike"), out), "exists: use overwrite = TRUE")
  # R expressions: calls to R (or functions in an r section when the language cannot read them directly)
  text = ilike::convert_ilike_file(conv_ext("test_rforms_gaussian_posterior.ilike"), out, overwrite = TRUE)
  expect_true(any(grepl("prior f1 = density\\(log = r::dnorm\\(parameters\\$", text)))
  expect_true(any(grepl("^ilike_converted_1 = function\\(variable, parameters\\) \\{$", text)))
  expect_error(ilike::convert_ilike_file(conv_ext("lang_builtins.ilike"), out, overwrite = TRUE), "already in the ilike language")
  # an algorithm file on its own: factor numbers become names, with a warning
  warnings = character(0)
  withCallingHandlers(ilike::convert_ilike_file(conv_ext("test_gibbs_m_proposal.ilike"), out, overwrite = TRUE),
                      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(length(warnings) == 2 && all(grepl("are not in these files", warnings)))
})

test_that("compiling an old-format recipe says the format is deprecated, once per session", {
  skip_on_cran()
  assign("old_format_message_shown", FALSE, envir = ilike:::.ilike_env)
  expect_message(ilike::compile(conv_ext("test_gaussian_posterior.ilike")), "deprecated: convert_ilike_file")
  expect_no_message(ilike::compile(conv_ext("test_gaussian_posterior.ilike")), message = "deprecated")
})
