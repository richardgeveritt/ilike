# Recipes in the ilike language (dev/ilike-language.md, I3b): compiled from
# the language, they give the same runs as the same models in the old format.

lang = function(f) system.file("extdata", f, package = "ilike")

lang_mcmc_values = function(recipe, name, ...)
{
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(recipe, results_name = name, initial_values = list(list(theta = 0)), seed = 3, ...)
  ilike::load_mcmc_output(name)$Value
}

lang_run = function(f, name, ...)
{
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  f(..., results_name = name)
}

test_that("a model in the ilike language gives the same MCMC, IS and SMC runs as in the old format", {
  skip_on_cran()
  old = ilike::compile(lang(c("test_gaussian_posterior.ilike", "test_m_proposal_norm_rw.ilike")))
  new = ilike::compile(lang("lang_gaussian_posterior.ilike"))
  expect_equal(names(new$algorithms), c("mcmc", "smc"))
  expect_identical(lang_mcmc_values(new, "lang_new", algorithm = "mcmc"), lang_mcmc_values(old, "lang_old"))
  expect_identical(lang_run(ilike::IS, "lang_is_new", new, number_of_importance_points = 2000, seed = 1, algorithm = "mcmc"),
                   lang_run(ilike::IS, "lang_is_old", old, number_of_importance_points = 2000, seed = 1))
  old_smc = ilike::compile(lang(c("test_gaussian_posterior.ilike", "test_smc_annealing.ilike")))
  expect_identical(lang_run(ilike::SMC_with_MCMC, "lang_smc_new", new, number_of_particles = 500,
                            smc_iterations_to_store = 100, seed = 1, algorithm = "smc"),
                   lang_run(ilike::SMC_with_MCMC, "lang_smc_old", old_smc, number_of_particles = 500,
                            smc_iterations_to_store = 100, seed = 1))
  expect_error(ilike::MCMC(new, results_name = "lang_x", initial_values = list(list(theta = 0))),
               "several algorithms \\(mcmc, smc\\)")
  expect_error(ilike::MCMC(new, results_name = "lang_x", initial_values = list(list(theta = 0)), algorithm = "rwm"),
               "no algorithm called rwm")
})

test_that("proposals written by the user and simulation-based factors match the old format", {
  skip_on_cran()
  new = ilike::compile(lang(c("lang_gaussian_posterior.ilike", "lang_proposals.ilike")))
  for (pair in list(c("test_mh_proposal.ilike", "mh"), c("test_imh_proposal.ilike", "imh")))
  {
    old = ilike::compile(lang(c("test_gaussian_posterior.ilike", pair[1])))
    expect_identical(lang_mcmc_values(new, paste0("lang_n_", pair[2]), algorithm = pair[2]),
                     lang_mcmc_values(old, paste0("lang_o_", pair[2])))
  }
  expect_identical(lang_run(ilike::IS, "lang_n_is", new, number_of_importance_points = 2000, seed = 1, algorithm = "is"),
                   lang_run(ilike::IS, "lang_o_is", ilike::compile(lang(c("test_gaussian_posterior.ilike", "test_is_custom_proposal.ilike"))),
                            number_of_importance_points = 2000, seed = 1))
  expect_identical(lang_run(ilike::IS, "lang_n_sbi", ilike::compile(lang("lang_sbi_sl.ilike")), number_of_importance_points = 2000, seed = 1),
                   lang_run(ilike::IS, "lang_o_sbi", ilike::compile(lang("test_sbi_sl.ilike")), number_of_importance_points = 2000, seed = 1))
})

test_that("built-in priors, literal data, named constants and adaptive moves work in the language", {
  skip_on_cran()
  recipe = ilike::compile(lang("lang_builtins.ilike"))
  values = lang_mcmc_values(recipe, "lang_builtins", model_parameter_list = list(step = 2))
  # prior N(1, 1), one observation x = 1 ~ N(theta, 1): the posterior is N(1, 0.5)
  expect_equal(mean(values), 1, tolerance = 0.1)
  expect_equal(var(values), 0.5, tolerance = 0.1)
  # step has no default
  expect_error(lang_mcmc_values(recipe, "lang_builtins_x"), "step is not the name of a run-time constant")
  expect_error(lang_mcmc_values(recipe, "lang_builtins_y", model_parameter_list = list(2)), "given by name")
})

test_that("mistakes in recipes in the language are reported with their place", {
  skip_on_cran()
  write = function(lines) { f = tempfile(fileext = ".ilike"); writeLines(c("ilike 1.0", lines), f); f }
  cpp = c("cpp {", "double loglik(const Parameters &parameters, const Data &data) { return 0.0; }",
          "double wrong(const Parameters &parameters) { return 0.0; }", "}")
  expect_error(ilike::compile(write(c("model { factor lik = cpp::missing }", cpp))),
               ":2:22: cpp::missing is not defined in a cpp section")
  expect_error(ilike::compile(write(c("model { factor lik = cpp::wrong }", cpp))),
               "cpp::wrong \\(used as evaluate_log_likelihood\\): it has arguments \\(const Parameters &parameters\\)")
  expect_error(ilike::compile(write(c("model { theta ~ norm(0, 1) }", "algorithm { move m = metropolis(norm_rw(theta, 1), factors = [lik]) }"))),
               "no factor called lik \\(factors: theta\\)")
  expect_error(ilike::compile(write(c("model { theta ~ norm(0, 1) }", "algorithm { move m = metropolis(gamma(theta, 1, 1)) }"))),
               "ilike::gamma cannot be used as a Metropolis proposal")
  expect_error(ilike::compile(write(c("model { theta ~ norm(0, 1) }", "algorithm { move m = gibbs(norm_rw(theta, 1)) }"))),
               "a move is metropolis\\(…\\)")
  expect_error(ilike::compile(write("model { theta ~ norm(0, 1)\ntheta ~ norm(0, 2) }")), "\"theta\" is defined more than once")
})

test_that("state-space models and algorithmic likelihoods in the language match the old format", {
  skip_on_cran()
  expect_identical(ilike::KF(ilike::compile(lang("lang_ssm_functions.ilike"))),
                   ilike::KF(ilike::compile(lang("test_kf_linear_gaussian.ilike"))))
  expect_identical(ilike::KF(ilike::compile(lang("lang_ssm_one_line.ilike"))),
                   ilike::KF(ilike::compile(lang("test_ssm_one_line.ilike"))))
  expect_identical(ilike::EnKF(ilike::compile(lang("lang_ssm_one_line.ilike")), number_of_ensemble_members = 500, seed = 1),
                   ilike::EnKF(ilike::compile(lang("test_ssm_one_line.ilike")), number_of_ensemble_members = 500, seed = 1))
  expect_identical(ilike::PF(ilike::compile(lang("lang_ssm_nonlinear.ilike")), number_of_particles = 500, seed = 1),
                   ilike::PF(ilike::compile(lang("test_ssm_nonlinear.ilike")), number_of_particles = 500, seed = 1))
  user = ilike::compile(lang("lang_pf_user.ilike"))
  expect_identical(ilike::PF(user, number_of_particles = 500, seed = 1, algorithm = "filter"),
                   ilike::PF(ilike::compile(lang("test_pf_user_transition.ilike")), number_of_particles = 500, seed = 1))
  expect_identical(ilike::PF(user, number_of_particles = 500, seed = 1, algorithm = "proposal_filter"),
                   ilike::PF(ilike::compile(lang("test_pf_transition_proposal.ilike")), number_of_particles = 500, seed = 1))
  # pf(inner): the inner file is imported; particles = 500 comes from its algorithm section
  expect_identical(lang_run(ilike::IS, "lang_alg_new", ilike::compile(lang("lang_algorithmic_pf.ilike")),
                            number_of_importance_points = 100, seed = 1),
                   lang_run(ilike::IS, "lang_alg_old", ilike::compile(lang("test_algorithmic_pf.ilike")),
                            number_of_importance_points = 100, seed = 1))
})

test_that("run-time settings in an algorithm section are used unless given in the call", {
  skip_on_cran()
  recipe = ilike::compile(lang("lang_ssm_settings.ilike"))
  from_file = ilike::PF(recipe, seed = 1, algorithm = "filter")
  expect_identical(from_file, ilike::PF(recipe, number_of_particles = 300, seed = 1, algorithm = "filter"))
  expect_false(identical(from_file, ilike::PF(recipe, number_of_particles = 50, seed = 1, algorithm = "filter")))
  expect_error(ilike::PF(recipe, seed = 1, algorithm = "bad"), "The setting particels in the algorithm section is not a setting")
})
