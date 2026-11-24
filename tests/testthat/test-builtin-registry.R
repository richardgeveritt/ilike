# The registry of built-ins (dev/plan.md 3.14, I2): one table gives every
# built-in's arguments; keyword arguments and named run-time constants work
# in today's format; mistakes give messages naming the built-in and argument.

registry_model = c("/***evaluate_log_prior***/",
                   "double f(const Parameters &p) { return ilike::dnorm(p[\"theta\"][0],0,1); }", "",
                   "/***evaluate_log_likelihood***/",
                   "double g(const Parameters &p, const Data &d) { return ilike::dnorm(1.0,p[\"theta\"][0],1); }", "")

registry_run = function(lines, constants = list(), iterations = 2000)
{
  f = tempfile(fileext = ".ilike")
  writeLines(c(registry_model, lines, paste0("/***mcmc_termination,ilike::iterations(", iterations, ")***/")), f)
  name = basename(tempfile("registry_"))
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, model_parameter_list = constants,
              initial_values = list(list(theta = 0)), seed = 1)
  ilike::load_mcmc_output(name)$Value
}

test_that("the registry lists every built-in with its arguments", {
  tab = ilike:::builtin_table()
  expect_true(all(c("builtin", "roles", "argument", "type", "required", "default") %in% names(tab)))
  expect_equal(tab$argument[tab$builtin == "mvnorm_rw"], c("variable", "cov", "scale"))
  expect_equal(tab$default[tab$builtin == "mvnorm_rw" & tab$argument == "scale"], "1")
  expect_true(all(c("norm", "gamma", "norm_rw", "mirror", "scaled_am", "gaussian", "iterations") %in% tab$builtin))
  expect_match(tab$roles[tab$builtin == "gamma"][1], "imh_proposal")
})

test_that("keyword arguments and named constants give the same chain as positional ones", {
  skip_on_cran()
  positional = registry_run("/***m_proposal,ilike::mvnorm_rw(theta,p1,2)***/", list(diag(1)))
  named = registry_run("/***m_proposal,ilike::mvnorm_rw(theta,scale=2,cov=Sigma)***/", list(Sigma = diag(1)))
  expect_identical(positional, named)
  adaptive = registry_run(c("/***m_proposal,ilike::mvnorm_rw(theta,p1)***/", "/***m_adaptor,ilike::scaled_am(theta,averaged)***/"),
                          list(diag(1)))
  adaptive_named = registry_run(c("/***m_proposal,ilike::mvnorm_rw(theta,Sigma)***/",
                                  "/***m_adaptor,ilike::scaled_am(theta,averaged=TRUE,kappa=0.6)***/"),
                                list(Sigma = diag(1)))
  expect_identical(adaptive, adaptive_named)
})

test_that("mistakes in built-in calls are reported with the built-in and argument", {
  skip_on_cran()
  constants = list(Sigma = diag(1))
  expect_error(registry_run("/***m_proposal,ilike::mvnorm_rw(theta,cov=Sigma,scal=2)***/", constants, 10),
               "no argument called scal \\(it takes ilike::mvnorm_rw\\(variable, cov, scale = 1\\)\\)")
  expect_error(registry_run("/***m_proposal,ilike::mvnorm_rw(theta)***/", constants, 10),
               "the argument cov is missing")
  expect_error(registry_run("/***m_proposal,ilike::mvnorm_rw(theta,Sgma)***/", constants, 10),
               "Sgma is not the name of a run-time constant \\(constants given: Sigma\\)")
  expect_error(registry_run("/***m_proposal,ilike::gamma(theta,1,1)***/", constants, 10),
               "ilike::gamma cannot be used as a Metropolis proposal")
  expect_error(registry_run("/***m_proposal,ilike::mvnrm_rw(theta,Sigma)***/", constants, 10),
               "Unknown Metropolis proposal ilike::mvnrm_rw")
  expect_error(registry_run("/***m_proposal,ilike::norm_rw(theta,2abc)***/", constants, 10),
               "\"2abc\" is not a real")
  expect_error(registry_run("/***m_proposal,ilike::norm_rw(theta,1,2)***/", constants, 10),
               "too many arguments")
  expect_error(registry_run("/***m_proposal,ilike::norm_rw(theta,sd=1,sd=2)***/", constants, 10),
               "given more than once")
  expect_error(registry_run("/***m_proposal,ilike::norm_rw(theta,p2)***/", constants, 10),
               "p2 refers to the 2nd run-time constant, but 1 was given")
  expect_error(registry_run(c("/***m_proposal,ilike::norm_rw(theta,1)***/",
                              "/***m_adaptor,ilike::scale(theta,target=1.5)***/"), constants, 10),
               "argument target: must be between 0 and 1")
})

# SMC settings (I2b): criteria, terminations and sequences through the registry
registry_smc = function(lines, constants = list())
{
  f = tempfile(fileext = ".ilike")
  writeLines(c("/***m_proposal,ilike::norm_rw(theta,1)***/", "", "/***mcmc_termination,ilike::iterations(5)***/", "",
               lines), f)
  name = basename(tempfile("registry_smc_"))
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  recipe = ilike::compile(c(system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike"), f))
  ilike::SMC_with_MCMC(recipe, results_name = name, number_of_particles = 500, smc_iterations_to_store = 100,
                       model_parameter_list = constants, seed = 1)
}

test_that("SMC sequences and criteria take keyword arguments and named constants", {
  skip_on_cran()
  positional = registry_smc(c("/***smc_sequence,ilike::annealing(power,c(0,0.25,0.5,1))***/",
                              "/***adaptive_resampling,ilike::ess(0.5)***/"))
  named = registry_smc(c("/***smc_sequence,ilike::annealing(schedule=s,variable=power)***/",
                         "/***adaptive_resampling,ilike::ess(proportion=half)***/"),
                       list(s = c(0, 0.25, 0.5, 1), half = 0.5))
  expect_identical(positional, named)
})

test_that("mistakes in SMC settings are reported", {
  skip_on_cran()
  # before the registry, an unknown sequence was taken to be a parameter sequence
  expect_error(registry_smc("/***smc_sequence,ilike::anealing(power,c(0,1))***/"),
               "Unknown SMC sequence ilike::anealing")
  expect_error(registry_smc(c("/***smc_sequence,ilike::annealing(power,c(0,1))***/",
                              "/***adaptive_resampling,ilike::cess(0.5)***/")),
               "ilike::cess cannot be used as an adaptive resampling method")
  expect_error(registry_smc(c("/***smc_sequence,ilike::annealing(power,c(0,1))***/",
                              "/***smc_termination,ilike::stable(0.1,0)***/")),
               "argument consecutive: must be at least 1")
  expect_error(registry_smc("/***smc_sequence,ilike::annealing(power,c(0,1),c(1),2)***/"),
               "too many arguments")
})

test_that("the SMC built-ins are in the table", {
  tab = ilike:::builtin_table()
  expect_true(all(c("ess", "cess", "positive", "stable", "always", "stochastic", "sqrt", "adjustment",
                    "annealing", "tempering", "parameter") %in% tab$builtin))
  expect_equal(tab$roles[tab$builtin == "ess"][1], "adaptive_resampling, adaptive_target")
})

# Factor kinds (I2c): simulation-based and algorithmic likelihoods
registry_is = function(fixture, from, to, points = 200)
{
  lines = readLines(system.file("extdata", fixture, package = "ilike"))
  f = tempfile(fileext = ".ilike")
  writeLines(sub(from, to, lines, fixed = TRUE), f)
  name = basename(tempfile("registry_is_"))
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  ilike::IS(ilike::compile(f), results_name = name, number_of_importance_points = points, seed = 1)
}

test_that("simulation-based and algorithmic likelihoods take keyword arguments", {
  skip_on_cran()
  expect_identical(registry_is("test_sbi_sl.ilike", "ilike::sl(100)", "ilike::sl(100)"),
                   registry_is("test_sbi_sl.ilike", "ilike::sl(100)", "ilike::sl(store_output=FALSE,number_of_points=100)"))
  inner = system.file("extdata", "test_algorithmic_inner_ssm.ilike", package = "ilike")
  pf_call = 'ilike::pf("test_algorithmic_inner_ssm.ilike",500)'
  expect_identical(registry_is("test_algorithmic_pf.ilike", pf_call, paste0('ilike::pf("', inner, '",100)'), 50),
                   registry_is("test_algorithmic_pf.ilike", pf_call, paste0('ilike::pf(particles=100,model="', inner, '")'), 50))
})

test_that("mistakes in simulation-based and algorithmic likelihoods are reported", {
  skip_on_cran()
  expect_error(registry_is("test_sbi_sl.ilike", "ilike::sl(100)", "ilike::sll(100)", 10),
               "Unknown simulation-based likelihood ilike::sll")
  expect_error(registry_is("test_sbi_sl.ilike", "ilike::sl(100)", "ilike::sl(100,unbiassed=TRUE)", 10),
               "no argument called unbiassed")
  inner = system.file("extdata", "test_algorithmic_inner_ssm.ilike", package = "ilike")
  expect_error(registry_is("test_algorithmic_pf.ilike", 'ilike::pf("test_algorithmic_inner_ssm.ilike",500)',
                           paste0('ilike::pf("', inner, '")'), 10),
               "the argument particles is missing")
})

test_that("every built-in is in the table", {
  tab = ilike:::builtin_table()
  expect_true(all(c("euclidean_uniform_abc", "gaussian_abc", "lp_uniform_abc", "ienki_abc", "sl", "kf", "enkf", "pf",
                    "ssm", "linear_gaussian_data_model", "linear_gaussian_transition_model") %in% tab$builtin))
  expect_equal(tab$argument[tab$builtin == "pf"], c("model", "particles", "iterations_to_store", "parallel", "grain_size"))
  expect_equal(length(unique(tab$builtin)), 66) # (hmc, nuts, esjd, normal_block, tnorm, tnorm_rw, dual_averaging, pretune, dirichlet_block, normal_invgamma_block, data, is, smc, pcn, adaptive_iterations, slice and pcn_mala added, C208, C210, C211, C216, C217, C219, C223, C228, C241, C249, C250, C252)
})

test_that("EnKI reads enk_sequence blocks, which allow only ilike::annealing", {
  skip_on_cran()
  model = c("/***prior,ilike::norm(theta,0,1)***/", "",
            "/***data***/", "Data set_data() { return Data(\"y\", 1.0); }", "",
            "/***linear_gaussian_data_model,ilike::linear_gaussian_data_model(y,theta,1,1)***/", "")
  run = function(lines) {
    f = tempfile(fileext = ".ilike"); writeLines(c(model, lines), f)
    name = basename(tempfile("registry_enki_"))
    on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
    z = ilike::IENKI(ilike::compile(f), results_name = name, number_of_ensemble_members = 200, seed = 1)
    list(z = z, targets = sort(unique(ilike::load_enk_output(name)$TargetParameters)))
  }
  # the default sequence: from the prior to the posterior in one step
  expect_equal(run(character(0))$targets, "power=1;")
  expect_equal(run("/***enk_sequence,ilike::annealing(power,c(0,0.25,0.5,1))***/")$targets,
               c("power=0.25;", "power=0.5;", "power=1;"))
  expect_error(run("/***enk_sequence,ilike::tempering(power,c(0,1))***/"),
               "ilike::tempering cannot be used as an EnK sequence")
})
