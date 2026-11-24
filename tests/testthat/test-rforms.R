# Models written with R expressions in the block headers (documented forms).
# These used to fail: "=" in named arguments (log=TRUE) or in "==" was taken as
# an assignment, strings and other non-call expressions were rejected, vector
# results were refused ("Not a matrix"), and the type of R-expression blocks
# (constant or depending on parameters) never reached the model.

test_that("R expressions for priors, likelihoods and gradients give the right posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("rforms_is", recursive=TRUE)
    unlink("rforms_langevin", recursive=TRUE)
  }, add=TRUE)

  z = ilike::IS(ilike::compile(rforms_files("test_rforms_gaussian_posterior.ilike")),
                results_name = "rforms_is", number_of_importance_points = 5000, seed = 1)
  expect_equal(z, dnorm(1,1,sqrt(2),log=TRUE), tolerance = 0.02)

  v = run_rforms_mcmc("test_rforms_langevin.ilike", "rforms_langevin")
  expect_equal(mean(v), 1, tolerance = 0.1)
  expect_equal(var(v), 0.5, tolerance = 0.1)
})

test_that("R expressions for MH and independent MH proposals give the right posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("rforms_mh", recursive=TRUE)
    unlink("rforms_imh", recursive=TRUE)
  }, add=TRUE)

  for (case in list(c("test_rforms_mh_proposal.ilike", "rforms_mh"),
                    c("test_rforms_imh_proposal.ilike", "rforms_imh")))
  {
    v = run_rforms_mcmc(case[1], case[2])
    expect_equal(mean(v), 1, tolerance = 0.1)
    expect_equal(var(v), 0.5, tolerance = 0.1)
  }
})

test_that("R expressions for state-space models work in KF, PF and EnKF", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    for (d in c("rforms_kf","rforms_pf_transition","rforms_pf_nonlinear","rforms_enkf_nonlinear"))
      unlink(d, recursive=TRUE)
  }, add=TRUE)

  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(ilike::KF(ilike::compile(rforms_files("test_rforms_ssm_linear.ilike")), results_name = "rforms_kf"),
               truth, tolerance = 1e-6)
  expect_equal(ilike::PF(ilike::compile(rforms_files("test_rforms_ssm_transition.ilike")), results_name = "rforms_pf_transition",
                         number_of_particles = 5000, seed = 1), truth, tolerance = 0.01)
  recipe = ilike::compile(rforms_files("test_rforms_ssm_nonlinear.ilike"))
  expect_equal(ilike::PF(recipe, results_name = "rforms_pf_nonlinear", number_of_particles = 5000, seed = 1),
               truth, tolerance = 0.01)
  expect_equal(ilike::EnKF(recipe, results_name = "rforms_enkf_nonlinear", number_of_ensemble_members = 5000, seed = 1),
               truth, tolerance = 0.01)
})

test_that("R expressions for a simulation-based model give the right evidence", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("rforms_sbi", recursive=TRUE), add=TRUE)

  # Synthetic likelihood of the mean of ten observations: evidence N(1; 1, 1.1).
  z = ilike::IS(ilike::compile(rforms_files("test_rforms_sbi.ilike")),
                results_name = "rforms_sbi", number_of_importance_points = 2000, seed = 1)
  expect_equal(z, dnorm(1,1,sqrt(1.1),log=TRUE), tolerance = 0.05)
})

test_that("a block header of exactly eight characters is accepted and an empty one is not", {
  # compile() used to reject every header whose name had eight characters
  # (such as data,x=1) as if it were the empty header /******/.
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("/***data,x=1***/", "",
               "/***evaluate_log_likelihood,dnorm(data$x,parameters$theta,1,log=TRUE)***/"), f)
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = list(theta = 0.5)),
               dnorm(1, 0.5, 1, log = TRUE), tolerance = 1e-8)
  writeLines(c("/******/", "/***data,x=1***/"), f)
  expect_error(ilike::compile(f), "needs a name")
})

test_that("variable names may be non-ASCII (in a UTF-8 locale)", {
  # (the other tests use ASCII names, so that they pass in any locale; this
  # one checks a Greek name, written with escapes to keep this file ASCII)
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]))
  theta_name = "θ"
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(enc2utf8(c("ilike 1.0", "data { y = [1.2] }", "model {", paste0("  ", theta_name, " ~ norm(0, 1)"),
                        paste0("  y ~ norm(", theta_name, ", 1)"), "}")), f, useBytes = TRUE)
  parameters = list(0.5)
  names(parameters) = theta_name
  expect_equal(ilike::evaluate_log(ilike::compile(f), parameters = parameters),
               dnorm(0.5, log = TRUE) + dnorm(1.2, 0.5, 1, log = TRUE), tolerance = 1e-10)
})
