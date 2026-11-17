# MFDS (Huang, Huang, Reich & Stuart, arXiv-2204.04386; dev/plan.md 3.16): for
# linear-Gaussian problems the ensemble converges to the posterior.

mfds_model = c("/***data***/",
               "Data set_data() { return Data(\"y\", 1.0); }", "",
               "/***simulate_data_model***/",
               "Data simulate_data_model(RandomNumberGenerator &rng, const Parameters &parameters)",
               "{ return Data(\"y\", rnorm(rng, parameters[\"theta\"][0], 1.0)); }", "",
               "/***data_variable***/",
               "std::string data_variable() { return \"y\"; }", "")

mfds_run = function(lines, iterations = 15, members = 2000, ...)
{
  f = tempfile(fileext = ".ilike")
  writeLines(c(lines, mfds_model), f)
  name = basename(tempfile("mfds_"))
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  ilike::MFDS(ilike::compile(f), results_name = name, number_of_ensemble_members = members,
              number_of_iterations = iterations, seed = 1, ...)
  o = ilike::load_enk_output(name)
  o[o$Target == max(o$Target), ]
}

test_that("MFDS converges to the posterior of a conjugate model, with either update", {
  skip_on_cran()
  # theta ~ N(0,1), y = 1 ~ N(theta,1): the posterior is N(0.5, 0.5)
  for (shifter in c("", "/***enk_shifter,ilike::adjustment()***/"))
  {
    final = mfds_run(c("/***prior,ilike::norm(theta,0,1)***/", "", shifter, ""))
    expect_equal(length(unique(final$Particle)), 2000)
    expect_equal(mean(final$Value), 0.5, tolerance = 0.1)
    expect_equal(var(final$Value), 0.5, tolerance = 0.1)
  }
})

test_that("MFDS refuses models and settings it cannot use", {
  skip_on_cran()
  expect_error(mfds_run("/***prior,ilike::gamma(theta,1,1)***/", iterations = 2, members = 10),
               "MFDS requires a Gaussian prior")
  expect_error(mfds_run(c("/***prior,ilike::norm(theta,0,1)***/", "/***enk_shifter,ilike::sqrt()***/"),
                        iterations = 2, members = 10),
               "ilike::sqrt is not available for MFDS")
  expect_error(mfds_run("/***prior,ilike::norm(theta,0,1)***/", iterations = 2, members = 10, delta_t = 1),
               "delta_t must be a number between 0 and 1")
})
