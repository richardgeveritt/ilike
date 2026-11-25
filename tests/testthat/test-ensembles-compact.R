# (C279-C281) The ensemble Kalman methods keep their members as rows of
# matrices when they can (compact ensembles), and as member objects
# otherwise; ILIKE_RICH_ENSEMBLES (a developer switch) forces member objects.
# The two give the same numbers: the same draws in the same order. Checked on
# short runs of the filter (two measurement factors, square-root shifter),
# SMC^2 with an ensemble Kalman filter (outer particles copied when
# resampled: their copies once shared their members' states), EnKI and MFDS
# (with a transform and given initial points).

both_ways = function(run)
{
  compact = run()
  Sys.setenv(ILIKE_RICH_ENSEMBLES = "1")
  on.exit(Sys.unsetenv("ILIKE_RICH_ENSEMBLES"), add = TRUE)
  rich = run()
  list(compact = compact, rich = rich)
}

test_that("compact ensembles give the same runs as member objects", {
  src = readLines("test-ssm.R")
  eval(parse(text = src[grep("^enkf_check_files", src):(grep("^test_that\\(\"EnKF runs a data covariance", src) - 1)]))
  d = enkf_check_files(TRUE, TRUE, 100, "sqrt")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  z = both_ways(function() enkf_check(d)[["enkf"]])
  expect_identical(z$compact, z$rich)

  src2 = readLines("test-smc2.R")
  eval(parse(text = src2[grep("^smc2_files", src2):(grep("^test_that", src2)[1]-1)]))
  d2 = smc2_files(6, 85, 10, pf = "enkf(inner, ensemble_members = 30, shifter = \"stochastic\")", inner_file = "inner_kf.ilike")
  on.exit(unlink(d2, recursive = TRUE), add = TRUE)
  r2 = suppressMessages(ilike::compile(file.path(d2, "smc2.ilike")))
  z = both_ways(function() suppressMessages(ilike::SMC_with_MCMC(r2, seed = 1, grain_size = 1000)))
  expect_identical(z$compact, z$rich)

  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("/***prior,ilike::mvnorm(theta,c(0,0),diag(2))***/", "", "/***data***/",
               "Data set_data() { arma::colvec y(2); y(0) = 1.0; y(1) = -0.5; return Data(\"y\", y); }", "",
               "/***linear_gaussian_data_model,ilike::linear_gaussian_data_model(y,theta,diag(2),diag(2))***/", "",
               "/***enk_sequence,ilike::annealing(power,c(0,0.25,0.5,1))***/"), f)
  r3 = suppressMessages(ilike::compile(f))
  z = both_ways(function() suppressMessages(ilike::IENKI(r3, number_of_ensemble_members = 200, seed = 1)))
  expect_identical(z$compact, z$rich)

  g = tempfile(fileext = ".ilike")
  on.exit(unlink(g), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.3, 1.4, 0.9] }",
               "model { theta ~ norm(0, 1); factor sim = sbi(simulate = cpp::simulate_y, data = y) }", "cpp {",
               "Data simulate_y(RandomNumberGenerator &rng, const Parameters &parameters)",
               "{ arma::colvec y(3); for (int i = 0; i < 3; ++i) y[i] = rnorm(rng, std::log(parameters[\"theta\"][0]), 1.0); return Data(\"y\", y); }",
               "Parameters to_positive(const Parameters &parameters) { return Parameters(\"theta\", arma::exp(parameters[\"theta\"])); }",
               "}", "algorithm mfds { transform = cpp::to_positive }"), g)
  r4 = suppressMessages(ilike::compile(g))
  mfds = function() {
    name = basename(tempfile("mfds_compact_"))
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    suppressMessages(ilike::MFDS(r4, results_name = name, number_of_ensemble_members = 100, number_of_iterations = 4, seed = 1))
    ilike::load_enk_output(name)$Value
  }
  z = both_ways(mfds)
  expect_identical(z$compact, z$rich)
})
