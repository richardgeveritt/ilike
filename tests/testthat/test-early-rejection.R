# Early rejection in Metropolis-Hastings with estimators that stop below a
# threshold (dev/plan.md 3.21; C245): ABC with a uniform kernel

early_rejection_abc = function(points, step, iterations)
{
  c("ilike 1.0", "data { y = [1.3, 0.2, 0.9, 1.7, 0.4] }", "model {", "  theta ~ norm(0, 3)",
    paste0("  factor sim = sbi(euclidean_uniform_abc(number_of_points = ", points, ", tolerance_variable = epsilon, tolerance = 0.5),"),
    "                   simulate = cpp::simulate_y, data = y, summaries = cpp::summaries)", "}", "cpp {",
    "static double early_rejection_count = 0.0;",
    "// [[Rcpp::export]]", "double early_rejection_simulations() { return early_rejection_count; }",
    "Data simulate_y(RandomNumberGenerator &rng, const Parameters &parameters)",
    "{ early_rejection_count += 1.0; arma::colvec y(5); for (arma::uword i = 0; i < 5; ++i) y[i] = rnorm(rng, parameters[\"theta\"][0], 1.0); return Data(\"y\", y); }",
    "Data summaries(const Data &data) { arma::colvec s(1); s[0] = arma::mean(arma::vectorise(data[\"y\"])); return Data(\"y\", s); }",
    "}", "algorithm mcmc {", paste0("  move m = metropolis(norm_rw(theta, ", step, "))"), paste0("  termination = iterations(", iterations, ")"), "}")
}

test_that("ABC-MCMC with a uniform kernel stops simulating for proposals that are certain to be rejected", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(early_rejection_abc(200, 3, 500), f)
  m = ilike::compile(f)
  before = early_rejection_simulations()
  ilike::MCMC(m, seed = 1, initial_values = list(list(theta = 0.9)))
  used = early_rejection_simulations() - before
  # (each iteration would use 200 simulations; about a third are saved here)
  expect_lt(used, 0.85 * 500 * 200)
})

test_that("ABC-MCMC with early stopping targets the ABC posterior", {
  skip_on_cran()
  y = c(1.3, 0.2, 0.9, 1.7, 0.4); eps = 0.5; sd = 1/sqrt(5)
  grid = seq(-5, 6, length.out = 20001)
  p = dnorm(grid, 0, 3) * (pnorm((mean(y) + eps - grid)/sd) - pnorm((mean(y) - eps - grid)/sd))
  p = p / sum(p)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "early_rejection_abc"), recursive = TRUE), add = TRUE)
  writeLines(early_rejection_abc(200, 3, 20000), f)
  ilike::MCMC(ilike::compile(f), results_name = "early_rejection_abc", seed = 1, initial_values = list(list(theta = 0.9)))
  o = ilike::load_mcmc_output("early_rejection_abc")
  b = o$Value[o$ParameterName == "theta"][-(1:1000)]
  expect_equal(mean(b), sum(p * grid), tolerance = 0.05)
  expect_equal(sd(b), sqrt(sum(p * grid^2) - sum(p * grid)^2), tolerance = 0.1)
})

test_that("a particle filter with a bounded measurement density stops early and targets the posterior", {
  # (C246) density(log = …, maximum = …) in the inner state-space model; the
  # posterior of a against the one from the Kalman filter's exact likelihood
  skip_on_cran()
  set.seed(4); n = 50; x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  d = tempfile(); dir.create(d)
  on.exit(unlink(c(d, "early_rejection_pf", "early_rejection_kf"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  x ~ norm(0, 1)",
               "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)",
               "  factor obs = density(log = cpp::meas, maximum = cpp::meas_max)",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }",
               "double meas(const Parameters &parameters, const Data &data) { return dnorm(data[\"y\"][0], parameters[\"x\"][0], 1.0); }",
               "double meas_max(const Parameters &parameters, const Data &data) { return -0.5*std::log(2.0*M_PI); }",
               "}", "algorithm filter {", "  particles = 200", "}"), file.path(d, "inner.ilike"))
  writeLines(c("ilike 1.0", "model {", "  x ~ norm(0, 1)",
               "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)",
               "  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }",
               "}", "algorithm kalman {", "}"), file.path(d, "inner_kf.ilike"))
  outer = function(import, factor) c("ilike 1.0", paste0("import \"", import, "\" as inner"),
                                     paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  a ~ unif(-1, 1)",
                                     paste0("  factor filter = ", factor), "}", "algorithm mcmc {",
                                     "  move m = metropolis(norm_rw(a, 0.3))", "  termination = iterations(4000)", "}")
  writeLines(outer("inner.ilike", "pf(inner)"), file.path(d, "pf.ilike"))
  writeLines(outer("inner_kf.ilike", "kf(inner)"), file.path(d, "kf.ilike"))
  ilike::MCMC(ilike::compile(file.path(d, "pf.ilike")), results_name = "early_rejection_pf", seed = 1, initial_values = list(list(a = 0.5)))
  ilike::MCMC(ilike::compile(file.path(d, "kf.ilike")), results_name = "early_rejection_kf", seed = 1, initial_values = list(list(a = 0.5)))
  a_pf = ilike::load_mcmc_output("early_rejection_pf"); a_pf = a_pf$Value[a_pf$ParameterName == "a"][-(1:200)]
  a_kf = ilike::load_mcmc_output("early_rejection_kf"); a_kf = a_kf$Value[a_kf$ParameterName == "a"][-(1:200)]
  expect_equal(mean(a_pf), mean(a_kf), tolerance = 0.05)
  expect_equal(sd(a_pf), sd(a_kf), tolerance = 0.15)
})
