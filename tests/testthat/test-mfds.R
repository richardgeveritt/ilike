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

# MFDS from the ilike language (dev/plan.md 3.16): a simulator alone (sbi
# without a method), likelihoods from model statements and Gaussian data
# models, Gaussian and log-normal priors, and a transformation of the user's.

mfds_language = function(lines, members = 2000, iterations = 15)
{
  f = tempfile(fileext = ".ilike")
  writeLines(c("ilike 1.0", lines), f)
  name = basename(tempfile("mfds_language_"))
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  ilike::MFDS(ilike::compile(f), results_name = name, number_of_ensemble_members = members,
              number_of_iterations = iterations, seed = 1)
  o = ilike::load_enk_output(name)
  o[o$Target == max(o$Target), ]
}

lognormal_simulator = c("cpp {",
  "Data simulate_y(RandomNumberGenerator &rng, const Parameters &parameters)",
  "{",
  "  arma::colvec y(3);",
  "  for (int i = 0; i < 3; ++i) y[i] = rnorm(rng, std::log(parameters[\"theta\"][0]), 1.0);",
  "  return Data(\"y\", y);",
  "}",
  "Parameters to_positive(const Parameters &parameters) { return Parameters(\"theta\", arma::exp(parameters[\"theta\"])); }",
  "}")

test_that("MFDS runs from the language: a simulator alone, statements, Gaussian data models", {
  skip_on_cran()
  y = c(0.3, 1.4, 0.9)
  # theta ~ N(0, 1), y_i ~ N(theta, 1): the posterior is N(sum(y)/4, 1/4)
  statements = mfds_language(c("data { y = [0.3, 1.4, 0.9] }",
                               "model { theta ~ norm(0, 1); for (i in 1:3) y[i] ~ norm(theta, 1) }"))
  expect_equal(mean(statements$Value), sum(y)/4, tolerance = 0.05)
  expect_equal(var(statements$Value), 1/4, tolerance = 0.1)
  # the same model with its likelihood as a simulator alone
  simulator = mfds_language(c("data { y = [0.3, 1.4, 0.9] }",
                              "model { theta ~ norm(0, 1); factor sim = sbi(simulate = cpp::simulate_y, data = y) }",
                              "cpp {", "Data simulate_y(RandomNumberGenerator &rng, const Parameters &parameters)",
                              "{ arma::colvec y(3); for (int i = 0; i < 3; ++i) y[i] = rnorm(rng, parameters[\"theta\"][0], 1.0);",
                              "  return Data(\"y\", y); }", "}"))
  expect_identical(simulator$Value, statements$Value)
  # a linear-Gaussian data model with a correlated prior
  linear = mfds_language(c("data { y = [1.0, -0.5] }", "model {",
                           "  theta ~ mvnorm([0, 0], [[1, 0.3], [0.3, 1]])",
                           "  factor obs = linear_gaussian(data = y, state = theta, matrix = [[1, 0.5], [0, 1]], covariance = [[0.5, 0], [0, 0.8]])",
                           "}"), members = 4000, iterations = 20)
  P0 = matrix(c(1, 0.3, 0.3, 1), 2); H = matrix(c(1, 0, 0.5, 1), 2); R = diag(c(0.5, 0.8))
  P = solve(solve(P0) + t(H) %*% solve(R) %*% H)
  m = P %*% t(H) %*% solve(R) %*% c(1, -0.5)
  expect_equal(as.numeric(tapply(linear$Value, linear$Dimension, mean)), as.numeric(m), tolerance = 0.1)
  # statements and a Gaussian data model together
  joint = mfds_language(c("data { y = [0.3, 1.4, 0.9]; z = 0.5 }", "model {", "  theta ~ norm(0, 1)",
                          "  for (i in 1:3) y[i] ~ norm(theta, 1)",
                          "  factor obs = linear_gaussian(data = z, state = theta, matrix = 2, covariance = 0.25)", "}"),
                        members = 4000, iterations = 20)
  precision = 1 + 3 + 4/0.25
  expect_equal(mean(joint$Value), (sum(y) + 2*0.5/0.25)/precision, tolerance = 0.05)
  expect_equal(var(joint$Value), 1/precision, tolerance = 0.1)
})

test_that("MFDS moves log-normal priors and the user's transformation on the Gaussian scale", {
  skip_on_cran()
  y = c(0.3, 1.4, 0.9)
  # theta ~ logN(0, 1), y_i ~ N(log theta, 1): log theta's posterior is N(sum(y)/4, 1/4)
  lognormal = mfds_language(c("data { y = [0.3, 1.4, 0.9] }",
                              "model { theta ~ lnorm(0, 1); for (i in 1:3) y[i] ~ norm(log(theta), 1) }"))
  expect_true(all(lognormal$Value > 0))
  expect_equal(mean(log(lognormal$Value)), sum(y)/4, tolerance = 0.05)
  expect_equal(var(log(lognormal$Value)), 1/4, tolerance = 0.1)
  # a simple log-normal prior with a simulator, and a normal prior with the
  # user's transformation exp: the same ensembles on the Gaussian scale
  simple = mfds_language(c("data { y = [0.3, 1.4, 0.9] }",
                           "model { theta ~ lnorm(0, 1); factor sim = sbi(simulate = cpp::simulate_y, data = y) }",
                           lognormal_simulator))
  expect_equal(simple$Value, lognormal$Value)
  transformed = mfds_language(c("data { y = [0.3, 1.4, 0.9] }",
                                "model { theta ~ norm(0, 1); factor sim = sbi(simulate = cpp::simulate_y, data = y) }",
                                lognormal_simulator,
                                "algorithm mfds { transform = cpp::to_positive }"))
  expect_equal(transformed$Value, lognormal$Value)
})

test_that("MFDS says why a model written in the language cannot be used", {
  skip_on_cran()
  expect_error(mfds_language(c("data { y = [0.3, 1.4, 0.9] }", "model {", "  mu ~ norm(0, 1)",
                               "  for (i in 1:3) theta[i] ~ norm(mu, 1)", "  for (i in 1:3) y[i] ~ norm(theta[i], 1)", "}"),
                             members = 10, iterations = 1),
               "the prior of theta depends on mu")
  expect_error(mfds_language(c("model { theta ~ norm(0, 1); factor pot += -theta^2 }"), members = 10, iterations = 1),
               "the factor pot is made of \\+= terms")
  expect_error(mfds_language(c("data { y = [0.3, 1.4, 0.9] }",
                               "model { theta ~ lnorm(0, 1); factor sim = sbi(simulate = cpp::simulate_y, data = y) }",
                               lognormal_simulator, "algorithm mfds { transform = cpp::to_positive }"),
                             members = 10, iterations = 1),
               "the prior gives the Gaussian on the scale the ensemble moves on") # (ASCII only)
  # other algorithms ask for a method
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [0.3, 1.4, 0.9] }",
               "model { theta ~ norm(0, 1); factor sim = sbi(simulate = cpp::simulate_y, data = y) }",
               lognormal_simulator, "algorithm mcmc { move m = metropolis(norm_rw(theta, 1)); termination = iterations(10) }"), f)
  expect_error(ilike::MCMC(ilike::compile(f), results_name = "mfds_mcmc_refused", initial_values = list(list(theta = 1))),
               "is a simulator without a method")
  unlink("mfds_mcmc_refused", recursive = TRUE)
})

test_that("MFDS uses deterministic nodes in the prior and the likelihood", {
  skip_on_cran()
  # lang_model_nodes.ilike: alpha[k] ~ N(s[k], 1), y[k] ~ N(alpha[k] - m[k], 1),
  # m = 2 w, s = m + min(1, w): alpha[k]'s posterior is N((s + y + m)/2, 1/2)
  name = basename(tempfile("mfds_nodes_"))
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MFDS(ilike::compile(system.file("extdata", "lang_model_nodes.ilike", package = "ilike")), results_name = name,
              number_of_ensemble_members = 4000, number_of_iterations = 20, seed = 1)
  o = ilike::load_enk_output(name)
  final = o[o$Target == max(o$Target), ]
  w = 1:3; m = 2 * w; s = m + pmin(1, w); y = c(0.5, 1.5, 2.5)
  expect_equal(as.numeric(tapply(final$Value, final$Dimension, mean)), (s + y + m)/2, tolerance = 0.05)
  expect_equal(as.numeric(tapply(final$Value, final$Dimension, var)), rep(0.5, 3), tolerance = 0.15)
})
