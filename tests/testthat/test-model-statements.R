# Model statements (dev/ilike-language.md 4.4, I4a): loops, indexed targets and
# deterministic nodes compiled to generated C++ factors. The log densities are
# compared with ones written out by hand in R; the prior is checked through its
# draws (IS with the prior as proposal stores draws from the prior).

model_file = function(f) system.file("extdata", f, package = "ilike")

model_is_draws = function(recipe, name)
{
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::IS(recipe, results_name = name, number_of_importance_points = 10000, seed = 1)
  output = ilike::load_smc_output(name)
  function(variable, dimension = 1)
    output$Value[output$ParameterName == variable & output$Dimension == dimension]
}

test_that("generated log densities equal the ones written by hand", {
  skip_on_cran()
  # latent AR(1) with a loop over time referring to the previous element
  y = c(0.5, 1.0, 0.2, -0.3, 0.8); phi = 0.6; s = 0.8; x = c(0.3, 0.5, 0.1, -0.2, 0.6)
  exact = dunif(phi, -1, 1, log = TRUE) + dgamma(s, 2, 2, log = TRUE) + dnorm(x[1], 0, 1, log = TRUE) +
    sum(dnorm(x[2:5], phi * x[1:4], s, log = TRUE)) + sum(dnorm(y, x, 1, log = TRUE))
  recipe = ilike::compile(model_file("lang_model_ar1.ilike"))
  expect_equal(ilike::evaluate_log(recipe, parameters = list(phi = phi, s = s, x = x)), exact, tolerance = 1e-10)

  # regression with a data-indexed group effect and a deterministic node
  yy = c(1.2, 0.4, 2.2, 1.9, 3.1, 0.2); xx = c(0.1, -0.5, 1.0, 0.7, 1.5, -1.0); g = c(1, 1, 2, 2, 3, 3)
  alpha = c(0.5, 1.5, 2.0); b = 0.7; sigma = 0.6
  exact = sum(dnorm(alpha, 0, 2, log = TRUE)) + dnorm(b, 0, 1, log = TRUE) + dlnorm(sigma, 0, 0.5, log = TRUE) +
    sum(dnorm(yy, alpha[g] + b * xx, sigma, log = TRUE))
  recipe = ilike::compile(model_file("lang_model_regression.ilike"))
  expect_equal(ilike::evaluate_log(recipe, parameters = list(alpha = alpha, b = b, sigma = sigma)), exact, tolerance = 1e-10)

  # beta and multivariate normal priors, Bernoulli data through inv_logit
  z = c(1, 0, 1, 1, 0, 1, 1, 1); p = 0.65; w = c(0.3, -1.2); S = matrix(c(1, 0.5, 0.5, 2), 2)
  q = plogis(qlogis(p) + 0.1 * w[1])
  log_mvnorm = -0.5 * log(det(2 * pi * S)) - 0.5 * drop(t(w) %*% solve(S) %*% w)
  exact = dbeta(p, 2, 3, log = TRUE) + log_mvnorm + sum(dbinom(z, 1, q, log = TRUE))
  recipe = ilike::compile(model_file("lang_model_beta_bern.ilike"))
  expect_equal(ilike::evaluate_log(recipe, parameters = list(p = p, w = w)), exact, tolerance = 1e-10)
})

test_that("the generated prior simulates the parameters in dependency order", {
  skip_on_cran()
  draws = model_is_draws(ilike::compile(model_file("lang_model_ar1.ilike")), "model_ar1_is")
  expect_equal(mean(draws("s")), 1, tolerance = 0.05)          # Gamma(2, rate 2)
  expect_equal(var(draws("x", 1)), 1, tolerance = 0.05)
  # x[2] = phi x[1] + noise: var(x[2]) = E[phi^2] var(x[1]) + E[s^2] = 1/3 + 1.5
  expect_equal(var(draws("x", 2)), 1 / 3 + 1.5, tolerance = 0.05)
  draws = model_is_draws(ilike::compile(model_file("lang_model_beta_bern.ilike")), "model_bb_is")
  expect_equal(mean(draws("p")), 0.4, tolerance = 0.05)        # Beta(2, 3)
  expect_equal(var(draws("w", 2)), 2, tolerance = 0.05)
  expect_equal(cov(draws("w", 1), draws("w", 2)), 0.5, tolerance = 0.1)
})

test_that("MCMC on a hierarchical model written as statements finds the exact posterior", {
  skip_on_cran()
  # eight schools with tau known: mu | y is normal with precision 1/100 + sum 1/(tau^2 + s^2).
  # The fixture runs 100,000 iterations (0.5 s): a random walk moving all eight
  # theta at once mixes slowly, and 10,000 gave sds 15% out depending on the
  # seed (element moves, I4c, will fix this).
  y = c(28, 8, -3, 7, -1, 1, 18, 12); s = c(15, 10, 16, 11, 9, 11, 10, 18); tau = 5
  precision = 1 / 100 + sum(1 / (tau^2 + s^2))
  exact_mean = sum(y / (tau^2 + s^2)) / precision
  name = "model_schools"
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(model_file("lang_model_schools.ilike")), results_name = name,
              initial_values = list(list(mu = 0, theta = rep(0, 8))), seed = 1)
  output = ilike::load_mcmc_output(name)
  mu = output$Value[output$ParameterName == "mu" & output$Iteration > 2000]
  expect_equal(mean(mu), exact_mean, tolerance = 0.1)
  expect_equal(sd(mu), sqrt(1 / precision), tolerance = 0.1)
})

test_that("model statements give clear errors", {
  bad = function(text)
  {
    f = tempfile(fileext = ".ilike")
    on.exit(unlink(f))
    writeLines(text, f)
    ilike::compile(f)
  }
  expect_error(bad(c("ilike 1.0", "model {", "  a ~ norm(0, 1)", "  y ~ wibble(a, 1)", "}")),
               "4:7: unknown distribution wibble in a model statement")
  expect_error(bad(c("ilike 1.0", "model {", "  a ~ norm(b, 1)", "  b ~ norm(a, 1)", "}")),
               "the prior's statements depend on each other in a cycle")
  expect_error(bad(c("ilike 1.0", "model { a ~ norm(0, 1) }", "algorithm mcmc {",
                     "  for (i in 1:2) m = metropolis(norm_rw(a, 1))", "}")),
               "4:3: loops \\(for\\) are only allowed in the model section, not in an algorithm section")
})
