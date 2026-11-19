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
               "a loop in an algorithm holds one move \\(move name = metropolis")
})

# Plates (I4b): each factor from model statements keeps its terms' values, and
# a move recomputes only the terms reading an element it changed; a move on
# known variables evaluates only the factors that may read them. The factor's
# value is always the sum of all its terms in order, so runs are identical to
# evaluating each factor whole (options(ilike.plate_factors = FALSE), I4a).

model_mcmc_values = function(file, name, initial_values, ...)
{
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(model_file(file)), results_name = name, initial_values = initial_values, seed = 2, ...)
  ilike::load_mcmc_output(name)$Value
}

without_plates = function(code)
{
  old = options(ilike.plate_factors = FALSE)
  on.exit(options(old), add = TRUE)
  code
}

test_that("plates and automatic move indices give the same runs as whole factors", {
  skip_on_cran()
  init = list(list(mu = 0, theta = rep(0, 8)))
  with_plates = model_mcmc_values("lang_model_schools.ilike", "plates_a", init)
  expect_identical(with_plates, without_plates(model_mcmc_values("lang_model_schools.ilike", "plates_b", init)))
  # the move on mu reads only the prior: its automatic index is factor 1
  # (model_prior), leaving out the likelihood of y
  f = model_file("lang_model_schools.ilike")
  lowered = ilike:::ilike_language_lower(paste(readLines(f), collapse = "\n"), f)
  expect_match(lowered$generated_code, "NumericVector::create\\(1\\);")
})

test_that("partial evaluations are right: elements, matrices and stochastic indices", {
  skip_on_cran()
  # a proposal changing one element of a matrix parameter (so only the terms
  # reading that element are recomputed), and an index that is a parameter
  # (alpha read whole); check_factor_index recomputes every term after each
  # partial evaluation and stops if one that was kept changed
  init = list(list(alpha = c(0, 0), z = 2, x = matrix(0, 3, 2)))
  checked = model_mcmc_values("lang_model_stress.ilike", "plates_c", init, check_factor_index = TRUE)
  expect_identical(checked, model_mcmc_values("lang_model_stress.ilike", "plates_d", init))
  expect_identical(checked, without_plates(model_mcmc_values("lang_model_stress.ilike", "plates_e", init)))
  y = c(0.3, -0.2, 1.1, 0.8, 2.0, 1.4); v = matrix(c(0.1, 0.3, 1.2, 0.4, -0.5, 0.9), 3)
  alpha = c(0.5, 1.5); z = 2.4; x = matrix(c(0.2, 0.1, 0.9, 1.1, -0.3, 0.8), 3)
  exact = sum(dnorm(alpha, 0, 2, log = TRUE)) + dunif(z, 1, 3, log = TRUE) +
    sum(dnorm(x, matrix(alpha, 3, 2, byrow = TRUE), 1, log = TRUE)) + sum(dnorm(v, x, 0.5, log = TRUE)) +
    sum(dnorm(y, alpha[2], 1, log = TRUE))
  recipe = ilike::compile(model_file("lang_model_stress.ilike"))
  expect_equal(ilike::evaluate_log(recipe, parameters = list(alpha = alpha, z = z, x = x)), exact, tolerance = 1e-10)
})

test_that("SMC with annealed plates gives the same run as whole factors", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  lines = readLines(model_file("lang_model_schools.ilike"))
  writeLines(c(lines[seq_len(grep("^algorithm", lines) - 1)],
               "algorithm smc {", "  move m1 = metropolis(norm_rw(mu, 3))",
               "  move m2 = metropolis(mvnorm_rw(theta, diag(8), scale = 4))", "  termination = iterations(3)",
               "  sequence = annealing(power, [0, 0.25, 0.5, 1])", "}"), f)
  run = function(name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = name, number_of_particles = 500, seed = 1)
  }
  # (plates add their terms in a summation tree, I4c-2: the evidence agrees to rounding)
  expect_equal(run("plates_smc_a"), without_plates(run("plates_smc_b")), tolerance = 1e-12)
})

test_that("deterministic nodes read literal data in the prior, and use other nodes", {
  skip_on_cran()
  recipe = ilike::compile(model_file("lang_model_nodes.ilike"))
  alpha = c(3, 5, 8); w = 1:3; y = c(0.5, 1.5, 2.5)
  exact = sum(dnorm(alpha, 2 * w + pmin(1, w), 1, log = TRUE)) + sum(dnorm(y, alpha - 2 * w, 1, log = TRUE))
  expect_equal(ilike::evaluate_log(recipe, parameters = list(alpha = alpha)), exact, tolerance = 1e-10)
  draws = model_is_draws(recipe, "model_nodes_is")
  expect_equal(mean(draws("alpha", 3)), 7, tolerance = 0.02)
})

# Moves on elements (I4c): for (j in 1:J) move m = metropolis(norm_rw(theta[j], s))
# is a scan of Metropolis moves on theta[1], ..., theta[J]; with plates each
# recomputes only the terms reading its element.

test_that("a scan of moves on elements finds the exact posterior with 10,000 iterations", {
  skip_on_cran()
  y = c(28, 8, -3, 7, -1, 1, 18, 12); s = c(15, 10, 16, 11, 9, 11, 10, 18); tau = 5
  precision = 1 / 100 + sum(1 / (tau^2 + s^2))
  init = list(list(mu = 0, theta = rep(0, 8)))
  values = model_mcmc_values("lang_model_schools_elements.ilike", "elements_a", init)
  name = "elements_b"
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(model_file("lang_model_schools_elements.ilike")), results_name = name,
              initial_values = init, seed = 1)
  output = ilike::load_mcmc_output(name)
  mu = output$Value[output$ParameterName == "mu" & output$Iteration > 1000]
  expect_equal(mean(mu), sum(y / (tau^2 + s^2)) / precision, tolerance = 0.1)
  expect_equal(sd(mu), sqrt(1 / precision), tolerance = 0.1)
  # the same run with whole factors, and with every term checked
  expect_identical(values, without_plates(model_mcmc_values("lang_model_schools_elements.ilike", "elements_c", init)))
  expect_identical(values, model_mcmc_values("lang_model_schools_elements.ilike", "elements_d", init,
                                             check_factor_index = TRUE))
})

test_that("moves on elements: matrices, constants given when compiling, and errors", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  write_model = function(algorithm)
    writeLines(c("ilike 1.0", "constants { int n = 2 }",
                 "data { v = [[0.1, 0.4], [0.3, -0.5], [1.2, 0.9]] }",
                 "model {", "  for (i in 1:3) for (j in 1:2) x[i, j] ~ norm(0, 1)",
                 "  for (i in 1:3) for (j in 1:2) v[i, j] ~ norm(x[i, j], 0.5)", "}",
                 "algorithm mcmc {", algorithm, "  termination = iterations(20000)", "}"), f)
  write_model("  for (i in 1:3) for (j in 1:n) move m = metropolis(unif_rw(x[i, j], 1))")
  name = "elements_matrix"
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, initial_values = list(list(x = matrix(0, 3, 2))), seed = 1)
  output = ilike::load_mcmc_output(name)
  # x[i, j] | v is normal with precision 1 + 4 and mean 4 v / 5
  last = output[output$ParameterName == "x" & output$Dimension == 6 & output$Iteration > 2000, "Value"]
  expect_equal(mean(last), 4 * 0.9 / 5, tolerance = 0.1)
  expect_equal(var(last), 1 / 5, tolerance = 0.1)

  # elements are fixed when compiling, from the constants' values then
  lowered = function(...) ilike:::ilike_language_lower(paste(readLines(f), collapse = "\n"), f, ...)
  expect_equal(lowered()$algorithms$mcmc$m_proposal[[1]]$m_elements$columns, c(0, 1, 0, 1, 0, 1))
  expect_equal(lowered(list(n = 1))$algorithms$mcmc$m_proposal[[1]]$m_elements$columns, c(0, 0, 0))

  write_model("  for (j in 1:2) {\n    move a = metropolis(norm_rw(x[1, j], 1))\n    move b = metropolis(norm_rw(x[2, j], 1))\n  }")
  expect_error(ilike::compile(f), "a loop in an algorithm holds one move")
  write_model("  for (j in 1:2) move a = metropolis_hastings(norm_rw(x[1, j], 1))")
  expect_error(ilike::compile(f), "a move on elements is a metropolis move")
  write_model("  for (j in 1:2) move a = metropolis(norm_rw(x, 1))")
  expect_error(ilike::compile(f), "the element is written theta\\[j\\]")
  write_model("  for (j in 1:2) move a = metropolis(norm_rw(v[1, j], 1))")
  expect_error(ilike::compile(f), "v is not a parameter")
})

# In-place moves on elements (I4c-2): the element is changed in the particle,
# the plates recompute the terms reading it, and a rejection restores them;
# the random numbers are those of the moves through proposed particles
# (options(ilike.in_place_moves = FALSE)), so the chains are the same.

through_proposed_particles = function(code)
{
  old = options(ilike.in_place_moves = FALSE)
  on.exit(options(old), add = TRUE)
  code
}

test_that("moves on elements in place give the runs of moves through proposed particles", {
  skip_on_cran()
  init = list(list(mu = 0, theta = rep(0, 8)))
  in_place = model_mcmc_values("lang_model_schools_elements.ilike", "in_place_a", init)
  expect_identical(in_place, through_proposed_particles(
    model_mcmc_values("lang_model_schools_elements.ilike", "in_place_b", init)))
  expect_identical(in_place, model_mcmc_values("lang_model_schools_elements.ilike", "in_place_c", init,
                                               check_factor_index = TRUE))

  # in an annealed SMC (the annealed factors get new outputs, as before)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  lines = readLines(model_file("lang_model_schools_elements.ilike"))
  writeLines(c(lines[seq_len(grep("^algorithm", lines) - 1)],
               "algorithm smc {", "  move m1 = metropolis(norm_rw(mu, 3))",
               "  for (j in 1:J) move m2 = metropolis(norm_rw(theta[j], 6))", "  termination = iterations(3)",
               "  sequence = annealing(power, [0, 0.25, 0.5, 1])", "}"), f)
  run = function(name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = name, number_of_particles = 300, seed = 1)
    ilike::load_smc_output(name)$Value
  }
  expect_identical(run("in_place_smc_a"), through_proposed_particles(run("in_place_smc_b")))
})

# Discrete elements (I4c-3): parameters with bern, binom and pois priors;
# metropolis(flip(s[j])) and gibbs(s[j]), checked against enumeration.

test_that("gibbs and flip moves sample a spin model as enumeration says", {
  skip_on_cran()
  edges = rbind(cbind(c(1, 2, 4, 5, 7, 8), c(2, 3, 5, 6, 8, 9)), cbind(c(1, 2, 3, 4, 5, 6), c(4, 5, 6, 7, 8, 9)))
  states = unname(as.matrix(expand.grid(rep(list(0:1), 9))))
  log_weight = apply(states, 1, function(s) sum(dbinom(s, 1, 0.3, log = TRUE)) +
                       sum(plogis(0.8 * (2 * s[edges[, 1]] - 1) * (2 * s[edges[, 2]] - 1), log.p = TRUE)))
  p = exp(log_weight - max(log_weight))
  p = p / sum(p)
  spins = function(file, name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::MCMC(ilike::compile(file), results_name = name, initial_values = list(list(s = rep(0, 9))), seed = 1)
    output = ilike::load_mcmc_output(name)
    sapply(1:9, function(k) output$Value[output$ParameterName == "s" & output$Dimension == k & output$Iteration > 500])
  }
  gibbs = spins(model_file("lang_model_spins.ilike"), "spins_gibbs")
  expect_true(all(gibbs %in% c(0, 1)))
  expect_equal(colMeans(gibbs), colSums(states * p), tolerance = 0.1)
  expect_equal(mean(gibbs[, 1] * gibbs[, 2]), sum(p * states[, 1] * states[, 2]), tolerance = 0.15)

  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(sub("move m = gibbs\\(s\\[k\\]\\)", "move m = metropolis(flip(s[k]))",
                 readLines(model_file("lang_model_spins.ilike"))), f)
  flip = spins(f, "spins_flip")
  expect_equal(colMeans(flip), colSums(states * p), tolerance = 0.1)

  # the prior of the spins simulates
  draws = model_is_draws(ilike::compile(model_file("lang_model_spins.ilike")), "spins_is")
  expect_equal(mean(draws("s", 1)), 0.3, tolerance = 0.05)
})

test_that("gibbs on a binomial element uses the states of its prior", {
  skip_on_cran()
  name = "binomial_gibbs"
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(model_file("lang_model_binomial.ilike")), results_name = name,
              initial_values = list(list(n = c(0, 0, 0))), seed = 1)
  output = ilike::load_mcmc_output(name)
  y = c(0, 2, 5)
  for (k in 1:3)
  {
    w = dbinom(0:4, 4, 0.3) * dpois(y[k], 0:4 + 0.5)
    v = output$Value[output$ParameterName == "n" & output$Dimension == k & output$Iteration > 500]
    expect_equal(mean(v), sum(w * (0:4)) / sum(w), tolerance = 0.1)
  }
})

test_that("discrete moves give clear errors", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  write_model = function(move)
    writeLines(c("ilike 1.0", "data { y = [0.5, 1.5] }", "model {", "  for (k in 1:2) theta[k] ~ norm(0, 1)",
                 "  for (k in 1:2) y[k] ~ norm(theta[k], 1)", "}", "algorithm mcmc {", move, "}"), f)
  # (a continuous element with a normal full conditional is a Gibbs move since
  # C206; the refusals for other continuous elements are in test-gibbs.R)
  write_model("  for (k in 1:2) move m = gibbs(theta[k])")
  expect_silent(ilike::compile(f))
  write_model("  for (k in 1:2) move m = metropolis(flip(theta[k], 1))")
  expect_error(ilike::compile(f), "flip\\(s\\[j\\]\\) takes only the element")
  write_model("  for (k in 1:2) move m = metropolis_hastings(flip(theta[k]))")
  expect_error(ilike::compile(f), "a move on elements is a metropolis move")
})

# Terms of a named factor (factor name += expression): log terms not tied to an
# observed variable, e.g. an Ising potential or a smoothing penalty.

test_that("an Ising potential written as terms samples as enumeration says", {
  skip_on_cran()
  edges = rbind(cbind(c(1, 2, 4, 5, 7, 8), c(2, 3, 5, 6, 8, 9)), cbind(c(1, 2, 3, 4, 5, 6), c(4, 5, 6, 7, 8, 9)))
  states = unname(as.matrix(expand.grid(rep(list(0:1), 9))))
  log_weight = apply(states, 1, function(s) sum(dbinom(s, 1, 0.3, log = TRUE)) +
                       sum(0.4 * (2 * s[edges[, 1]] - 1) * (2 * s[edges[, 2]] - 1)))
  p = exp(log_weight - max(log_weight))
  p = p / sum(p)
  spins = function(name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::MCMC(ilike::compile(model_file("lang_model_ising.ilike")), results_name = name,
                initial_values = list(list(s = rep(0, 9))), seed = 1)
    output = ilike::load_mcmc_output(name)
    sapply(1:9, function(k) output$Value[output$ParameterName == "s" & output$Dimension == k & output$Iteration > 500])
  }
  gibbs = spins("ising_gibbs")
  expect_equal(colMeans(gibbs), colSums(states * p), tolerance = 0.1)
  expect_equal(mean(gibbs[, 1] * gibbs[, 2]), sum(p * states[, 1] * states[, 2]), tolerance = 0.15)
  # the same chain when the factor is evaluated whole
  old = options(ilike.plate_factors = FALSE)
  on.exit(options(old), add = TRUE)
  expect_identical(spins("ising_whole"), gibbs)
})

test_that("a smoothing term gives the exact Gaussian posterior, by MCMC and SMC", {
  skip_on_cran()
  y = c(1.2, 0.4, -0.3, 0.8, 2.1, 1.7)
  K = 6
  L = matrix(0, K, K)
  for (k in 1:(K - 1)) L[c(k, k + 1), c(k, k + 1)] = L[c(k, k + 1), c(k, k + 1)] + matrix(c(1, -1, -1, 1), 2)
  Q = diag(K) / 100 + 2 * L + diag(K)
  run = function(name)
  {
    on.exit(unlink(name, recursive = TRUE), add = TRUE)
    ilike::MCMC(ilike::compile(model_file("lang_model_smooth.ilike")), results_name = name,
                initial_values = list(list(theta = rep(0, K))), seed = 1)
    output = ilike::load_mcmc_output(name)
    sapply(1:K, function(k) output$Value[output$ParameterName == "theta" & output$Dimension == k & output$Iteration > 500])
  }
  theta = run("smooth_in_place")
  expect_equal(colMeans(theta), solve(Q, y), tolerance = 0.1)
  expect_equal(apply(theta, 2, sd), sqrt(diag(solve(Q))), tolerance = 0.1)
  # identical chains through proposed particles and with whole factors
  old = options(ilike.in_place_moves = FALSE)
  expect_identical(run("smooth_copies"), theta)
  options(old)
  old = options(ilike.plate_factors = FALSE)
  expect_identical(run("smooth_whole"), theta)
  options(old)

  # the terms are annealed with the likelihood
  log_z = -K / 2 * log(2 * pi) - K / 2 * log(100) - 0.5 * as.numeric(determinant(Q)$modulus) -
    0.5 * (sum(y^2) - sum(y * solve(Q, y)))
  lines = readLines(model_file("lang_model_smooth.ilike"))
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(lines[seq_len(grep("^algorithm", lines) - 1)], "algorithm smc {",
               "  for (k in 1:K) move m = metropolis(norm_rw(theta[k], 0.8))", "  termination = iterations(5)",
               "  sequence = annealing(power, [0, 0.01, 0.03, 0.1, 0.3, 1])", "}"), f)
  on.exit(unlink("smooth_smc", recursive = TRUE), add = TRUE)
  z = ilike::SMC_with_MCMC(ilike::compile(f), results_name = "smooth_smc", number_of_particles = 2000, seed = 1)
  expect_equal(z, log_z, tolerance = 0.05)
})

test_that("terms of a named factor give clear errors", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  write_model = function(...)
    writeLines(c("ilike 1.0", "data { y = [0.5, 1.5] }", "model {", "  for (k in 1:2) theta[k] ~ norm(0, 1)", ..., "}"), f)
  write_model("  for (k in 1:2) pot += theta[k]")
  expect_error(ilike::compile(f), "adds a term to a named factor: factor name \\+= expression")
  write_model("  for (k in 1:2) prior pot += theta[k]")
  expect_error(ilike::compile(f), "adds a term to a named factor: factor name \\+= expression")
  write_model("  for (k in 1:2) factor pot[k] += theta[k]")
  expect_error(ilike::compile(f), "a factor's name has no index")
  write_model("  for (k in 1:2) factor y += theta[k]")
  expect_error(ilike::compile(f), "y is data, so it cannot name a factor of terms")
  write_model("  for (k in 1:2) factor theta += theta[k]")
  expect_error(ilike::compile(f), "theta is a variable in the model statements, so it cannot name a factor of terms")
  write_model("  factor pot += theta[1]", "  factor pot = cpp::f", "}",
              "cpp { double f(const Parameters &parameters, const Data &data) { return 0.0; }")
  expect_error(ilike::compile(f), "\"pot\" is defined more than once")
  writeLines(c("ilike 1.0", "data { y += [0.5, 1.5] }", "model { theta ~ norm(0, 1) }"), f)
  expect_error(ilike::compile(f), "\"\\+=\" is only allowed in the model section, not in data")
})

test_that("a constant fixed when compiling cannot be given another value when running", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [1, 2, 3] }", "constants { real s; real step = 1 }",
               "model {", "  mu ~ norm(0, 10)", "  for (i in 1:3) y[i] ~ norm(mu, s)", "}",
               "algorithm mcmc { move m = metropolis(norm_rw(mu, step)); termination = iterations(10) }"), f)
  recipe = ilike::compile(f, model_parameter_list = list(s = 1))
  run = function(...)
  {
    on.exit(unlink("fixed_constant", recursive = TRUE), add = TRUE)
    ilike::MCMC(recipe, results_name = "fixed_constant", initial_values = list(list(mu = 0)), seed = 1, ...)
  }
  # the compiled value again, and a constant not fixed when compiling (step), are fine
  expect_no_error(run(model_parameter_list = list(s = 1, step = 0.5)))
  expect_error(run(model_parameter_list = list(s = 2)),
               "The constant s was fixed when the recipe was compiled")
})

test_that("an index outside its variable: -Inf when it reads parameters, an error otherwise", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  # a data index outside the variable is an error in the model
  writeLines(c("ilike 1.0", "data { y = [1, 2, 3]; g = [1, 2, 5] }", "model {",
               "  for (k in 1:2) a[k] ~ norm(0, 1)", "  for (i in 1:3) y[i] ~ norm(a[g[i]], 1)", "}"), f)
  recipe = ilike::compile(f)
  expect_error(ilike::evaluate_log(recipe, parameters = list(a = c(0, 0))),
               "a\\[...\\] \\(.*:5:30\\): the index \\(5\\) is outside the variable, which has 2 row\\(s\\) and 1 column\\(s\\)\\.")
  # an index that reads a parameter, proposed outside the variable: density 0
  writeLines(c("ilike 1.0", "data { y = 2.5 }", "model {", "  for (j in 1:3) mu[j] ~ norm(j, 0.1)",
               "  k ~ unif(1, 4)", "  y ~ norm(mu[k], 1)", "}",
               "algorithm mcmc { move m = metropolis(norm_rw(k, 2)); termination = iterations(2000) }"), f)
  recipe = ilike::compile(f)
  expect_equal(ilike::evaluate_log(recipe, parameters = list(mu = c(1, 2, 3), k = 4.5)), -Inf)
  expect_true(is.finite(ilike::evaluate_log(recipe, parameters = list(mu = c(1, 2, 3), k = 3.5))))
  on.exit(unlink("outside_index", recursive = TRUE), add = TRUE)
  expect_no_error(ilike::MCMC(recipe, results_name = "outside_index", seed = 1,
                              initial_values = list(list(mu = c(1, 2, 3), k = 2))))
  k = ilike::load_mcmc_output("outside_index")
  k = k$Value[k$ParameterName == "k"]
  expect_true(all(k >= 1 & k < 4))
})
