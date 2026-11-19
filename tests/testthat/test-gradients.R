# Gradients of model statements (dev/plan.md I4e): code generated from each
# term, checked against finite differences, and used by gradient-based moves.

gradient_parameters = list(mu = 0.3, sigma = 0.8, a = 1.2, r = 0.7, q = 1.1, t = 0.4, p = 0.35,
                           v = c(0.2, -0.5), th = c(0.1, -0.4, 0.6))

gradient_errors = function(recipe, P)
{
  errors = c()
  for (name in names(P))
  {
    g = ilike:::evaluate_gradient_log(recipe, P, name)
    for (k in seq_along(P[[name]]))
    {
      h = 1e-5; Pp = P; Pm = P
      Pp[[name]][k] = Pp[[name]][k] + h
      Pm[[name]][k] = Pm[[name]][k] - h
      fd = (ilike::evaluate_log(recipe, parameters = Pp) - ilike::evaluate_log(recipe, parameters = Pm))/(2*h)
      errors = c(errors, abs(g[k] - fd)/(1 + abs(fd)))
    }
  }
  errors
}

test_that("gradients generated from model statements agree with finite differences", {
  skip_on_cran()
  file = system.file("extdata", "lang_model_gradients.ilike", package = "ilike")
  # every distribution (norm, lnorm, gamma, exp, unif, tnorm, beta, pois,
  # bern, binom, mvnorm) and function, nodes and += terms; plates
  expect_lt(max(gradient_errors(ilike::compile(file), gradient_parameters)), 1e-6)
  # and the factors evaluated whole
  old = options(ilike.plate_factors = FALSE)
  on.exit(options(old), add = TRUE)
  expect_lt(max(gradient_errors(ilike::compile(file), gradient_parameters)), 1e-6)
})

test_that("the gradients of several variables in one pass equal those taken one at a time", {
  skip_on_cran()
  # (dev/plan.md C204) a stacked name "a;b;c" asks each factor once for all
  # of its variables; nine variables use the set's hash index, three its
  # linear search; with plates off, the factors are asked one variable at a time
  file = system.file("extdata", "lang_model_gradients.ilike", package = "ilike")
  for (plates in c(TRUE, FALSE))
  {
    old = options(ilike.plate_factors = plates)
    recipe = ilike::compile(file)
    for (names in list(names(gradient_parameters), c("th", "mu", "v")))
    {
      one_at_a_time = unlist(lapply(names, function(name)
        as.vector(ilike:::evaluate_gradient_log(recipe, gradient_parameters, name))))
      stacked = as.vector(ilike:::evaluate_gradient_log(recipe, gradient_parameters, paste(names, collapse = ";")))
      expect_identical(stacked, one_at_a_time)
    }
    options(old)
  }
})

test_that("Langevin and Barker moves use the generated gradients", {
  skip_on_cran()
  y = c(1.2, 0.4, -0.3, 0.8, 2.1, 1.7); K = 6
  L = matrix(0, K, K)
  for (k in 1:(K - 1)) L[c(k, k + 1), c(k, k + 1)] = L[c(k, k + 1), c(k, k + 1)] + matrix(c(1, -1, -1, 1), 2)
  Q = diag(K)/100 + 2*L + diag(K)
  S = solve(Q)
  covariance = paste0("[", paste(apply(S, 1, function(r) paste0("[", paste(signif(r, 8), collapse = ", "), "]")), collapse = ", "), "]")
  lines = readLines(system.file("extdata", "lang_model_smooth.ilike", package = "ilike"))
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  for (proposal in c("langevin", "barker"))
  {
    writeLines(c(lines[seq_len(grep("^algorithm", lines) - 1)], "algorithm mcmc {",
                 paste0("  move m = metropolis_hastings(", proposal, "(theta, ", covariance, ", 0.8))"),
                 "  termination = iterations(20000)", "}"), f)
    name = paste0("gradient_", proposal)
    ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(theta = rep(0, K))))
    o = ilike::load_mcmc_output(name)
    unlink(name, recursive = TRUE)
    theta = sapply(1:K, function(k) o$Value[o$Dimension == k & o$Iteration > 2000])
    expect_lt(max(abs(colMeans(theta) - solve(Q, y))/sqrt(diag(S))), 0.1)
    expect_equal(apply(theta, 2, sd), sqrt(diag(S)), tolerance = 0.05)
  }
})

test_that("generated and user-written gradients are added", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "gradient_mixed")), add = TRUE)
  writeLines(c("ilike 1.0", "data { x = [1.3, 0.2, 0.9] }", "model {", "  theta ~ norm(0, 10)",
               "  factor lik = density(log = cpp::ll, gradient = cpp::gll)", "}", "cpp {",
               "double ll(const Parameters &p, const Data &d) { double s = 0; for (int i = 0; i < 3; ++i) s += dnorm(d[\"x\"][i], p[\"theta\"][0], 1.0); return s; }",
               "arma::mat gll(const std::string &variable, const Parameters &parameters, const Data &data) { arma::mat g(1, 1); g(0, 0) = arma::accu(data[\"x\"] - parameters[\"theta\"][0]); return g; }",
               "}", "algorithm mcmc { move m = metropolis_hastings(langevin(theta, 0.3, 1)); termination = iterations(20000) }"), f)
  recipe = ilike::compile(f)
  # the prior's generated gradient and the likelihood's own, at theta = 0.5
  expect_equal(as.numeric(ilike:::evaluate_gradient_log(recipe, list(theta = 0.5), "theta")),
               -0.5/100 + sum(c(1.3, 0.2, 0.9) - 0.5), tolerance = 1e-12)
  on.exit(unlink("gradient_mixed", recursive = TRUE), add = TRUE)
  ilike::MCMC(recipe, results_name = "gradient_mixed", seed = 1, initial_values = list(list(theta = 0)))
  v = ilike::load_mcmc_output("gradient_mixed")$Value
  v = v[2001:length(v)]
  variance = 1/(3 + 1/100)
  expect_equal(mean(v), variance*2.4, tolerance = 0.05)
  expect_equal(sd(v), sqrt(variance), tolerance = 0.05)
})

test_that("gradients that are not available say so when asked for", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  # a covariance that is a parameter
  writeLines(c("ilike 1.0", "data { obs = [0.4, -0.3] }", "model {", "  v ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
               "  C ~ mvnorm([1, 0, 0, 1], [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])",
               "  obs ~ mvnorm(v, C)", "}"), f)
  recipe = ilike::compile(f)
  P = list(v = c(0.1, 0.2), C = matrix(c(1, 0.2, 0.2, 1), 2))
  expect_true(all(is.finite(ilike:::evaluate_gradient_log(recipe, P, "v"))))
  expect_error(ilike:::evaluate_gradient_log(recipe, P, "C"), "covariance of mvnorm \\(C\\) is not available")
  # an index that reads a parameter has no derivative through the index
  writeLines(c("ilike 1.0", "data { y = 2.5 }", "model {", "  for (j in 1:3) mu[j] ~ norm(j, 0.1)",
               "  k ~ unif(1, 4)", "  y ~ norm(mu[k], 1)", "}"), f)
  recipe = ilike::compile(f)
  expect_equal(as.numeric(ilike:::evaluate_gradient_log(recipe, list(mu = c(1, 2, 3), k = 2.5), "k")), 0)
  expect_equal(as.numeric(ilike:::evaluate_gradient_log(recipe, list(mu = c(1, 2, 3), k = 2.5), "mu")),
               c(0, -(2 - 2)/0.01 + (2.5 - 2), -(3 - 3)/0.01), tolerance = 1e-12)
})
