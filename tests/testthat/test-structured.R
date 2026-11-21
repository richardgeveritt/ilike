# Structured models (dev/plan.md 3.6, C233): sums over index sets in
# statements, sum(k in a:b, expression), and the reference problems: a CAR
# model, an Ising model with neighbour sums, a stochastic block model.

test_that("a CAR model with neighbour sums: its density, gradient and element moves", {
  g = lattice(3)
  W = matrix(0, g$n, g$n); for (i in 1:g$n) W[i, g$nb_list[[i]]] = 1
  Q = diag(g$nn) - 0.8*W
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "structured_car"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", "  for (i in 1:9) {",
               "    factor car += -0.5 / pow(1.5, 2) * (nn[i] * pow(x[i], 2) - 0.8 * x[i] * sum(k in 1:nn[i], x[nb[i, k]]))",
               "  }", "}", "algorithm mcmc {", "  for (j in 1:9) move m = metropolis(norm_rw(x[j], 1.5))",
               "  termination = iterations(20000)", "}"), f)
  recipe = ilike::compile(f)
  x = seq(-1, 1, length = 9)
  expect_equal(ilike::evaluate_log(recipe, parameters = list(x = x)), -0.5*sum(x*(Q %*% x))/1.5^2, tolerance = 1e-12)
  expect_equal(as.vector(ilike:::evaluate_gradient_log(recipe, list(x = x), "x")), as.vector(-(Q %*% x)/1.5^2), tolerance = 1e-12)
  skip_on_cran()
  unlink("structured_car", recursive = TRUE)
  ilike::MCMC(recipe, results_name = "structured_car", seed = 1, initial_values = list(list(x = rep(0, 9))))
  o = ilike::load_mcmc_output("structured_car", wide = TRUE)
  X = as.matrix(o[o$Iteration > 2000, grep("^x_", names(o))])
  expect_equal(as.vector(apply(X, 2, var)), diag(solve(Q))*1.5^2, tolerance = 0.06)
})

test_that("an Ising model with neighbour sums, by gibbs on the spins, against enumeration", {
  skip_on_cran()
  g = lattice(3)
  beta = 0.4
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "structured_ising"), recursive = TRUE), add = TRUE)
  # each edge is counted from both ends: beta/2 per end
  writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", "  for (i in 1:9) s[i] ~ bern(0.5)",
               "  for (i in 1:9) {", "    factor ising += 0.2 * (2 * s[i] - 1) * sum(k in 1:nn[i], 2 * s[nb[i, k]] - 1)", "  }", "}",
               "algorithm mcmc {", "  for (j in 1:9) move m = gibbs(s[j])", "  termination = iterations(20000)", "}"), f)
  unlink("structured_ising", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "structured_ising", seed = 1, initial_values = list(list(s = rep(0, 9))))
  o = ilike::load_mcmc_output("structured_ising", wide = TRUE)
  S = 2*as.matrix(o[, grep("^s_", names(o))]) - 1
  # exact E[s1 s2] (neighbours) and E[s1 s9] (opposite corners) by enumeration
  configurations = as.matrix(expand.grid(rep(list(c(-1, 1)), 9)))
  energy = apply(configurations, 1, function(v) sum(sapply(1:9, function(i) v[i]*sum(v[g$nb_list[[i]]])))/2)
  w = exp(beta*energy); w = w/sum(w)
  expect_equal(mean(S[, 1]*S[, 2]), sum(w*configurations[, 1]*configurations[, 2]), tolerance = 0.03)
  expect_equal(mean(S[, 1]*S[, 9]), sum(w*configurations[, 1]*configurations[, 9]), tolerance = 0.08)
})

test_that("a stochastic block model: gibbs on the block labels against enumeration", {
  skip_on_cran()
  set.seed(4)
  n = 6; pairs = t(combn(n, 2)); B = matrix(c(0.6, 0.3, 0.3, 0.5), 2); prior = c(0.6, 0.4)
  A = rbinom(nrow(pairs), 1, 0.5)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "structured_sbm"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { A = [", paste(A, collapse = ", "), "]; e1 = [", paste(pairs[, 1], collapse = ", "), "]; e2 = [",
                                    paste(pairs[, 2], collapse = ", "), "]; B = [[0.6, 0.3], [0.3, 0.5]] }"), "model {",
               "  for (i in 1:6) z[i] ~ categorical([0.6, 0.4])", "  for (m in 1:15) A[m] ~ bern(B[z[e1[m]], z[e2[m]]])", "}",
               "algorithm mcmc {", "  for (j in 1:6) move m = gibbs(z[j])", "  termination = iterations(20000)", "}"), f)
  unlink("structured_sbm", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "structured_sbm", seed = 1, initial_values = list(list(z = rep(1, 6))))
  o = ilike::load_mcmc_output("structured_sbm", wide = TRUE)
  Z = as.matrix(o[, grep("^z_", names(o))])
  labels = as.matrix(expand.grid(rep(list(1:2), n)))
  logw = apply(labels, 1, function(z) sum(log(prior[z])) +
                 sum(dbinom(A, 1, B[cbind(z[pairs[, 1]], z[pairs[, 2]])], log = TRUE)))
  w = exp(logw - max(logw)); w = w/sum(w)
  expect_equal(as.vector(colMeans(Z == 1)), as.vector(colSums(w*(labels == 1))), tolerance = 0.04)
})

test_that("sum is refused where its range or index cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  check = function(statement, message)
  {
    writeLines(c("ilike 1.0", "data { y = [1, 2]; n = 2 }", "model {", "  for (i in 1:2) x[i] ~ norm(0, 1)", "  m ~ unif(1, 2)", statement, "}"), f)
    expect_error(ilike::compile(f), message, fixed = TRUE)
  }
  check("  for (i in 1:2) y[i] ~ norm(sum(k in 1:m, x[k]), 1)", "sum's range depends on the parameter m")
  check("  for (i in 1:2) y[i] ~ norm(sum(i in 1:n, x[i]), 1)", "sum's index i is already the name of a loop index")
})

test_that("a CAR model with random rho and sigma: its normalising constant from the eigenvalues", {
  # (C234) det(D - rho W) = prod(nn) prod(1 - rho lambda_k), lambda the
  # eigenvalues of D^-1/2 W D^-1/2 (data)
  g = lattice(3)
  W = matrix(0, g$n, g$n); for (i in 1:g$n) W[i, g$nb_list[[i]]] = 1
  D = diag(g$nn)
  lambda = eigen(diag(1/sqrt(g$nn)) %*% W %*% diag(1/sqrt(g$nn)), symmetric = TRUE)$values
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, "; lambda = [", paste(signif(lambda, 17), collapse = ", "), "] }"), "model {",
               "  rho ~ unif(-1, 1)", "  sigma ~ halfnorm(2)",
               "  factor car_constant += -4.5 * log(2 * 3.14159265358979) + 0.5 * sum(k in 1:9, log(nn[k])) - 9 * log(sigma) + 0.5 * sum(k in 1:9, log(1 - rho * lambda[k]))",
               "  for (i in 1:9) {",
               "    factor car += -0.5 / pow(sigma, 2) * (nn[i] * pow(x[i], 2) - rho * x[i] * sum(k in 1:nn[i], x[nb[i, k]]))",
               "  }", "}"), f)
  recipe = ilike::compile(f)
  P = list(rho = 0.6, sigma = 1.3, x = seq(-1, 1, length = 9))
  Sigma = P$sigma^2*solve(D - P$rho*W)
  log_mvnorm = -4.5*log(2*pi) - 0.5*determinant(Sigma)$modulus[1] - 0.5*sum(P$x*solve(Sigma, P$x))
  expect_equal(ilike::evaluate_log(recipe, parameters = P), dunif(0.6, -1, 1, log = TRUE) + log(2) + dnorm(1.3, 0, 2, log = TRUE) + log_mvnorm,
               tolerance = 1e-10)
  for (name in c("rho", "sigma"))
  {
    h = 1e-6
    fd = (ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] + h)) -
          ilike::evaluate_log(recipe, parameters = replace(P, name, P[[name]] - h)))/(2*h)
    expect_equal(ilike:::evaluate_gradient_log(recipe, P, name)[1], fd, tolerance = 1e-5, label = name)
  }
})

test_that("a Potts model with equals(), by gibbs on the colours, against enumeration", {
  # (C234) three colours on a 3 x 3 lattice, beta = 0.5 per edge (beta/2 per end)
  skip_on_cran()
  g = lattice(3)
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "structured_potts"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", paste0("data { ", g$data, " }"), "model {", "  for (i in 1:9) x[i] ~ categorical([1, 1, 1] / 3)",
               "  for (i in 1:9) {", "    factor potts += 0.25 * sum(k in 1:nn[i], equals(x[i], x[nb[i, k]]))", "  }", "}",
               "algorithm mcmc {", "  for (j in 1:9) move m = gibbs(x[j])", "  termination = iterations(20000)", "}"), f)
  unlink("structured_potts", recursive = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = "structured_potts", seed = 1, initial_values = list(list(x = rep(1, 9))))
  o = ilike::load_mcmc_output("structured_potts", wide = TRUE)
  X = as.matrix(o[, grep("^x_", names(o))])
  configurations = as.matrix(expand.grid(rep(list(1:3), 9)))
  agreements = apply(configurations, 1, function(v) sum(sapply(1:9, function(i) sum(v[i] == v[g$nb_list[[i]]])))/2)
  w = exp(0.5*agreements); w = w/sum(w)
  expect_equal(mean(X[, 1] == X[, 2]), sum(w*(configurations[, 1] == configurations[, 2])), tolerance = 0.03)
  expect_equal(mean(X[, 1] == X[, 9]), sum(w*(configurations[, 1] == configurations[, 9])), tolerance = 0.05)
})

test_that("terms written in C++ with declared reads behave as model statements", {
  # (C235) for (i in 1:n) factor lik += cpp::term(i, reads = [...])
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "structured_cpp_a", "structured_cpp_b"), recursive = TRUE), add = TRUE)
  head = c("ilike 1.0", "data { y = [1.2, 0.4, 2.1, 1.7, 0.9, 1.4]; group = [1, 1, 2, 2, 3, 3] }", "model {",
           "  for (j in 1:3) mu[j] ~ norm(0, 2)", "  sigma ~ halfnorm(1)")
  algorithm = c("algorithm mcmc {", "  for (j in 1:3) move a = metropolis(norm_rw(mu[j], 0.5))",
                "  move b = metropolis(norm_rw(sigma, 0.3))", "  termination = iterations(500)", "}")
  cpp = c("cpp {", "double term(const Parameters &parameters, const Data &data, int i)", "{",
          "  int g = (int)data[\"group\"][i - 1];",
          "  return dnorm(data[\"y\"][i - 1], parameters[\"mu\"][g - 1], parameters[\"sigma\"][0]);", "}", "}")
  run = function(likelihood, name, extra = character(0))
  {
    writeLines(c(head, likelihood, "}", sub("termination", paste0(extra, "termination"), algorithm), cpp), f)
    unlink(name, recursive = TRUE)
    ilike::MCMC(ilike::compile(f), results_name = name, seed = 1, initial_values = list(list(mu = c(0, 0, 0), sigma = 1)))
    ilike::load_mcmc_output(name)
  }
  a = run("  for (i in 1:6) y[i] ~ norm(mu[group[i]], sigma)", "structured_cpp_a")
  b = run("  for (i in 1:6) factor lik += cpp::term(i, reads = [mu[group[i]], sigma])", "structured_cpp_b")
  # the same chain: the same terms, recomputed for the same elements
  expect_identical(a$Value, b$Value)
  # a wrong declaration is caught by check_factor_index (term i reads mu[group[i]], not mu[1])
  expect_error(run("  for (i in 1:6) factor lik += cpp::term(i, reads = [mu[1], sigma])", "structured_cpp_b", "check_factor_index = true\n  "),
               "check_factor_index")
  # and a missing one is refused
  writeLines(c(head, "  for (i in 1:6) factor lik += cpp::term(i)", "}", algorithm, cpp), f)
  expect_error(ilike::compile(f), "a term written in C++ declares what it reads", fixed = TRUE)
})
