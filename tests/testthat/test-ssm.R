test_that("KF gives the exact log-likelihood of a linear-Gaussian SSM", {

  on.exit(unlink("kf_linear_gaussian", recursive=TRUE), add=TRUE)

  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_kf_fixture("test_kf_linear_gaussian.ilike", "kf_linear_gaussian"), truth, tolerance = 1e-6)
})

test_that("KF explains that the prior must be a built-in Gaussian", {

  on.exit(unlink("kf_user_prior", recursive=TRUE), add=TRUE)

  # Used to stop with "Invalid factor."
  expect_error(run_kf_fixture("test_kf_user_prior.ilike", "kf_user_prior"), "built-in Gaussian")
})

test_that("PF with a user-defined transition model (simulate and evaluate) gives the log-likelihood", {

  on.exit(unlink("pf_user_transition", recursive=TRUE), add=TRUE)

  # evaluate_log_transition_model used to be generated with one argument
  # instead of two, so the model did not compile (B15).
  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_pf_fixture("test_pf_user_transition.ilike", "pf_user_transition"), truth, tolerance = 0.01)
})

test_that("PF with a user-defined transition proposal taking the data gives the log-likelihood", {

  on.exit(unlink("pf_transition_proposal", recursive=TRUE), add=TRUE)

  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_pf_fixture("test_pf_transition_proposal.ilike", "pf_transition_proposal"), truth, tolerance = 0.01)
})

test_that("the one-line linear-Gaussian transition and data models work in KF and PF", {

  on.exit({
    unlink("kf_one_line", recursive=TRUE)
    unlink("pf_one_line", recursive=TRUE)
  }, add=TRUE)

  # These built-ins were parsed but nothing used them (B17).
  truth = kf_log_likelihood(c(0.5,1.0,0.2,-0.3,0.8))
  expect_equal(run_kf_fixture("test_ssm_one_line.ilike", "kf_one_line"), truth, tolerance = 1e-6)
  expect_equal(run_pf_fixture("test_ssm_one_line.ilike", "pf_one_line"), truth, tolerance = 0.01)
})

test_that("PF and EnKF on a linear-Gaussian SSM give the log-likelihood and the final filtering distribution", {

  on.exit({
    unlink("pf_linear_gaussian", recursive=TRUE)
    unlink("enkf_linear_gaussian", recursive=TRUE)
  }, add=TRUE)

  truth = kf_log_likelihood(ssm_y)
  filter = kf_final_filter(ssm_y)
  recipe = ilike::compile(system.file("extdata", "test_kf_linear_gaussian.ilike", package = "ilike"))

  z = ilike::PF(recipe, results_name = "pf_linear_gaussian", number_of_particles = 5000, seed = 1,
                smc_iterations_to_store = 10)
  expect_equal(z, truth, tolerance = 0.01)
  o = ilike::load_smc_output("pf_linear_gaussian")
  o = o[o$Target == max(o$Target), ]
  w = exp(o$LogWeight - max(o$LogWeight)); w = w/sum(w)
  m = sum(w*o$Value)
  expect_equal(m, unname(filter["mean"]), tolerance = 0.1)
  expect_equal(sum(w*(o$Value - m)^2), unname(filter["var"]), tolerance = 0.1)

  z = ilike::EnKF(recipe, results_name = "enkf_linear_gaussian", number_of_ensemble_members = 5000, seed = 1)
  expect_equal(z, truth, tolerance = 0.01)
})

test_that("PF and EnKF work with nonlinear-Gaussian transition and data models given in C++", {

  on.exit({
    unlink("pf_nonlinear", recursive=TRUE)
    unlink("enkf_nonlinear", recursive=TRUE)
  }, add=TRUE)

  truth = kf_log_likelihood(ssm_y)
  recipe = ilike::compile(system.file("extdata", "test_ssm_nonlinear.ilike", package = "ilike"))
  expect_equal(ilike::PF(recipe, results_name = "pf_nonlinear", number_of_particles = 5000, seed = 1), truth, tolerance = 0.01)
  expect_equal(ilike::EnKF(recipe, results_name = "enkf_nonlinear", number_of_ensemble_members = 5000, seed = 1), truth, tolerance = 0.01)
})

test_that("a user function with the wrong arguments gives a clear message", {

  # Used to fail inside the compiler ("different number of parameters (3 vs 2)").
  filename = system.file("extdata", "test_pf_transition_proposal_wrong_arguments.ilike", package = "ilike")
  expect_error(ilike::compile(filename), "simulate_transition_proposal function must have arguments")
})

test_that("KF and the stochastic and square-root EnKF handle a 2-D measurement with correlated noise", {
  # Until C82-C84 the KF stopped on vector measurements (a 1 x 2 data row was
  # assigned to a column vector), and the square-root EnKF used Cholesky
  # factors where the update needs symmetric square roots (about 0.45 too low
  # on this model). Both were correct only for scalar measurements.
  skip_on_cran()
  pl = list(c(0, 0), diag(2))
  model = system.file("extdata", "test_enkf_bivariate.ilike", package = "ilike")
  truth = kf_bivariate_log_likelihood()
  expect_equal(ilike::KF(ilike::compile(model, model_parameter_list = pl), model_parameter_list = pl),
               truth, tolerance = 1e-6)
  enkf = ilike::compile(model, model_parameter_list = pl)
  sqrt_enkf = ilike::compile(c(model, system.file("extdata", "test_enkf_sqrt.ilike", package = "ilike")),
                             model_parameter_list = pl)
  expect_lt(abs(ilike::EnKF(enkf, model_parameter_list = pl, number_of_ensemble_members = 20000, seed = 1) - truth), 0.15)
  expect_lt(abs(ilike::EnKF(sqrt_enkf, model_parameter_list = pl, number_of_ensemble_members = 20000, seed = 1) - truth), 0.15)
  # the ensemble adjustment update (before C183 it used the full SVD with V
  # in place of V^T, and a J x J inverse, too slow at this ensemble size)
  adjustment = tempfile(fileext = ".ilike")
  writeLines("/***enk_shifter,ilike::adjustment()***/", adjustment)
  on.exit(unlink(adjustment), add = TRUE)
  adjustment_enkf = ilike::compile(c(model, adjustment), model_parameter_list = pl)
  expect_lt(abs(ilike::EnKF(adjustment_enkf, model_parameter_list = pl, number_of_ensemble_members = 20000, seed = 1) - truth), 0.15)
})

# (C260) a measurement written as a model statement in a state-space model:
# with its data declared (data { vector y }) the filter uses it; undeclared,
# it would be a prior on a parameter y depending on the state (a filter that
# ignored the measurement, silently), which is refused
test_that("a state-space model's measurement can be a model statement on declared data", {
  set.seed(4); n = 20; x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  ssm = function(meas, first, last) c(first, "model {", "  x ~ norm(0, 1)",
    "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)", meas,
    paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
    "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
    "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }", "}", last)
  outer = function(factor, inner) c("ilike 1.0", paste0("import \"", inner, "\" as inner"), paste0("data { y = [", paste(y, collapse = ", "), "] }"),
    "model {", "  a ~ unif(-1, 1)", paste0("  factor filter = ", factor), "}")
  d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines(ssm("  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)", "ilike 1.0", c("algorithm kalman {", "}")),
             file.path(d, "inner_kf.ilike"))
  writeLines(outer("kf(inner)", "inner_kf.ilike"), file.path(d, "kf.ilike"))
  writeLines(ssm("  y ~ norm(x, 1)", c("ilike 1.0", "data { vector y }"), c("algorithm filter {", "  particles = 2000", "}")),
             file.path(d, "inner.ilike"))
  writeLines(outer("pf(inner)", "inner.ilike"), file.path(d, "pf.ilike"))
  exact = suppressMessages(ilike::evaluate_log(ilike::compile(file.path(d, "kf.ilike")), parameters = list(a = 0.7)))
  estimate = suppressMessages(ilike::evaluate_log(ilike::compile(file.path(d, "pf.ilike")), parameters = list(a = 0.7), seed = 1))
  expect_lt(abs(estimate - exact), 0.3)
  # undeclared: refused
  writeLines(ssm("  y ~ norm(x, 1)", "ilike 1.0", c("algorithm filter {", "  particles = 20", "}")), file.path(d, "inner.ilike"))
  expect_error(ilike::compile(file.path(d, "pf.ilike")), "y is not data, so this would be a prior on a parameter y that depends on the state x",
               fixed = TRUE)
})

# (C264) a Kalman or ensemble Kalman filter as the likelihood in MCMC stores
# each iteration compactly (nothing of the filter is stored by default): the
# chain and the filters' output files are as before
test_that("KF and EnKF likelihoods in MCMC write their output and the chain loads", {
  n = 10
  set.seed(4); x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  x ~ norm(0, 1)", "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)",
               "  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }", "}",
               "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
  for (factor in c("kf(inner)", "enkf(inner, ensemble_members = 20, shifter = \"stochastic\")")) {
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {",
                 "  a ~ unif(-1, 1)", paste0("  factor filter = ", factor), "}",
                 "algorithm mcmc {", "  move m = metropolis(norm_rw(a, 0.2))", "  termination = iterations(50)", "}"), file.path(d, "outer.ilike"))
    out = "kf_enkf_mcmc_output"; unlink(out, recursive = TRUE)
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    suppressMessages(ilike::MCMC(ilike::compile(file.path(d, "outer.ilike")), results_name = out, seed = 1, initial_values = list(list(a = 0.5))))
    o = ilike::load_mcmc_output(out)
    expect_equal(length(unique(o$Iteration)), 50)
    expect_true(any(grepl(if (grepl("^kf", factor)) "factor2_kf/output.h5" else "factor2_enk/output.h5",
                          list.files(out, recursive = TRUE), fixed = TRUE)))
  }
})

# (C269, PMMH P1) a transition's matrix and covariance written as C++
# functions, with reads = [...] declaring the parameters they read: computed
# once per filter run (the same chain as without the declaration); the
# declaration is checked (a function reading another variable stops), and the
# state or the time index cannot be declared
test_that("a transition's reads = [...] gives the same chain, and is checked", {
  n = 5
  set.seed(3); y = matrix(round(rnorm(2*n), 3), n, 2)
  d = tempfile(); dir.create(d); on.exit(unlink(c(d, "transition_reads_output"), recursive = TRUE), add = TRUE)
  inner = function(reads, q = "parameters[\"s\"][0]*parameters[\"s\"][0]") c("ilike 1.0", "data { vector y }", "model {",
    "  x ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
    paste0("  transition x = linear_gaussian(matrix = cpp::A, covariance = cpp::Q", reads, ")"),
    "  factor obs = density(log = cpp::meas)",
    paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
    "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
    "arma::mat A(const Parameters &parameters) { return parameters[\"a\"][0]*arma::eye(2, 2); }",
    paste0("arma::mat Q(const Parameters &parameters) { return (", q, ")*arma::eye(2, 2); }"),
    "double meas(const Parameters &parameters, const Data &data) { arma::colvec yy = arma::vectorise(data[\"y\"]); return -0.5*arma::accu(arma::square(yy - arma::vectorise(parameters[\"x\"]))); }",
    "}", "algorithm filter {", "  particles = 20", "}")
  outer = c("ilike 1.0", "import \"inner.ilike\" as inner",
            paste0("data { y = [", paste(apply(y, 1, function(r) paste0("[", paste(r, collapse = ", "), "]")), collapse = ", "), "] }"),
            "model {", "  a ~ unif(-1, 1)", "  s ~ unif(0.2, 2)", "  factor filter = pf(inner)", "}",
            "algorithm mcmc {", "  move m = metropolis(mvnorm_rw([a, s], [[0.01, 0], [0, 0.01]]))", "  termination = iterations(30)", "}")
  writeLines(outer, file.path(d, "outer.ilike"))
  run = function(reads, ...) {
    writeLines(inner(reads, ...), file.path(d, "inner.ilike"))
    unlink("transition_reads_output", recursive = TRUE)
    suppressMessages(ilike::MCMC(ilike::compile(file.path(d, "outer.ilike")), results_name = "transition_reads_output", seed = 1,
                                 initial_values = list(list(a = 0.5, s = 0.8))))
    ilike::load_mcmc_output("transition_reads_output")$Value
  }
  expect_identical(run(", reads = [a, s]"), run(""))
  # Q reads s, but only a is declared
  expect_error(run(", reads = [a]"), "do they read another variable?", fixed = TRUE)
  # the state, and the time index, cannot be declared
  writeLines(inner(", reads = [a, x]"), file.path(d, "inner.ilike"))
  expect_error(ilike::compile(file.path(d, "outer.ilike")), "reads may not name the state", fixed = TRUE)
  expect_error(run(", reads = [a, s, i]"), "may name only parameters, not the state", fixed = TRUE)
})

# (C271) the Kalman filter chains its predictions (predictions_per_update > 1)
# and its updates (several measurement factors): each used to start again
# from the stored posterior or prediction. Checked against a Kalman filter
# written here, on a 2-d state with a 2-d and a 1-d measurement factor
test_that("KF chains several predictions per update and several measurement factors", {
  m = 2; n = 8; H1 = matrix(c(1, 0.3, -0.5, 1), 2, 2); H2 = matrix(c(0.7, 0.2), 1, 2)
  set.seed(3); x = rnorm(m); Y1 = matrix(0, n, 2); Y2 = matrix(0, n, 1)
  for (t in 1:n) { x = 0.8*x + rnorm(m, 0, sqrt(0.5)); Y1[t, ] = round(H1 %*% x + rnorm(2), 3); Y2[t, ] = round(H2 %*% x + rnorm(1, 0, 2), 3) }
  mat = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(r, collapse = ", "), "]")), collapse = ", "), "]")
  kf_here = function(a, k) {
    H = rbind(H1, H2); R = diag(c(1, 1, 4)); Y = cbind(Y1, Y2); mu = rep(0, m); P = diag(m); ll = 0
    for (t in 1:n) {
      if (t > 1) for (j in 1:k) { mu = a*mu; P = a^2*P + 0.5*diag(m) }
      S = H %*% P %*% t(H) + R; v = Y[t, ] - H %*% mu
      ll = ll - 0.5*(3*log(2*pi) + determinant(S)$modulus[1] + sum(v*solve(S, v)))
      K = P %*% t(H) %*% solve(S); mu = mu + K %*% v; P = P - K %*% H %*% P
    }
    ll
  }
  d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  for (k in 1:2) {
    writeLines(c("ilike 1.0", "model {", "  x ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
                 "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = [[0.5, 0], [0, 0.5]])",
                 paste0("  factor obs1 = linear_gaussian(data = y1, state = x, matrix = ", mat(H1), ", covariance = [[1, 0], [0, 1]])"),
                 paste0("  factor obs2 = linear_gaussian(data = y2, state = x, matrix = ", mat(H2), ", covariance = [[4]])"),
                 paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
                 paste0("             time_difference = dt, time_step = 1, predictions_per_update = ", k, ")"), "}", "cpp {",
                 "arma::mat transition_matrix(const Parameters &parameters) { return parameters[\"a\"][0]*arma::eye(2, 2); }",
                 "}", "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", "data {", paste0("  y1 = ", mat(Y1)), paste0("  y2 = ", mat(Y2)), "}",
                 "model {", "  a ~ unif(-1, 1)", "  factor filter = kf(inner)", "}"), file.path(d, "outer.ilike"))
    r = suppressMessages(ilike::compile(file.path(d, "outer.ilike")))
    expect_equal(suppressMessages(ilike::evaluate_log(r, parameters = list(a = 0.6))) - log(0.5), kf_here(0.6, k), tolerance = 1e-10)
  }
})

# (C273) matrices that depend on time: each step's matrix functions are given
# the step's index i, time t and time difference dt, as the particle filter's
# are (the k predictions of a step all see the step's values). (C276) With
# steady_state the covariances are frozen once they stop changing, and a
# change of the matrices ends it: a model whose matrices change at every step
# is filtered exactly, and one whose data covariance jumps at step 30 within
# about the tolerance. Both against a Kalman filter written here.
test_that("KF handles time-varying matrices, and steady_state freezes and unfreezes", {
  mat = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(r, collapse = ", "), "]")), collapse = ", "), "]")
  H = matrix(c(1, 0.5), 1, 2); t0 = 0.5; step = 0.25
  run = function(n, k, A_of, R_of, A_cpp, R_cpp, steady) {
    set.seed(5); Y = matrix(round(rnorm(n, 0, 2), 3), n, 1)
    d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    writeLines(c("ilike 1.0", "model {", "  x ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
                 "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = [[0.5, 0], [0, 0.5]])",
                 paste0("  factor obs = linear_gaussian(data = y, state = x, matrix = ", mat(H), ", covariance = cpp::data_cov)"),
                 paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = ", t0, ","),
                 paste0("             time_difference = dt, time_step = ", step, ", predictions_per_update = ", k, ")"), "}", "cpp {",
                 paste0("arma::mat transition_matrix(const Parameters &parameters) { double a = parameters[\"a\"][0], t = parameters[\"t\"][0]; return (", A_cpp, ")*arma::eye(2, 2); }"),
                 paste0("arma::mat data_cov(const Parameters &parameters) { double i = parameters[\"i\"][0]; arma::mat R(1, 1); R(0, 0) = ", R_cpp, "; return R; }"),
                 "}", "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", "data {", paste0("  y = ", mat(Y)), "}", "model {", "  a ~ unif(-1, 1)",
                 paste0("  factor filter = kf(inner", if (steady > 0) paste0(", steady_state = ", steady), ")"), "}"), file.path(d, "outer.ilike"))
    r = suppressMessages(ilike::compile(file.path(d, "outer.ilike")))
    mu = c(0, 0); P = diag(2); ll = 0; a = 0.7
    for (j in 1:n) { i = j - 1; t = t0 + i*step
      if (j > 1) for (s in 1:k) { A = A_of(a, t); mu = A*mu; P = A^2*P + 0.5*diag(2) }
      S = c(H %*% P %*% t(H)) + R_of(i); v = Y[j] - c(H %*% mu)
      ll = ll - 0.5*(log(2*pi) + log(S) + v^2/S); K = P %*% t(H)/S; mu = mu + K*v; P = P - K %*% H %*% P }
    c(ilike = suppressMessages(ilike::evaluate_log(r, parameters = list(a = a))) - log(0.5), exact = ll)
  }
  z = run(10, 2, function(a, t) a*(1 + 0.1*t), function(i) 1 + 0.2*i, "a*(1.0 + 0.1*t)", "1.0 + 0.2*i", 1e-6)
  expect_equal(z[["ilike"]], z[["exact"]], tolerance = 1e-10)
  z = run(80, 1, function(a, t) a, function(i) if (i >= 30) 3 else 1, "a", "(i >= 30) ? 3.0 : 1.0", 1e-9)
  expect_lt(abs(z[["ilike"]] - z[["exact"]]), 1e-6)
})

# (C277, C278) The ensemble Kalman filter against the Kalman filter on models
# it used to get wrong: two measurement factors (it took them from the same
# ensemble: the product of their marginal densities; now one after another),
# and a data covariance depending on the index (its members kept the first
# step's index). The square-root shifter with fewer measurement components
# than states stopped with a dimension error (an identity of the wrong size).
enkf_check_files = function(two_factors, time_varying, members, shifter)
{
  mat = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(r, collapse = ", "), "]")), collapse = ", "), "]")
  n = 10; set.seed(7); H1 = matrix(c(1, 0.3, -0.5, 1), 2, 2); H2 = matrix(c(0.7, 0.2), 1, 2)
  R1 = function(i) if (time_varying) 1 + 0.2*i else 1
  x = rnorm(2); Y1 = matrix(0, n, 2); Y2 = matrix(0, n, 1)
  for (j in 1:n) { if (j > 1) x = 0.7*x + rnorm(2, 0, sqrt(0.5))
    Y1[j, ] = round(H1 %*% x + rnorm(2, 0, sqrt(R1(j - 1))), 3); Y2[j, ] = round(H2 %*% x + rnorm(1, 0, 2), 3) }
  d = tempfile(); dir.create(d)
  inner = c("ilike 1.0", "model {", "  x ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
            "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = [[0.5, 0], [0, 0.5]])",
            paste0("  factor obs1 = linear_gaussian(data = y1, state = x, matrix = ", mat(H1), ", covariance = cpp::data_cov)"),
            if (two_factors) paste0("  factor obs2 = linear_gaussian(data = y2, state = x, matrix = ", mat(H2), ", covariance = [[4]])"),
            paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
            "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
            "arma::mat transition_matrix(const Parameters &parameters) { return parameters[\"a\"][0]*arma::eye(2, 2); }",
            paste0("arma::mat data_cov(const Parameters &parameters) { return ", if (time_varying) "(1.0 + 0.2*parameters[\"i\"][0])" else "1.0", "*arma::eye(2, 2); }"),
            "}", "algorithm kalman {", "}")
  writeLines(inner, file.path(d, "inner.ilike"))
  outer = function(f) c("ilike 1.0", "import \"inner.ilike\" as inner", "data {", paste0("  y1 = ", mat(Y1)), paste0("  y2 = ", mat(Y2)), "}",
                        "model {", "  a ~ unif(-1, 1)", paste0("  factor filter = ", f), "}")
  writeLines(outer("kf(inner)"), file.path(d, "kf.ilike"))
  writeLines(outer(paste0("enkf(inner, ensemble_members = ", members, ", shifter = \"", shifter, "\")")), file.path(d, "enkf.ilike"))
  d
}

enkf_check = function(d, seed = 1)
  sapply(c("kf", "enkf"), function(f) suppressMessages(ilike::evaluate_log(ilike::compile(file.path(d, paste0(f, ".ilike"))),
                                                                           parameters = list(a = 0.7), seed = seed)))

test_that("EnKF runs a data covariance that reads the index, and the square-root shifter with a 1-d measurement of a 2-d state", {
  d = enkf_check_files(TRUE, TRUE, 50, "sqrt")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_true(is.finite(enkf_check(d)[["enkf"]]))
})

test_that("EnKF agrees with the KF with two measurement factors and with a time-varying data covariance", {
  skip_on_cran()
  for (case in list(list(TRUE, FALSE, "stochastic"), list(TRUE, FALSE, "sqrt"), list(FALSE, TRUE, "stochastic"))) {
    d = enkf_check_files(case[[1]], case[[2]], 5000, case[[3]])
    z = enkf_check(d)
    unlink(d, recursive = TRUE)
    expect_lt(abs(z[["enkf"]] - z[["kf"]]), 0.15)
  }
})

# (C305, dev/plan.md 3.22 U9 R3b) correlated pseudo-marginal MCMC with the
# EnKF likelihood: each member's random numbers from its own seed, and a
# proposal keeping all members' but `refresh` of them
correlated_enkf_files = function(n = 50)
{
  set.seed(4); x = numeric(n); x[1] = rnorm(1)
  for (t in 2:n) x[t] = 0.7*x[t-1] + rnorm(1, 0, sqrt(0.5))
  y = round(x + rnorm(n), 3)
  d = tempfile(); dir.create(d)
  writeLines(c("ilike 1.0", "model {", "  x ~ norm(0, 1)", "  transition x = linear_gaussian(matrix = cpp::transition_matrix, covariance = 0.5)",
               "  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "arma::mat transition_matrix(const Parameters &parameters) { arma::mat A(1, 1); A(0, 0) = parameters[\"a\"][0]; return A; }", "}",
               "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
  list(d = d, outer = function(factor, iterations) {
    f = file.path(d, "outer.ilike")
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {",
                 "  a ~ unif(-1, 1)", paste0("  factor filter = ", factor), "}",
                 "algorithm mcmc {", "  move m = metropolis(norm_rw(a, 0.1))", paste0("  termination = iterations(", iterations, ")"), "}"), f)
    out = "correlated_enkf_output"; unlink(out, recursive = TRUE)
    suppressMessages(ilike::MCMC(ilike::compile(f), results_name = out, seed = 1, initial_values = list(list(a = 0.5))))
    o = ilike::load_mcmc_output(out); o$Value[o$ParameterName == "a"]
  })
}

test_that("correlated EnKF estimates in MCMC: the default unchanged, and more proposals accepted", {
  m = correlated_enkf_files()
  on.exit(unlink(c(m$d, "correlated_enkf_output"), recursive = TRUE), add = TRUE)
  enkf = function(extra = "") paste0("enkf(inner, ensemble_members = 10, shifter = \"stochastic\"", extra, ")")
  expect_identical(m$outer(enkf(), 200), m$outer(enkf(", correlated = false"), 200))
  fresh = m$outer(enkf(), 2000)
  correlated = m$outer(enkf(", correlated = true"), 2000)
  # (10 members, 50 observations: about 0.22 against 0.59)
  expect_gt(mean(diff(correlated) != 0), mean(diff(fresh) != 0) + 0.15)
  # (C306) Crank-Nicolson on the members' standard normals (about 0.66)
  cn = m$outer(enkf(", correlated = crank_nicolson(sd = 0.2)"), 2000)
  expect_gt(mean(diff(cn) != 0), mean(diff(fresh) != 0) + 0.15)
  expect_error(m$outer(enkf(", correlated = crank_nicolson(sd = 1.5)"), 10), "between 0 and 1", fixed = TRUE)
})

test_that("Crank-Nicolson EnKF estimates refuse a model with draws other than built-in Gaussians", {
  n = 10; set.seed(4); y = round(rnorm(n), 3)
  d = tempfile(); dir.create(d)
  on.exit(unlink(c(d, "correlated_enkf_output"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "model {", "  x ~ norm(0, 1)", "  transition x = density(simulate = cpp::sim, log = cpp::lt)",
               "  factor obs = linear_gaussian(data = y, state = x, matrix = 1, covariance = 1)",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "Parameters sim(RandomNumberGenerator &rng, const Parameters &p) { Parameters o; o[\"x\"] = p[\"a\"][0]*p[\"x\"][0] + rnorm(rng, 0.0, std::sqrt(0.5)); return o; }",
               "double lt(const Parameters &p0, const Parameters &p1) { return dnorm(p1[\"x\"][0], p0[\"a\"][0]*p0[\"x\"][0], std::sqrt(0.5)); }", "}",
               "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
  outer = function(correlated, filter = "enkf(inner, ensemble_members = 10, shifter = \"stochastic\"") {
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {",
                 "  a ~ unif(-1, 1)", paste0("  factor filter = ", filter, ", correlated = ", correlated, ")"), "}",
                 "algorithm mcmc {", "  move m = metropolis(norm_rw(a, 0.1))", "  termination = iterations(20)", "}"), file.path(d, "outer.ilike"))
    unlink("correlated_enkf_output", recursive = TRUE)
    suppressMessages(ilike::MCMC(ilike::compile(file.path(d, "outer.ilike")), results_name = "correlated_enkf_output", seed = 1,
                                 initial_values = list(list(a = 0.5))))
  }
  # member refresh works with any model; Crank-Nicolson needs stored normals
  expect_no_error(outer("true"))
  expect_error(outer("crank_nicolson(sd = 0.2)"), "needs every draw of the inner model to be a built-in Gaussian", fixed = TRUE)
  # (C307) and so does the correlated particle filter
  expect_error(outer("crank_nicolson(sd = 0.2)", "pf(inner, particles = 10"), "needs every draw of the particles to be a built-in Gaussian", fixed = TRUE)
})

test_that("correlated EnKF estimates in MCMC: the posterior agrees with the Kalman filter's", {
  skip_on_cran()
  m = correlated_enkf_files()
  on.exit(unlink(c(m$d, "correlated_enkf_output"), recursive = TRUE), add = TRUE)
  kf = m$outer("kf(inner)", 10000)[-(1:500)]
  enkf = m$outer("enkf(inner, ensemble_members = 10, shifter = \"stochastic\", correlated = true, refresh = 2)", 10000)[-(1:500)]
  expect_equal(mean(enkf), mean(kf), tolerance = 0.03)
  expect_equal(sd(enkf), sd(kf), tolerance = 0.15)
  cn = m$outer("enkf(inner, ensemble_members = 10, shifter = \"stochastic\", correlated = crank_nicolson(sd = 0.2))", 10000)[-(1:500)]
  expect_equal(mean(cn), mean(kf), tolerance = 0.03)
  expect_equal(sd(cn), sd(kf), tolerance = 0.15)
})

# (C307, dev/plan.md 3.22 U9 R3c) the correlated particle filter of
# Deligiannidis, Doucet & Pitt: the filter's standard normals moved by
# Crank-Nicolson, and systematic resampling of sorted particles
test_that("the correlated particle filter in MCMC: more proposals accepted, and the refusals", {
  m = correlated_enkf_files()
  on.exit(unlink(c(m$d, "correlated_enkf_output"), recursive = TRUE), add = TRUE)
  fresh = m$outer("pf(inner, particles = 20)", 2000)
  cpf = m$outer("pf(inner, particles = 20, correlated = crank_nicolson(sd = 0.05))", 2000)
  # (20 particles, 50 observations: about 0.32 against 0.68)
  expect_gt(mean(diff(cpf) != 0), mean(diff(fresh) != 0) + 0.2)
  expect_error(m$outer("pf(inner, particles = 20, correlated = true)", 10), "use correlated = crank_nicolson(sd = s)", fixed = TRUE)
  expect_error(m$outer("pf(inner, particles = 20, correlated = crank_nicolson(sd = 0.1), parallel = TRUE)", 10),
               "parallel = TRUE is not available", fixed = TRUE)
})

test_that("the correlated particle filter in MCMC: the posterior is the Kalman filter's, for a scalar and a 2-d state", {
  skip_on_cran()
  m = correlated_enkf_files()
  on.exit(unlink(c(m$d, "correlated_enkf_output"), recursive = TRUE), add = TRUE)
  kf = m$outer("kf(inner)", 10000)[-(1:500)]
  # (s = 0.2: with only 50 observations a smaller step leaves the stored
  # normals nearly fixed, and the chain mixes slowly over them)
  cpf = m$outer("pf(inner, particles = 20, correlated = crank_nicolson(sd = 0.2))", 10000)[-(1:500)]
  expect_equal(mean(cpf), mean(kf), tolerance = 0.03)
  expect_equal(sd(cpf), sd(kf), tolerance = 0.15)
  # a 2-d state: the particles are sorted along the Hilbert curve
  n = 50; set.seed(7); x = rnorm(2); Y = matrix(0, n, 2)
  for (j in 1:n) { if (j > 1) x = 0.7*x + rnorm(2, 0, sqrt(0.5)); Y[j, ] = round(x + rnorm(2), 3) }
  d = tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  mat = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(r, collapse = ", "), "]")), collapse = ", "), "]")
  writeLines(c("ilike 1.0", "model {", "  x ~ mvnorm([0, 0], [[1, 0], [0, 1]])",
               "  transition x = linear_gaussian(matrix = cpp::A, covariance = [[0.5, 0], [0, 0.5]])",
               "  factor obs = linear_gaussian(data = y, state = x, matrix = [[1, 0], [0, 1]], covariance = [[1, 0], [0, 1]])",
               paste0("  time = ssm(index = i, first_index = 0, last_index = ", n - 1, ", time = t, initial_time = 0,"),
               "             time_difference = dt, time_step = 1, predictions_per_update = 1)", "}", "cpp {",
               "arma::mat A(const Parameters &parameters) { return parameters[\"a\"][0]*arma::eye(2, 2); }", "}",
               "algorithm kalman {", "}"), file.path(d, "inner.ilike"))
  run = function(factor) {
    writeLines(c("ilike 1.0", "import \"inner.ilike\" as inner", "data {", paste0("  y = ", mat(Y)), "}", "model {",
                 "  a ~ unif(-1, 1)", paste0("  factor filter = ", factor), "}",
                 "algorithm mcmc {", "  move m = metropolis(norm_rw(a, 0.1))", "  termination = iterations(10000)", "}"), file.path(d, "outer.ilike"))
    unlink("correlated_enkf_output", recursive = TRUE)
    suppressMessages(ilike::MCMC(ilike::compile(file.path(d, "outer.ilike")), results_name = "correlated_enkf_output", seed = 1,
                                 initial_values = list(list(a = 0.5))))
    o = ilike::load_mcmc_output("correlated_enkf_output"); o$Value[o$ParameterName == "a"][-(1:500)]
  }
  kf2 = run("kf(inner)")
  cpf2 = run("pf(inner, particles = 30, correlated = crank_nicolson(sd = 0.2))")
  expect_equal(mean(cpf2), mean(kf2), tolerance = 0.03)
  expect_equal(sd(cpf2), sd(kf2), tolerance = 0.15)
})
