# (C313, dev/plan.md 3.22 U6) reversible jump: move bd = birth_death(k, [mu],
# points = N), on a model indicator k ~ categorical(…) whose variables' loops
# run to k (padded to k's largest value); points = N is RmJ-MCMC (Andrieu,
# Yildirim, Doucet & Chopin, arXiv 2101.01253, Algorithm 3)

rj_model = function(move, prior = NULL, mu_statement = "  for (j in 1:k) mu[j] ~ norm(0, 1)", iterations = 20000)
{
  n = 4; set.seed(3); Y = round(cbind(rnorm(n, 0.8), rnorm(n, 0.6), rnorm(n, 0.2)), 3)
  logm = function(y) { n = length(y); -n/2*log(2*pi) - 0.5*log(1 + n) - 0.5*(sum(y^2) - sum(y)^2/(1 + n)) }
  lp = cumsum(apply(Y, 2, logm))
  # the prior on k chosen so that the posterior is (0.5, 0.3, 0.2)
  if (is.null(prior)) { prior = c(0.5, 0.3, 0.2)*exp(-(lp - max(lp))); prior = prior/sum(prior) }
  posterior = prior*exp(lp - max(lp)); posterior = posterior/sum(posterior)
  mat = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(r, collapse = ", "), "]")), collapse = ", "), "]")
  list(lines = c("ilike 1.0", paste0("data { y = ", mat(Y), " }"), "model {",
                 paste0("  k ~ categorical([", paste(signif(prior, 12), collapse = ", "), "])"), mu_statement,
                 paste0("  for (j in 1:k) for (i in 1:", n, ") y[i, j] ~ norm(mu[j], 1)"), "}",
                 "algorithm mcmc {", "  for (j in 1:3) move w = metropolis(norm_rw(mu[j], 0.7))", paste0("  move bd = ", move),
                 paste0("  termination = iterations(", iterations, ")"), "}"),
       posterior = posterior)
}

test_that("birth_death: the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(rj_model("birth_death(k, [mu], points = 0)")$lines, f)
  expect_error(ilike::compile(f), "a whole number of at least 1", fixed = TRUE)
  writeLines(rj_model("birth_death(mu, [k])")$lines, f)
  expect_error(ilike::compile(f), "is a model indicator", fixed = TRUE)
  writeLines(rj_model("birth_death(k, [mu])", mu_statement = "  for (j in 1:k) mu[j] ~ norm(j, 1)")$lines, f)
  expect_error(ilike::compile(f), "reads the entry's position", fixed = TRUE)
})

# (C314) the padded entries after k's value are missing in the loaded results
test_that("birth_death: the entries after the indicator's value are NA in load_mcmc_output", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rj_output"), recursive = TRUE), add = TRUE)
  writeLines(rj_model("birth_death(k, [mu], points = 2)", iterations = 300)$lines, f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "rj_output", seed = 1, initial_values = list(list(k = 1, mu = c(0, 0, 0))))
  o = ilike::load_mcmc_output("rj_output")
  mu = o[o$ParameterName == "mu", ]; k = o[o$ParameterName == "k", ]
  kk = k$Value[match(paste(mu$Chain, mu$Iteration), paste(k$Chain, k$Iteration))]
  expect_true(any(kk < 3))
  expect_true(all(is.na(mu$Value[mu$Dimension > kk])))
  expect_true(all(!is.na(mu$Value[mu$Dimension <= kk])))
  w = ilike::load_mcmc_output("rj_output", wide = TRUE)
  expect_true(all(is.na(w$mu_3[w$k_1 < 3])) && all(!is.na(w$mu_1)))
})

test_that("birth_death: the posterior model probabilities, with and without averaging (RmJ)", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rj_output"), recursive = TRUE), add = TRUE)
  for (N in c(1, 5)) {
    model = rj_model(paste0("birth_death(k, [mu], points = ", N, ")"), iterations = 40000)
    writeLines(model$lines, f)
    unlink("rj_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "rj_output", seed = 1, initial_values = list(list(k = 1, mu = c(0, 0, 0))))
    o = ilike::load_mcmc_output("rj_output")
    k = o$Value[(o$ParameterName == "k") & (o$Iteration > 2000)]
    expect_lt(max(abs(tabulate(k, 3)/length(k) - model$posterior)), 0.02)
  }
})

# (C315) rj(k, [mu], auxiliary, log_auxiliary, map): reversible jump with the
# user's dimension-matching map, here birth and death of the last entry
rj_cpp = function(shift = "")
  c("cpp {",
    "Parameters draw_u(RandomNumberGenerator &rng, const Parameters &p, const Data &d) { Parameters u; if (p[\"k_proposed\"][0] > p[\"k\"][0]) u[\"u\"] = arma::mat(1, 1, arma::fill::value(rnorm(rng, 0.0, 1.0))); return u; }",
    "double log_u(const Parameters &p, const Data &d) { if (p[\"k_proposed\"][0] > p[\"k\"][0]) return dnorm(p[\"u\"][0], 0.0, 1.0); return 0.0; }",
    "Parameters jump(const Parameters &p, const Data &d) { Parameters out; arma::mat mu = p[\"mu\"]; int k = (int)p[\"k\"][0];",
    paste0("  if (p[\"k_proposed\"][0] > k) mu(k) = p[\"u\"][0]", shift, "; else out[\"u\"] = arma::mat(1, 1, arma::fill::value(mu(k - 1)));"),
    "  out[\"mu\"] = mu; return out; }", "}")
rj_map_model = function(move, shift = "", iterations = 20000)
{
  model = rj_model(move, iterations = iterations)
  model$lines = append(model$lines, rj_cpp(shift), after = which(model$lines == "}")[1])
  model
}

test_that("rj with the user's map: the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rj_output"), recursive = TRUE), add = TRUE)
  writeLines(rj_map_model("rj(k, [mu], auxiliary = cpp::draw_u, log_auxiliary = cpp::log_u)")$lines, f)
  expect_error(ilike::compile(f), "map = cpp::… is needed", fixed = TRUE)
  # a map that the jump back does not undo
  writeLines(rj_map_model("rj(k, [mu], auxiliary = cpp::draw_u, log_auxiliary = cpp::log_u, map = cpp::jump)", shift = " + 1.0", iterations = 50)$lines, f)
  expect_error(ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "rj_output", seed = 1, initial_values = list(list(k = 1, mu = c(0, 0, 0)))),
               "does not undo the map", fixed = TRUE)
})

test_that("rj with the user's map: the posterior model probabilities, with and without averaging (RmJ)", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rj_output"), recursive = TRUE), add = TRUE)
  for (N in c(1, 5)) {
    model = rj_map_model(paste0("rj(k, [mu], auxiliary = cpp::draw_u, log_auxiliary = cpp::log_u, map = cpp::jump, points = ", N, ")"), iterations = 40000)
    writeLines(model$lines, f)
    unlink("rj_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "rj_output", seed = 1, initial_values = list(list(k = 1, mu = c(0, 0, 0))))
    o = ilike::load_mcmc_output("rj_output")
    k = o$Value[(o$ParameterName == "k") & (o$Iteration > 2000)]
    expect_lt(max(abs(tabulate(k, 3)/length(k) - model$posterior)), 0.02)
  }
})

# (stage B) annealed reversible jump (Karagiannis & Andrieu 2013): bridges = K
# between the two models' extended targets, with the named moves as the kernel
test_that("annealed reversible jump: no bridges is the plain move, and the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rj_output"), recursive = TRUE), add = TRUE)
  chain = function(move) {
    writeLines(rj_model(move, iterations = 300)$lines, f)
    unlink("rj_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "rj_output", seed = 1, initial_values = list(list(k = 1, mu = c(0, 0, 0))))
    ilike::load_mcmc_output("rj_output")$Value
  }
  expect_identical(chain("birth_death(k, [mu], bridges = 0, kernel = w)"), chain("birth_death(k, [mu])"))
  writeLines(rj_model("birth_death(k, [mu], bridges = 2)")$lines, f)
  expect_error(ilike::compile(f), "the moves run at each bridge, by name", fixed = TRUE)
  model = rj_model("birth_death(k, [mu], bridges = 2, kernel = [w, j])")
  model$lines = append(model$lines, "  move j = metropolis_hastings(dimension_jump(k, 1, 3))", after = grep("move w = ", model$lines))
  writeLines(model$lines, f)
  expect_error(ilike::compile(f), "move j changes k, the model indicator the jump changes", fixed = TRUE)
})

test_that("annealed reversible jump: the posterior model probabilities, and acceptance rising with the bridges", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "rj_output"), recursive = TRUE), add = TRUE)
  probabilities = function(model, iterations) {
    writeLines(model$lines, f)
    unlink("rj_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "rj_output", seed = 1, initial_values = list(list(k = 1, mu = c(0, 0, 0))))
    o = ilike::load_mcmc_output("rj_output")
    k = o$Value[(o$ParameterName == "k") & (o$Iteration > iterations/10)]
    list(p = tabulate(k, 3)/length(k), changes = mean(diff(k) != 0))
  }
  for (move in c("birth_death(k, [mu], bridges = 5, kernel = w)", "birth_death(k, [mu], points = 3, bridges = 5, kernel = w)")) {
    model = rj_model(move, iterations = 40000)
    expect_lt(max(abs(probabilities(model, 40000)$p - model$posterior)), 0.02)
  }
  model = rj_map_model("rj(k, [mu], auxiliary = cpp::draw_u, log_auxiliary = cpp::log_u, map = cpp::jump, bridges = 5, kernel = w)", iterations = 40000)
  expect_lt(max(abs(probabilities(model, 40000)$p - model$posterior)), 0.02)
  # births from a wide prior (sd 5), poorly matched to the likelihood: the
  # bridges carry the new entry towards it (k changes in 0.17 of iterations
  # with no bridges, 0.34 with 20)
  n = 4; set.seed(3); Y = round(cbind(rnorm(n, 0.8), rnorm(n, 0.6), rnorm(n, 0.2)), 3)
  logm5 = function(y) { n = length(y); -n/2*log(2*pi) - 0.5*log(1 + 25*n) - 0.5*(sum(y^2) - 25*sum(y)^2/(1 + 25*n)) }
  lp = cumsum(apply(Y, 2, logm5)); prior = c(0.5, 0.3, 0.2)*exp(-(lp - max(lp))); prior = prior/sum(prior)
  changes = sapply(c(0, 20), function(K) probabilities(rj_model(paste0("birth_death(k, [mu], bridges = ", K, ", kernel = w)"), prior = prior,
                                                                mu_statement = "  for (j in 1:k) mu[j] ~ norm(0, 5)", iterations = 20000), 20000)$changes)
  expect_gt(changes[2], changes[1] + 0.1)
})
