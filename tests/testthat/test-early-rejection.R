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
