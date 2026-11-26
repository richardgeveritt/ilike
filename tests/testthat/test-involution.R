# Moves written as involutions by the user (dev/plan.md 3.22 E2, C288):
# move m = involution(x, auxiliary = cpp::draw, log_auxiliary = cpp::log_q,
# map = cpp::phi, log_jacobian = cpp::log_jacobian, check = k). The move draws
# v ~ q(. | x), maps (x, v) to (x', v') = phi(x, v), and accepts with
# pi(x') q(v' | x') |J| / [pi(x) q(v | x)].

involution_file = function(model, cpp, move, iterations)
{
  f = tempfile(fileext = ".ilike")
  writeLines(c("ilike 1.0", "model {", model, "}", "cpp {", cpp, "}", "algorithm mcmc {", paste0("  move m = ", move),
               paste0("  termination = iterations(", iterations, ")"), "}"), f)
  f
}

involution_chain = function(f, initial, seed = 1)
{
  name = basename(tempfile("involution_"))
  on.exit(unlink(name, recursive = TRUE), add = TRUE)
  ilike::MCMC(ilike::compile(f), results_name = name, seed = seed, initial_values = list(initial))
  ilike::load_mcmc_output(name)
}

reflection_cpp = c(
  "Parameters draw_v(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data) { return Parameters(\"v\", arma::mat(1, 1, arma::fill::value(runif(rng) < 0.5 ? -1.0 : 1.0))); }",
  "double log_v(const Parameters &parameters, const Data &data) { return std::log(0.5); }",
  "Parameters reflect(const Parameters &parameters, const Data &data) { double x = parameters[\"x\"][0], v = parameters[\"v\"][0];",
  "  Parameters out; out[\"x\"] = arma::mat(1, 1, arma::fill::value(x + v)); out[\"v\"] = arma::mat(1, 1, arma::fill::value(-v)); return out; }",
  "Parameters not_involution(const Parameters &parameters, const Data &data) { double x = parameters[\"x\"][0], v = parameters[\"v\"][0];",
  "  Parameters out; out[\"x\"] = arma::mat(1, 1, arma::fill::value(x + v)); out[\"v\"] = arma::mat(1, 1, arma::fill::value(v)); return out; }")

test_that("a discrete reflection written as an involution samples a categorical distribution", {
  skip_on_cran()
  p = c(0.05, 0.1, 0.3, 0.25, 0.2, 0.1)
  f = involution_file(paste0("  x ~ categorical([", paste(p, collapse = ", "), "])"), reflection_cpp,
                      "involution(x, auxiliary = cpp::draw_v, log_auxiliary = cpp::log_v, map = cpp::reflect)", 50000)
  on.exit(unlink(f), add = TRUE)
  o = involution_chain(f, list(x = 3), seed = 2)
  x = o$Value[o$ParameterName == "x"]
  expect_equal(tabulate(x, 6) / length(x), p, tolerance = 0.06)
})

test_that("HMC written as an involution (leapfrog, then the momentum flipped) samples a correlated normal", {
  skip_on_cran()
  cpp = c(
    "arma::colvec grad_x(const arma::colvec &x) { arma::mat S = {{1, 0.8}, {0.8, 1}}; arma::colvec mu = {1, -1}; return -arma::solve(S, x - mu); }",
    "Parameters draw_p(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data) { arma::colvec p(2); p[0] = rnorm(rng, 0.0, 1.0); p[1] = rnorm(rng, 0.0, 1.0); return Parameters(\"p\", p); }",
    "double log_p(const Parameters &parameters, const Data &data) { arma::colvec p = parameters[\"p\"]; return -0.5 * arma::dot(p, p); }",
    "Parameters leapfrog(const Parameters &parameters, const Data &data) { arma::colvec x = parameters[\"x\"], p = parameters[\"p\"]; double eps = 0.15;",
    "  for (int l = 0; l < 10; ++l) { p += 0.5 * eps * grad_x(x); x += eps * p; p += 0.5 * eps * grad_x(x); }",
    "  Parameters out; out[\"x\"] = x; out[\"p\"] = -p; return out; }")
  f = involution_file("  x ~ mvnorm([1, -1], [[1, 0.8], [0.8, 1]])", cpp,
                      "involution(x, auxiliary = cpp::draw_p, log_auxiliary = cpp::log_p, map = cpp::leapfrog)", 20000)
  on.exit(unlink(f), add = TRUE)
  o = involution_chain(f, list(x = c(0, 0)))
  X = matrix(o$Value[o$ParameterName == "x"], ncol = 2, byrow = TRUE)[-(1:500), ]
  expect_equal(colMeans(X), c(1, -1), tolerance = 0.05)
  expect_equal(diag(var(X)), c(1, 1), tolerance = 0.05)
  expect_equal(cor(X)[1, 2], 0.8, tolerance = 0.03)
})

test_that("a multiplicative move with its Jacobian, (x, u) -> (x u, 1/u), samples a gamma distribution", {
  skip_on_cran()
  # u ~ U(1/2, 2), whose density is the same at 1/u; |det d(xu, 1/u)/d(x, u)| = 1/u
  cpp = c(
    "Parameters draw_u(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data) { return Parameters(\"u\", arma::mat(1, 1, arma::fill::value(0.5 + 1.5 * runif(rng)))); }",
    "double log_u(const Parameters &parameters, const Data &data) { double u = parameters[\"u\"][0]; return ((u >= 0.5) && (u <= 2.0)) ? std::log(1.0 / 1.5) : -arma::datum::inf; }",
    "Parameters scale(const Parameters &parameters, const Data &data) { double x = parameters[\"x\"][0], u = parameters[\"u\"][0];",
    "  Parameters out; out[\"x\"] = arma::mat(1, 1, arma::fill::value(x * u)); out[\"u\"] = arma::mat(1, 1, arma::fill::value(1.0 / u)); return out; }",
    "double log_jacobian_scale(const Parameters &parameters, const Data &data) { return -std::log(parameters[\"u\"][0]); }")
  f = involution_file("  x ~ gamma(3, 1)", cpp,
                      "involution(x, auxiliary = cpp::draw_u, log_auxiliary = cpp::log_u, map = cpp::scale, log_jacobian = cpp::log_jacobian_scale)", 40000)
  on.exit(unlink(f), add = TRUE)
  o = involution_chain(f, list(x = 2))
  x = o$Value[o$ParameterName == "x"][-(1:1000)]
  expect_equal(mean(x), 3, tolerance = 0.05)
  expect_equal(var(x), 3, tolerance = 0.1)
})

test_that("involution(...) checks that its map is an involution, and refuses what does not apply", {
  model = "  x ~ categorical([0.2, 0.3, 0.5])"
  f = involution_file(model, reflection_cpp, "involution(x, auxiliary = cpp::draw_v, log_auxiliary = cpp::log_v, map = cpp::not_involution)", 50)
  on.exit(unlink(f), add = TRUE)
  expect_error(involution_chain(f, list(x = 2)), "the map is not an involution", fixed = TRUE)
  g = involution_file(model, reflection_cpp, "involution(x, auxiliary = cpp::draw_v, log_auxiliary = cpp::log_v)", 50)
  on.exit(unlink(g), add = TRUE)
  expect_error(ilike::compile(g), "map = cpp::… is needed", fixed = TRUE)
  h = involution_file(model, reflection_cpp, "involution(x, auxiliary = cpp::draw_v, log_auxiliary = cpp::log_v, map = 3)", 50)
  on.exit(unlink(h), add = TRUE)
  expect_error(ilike::compile(h), "is a function written in C++", fixed = TRUE)
})
