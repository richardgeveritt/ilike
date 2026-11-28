# (C310, dev/plan.md 3.22 U5a) AIS MCMC: move m = annealed(proposal, bridges = K,
# kernel = g), a joint move on the parameters and the latent variables, the
# latter taken along by the moves g at K values between the current and
# proposed parameters (Karagiannis & Andrieu 2013; arXiv 1809.02527)

annealed_model = function(move, others = character(0), T = 50, iterations = 5000)
{
  set.seed(1); z = rnorm(T, 1.5, 1); y = round(rnorm(T, z, 1), 3)
  prec = 1/100 + T/2
  list(lines = c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  theta ~ norm(0, 10)",
                 paste0("  for (t in 1:", T, ") z[t] ~ norm(theta, 1)"), paste0("  for (t in 1:", T, ") y[t] ~ norm(z[t], 1)"), "}",
                 "algorithm mcmc {", paste0("  for (t in 1:", T, ") move g = gibbs(z[t])"), others, paste0("  move m = ", move),
                 paste0("  termination = iterations(", iterations, ")"), "}"),
       mean = (sum(y)/2)/prec, sd = sqrt(1/prec), T = T)
}

test_that("annealed moves: the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  run = function(model) {
    writeLines(model$lines, f)
    unlink("annealed_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
                initial_values = list(list(theta = 0, z = rep(0, model$T))))
  }
  writeLines(annealed_model("annealed(norm_rw(theta, 0.3), bridges = 2, kernel = h)")$lines, f)
  expect_error(ilike::compile(f), "there is no move called h", fixed = TRUE)
  writeLines(annealed_model("annealed(norm_rw(theta, 0.3), bridges = 2)")$lines, f)
  expect_error(ilike::compile(f), "the moves on the latent variables, by name", fixed = TRUE)
  # an annealed move inside another's kernel, and a kernel move that changes theta
  writeLines(annealed_model("annealed(norm_rw(theta, 0.3), bridges = 2, kernel = a)",
                            "  move a = annealed(norm_rw(theta, 0.1), bridges = 1, kernel = g)")$lines, f)
  expect_error(ilike::compile(f), "cannot be part of another's kernel", fixed = TRUE)
  # (stage A: refused when compiling, from the variables each move changes)
  writeLines(annealed_model("annealed(norm_rw(theta, 0.3), bridges = 2, kernel = w)", "  move w = metropolis(norm_rw(theta, 0.1))")$lines, f)
  expect_error(ilike::compile(f), "move w changes theta, which the annealed move's proposal changes", fixed = TRUE)
})

# (C311) with no bridges an annealed move is the plain move on theta: the same chain
test_that("annealed moves with no bridges give the plain move's chain, for symmetric and Langevin proposals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  chain = function(move) {
    model = annealed_model(move, iterations = 200)
    writeLines(model$lines, f)
    unlink("annealed_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
                initial_values = list(list(theta = 1.5, z = rep(1.5, model$T))))
    ilike::load_mcmc_output("annealed_output")$Value
  }
  expect_identical(chain("annealed(norm_rw(theta, 0.3), bridges = 0, kernel = g)"), chain("metropolis(norm_rw(theta, 0.3))"))
  expect_identical(chain("annealed(langevin(theta, [[1]], 0.15), bridges = 0, kernel = g)"),
                   chain("metropolis_hastings(langevin(theta, [[1]], 0.15))"))
})

# (stage A) the path is a run of the package's SMC, whose moves are made from
# the kernel's statements: a kernel of random-walk moves (proposal kernels:
# the SMC's copies once pointed at the original's) and the cached factors of
# the new state (the path changes the latent variables, read by factors that
# do not read theta)
test_that("annealed moves: a random-walk kernel, and the new state's factors", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  model = annealed_model("annealed(norm_rw(theta, 1), bridges = 2, kernel = g)", T = 5, iterations = 100)
  model$lines = sub("gibbs(z[t])", "metropolis(norm_rw(z[t], 0.5))", model$lines, fixed = TRUE)
  writeLines(model$lines, f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
              initial_values = list(list(theta = 1.5, z = rep(1.5, model$T))))
  theta = ilike::load_mcmc_output("annealed_output")
  theta = theta$Value[theta$ParameterName == "theta"]
  expect_true(all(is.finite(theta)) && (length(unique(theta)) > 10))
  # gibbs moves on z after an annealed move read the factors of y at the new z
  model = annealed_model("annealed(norm_rw(theta, 1), bridges = 2, kernel = g)", T = 5, iterations = 100)
  writeLines(model$lines, f)
  unlink("annealed_output", recursive = TRUE)
  expect_no_error(ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
                              initial_values = list(list(theta = 1.5, z = rep(1.5, model$T)))))
})

# (A2) the path runs the chain's own kernel moves: those a chain adapts, with
# the values reached so far (fixed during each path). A random-walk kernel
# started at a poor scale and adapted gives the annealed move the acceptance
# rate of a good fixed scale (with copies of the moves at their starting
# scale, as before A2, it stayed near the poor one's: 0.29 against 0.27)
test_that("annealed moves: the kernel's moves are the chain's, with their adapted scales", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  acceptance = function(kernel_move) {
    model = annealed_model("annealed(norm_rw(theta, 1), bridges = 5, kernel = g)", T = 20, iterations = 3000)
    model$lines = sub("gibbs(z[t])", kernel_move, model$lines, fixed = TRUE)
    writeLines(model$lines, f)
    unlink("annealed_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
                initial_values = list(list(theta = 1.5, z = rep(1.5, model$T))))
    theta = ilike::load_mcmc_output("annealed_output")
    theta = theta$Value[theta$ParameterName == "theta"][-(1:1000)]
    mean(diff(theta) != 0)
  }
  poor = acceptance("metropolis(norm_rw(z[t], 20))")
  adapted = acceptance("metropolis(norm_rw(z[t], 20), adapt = scale())")
  good = acceptance("metropolis(norm_rw(z[t], 1.5))")
  expect_gt(adapted, poor + 0.03)
  expect_lt(abs(adapted - good), 0.02)
})

# (A2) in SMC_with_MCMC, a kernel tuned by esjd: the path's particle carries
# the outer particle's tuning (each move in the path reads that particle's
# own scale), and the path's SMC tunes nothing itself
test_that("annealed moves in SMC with a kernel tuned by esjd", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  set.seed(1); T = 10; z = rnorm(T, 1.5, 1); y = round(rnorm(T, z, 1), 3)
  writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  theta ~ norm(0, 10)",
               paste0("  for (t in 1:", T, ") z[t] ~ norm(theta, 1)"), paste0("  for (t in 1:", T, ") y[t] ~ norm(z[t], 1)"), "}",
               "algorithm smc {", paste0("  for (t in 1:", T, ") move g = metropolis(norm_rw(z[t], 1), smc_adapt = esjd(scale = [0.1, 5]))"),
               "  move m = annealed(norm_rw(theta, 1), bridges = 3, kernel = g)", "  sequence = annealing(power, [0, 0.25, 0.5, 0.75, 1])",
               "  particles = 200", "  termination = iterations(2)", "}"), f)
  ilike::SMC_with_MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1)
  o = ilike::load_smc_output("annealed_output")
  final = o[(o$Target == max(o$Target)) & (o$ParameterName == "theta"), ]
  w = exp(final$LogWeight - max(final$LogWeight))
  prec = 1/100 + T/2
  expect_equal(sum(w*final$Value)/sum(w), (sum(y)/2)/prec, tolerance = 0.05)
})

# (A2) an estimated factor among those that read theta is refused until the
# state-space form (AIS with conditional SMC, U8)
test_that("annealed moves: refused with an estimated factor", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = [1.3, 0.2, 0.9, 1.7, 0.4] }", "model {", "  theta ~ norm(0, 3)", "  z ~ norm(0, 1)",
               "  factor sim = sbi(euclidean_uniform_abc(number_of_points = 20, tolerance_variable = epsilon, tolerance = 0.5),",
               "                   simulate = cpp::simulate_y, data = y, summaries = cpp::summaries)", "}", "cpp {",
               "Data simulate_y(RandomNumberGenerator &rng, const Parameters &parameters)",
               "{ arma::colvec y(5); for (arma::uword i = 0; i < 5; ++i) y[i] = rnorm(rng, parameters[\"theta\"][0], 1.0); return Data(\"y\", y); }",
               "Data summaries(const Data &data) { arma::colvec s(1); s[0] = arma::mean(arma::vectorise(data[\"y\"])); return Data(\"y\", s); }",
               "}", "algorithm mcmc {", "  move g = metropolis(norm_rw(z, 1))",
               "  move m = annealed(norm_rw(theta, 1), bridges = 2, kernel = g)", "  termination = iterations(5)", "}"), f)
  expect_error(ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
                           initial_values = list(list(theta = 0.9, z = 0))),
               "annealed moves with estimated factors are not available yet", fixed = TRUE)
})

test_that("annealed moves: the posterior of the Gaussian latent model, and acceptance rising with the bridges", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  accept = numeric(0)
  for (K in c(0, 5)) {
    model = annealed_model(paste0("annealed(norm_rw(theta, 1), bridges = ", K, ", kernel = g)"))
    writeLines(model$lines, f)
    unlink("annealed_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
                initial_values = list(list(theta = 0, z = rep(0, model$T))))
    o = ilike::load_mcmc_output("annealed_output")
    a = o$Value[o$ParameterName == "theta"][-(1:500)]
    expect_equal(mean(a), model$mean, tolerance = 0.03)
    expect_equal(sd(a), model$sd, tolerance = 0.1)
    accept = c(accept, mean(diff(a) != 0))
  }
  # (about 0.18 with z fixed, 0.23 with five bridges)
  expect_gt(accept[2], accept[1] + 0.03)
})

test_that("annealed moves with a Langevin proposal (the Metropolis-Hastings form): the posterior", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "annealed_output"), recursive = TRUE), add = TRUE)
  model = annealed_model("annealed(langevin(theta, [[1]], 0.15), bridges = 5, kernel = g)")
  writeLines(model$lines, f)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "annealed_output", seed = 1,
              initial_values = list(list(theta = 0, z = rep(0, model$T))))
  o = ilike::load_mcmc_output("annealed_output")
  a = o$Value[o$ParameterName == "theta"][-(1:500)]
  expect_equal(mean(a), model$mean, tolerance = 0.03)
  expect_equal(sd(a), model$sd, tolerance = 0.1)
})
