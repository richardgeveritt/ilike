# Tuning MCMC kernels in SMC from the particles' squared jumps (dev/plan.md
# 3.19, C211): smc_adapt = esjd(…), Fearnhead & Taylor (2013) and Buchholz,
# Chopin & Jacob (2021).

esjd_problem = function()
{
  # tempering N(0, I) to a shifted, correlated normal in dimension 10
  d = 10
  set.seed(3)
  m = seq(-1, 1, length = d)
  A = matrix(rnorm(d*d), d)/sqrt(d)
  C = crossprod(A) + diag(d)*0.3
  text = function(M) paste0("[", paste(apply(M, 1, function(r) paste0("[", paste(signif(r, 8), collapse = ", "), "]")), collapse = ", "), "]")
  vector = function(v) paste0("[", paste(signif(v, 8), collapse = ", "), "]")
  S = C + diag(d)
  list(identity = text(diag(d)),
       exact = -0.5*d*log(2*pi) - 0.5*determinant(S)$modulus[1] - 0.5*sum(m*solve(S, m)),
       lines = function(move) c("ilike 1.0", paste0("data { y = ", vector(m), " }"), "model {",
                                paste0("  theta ~ mvnorm(", vector(rep(0, d)), ", ", text(diag(d)), ")"),
                                paste0("  y ~ mvnorm(theta, ", text(C), ")"), "}", "algorithm smc {", paste0("  move mv = ", move),
                                "  termination = iterations(10)",
                                "  sequence = annealing(power, [0, 0.0005, 0.002, 0.008, 0.03, 0.07, 0.15, 0.3, 0.5, 0.75, 1])", "}"))
}

esjd_run = function(lines, name, particles = 1000)
{
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, name), recursive = TRUE), add = TRUE)
  writeLines(lines, f)
  unlink(name, recursive = TRUE)
  z = ilike::SMC_with_MCMC(ilike::compile(f), results_name = name, number_of_particles = particles, seed = 1)
  tuning = tryCatch(ilike::load_smc_tuning(name), error = function(e) NULL)
  list(z = z, tuning = tuning)
}

test_that("a fixed esjd value is the kernel with that setting", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "esjd_fixed"), recursive = TRUE), add = TRUE)
  model = function(move) c("ilike 1.0", "data { y = 1 }", "model {", "  theta ~ norm(0, 1)", "  y ~ norm(theta, 0.5)", "}",
                           "algorithm smc {", paste0("  move m = ", move), "  termination = iterations(3)",
                           "  sequence = annealing(power, [0, 0.1, 0.5, 1])", "}")
  run = function(move) { writeLines(model(move), f); unlink("esjd_fixed", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "esjd_fixed", number_of_particles = 200, seed = 2) }
  expect_identical(run("metropolis(norm_rw(theta, 0.5))"),
                   run("metropolis(norm_rw(theta, 1), smc_adapt = esjd(scale = [0.5, 0.5]))"))
  expect_error(ilike::load_smc_tuning("esjd_fixed_none"), "no tuning values")
})

test_that("esjd tunes a random walk's step towards the optimum and keeps the evidence", {
  skip_on_cran()
  problem = esjd_problem()
  r = esjd_run(problem$lines(paste0("metropolis(mvnorm_rw(theta, ", problem$identity, ", 1), ",
                                    "smc_adapt = [gaussian(), esjd(scale = [0.1, 3])])")), "esjd_rw")
  expect_lt(abs(r$z - problem$exact), 0.2)
  # the multiplier of the particles' covariance settles near 2.38/sqrt(d) = 0.75
  final = tail(r$tuning$scale.q50, 1)
  expect_gt(final, 0.5)
  expect_lt(final, 1.1)
  expect_equal(r$tuning$step, seq_len(nrow(r$tuning)))
})

test_that("esjd tunes hmc and nuts and keeps the evidence", {
  skip_on_cran()
  problem = esjd_problem()
  h = esjd_run(problem$lines("metropolis_hastings(hmc(theta, step = 0.1, steps = 5), smc_adapt = esjd(step = [0.01, 1], steps = [1, 20]))"),
               "esjd_hmc")
  expect_lt(abs(h$z - problem$exact), 0.2)
  expect_true(all(c("step.q50", "steps.q50") %in% names(h$tuning)))
  n = esjd_run(problem$lines("nuts(theta, step = 0.1, smc_adapt = esjd(step = [0.01, 1]))"), "esjd_nuts")
  expect_lt(abs(n$z - problem$exact), 0.2)
  # (nuts has no number of steps to tune: it stays at 1)
  expect_true(all(n$tuning$steps.mean == 1))
})

test_that("esjd is refused where it cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "esjd_refused"), recursive = TRUE), add = TRUE)
  run = function(move, extra = character(0))
  {
    writeLines(c("ilike 1.0", "data { y = 1 }", "model {", "  theta ~ norm(0, 1)", "  y ~ norm(theta, 1)", "}", "algorithm smc {",
                 paste0("  move m = ", move), "  termination = iterations(2)", "  sequence = annealing(power, [0, 0.5, 1])", extra, "}"), f)
    unlink("esjd_refused", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "esjd_refused", number_of_particles = 50, seed = 1)
  }
  expect_error(run("metropolis(unif_rw(theta, 1), smc_adapt = esjd(scale = [0.5, 2]))"), "esjd is available for", fixed = TRUE)
  expect_error(run("metropolis_hastings(hmc(theta, step = 0.1, steps = 2), smc_adapt = esjd(scale = [0.5, 2]))"),
               "hmc and nuts are tuned by step", fixed = TRUE)
  expect_error(run("metropolis(norm_rw(theta, 1), smc_adapt = esjd(step = [0.5, 2]))"), "tuned by scale", fixed = TRUE)
  expect_error(run("metropolis(norm_rw(theta, 1), smc_adapt = esjd(scale = [2, 1]))"), "0 < lower <= upper", fixed = TRUE)
})

test_that("pretune tunes hmc from trial moves and keeps the evidence", {
  # (dev/plan.md 3.19 stage S3, C223)
  skip_on_cran()
  problem = esjd_problem()
  r = esjd_run(problem$lines("metropolis_hastings(hmc(theta, step = 0.1, steps = 5), smc_adapt = pretune(max_steps = 20))"),
               "esjd_pretune")
  expect_lt(abs(r$z - problem$exact), 0.2)
  expect_true(all(c("step.q50", "steps.q50", "step_max", "steps_max") %in% names(r$tuning)))
  # eps* follows the targets as they narrow along the tempering
  expect_lt(tail(r$tuning$step_max, 1), r$tuning$step_max[1])
  expect_true(all(r$tuning$steps.q50 >= 1 & r$tuning$steps.q50 <= 25))
  # the rule for L_max can be changed: with steps_change = 0 it stays as given
  r = esjd_run(problem$lines("metropolis_hastings(hmc(theta, step = 0.1, steps = 5), smc_adapt = pretune(max_steps = 12, steps_change = 0))"),
               "esjd_pretune", particles = 300)
  expect_true(all(r$tuning$steps_max == 12))
  expect_true(all(r$tuning$steps.q90 <= 12))
})

test_that("pretune is refused where it cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "esjd_pretune_refused"), recursive = TRUE), add = TRUE)
  run = function(move)
  {
    writeLines(c("ilike 1.0", "data { y = 1 }", "model {", "  theta ~ norm(0, 1)", "  y ~ norm(theta, 1)", "}", "algorithm smc {",
                 paste0("  move m = ", move), "  termination = iterations(2)", "  sequence = annealing(power, [0, 0.5, 1])", "}"), f)
    unlink("esjd_pretune_refused", recursive = TRUE)
    ilike::SMC_with_MCMC(ilike::compile(f), results_name = "esjd_pretune_refused", number_of_particles = 50, seed = 1)
  }
  expect_error(run("nuts(theta, step = 0.1, smc_adapt = pretune())"), "for nuts, use smc_adapt = esjd", fixed = TRUE)
  expect_error(run("metropolis(norm_rw(theta, 1), smc_adapt = pretune())"), "pretune is available for ilike::hmc", fixed = TRUE)
})

test_that("choose learns which proposal moves the particles best, and keeps the evidence", {
  # (dev/plan.md 3.19 stage S4, C224)
  skip_on_cran()
  problem = esjd_problem()
  # a random walk that barely moves against one that moves well: the second takes over
  # (the share left on the first is about change = 0.1)
  r = esjd_run(problem$lines(paste0("metropolis_hastings(choose([mvnorm_rw(theta, ", problem$identity, ", 0.000001), ",
                                    "mvnorm_rw(theta, ", problem$identity, ", 0.1)]))")), "esjd_choose")
  expect_lt(abs(r$z - problem$exact), 0.2)
  expect_gt(tail(r$tuning$kernel2.share, 1), 0.8)
  # a random walk tuned by esjd against hmc (cost 11): each proposal's own tuning alongside the choice
  r = esjd_run(problem$lines(paste0("metropolis_hastings(choose([mvnorm_rw(theta, ", problem$identity, ", 0.1), ",
                                    "hmc(theta, step = 0.15, steps = 10)], smc_adapt = [esjd(scale = [0.1, 3]), none]))")),
               "esjd_choose")
  expect_lt(abs(r$z - problem$exact), 0.2)
  expect_equal(sort(unique(r$tuning$tuner)), 1:2)
  expect_true(all(c("kernel1.share", "kernel2.share", "scale.q50") %in% names(r$tuning)))
})

test_that("choose is refused where it cannot be used", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "esjd_choose_refused"), recursive = TRUE), add = TRUE)
  model = function(move, algorithm = "smc")
    c("ilike 1.0", "data { y = 1 }", "model {", "  theta ~ norm(0, 1)", "  b ~ norm(0, 1)", "  y ~ norm(theta + b, 1)", "}",
      paste0("algorithm ", algorithm, " {"), paste0("  move m = ", move), "  termination = iterations(2)",
      if (algorithm == "smc") "  sequence = annealing(power, [0, 0.5, 1])", "}")
  compile_only = function(move) { writeLines(model(move), f); ilike::compile(f) }
  expect_error(compile_only("metropolis(choose([norm_rw(theta, 1), norm_rw(theta, 0.1)]))"), "goes in metropolis_hastings", fixed = TRUE)
  expect_error(compile_only("metropolis_hastings(choose([norm_rw(theta, 1), norm_rw(b, 0.1)]))"), "move the same variables", fixed = TRUE)
  expect_error(compile_only("metropolis_hastings(choose([norm_rw(theta, 1), norm_rw(theta, 0.1)]), smc_adapt = esjd(scale = [0.5, 2]))"),
               "give each proposal's smc_adapt inside it", fixed = TRUE)
  expect_error(compile_only("metropolis_hastings(choose([norm_rw(theta, 1), norm_rw(theta, 0.1)], cost = [1]))"),
               "a positive number per proposal", fixed = TRUE)
  # outside SMC no particle carries a choice
  writeLines(model("metropolis_hastings(choose([norm_rw(theta, 1), norm_rw(theta, 0.1)]))", "mcmc"), f)
  unlink("esjd_choose_refused", recursive = TRUE)
  expect_error(ilike::MCMC(ilike::compile(f), results_name = "esjd_choose_refused", initial_values = list(list(theta = 0, b = 0))),
               "used in SMC_with_MCMC", fixed = TRUE)
})

test_that("esjd tunes moves on elements, one multiplier of every element's step per particle", {
  # (C226)
  skip_on_cran()
  problem = esjd_problem()
  lines = problem$lines("X")
  lines[grep("move mv = X", lines)] = "  for (k in 1:10) move mv = metropolis(norm_rw(theta[k], 0.3), smc_adapt = esjd(scale = [0.1, 10]))"
  r = esjd_run(lines, "esjd_elements")
  expect_lt(abs(r$z - problem$exact), 0.2)
  # the step given (0.3) is too small: the multiplier grows
  expect_gt(tail(r$tuning$scale.q50, 1), 2)
})
