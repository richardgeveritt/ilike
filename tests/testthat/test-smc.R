test_that("SMC with MCMC moves and fixed annealing gives the right evidence and weighted particles", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("smc_annealing", recursive=TRUE), add=TRUE)

  expect_equal(test_smc_annealing_z("smc_annealing"), annealed_log_evidence(1), tolerance = 0.05)

  # Each target after the first is reached by reweighting the particles from
  # the previous one and storing several MCMC iterations per particle. The
  # weights must be attached to every row of their particle (B9).
  output = ilike::load_smc_output("smc_annealing")
  moments = smc_weighted_moments_by_target(output)
  # The power-0 target (the prior) is not written out: target 1 is the prior
  # sample reweighted to power 0.25.
  powers = c(0.25,0.5,0.75,1)
  target_parameters = unique(output[,c("Target","TargetParameters")])
  expect_equal(target_parameters$TargetParameters, paste("power=",powers,";",sep=""))
  expect_equal(nrow(moments), length(powers))
  expect_equal(moments[,"mean"], rep(1,length(powers)), tolerance = 0.1)
  expect_equal(moments[,"var"], 1/(1+powers), tolerance = 0.1)
})

test_that("SMC with MCMC moves and no smc_sequence defaults to importance sampling", {

  on.exit(unlink("smc_default_is", recursive=TRUE), add=TRUE)

  # Used to segfault when writing the output (B5).
  expect_equal(test_smc_z("test_smc_default_is.ilike", "smc_default_is"), annealed_log_evidence(1), tolerance = 0.05)
  output = ilike::load_smc_output("smc_default_is")
  moments = smc_weighted_moments_by_target(output)
  expect_equal(unname(moments[nrow(moments),"mean"]), 1, tolerance = 0.1)
  expect_equal(unname(moments[nrow(moments),"var"]), 0.5, tolerance = 0.1)
})

test_that("ilike::stable stops SMC once the target parameter changes little for several iterations", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit(unlink("smc_stable", recursive=TRUE), add=TRUE)

  # Used to stop after the first step (B12): stops at power 0.54.
  expect_equal(test_smc_z("test_smc_stable.ilike", "smc_stable"), annealed_log_evidence(0.54), tolerance = 0.05)
  output = ilike::load_smc_output("smc_stable")
  expect_equal(unique(output$TargetParameters), paste("power=",c(0.1,0.5,0.52,0.54),";",sep=""))
})

test_that("ilike::annealing with a factor argument anneals only those factors", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("smc_annealing_factor", recursive=TRUE)
    unlink("smc_annealing_factor_is", recursive=TRUE)
    unlink("smc_annealing_prior", recursive=TRUE)
  }, add=TRUE)

  # Was read as the schedule, annealing to power 2 (B11). Annealing only the
  # likelihood gives the targets prior x likelihood^power.
  powers = c(0.25,0.5,0.75,1)
  for (case in list(c("test_smc_annealing_factor.ilike", "smc_annealing_factor"),
                    c("test_smc_annealing_factor_is_proposal.ilike", "smc_annealing_factor_is")))
  {
    expect_equal(test_smc_z(case[1], case[2]), annealed_log_evidence(1), tolerance = 0.05)
    output = ilike::load_smc_output(case[2])
    moments = smc_weighted_moments_by_target(output)
    expect_equal(unique(output$TargetParameters), paste("power=",powers,";",sep=""))
    expect_equal(moments[,"mean"], rep(1,length(powers)), tolerance = 0.1)
    expect_equal(moments[,"var"], 1/(1+powers), tolerance = 0.1)
  }

  # The prior is the proposal here, so it cancels and cannot be annealed.
  expect_error(test_smc_z("test_smc_annealing_prior_as_proposal.ilike", "smc_annealing_prior"), "cannot be annealed")
})

test_that("adaptive_target ilike::ess chooses each power to give the requested effective sample size", {

  on.exit(unlink("smc_adaptive_ess", recursive=TRUE), add=TRUE)

  # Used to run indefinitely, the bisection reading a stale ESS (B10).
  expect_equal(test_smc_z("test_smc_adaptive_ess.ilike", "smc_adaptive_ess"), annealed_log_evidence(1), tolerance = 0.05)
  output = ilike::load_smc_output("smc_adaptive_ess")
  targets = unique(output[,c("Target","TargetParameters","ISESS")])
  expect_true(nrow(targets) > 1)
  expect_equal(targets$TargetParameters[nrow(targets)], "power=1;")
  # Every intermediate target has ESS 0.95 x 2000 particles; the last may be higher.
  expect_equal(targets$ISESS[-nrow(targets)], rep(1900, nrow(targets)-1), tolerance = 0.01)
  expect_true(targets$ISESS[nrow(targets)] >= 1880)
})

test_that("ilike::parameter sequences work with and without the factor argument", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  on.exit({
    unlink("smc_parameter", recursive=TRUE)
    unlink("smc_parameter_factors", recursive=TRUE)
  }, add=TRUE)

  # The factor argument used to be required; without it every factor is
  # re-evaluated when the sequence variable changes.
  for (n in c("smc_parameter", "smc_parameter_factors"))
  {
    z = test_smc_z(paste("test_", n, ".ilike", sep=""), n, model_file = "test_gaussian_posterior_tau.ilike")
    expect_equal(z, annealed_log_evidence(1), tolerance = 0.05)
    moments = smc_weighted_moments_by_target(ilike::load_smc_output(n))
    expect_equal(unname(moments[nrow(moments),"mean"]), 1, tolerance = 0.1)
    expect_equal(unname(moments[nrow(moments),"var"]), 0.5, tolerance = 0.1)
  }
})

test_that("tempering, an adaptive CESS target and ESS-based resampling give the evidence and posterior", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally

  cases = c("test_smc_tempering.ilike", "test_smc_adaptive_cess.ilike", "test_smc_resample_ess.ilike")
  on.exit(for (k in seq_along(cases)) unlink(paste("smc_option_", k, sep=""), recursive=TRUE), add=TRUE)

  for (k in seq_along(cases))
  {
    name = paste("smc_option_", k, sep = "")
    expect_equal(test_smc_z(cases[k], name), annealed_log_evidence(1), tolerance = 0.05)
    moments = smc_weighted_moments_by_target(ilike::load_smc_output(name))
    expect_equal(unname(moments[nrow(moments), "mean"]), 1, tolerance = 0.1)
    expect_equal(unname(moments[nrow(moments), "var"]), 0.5, tolerance = 0.1)
  }
})

test_that("SMC with MCMC moves runs without writing results", {

  recipe = ilike::compile(system.file("extdata", c("test_gaussian_posterior.ilike", "test_smc_annealing.ilike"),
                                      package = "ilike"))
  z = ilike::SMC_with_MCMC(recipe, number_of_particles = 2000, seed = 1)
  expect_equal(z, annealed_log_evidence(1), tolerance = 0.05)
})

# (C241, dev/plan.md 3.20 S2) data tempering: the data points of a factor
# written as model statements enter along the sequence, each annealed in
data_tempering_file = function(sequence, target = "  target = cess(0.9)", move = "metropolis(norm_rw(mu, 0.3))",
                               likelihood = NULL)
{
  set.seed(3)
  y = round(rnorm(20, 1.5, 1), 3)
  f = tempfile(fileext = ".ilike")
  writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  mu ~ norm(0, 10)",
               if (is.null(likelihood)) "  for (i in 1:20) y[i] ~ norm(mu, 1)" else likelihood, "}",
               "algorithm smc {", paste0("  move = ", move), "  termination = iterations(5)",
               paste0("  sequence = ", sequence), target, "  particles = 1000", "}"), f)
  list(file = f, y = y)
}

test_that("data tempering gives the evidence and the posterior", {
  skip_on_cran()
  d = data_tempering_file("data(fraction, [0, 1])")
  on.exit(unlink(c(d$file, "smc_data_tempering"), recursive = TRUE), add = TRUE)
  z = ilike::SMC_with_MCMC(ilike::compile(d$file), results_name = "smc_data_tempering", seed = 1)
  n = length(d$y); tau = 10
  exact = -n/2 * log(2*pi) - 0.5 * log(1 + n * tau^2) - 0.5 * (sum(d$y^2) - tau^2 * sum(d$y)^2 / (1 + n * tau^2))
  expect_lt(abs(z - exact), 0.2)
  moments = smc_weighted_moments_by_target(ilike::load_smc_output("smc_data_tempering"))
  expect_equal(unname(moments[nrow(moments), "mean"]), sum(d$y) / (n + 1/tau^2), tolerance = 0.05)
  expect_equal(unname(moments[nrow(moments), "var"]), 1 / (n + 1/tau^2), tolerance = 0.2)
})

test_that("data tempering refuses what it cannot do", {
  # whole data points with an adaptive target (the criterion is flat between data points)
  d = data_tempering_file("data(fraction, [0, 1], anneal = false)")
  expect_error(ilike::SMC_with_MCMC(ilike::compile(d$file), seed = 1), "anneal = false", fixed = TRUE)
  # a factor written as a function
  d = data_tempering_file("data(fraction, [0, 1])", likelihood = "  factor lik = r::sum(dnorm(y, parameters$mu, 1, log = TRUE))")
  expect_error(ilike::SMC_with_MCMC(ilike::compile(d$file), seed = 1), "is not written as model statements", fixed = TRUE)
  # gradients include every data point
  d = data_tempering_file("data(fraction, [0, 1])", move = "metropolis_hastings(langevin(mu, 0.1, 1))")
  expect_error(ilike::SMC_with_MCMC(ilike::compile(d$file), seed = 1), "under data tempering", fixed = TRUE)
  unlink(d$file)
})

# (dev/plan.md 3.21 (c6), step 2) with few stored generations, the dropped
# generations' chain outputs are reused by the next moves (Particles::spare,
# SMCMCMCMove::move_reusing): the same random numbers and results as building
# them anew (all generations stored, none dropped)
test_that("SMC moves that reuse dropped generations' outputs give the same run", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = 1 }", "model {", "  theta ~ norm(0, 1)", "  y ~ norm(theta, 0.5)", "}",
               "algorithm smc {", "  move m = metropolis(norm_rw(theta, 0.5))", "  termination = iterations(3)",
               "  sequence = annealing(power, [0, 0.1, 0.3, 0.6, 1])", "}"), f)
  r = suppressMessages(ilike::compile(f))
  z = sapply(c(2, 100), function(k) suppressMessages(ilike::SMC_with_MCMC(r, number_of_particles = 200, smc_iterations_to_store = k, seed = 3)))
  expect_identical(z[1], z[2])
})

# (2026-11-27) annealing from an importance proposal adds q^(1 - p) after the
# model's factors; a move on one of two independent parameters evaluates only
# some of the factors, and its factor index must include q^(1 - p) too (before,
# the particles followed (prior x likelihood)^p alone: at p = 0.5, theta had
# mean 0.48 and variance 1.15 against 0.667 and 0.667, and the evidence was
# -1.62 against -3.03)
test_that("annealing from an importance proposal, with moves on some of the parameters", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "smc_proposal_moves"), recursive = TRUE), add = TRUE)
  writeLines(c("ilike 1.0", "data { y = 1", "  w = 1 }", "model {", "  theta ~ norm(0, 1)", "  phi ~ norm(0, 1)",
               "  y ~ norm(theta, 1)", "  w ~ norm(phi, 1)", "}", "algorithm smc {",
               "  move a = metropolis(norm_rw(theta, 0.5))", "  move b = metropolis(norm_rw(phi, 0.5))", "  termination = iterations(20)",
               "  sequence = annealing(power, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])",
               "  proposal = density(simulate = cpp::draw, log = cpp::log_draw)", "  particles = 4000", "}", "cpp {",
               "Parameters draw(RandomNumberGenerator &rng) { Parameters p; p[\"theta\"] = arma::mat(1, 1, arma::fill::value(rnorm(rng, 1.0, 1.0))); p[\"phi\"] = arma::mat(1, 1, arma::fill::value(rnorm(rng, 1.0, 1.0))); return p; }",
               "double log_draw(const Parameters &p) { return dnorm(p[\"theta\"][0], 1.0, 1.0) + dnorm(p[\"phi\"][0], 1.0, 1.0); }", "}"), f)
  z = suppressMessages(ilike::SMC_with_MCMC(suppressMessages(ilike::compile(f)), results_name = "smc_proposal_moves", seed = 1))
  expect_equal(z, 2*dnorm(1, 0, sqrt(2), log = TRUE), tolerance = 0.02)
  o = ilike::load_smc_output("smc_proposal_moves")
  o = o[o$ParameterName == "theta", ]
  half = o[o$TargetParameters == "power=0.5;", ]
  w = exp(half$LogWeight - max(half$LogWeight)); w = w/sum(w)
  m = sum(w*half$Value)
  # the bridge q^(1/2) (prior x likelihood)^(1/2): N(2/3, 2/3)
  expect_equal(m, 2/3, tolerance = 0.03)
  expect_equal(sum(w*(half$Value - m)^2), 2/3, tolerance = 0.05)
})
