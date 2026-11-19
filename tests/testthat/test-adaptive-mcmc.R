# Adaptive MCMC (dev/plan.md 3.2; Andrieu & Thoms 2008, Section 5.1): an
# adaptor block after a proposal adapts that proposal along each chain.
#   scaled_am: Algorithm 4 (scale by Robbins-Monro towards a target acceptance
#              rate, mean and covariance from the chain)
#   am, rb_am: Algorithms 2 and 3 (covariance from the chain, scale 2.38^2/d)
#   scale:     the scale recursion alone

adaptive_run <- function(model, algorithm_lines, model_parameter_list, initial_values,
                         seed = 1, results_name = basename(tempfile("adaptive")), ...)
{
  algorithm = tempfile(fileext = ".ilike")
  writeLines(algorithm_lines, algorithm)
  recipe = suppressMessages(ilike::compile(c(system.file("extdata", model, package = "ilike"), algorithm),
                                           model_parameter_list = model_parameter_list))
  suppressMessages(ilike::MCMC(recipe, model_parameter_list = model_parameter_list,
                               results_name = results_name, seed = seed,
                               initial_values = initial_values, ...))
  o = ilike::load_mcmc_output(results_name)
  attr(o, "adaptation") = tryCatch(ilike::load_mcmc_adaptation(results_name), error = function(e) NULL)
  unlink(results_name, recursive = TRUE)
  o
}

# the second half of chain k of a 2-D parameter, one row per iteration
second_half <- function(o, chain = 1, parameter = "theta")
{
  if ("Chain" %in% names(o))
    o = o[o$Chain == chain, ]
  o = o[o$ParameterName == parameter, ]
  x = cbind(o$Value[o$Dimension == 1], o$Value[o$Dimension == 2])
  x[(nrow(x)/2 + 1):nrow(x), , drop = FALSE]
}

acceptance_rate <- function(x) mean(rowSums(abs(diff(x))) > 0)

# correlated 2-D Gaussian posterior (test_mvnorm_correlated_prior.ilike)
m = c(0.5, -0.5); S = 4*matrix(c(1, 0.95, 0.95, 1), 2)
V = solve(solve(S) + diag(2)); mu = drop(V %*% (solve(S, m) + c(1, 1)))
correlated_lines <- function(adaptor, n = 40000)
  c("/***m_proposal,ilike::mvnorm_rw(theta,p3,p4)***/", "", adaptor, "",
    sprintf("/***mcmc_termination,ilike::iterations(%d)***/", n))

test_that("scaled_am (Algorithm 4) reaches the target acceptance rate and the posterior from badly scaled starts", {
  skip_on_cran()
  for (start in c(100, 0.01))
  {
    pl = list(m, S, start*diag(2), 1)
    o = adaptive_run("test_mvnorm_correlated_prior.ilike",
                     correlated_lines("/***m_adaptor,ilike::scaled_am(theta)***/"),
                     pl, list(list(theta = c(0, 0))))
    x = second_half(o)
    expect_equal(acceptance_rate(x), 0.234, tolerance = 0.03/0.234)
    expect_lt(max(abs(colMeans(x) - mu)), 0.1)
    expect_equal(as.vector(cov(x)), as.vector(V), tolerance = 0.15)
  }
})

test_that("am and rb_am (Algorithms 2 and 3) learn the posterior covariance", {
  skip_on_cran()
  for (adaptor in c("am", "rb_am"))
  {
    pl = list(m, S, 100*diag(2), 1)
    o = adaptive_run("test_mvnorm_correlated_prior.ilike",
                     correlated_lines(paste0("/***m_adaptor,ilike::", adaptor, "(theta)***/")),
                     pl, list(list(theta = c(0, 0))))
    x = second_half(o)
    expect_lt(max(abs(colMeans(x) - mu)), 0.1)
    expect_equal(as.vector(cov(x)), as.vector(V), tolerance = 0.15)
    # with the covariance learnt and scale 2.38^2/d, a 2-D Gaussian accepts
    # about 35% of proposals
    expect_gt(acceptance_rate(x), 0.25)
    expect_lt(acceptance_rate(x), 0.45)
  }
})

test_that("the scale adaptor alone tunes random walk, Langevin and Barker proposals", {
  skip_on_cran()
  pl = list(m, S, 0.01*diag(2), 1)
  o = adaptive_run("test_mvnorm_correlated_prior.ilike",
                   correlated_lines("/***m_adaptor,ilike::scale(theta,0.3)***/"),
                   pl, list(list(theta = c(0, 0))))
  expect_equal(acceptance_rate(second_half(o)), 0.3, tolerance = 0.03/0.3)

  # test_gaussian_posterior_vector.ilike: each component N(1, 1/2); default
  # targets 0.574 (Langevin) and 0.4 (Barker)
  for (case in list(c("langevin", 0.574), c("barker", 0.4)))
  {
    lines = c(sprintf("/***mh_proposal,ilike::%s(theta,p1,p2)***/", case[1]), "",
              "/***mh_adaptor,ilike::scale(theta)***/", "",
              "/***mcmc_termination,ilike::iterations(20000)***/")
    o = adaptive_run("test_gaussian_posterior_vector.ilike", lines,
                     list(diag(2), 25), list(list(theta = c(1, 1))))
    x = second_half(o)
    expect_equal(acceptance_rate(x), as.numeric(case[2]), tolerance = 0.05/as.numeric(case[2]))
    expect_equal(colMeans(x), c(1, 1), tolerance = 0.1)
    expect_equal(apply(x, 2, var), c(0.5, 0.5), tolerance = 0.15)
  }
})

test_that("each move of a Gibbs sweep adapts its own proposal", {
  skip_on_cran()
  lines = c("/***m_proposal,ilike::norm_rw(theta1,10)***/", "",
            "/***m_adaptor,ilike::scaled_am(theta1,0.44)***/", "",
            "/***m_factor_index,c(1,2)***/", "",
            "/***m_proposal,ilike::norm_rw(theta2,0.01)***/", "",
            "/***m_adaptor,ilike::scaled_am(theta2,0.44)***/", "",
            "/***m_factor_index,c(3,4)***/", "",
            "/***mcmc_termination,ilike::iterations(20000)***/")
  o = adaptive_run("test_gaussian_posterior_two_params.ilike", lines, list(),
                   list(list(theta1 = 1, theta2 = 1)))
  for (p in c("theta1", "theta2"))
  {
    v = o$Value[o$ParameterName == p]
    v = v[(length(v)/2 + 1):length(v)]
    expect_equal(mean(v != c(NA, v[-length(v)]), na.rm = TRUE), 0.44, tolerance = 0.04/0.44)
    expect_equal(mean(v), 1, tolerance = 0.1)
    expect_equal(var(v), 0.5, tolerance = 0.15)
  }
})

test_that("chains adapt separately, and an adaptive run is reproducible from its seed", {
  skip_on_cran()
  pl = list(m, S, 100*diag(2), 1)
  lines = correlated_lines("/***m_adaptor,ilike::scaled_am(theta)***/", n = 10000)
  starts = list(list(theta = c(0, 0)), list(theta = c(3, 3)))
  o1 = adaptive_run("test_mvnorm_correlated_prior.ilike", lines, pl, starts, seed = 5)
  o2 = adaptive_run("test_mvnorm_correlated_prior.ilike", lines, pl, starts, seed = 5)
  expect_identical(o1$Value, o2$Value)
  x1 = second_half(o1, 1); x2 = second_half(o1, 2)
  expect_false(isTRUE(all.equal(x1, x2)))
  for (x in list(x1, x2))
    expect_equal(acceptance_rate(x), 0.234, tolerance = 0.04/0.234)
})

test_that("componentwise_am (Algorithm 5) gives each component the target acceptance rate", {
  skip_on_cran()
  # posterior variances 0.0099 and 0.99: one proposal variance cannot suit
  # both components; start with the same variance (1) for both
  Sd = diag(c(0.01, 100))
  Vd = solve(solve(Sd) + diag(2)); mud = drop(Vd %*% (solve(Sd, m) + c(1, 1)))
  lines = c("/***m_proposal,ilike::componentwise_rw(theta,p3)***/", "",
            "/***m_adaptor,ilike::componentwise_am(theta)***/", "",
            "/***mcmc_termination,ilike::iterations(40000)***/")
  o = adaptive_run("test_mvnorm_correlated_prior.ilike", lines, list(m, Sd, diag(2)),
                   list(list(theta = c(0, 0))))
  x = second_half(o)
  # each component is chosen half the time and accepted 44% of those times
  moved = colMeans(abs(diff(x)) > 0)
  expect_equal(moved, c(0.22, 0.22), tolerance = 0.025/0.22)
  expect_lt(max(abs(colMeans(x) - mud)), 0.1)
  expect_equal(apply(x, 2, var), diag(Vd), tolerance = 0.15)
})

test_that("the adaptation is written to the output: a trace per stored iteration and the final state per chain", {
  skip_on_cran()
  pl = list(m, S, 100*diag(2), 1)
  starts = list(list(theta = c(0, 0)), list(theta = c(3, 3)))
  o = adaptive_run("test_mvnorm_correlated_prior.ilike",
                   correlated_lines("/***m_adaptor,ilike::scaled_am(theta)***/", n = 10000),
                   pl, starts, thin = 10)
  a = attr(o, "adaptation")
  expect_setequal(unique(a$trace$Parameter), c("theta.scale", "theta.mean[1]", "theta.mean[2]"))
  # one trace row per stored point (1000 per chain with thin = 10)
  expect_equal(nrow(a$trace), 2*1000*3)
  expect_equal(max(a$trace$Iteration), 1000)
  expect_setequal(unique(a$final$Parameter),
                  c("theta.scale", "theta.mean[1]", "theta.mean[2]",
                    "theta.covariance[1,1]", "theta.covariance[1,2]", "theta.covariance[2,2]"))
  expect_equal(sort(unique(a$final$Chain)), 1:2)
  # the final scale is the trace's last value, and it has moved away from 1
  last = a$trace[a$trace$Parameter == "theta.scale" & a$trace$Iteration == 1000, ]
  final = a$final[a$final$Parameter == "theta.scale", ]
  expect_equal(last$Value[order(last$Chain)], final$Value[order(final$Chain)])
  expect_true(all(final$Value > 2))

  # several moves: names prefixed by the move
  lines = c("/***m_proposal,ilike::norm_rw(theta1,10)***/", "",
            "/***m_adaptor,ilike::scale(theta1,0.44)***/", "",
            "/***m_factor_index,c(1,2)***/", "",
            "/***m_proposal,ilike::norm_rw(theta2,10)***/", "",
            "/***m_factor_index,c(3,4)***/", "",
            "/***mcmc_termination,ilike::iterations(100)***/")
  o = adaptive_run("test_gaussian_posterior_two_params.ilike", lines, list(),
                   list(list(theta1 = 1, theta2 = 1)))
  expect_equal(unique(attr(o, "adaptation")$trace$Parameter), "move1.theta1.scale")
})

test_that("the mirror proposal adapts its mean, covariance and scale", {
  skip_on_cran()
  # a poor start: mean (5, 5), covariance 10 I
  lines = c("/***m_proposal,ilike::mirror(theta,p3,p4,p5)***/", "",
            "/***m_adaptor,ilike::scaled_am(theta)***/", "",
            "/***mcmc_termination,ilike::iterations(20000)***/")
  o = adaptive_run("test_mvnorm_correlated_prior.ilike", lines,
                   list(m, S, c(5, 5), 10*diag(2), 1), list(list(theta = c(0, 0))))
  x = second_half(o)
  # default target for a 2-D mirror proposal: 0.42 (found by experiment, C171)
  expect_equal(acceptance_rate(x), 0.42, tolerance = 0.03/0.42)
  expect_lt(max(abs(colMeans(x) - mu)), 0.1)
  expect_equal(as.vector(cov(x)), as.vector(V), tolerance = 0.15)
  final = attr(o, "adaptation")$final
  expect_lt(max(abs(final$Value[final$Parameter %in% c("theta.mean[1]", "theta.mean[2]")] - mu)), 0.2)
})

test_that("pc_am and online_pc_am (Algorithm 8) learn the principal directions", {
  skip_on_cran()
  # eigenvectors of the posterior covariance: (1, 1) and (1, -1) over sqrt(2)
  for (adaptor in c("pc_am", "online_pc_am"))
  {
    lines = c("/***m_proposal,ilike::pc_rw(theta,p3)***/", "",
              paste0("/***m_adaptor,ilike::", adaptor, "(theta)***/"), "",
              "/***mcmc_termination,ilike::iterations(40000)***/")
    o = adaptive_run("test_mvnorm_correlated_prior.ilike", lines, list(m, S, diag(2)),
                     list(list(theta = c(0, 0))))
    x = second_half(o)
    expect_equal(acceptance_rate(x), 0.44, tolerance = 0.03/0.44)
    expect_lt(max(abs(colMeans(x) - mu)), 0.1)
    expect_equal(as.vector(cov(x)), as.vector(V), tolerance = 0.15)
    final = attr(o, "adaptation")$final
    w1 = final$Value[final$Parameter %in% c("theta.direction[1,1]", "theta.direction[2,1]")]
    expect_gt(abs(sum(w1*c(1, 1)/sqrt(2))), 0.97)
  }
})

test_that("named proposal parameters adapt, by the acceptance rate or by the user's function", {
  skip_on_cran()
  # test_gaussian_posterior.ilike: posterior N(1, 1/2); the proposal's sd is
  # the proposal parameter "sd", given far too large
  proposal = readLines(system.file("extdata", "test_m_proposal_proposal_parameters.ilike", package = "ilike"))
  adaptations = list(
    "/***m_adaptor,ilike::acceptance_rate(sd,0.44)***/",
    c("/***adapt_m_proposal***/",
      "Parameters adapt(const Parameters &proposal_parameters, const Parameters &parameters,",
      "                 const Parameters &proposed_parameters, const Parameters &new_parameters,",
      "                 double acceptance_probability) {",
      "  return Parameters(\"sd\", proposal_parameters[\"sd\"]*(acceptance_probability - 0.44));",
      "}"))
  for (adaptation in adaptations)
  {
    o = adaptive_run("test_gaussian_posterior.ilike", c(proposal, adaptation), list(),
                     list(list(theta = 1)), algorithm_parameter_list = list(sd = 50))
    v = o$Value[(nrow(o)/2 + 1):nrow(o)]
    expect_equal(mean(v != c(NA, v[-length(v)]), na.rm = TRUE), 0.44, tolerance = 0.03/0.44)
    expect_equal(mean(v), 1, tolerance = 0.1)
    expect_equal(var(v), 0.5, tolerance = 0.15)
    final = attr(o, "adaptation")$final
    expect_equal(final$Parameter, "sd[1]")
    expect_lt(final$Value, 10)
  }
})

test_that("Polyak-Juditsky averaging (option averaged) gives a much less noisy adapted covariance", {
  skip_on_cran()
  # with kappa = 0.6 the recursion's own covariance is noisy (up to 30-40 %
  # off after 20,000 iterations); the average of its iterates is not
  pl = list(m, S, 100*diag(2), 1)
  for (seed in 1:2)
  {
    o = adaptive_run("test_mvnorm_correlated_prior.ilike",
                     correlated_lines("/***m_adaptor,ilike::scaled_am(theta,averaged)***/", n = 20000),
                     pl, list(list(theta = c(0, 0))), seed = seed)
    x = second_half(o)
    expect_equal(acceptance_rate(x), 0.234, tolerance = 0.03/0.234)
    expect_lt(max(abs(colMeans(x) - mu)), 0.1)
    final = attr(o, "adaptation")$final
    adapted = final$Value[match(c("theta.covariance[1,1]", "theta.covariance[1,2]", "theta.covariance[2,2]"), final$Parameter)]
    expect_equal(adapted, V[c(1, 3, 4)], tolerance = 0.1)
  }
  # the other adaptors take the option too
  o = adaptive_run("test_mvnorm_correlated_prior.ilike",
                   c("/***m_proposal,ilike::pc_rw(theta,p3)***/", "",
                     "/***m_adaptor,ilike::pc_am(theta,averaged)***/", "",
                     "/***mcmc_termination,ilike::iterations(20000)***/"),
                   list(m, S, diag(2)), list(list(theta = c(0, 0))))
  expect_equal(acceptance_rate(second_half(o)), 0.44, tolerance = 0.04/0.44)
})

test_that("named proposal parameters adapt separately in each chain, with or without averaging", {
  skip_on_cran()
  proposal = readLines(system.file("extdata", "test_m_proposal_proposal_parameters.ilike", package = "ilike"))
  H = c("/***adapt_m_proposal***/",
        "Parameters adapt(const Parameters &proposal_parameters, const Parameters &parameters,",
        "                 const Parameters &proposed_parameters, const Parameters &new_parameters,",
        "                 double acceptance_probability) {",
        "  return Parameters(\"sd\", proposal_parameters[\"sd\"]*(acceptance_probability - 0.44));",
        "}")
  cases = list("/***m_adaptor,ilike::acceptance_rate(sd,0.44)***/",
               "/***m_adaptor,ilike::acceptance_rate(sd,0.44,averaged)***/",
               H, c(H, "/***m_adaptor,ilike::averaged()***/"))
  finals = list()
  for (k in seq_along(cases))
  {
    o = adaptive_run("test_gaussian_posterior.ilike", c(proposal, cases[[k]]), list(),
                     list(list(theta = 1), list(theta = 3)), algorithm_parameter_list = list(sd = 50))
    final = attr(o, "adaptation")$final
    final = final$Value[order(final$Chain)]
    # each chain has its own sd (they used to share the algorithm's matrix)
    expect_false(isTRUE(all.equal(final[1], final[2])))
    expect_true(all(final > 0.5 & final < 5))
    for (chain in 1:2)
    {
      v = o$Value[o$Chain == chain]
      v = v[(length(v)/2 + 1):length(v)]
      expect_equal(mean(v[-1] != v[-length(v)]), 0.44, tolerance = 0.05/0.44)
    }
    finals[[k]] = final
  }
  # averaging changes the result
  expect_false(isTRUE(all.equal(finals[[3]], finals[[4]])))
})

test_that("adaptors that do not fit the proposal are refused with a clear message", {
  skip_on_cran()
  lines = c("/***mh_proposal,ilike::langevin(theta,p1,p2)***/", "",
            "/***mh_adaptor,ilike::am(theta)***/", "",
            "/***mcmc_termination,ilike::iterations(10)***/")
  expect_error(adaptive_run("test_gaussian_posterior_vector.ilike", lines, list(diag(2), 1),
                            list(list(theta = c(1, 1)))),
               "applies to the proposals ilike::norm_rw, ilike::mvnorm_rw, ilike::mirror and ilike::hmc")
  lines[3] = "/***mh_adaptor,ilike::nonsense(theta)***/"
  expect_error(adaptive_run("test_gaussian_posterior_vector.ilike", lines, list(diag(2), 1),
                            list(list(theta = c(1, 1)))),
               "Unknown MCMC adaptor")
  expect_error(adaptive_run("test_mvnorm_correlated_prior.ilike",
                            correlated_lines("/***m_adaptor,ilike::averaged()***/", n = 10),
                            list(m, S, diag(2), 1), list(list(theta = c(0, 0)))),
               "goes with an adapt_m_proposal function")
  expect_error(adaptive_run("test_mvnorm_correlated_prior.ilike",
                            correlated_lines("/***m_adaptor,ilike::componentwise_am(theta)***/", n = 10),
                            list(m, S, diag(2), 1), list(list(theta = c(0, 0)))),
               "applies to the proposal ilike::componentwise_rw")
})
