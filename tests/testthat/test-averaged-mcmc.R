# (C312, dev/plan.md 3.22 U5b) MHAAR-RB: move m = averaged(proposal, latent = z,
# points = M, refresh = true), a joint move on the parameters and a latent
# variable whose elements are independent given them (Andrieu, Yildirim,
# Doucet & Chopin, arXiv 2101.01253, Algorithm 4)

averaged_model = function(move, T = 50, latent = NULL, data = NULL, iterations = 5000)
{
  set.seed(1); z = rnorm(T, 1.5, 1); y = round(rnorm(T, z, 1), 3)
  prec = 1/100 + T/2
  list(lines = c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  theta ~ norm(0, 10)",
                 if (is.null(latent)) paste0("  for (t in 1:", T, ") z[t] ~ norm(theta, 1)") else latent,
                 if (is.null(data)) paste0("  for (t in 1:", T, ") y[t] ~ norm(z[t], 1)") else data, "}",
                 "algorithm mcmc {", paste0("  move m = ", move), paste0("  termination = iterations(", iterations, ")"), "}"),
       y = y, mean = (sum(y)/2)/prec, var = 1/prec, T = T)
}

test_that("averaged moves: the refusals", {
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "averaged_output"), recursive = TRUE), add = TRUE)
  move = "averaged(norm_rw(theta, 0.3), latent = z, points = 3)"
  writeLines(averaged_model(move, T = 10, latent = c("  z[1] ~ norm(theta, 1)", "  for (t in 2:10) z[t] ~ norm(z[t - 1], 1)"))$lines, f)
  expect_error(ilike::compile(f), "reads z itself", fixed = TRUE)
  writeLines(averaged_model("averaged(norm_rw(theta, 0.3), latent = z, points = 1)", T = 10)$lines, f)
  expect_error(ilike::compile(f), "a whole number of at least 2", fixed = TRUE)
  writeLines(averaged_model("averaged(norm_rw(z, 0.3), latent = z, points = 3)", T = 10)$lines, f)
  expect_error(ilike::compile(f), "the proposal changes z", fixed = TRUE)
  # a data term reading two elements of z
  model = averaged_model(move, T = 10, data = "  for (t in 1:9) y[t] ~ norm(z[t] + z[t + 1], 1)", iterations = 5)
  model$lines[2] = paste0("data { y = [", paste(model$y[1:9], collapse = ", "), "] }")
  writeLines(model$lines, f)
  expect_error(ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "averaged_output", seed = 1,
                           initial_values = list(list(theta = 0, z = rep(0, 10)))),
               "has a term reading more than one element of z", fixed = TRUE)
})

test_that("averaged moves: the posterior of theta and z in the Gaussian latent model, and of a Bernoulli latent model", {
  skip_on_cran()
  f = tempfile(fileext = ".ilike")
  on.exit(unlink(c(f, "averaged_output"), recursive = TRUE), add = TRUE)
  for (refresh in c("true", "false")) {
    model = averaged_model(paste0("averaged(norm_rw(theta, 0.3), latent = z, points = 5, refresh = ", refresh, ")"))
    writeLines(model$lines, f)
    unlink("averaged_output", recursive = TRUE)
    ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "averaged_output", seed = 1,
                initial_values = list(list(theta = 0, z = rep(0, model$T))))
    o = ilike::load_mcmc_output("averaged_output")
    a = o$Value[(o$ParameterName == "theta") & (o$Iteration > 500)]
    expect_equal(mean(a), model$mean, tolerance = 0.03)
    expect_equal(sd(a), sqrt(model$var), tolerance = 0.1)
    # z_t given y: mean (E[theta | y] + y_t)/2, sd sqrt(1/2 + var(theta | y)/4)
    zz = o[(o$ParameterName == "z") & (o$Iteration > 500), ]
    expect_lt(max(abs(tapply(zz$Value, zz$Dimension, mean) - (model$mean + model$y)/2)), 0.15)
    expect_equal(mean(tapply(zz$Value, zz$Dimension, sd)), sqrt(0.5 + model$var/4), tolerance = 0.05)
  }
  # a discrete latent: z_t ~ bern(theta), y_t ~ norm(z_t, 1)
  T = 30; set.seed(2); z = rbinom(T, 1, 0.3); y = round(rnorm(T, z, 1), 3)
  grid = seq(0.0005, 0.9995, by = 0.001)
  lp = sapply(grid, function(th) sum(log(th * dnorm(y, 1, 1) + (1 - th) * dnorm(y, 0, 1)))); p = exp(lp - max(lp)); p = p / sum(p)
  writeLines(c("ilike 1.0", paste0("data { y = [", paste(y, collapse = ", "), "] }"), "model {", "  theta ~ unif(0, 1)",
               paste0("  for (t in 1:", T, ") z[t] ~ bern(theta)"), paste0("  for (t in 1:", T, ") y[t] ~ norm(z[t], 1)"), "}",
               "algorithm mcmc {", "  move m = averaged(norm_rw(theta, 0.2), latent = z, points = 4)", "  termination = iterations(20000)", "}"), f)
  unlink("averaged_output", recursive = TRUE)
  ilike::MCMC(suppressMessages(ilike::compile(f)), results_name = "averaged_output", seed = 1,
              initial_values = list(list(theta = 0.5, z = rep(0, T))))
  o = ilike::load_mcmc_output("averaged_output")
  a = o$Value[(o$ParameterName == "theta") & (o$Iteration > 1000)]
  expect_equal(mean(a), sum(p * grid), tolerance = 0.03)
  expect_equal(sd(a), sqrt(sum(p * grid^2) - sum(p * grid)^2), tolerance = 0.1)
})
