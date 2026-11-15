# Multivariate normal simulation with a correlated covariance. Until C80
# rmvnorm_using_chol returned mu + R z with R = chol(Sigma) upper triangular,
# whose covariance R R^T is not Sigma unless Sigma is diagonal: every Gaussian
# proposal (random walk, Langevin, mirror, independent, log-normal, noise
# kernels) and rmvnorm() drew from the wrong distribution. The earlier tests
# all used diagonal covariances.

test_that("IS with a correlated multivariate normal prior gives the right posterior", {
  skip_on_cran()
  on.exit(unlink("is_mvnorm_correlated", recursive = TRUE), add = TRUE)
  m = c(0.5, -0.5); S = matrix(c(1, 0.8, 0.8, 1), 2)
  recipe = ilike::compile(system.file("extdata", "test_mvnorm_correlated_prior.ilike", package = "ilike"),
                          model_parameter_list = list(m, S))
  ilike::IS(recipe, model_parameter_list = list(m, S), results_name = "is_mvnorm_correlated",
            number_of_importance_points = 20000, seed = 1)
  o = ilike::load_smc_output("is_mvnorm_correlated")
  x = cbind(o$Value[o$Dimension == 1], o$Value[o$Dimension == 2])
  w = exp(o$LogWeight[o$Dimension == 1] - max(o$LogWeight)); w = w/sum(w)
  P = solve(S) + diag(2); V = solve(P); mu = drop(V %*% (solve(S, m) + c(1, 1)))
  est_mean = colSums(w*x); r = sweep(x, 2, est_mean); est_cov = t(r) %*% (w*r)
  expect_equal(est_mean, mu, tolerance = 0.05)
  expect_equal(as.vector(est_cov), as.vector(V), tolerance = 0.1)
})

test_that("independent MH with a correlated multivariate normal proposal targets the posterior", {
  skip_on_cran()
  on.exit(unlink("imh_mvnorm_correlated", recursive = TRUE), add = TRUE)
  o = run_mcmc_fixture("test_gaussian_posterior_vector.ilike", "test_imh_proposal_mvnorm_2d.ilike",
                       "imh_mvnorm_correlated", list(theta = c(1, 1)),
                       model_parameter_list = list(c(1, 1), matrix(c(1, 0.5, 0.5, 1), 2)))
  x = cbind(o$Value[o$Dimension == 1], o$Value[o$Dimension == 2])
  expect_equal(colMeans(x), c(1, 1), tolerance = 0.1)
  expect_equal(apply(x, 2, var), c(0.5, 0.5), tolerance = 0.1)
  expect_lt(abs(cor(x)[1, 2]), 0.1)
})
