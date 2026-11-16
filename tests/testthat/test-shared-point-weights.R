# The parameter estimators fitted to MoveOutputs (used by the Gaussian SMC
# adaptor) share each particle's weight equally between its points
# (src/shared_point_weights.h). Not reachable from R, so checked directly.

test_that("a particle's weight is shared equally between its points", {
  set.seed(1)
  points = list(matrix(rnorm(2), 1, 2), matrix(rnorm(4), 2, 2), matrix(rnorm(6), 3, 2))
  w = c(0.5, 0.3, 0.2)
  out = ilike:::shared_point_weights_check(points, log(w))
  expect_equal(as.vector(exp(out$log_weights)), c(0.5, 0.15, 0.15, rep(0.2/3, 3)))
  expect_equal(sum(exp(out$log_weights)), 1)
  expect_equal(out$points, do.call(rbind, points))
  # the weighted mean of the points is the weighted mean of the particles'
  # own means
  expect_equal(colSums(exp(as.vector(out$log_weights)) * out$points),
               colSums(w * t(sapply(points, colMeans))))
})
