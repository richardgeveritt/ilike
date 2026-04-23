test_that("evaluate_log works for various priors", {

  x = 1
  x2 = matrix(x,2,1)
  truth = 0
  truth = truth + dunif(x,log=TRUE)
  truth = truth + dunif(x,-2,2,log=TRUE)
  truth = truth + sum(dunif(x2,matrix(-2,2,1),matrix(2,2,1),log=TRUE))
  truth = truth + dexp(x,2,log=TRUE)
  truth = truth + dexp(x-0.2,2,log=TRUE)
  truth = truth + dpois(x,2,log=TRUE)
  truth = truth + dnorm(x,log=TRUE)
  truth = truth + dnorm(x,3,2,log=TRUE)
  truth = truth + sum(dnorm(x2,3,2,log=TRUE))
  truth = truth + sum(dnorm(x2,matrix(3,2,1),2,log=TRUE))
  truth = truth + sum(dnorm(x2,3,matrix(2,2,1),log=TRUE))
  truth = truth + sum(dnorm(x2,matrix(3,2,1),matrix(2,2,1),log=TRUE))
  truth = truth - 0.1664315 # The result of dtnorm(1,0,1,0.5,2,log=TRUE) from the msm package
  truth = truth - 0.5169834 # The result of dtnorm(1,3,2,0.5,2,log=TRUE) from the msm package
  truth = truth - 1.033967 # The result of sum(dtnorm(x2,3,2,0.5,2,log=TRUE)) from the msm package
  mu = matrix(3,2,1)
  Sigma = 2*diag(2)
  truth = truth - 4.531024 # The result of dmvnorm(as.vector(x2),as.vector(mu),Sigma,log=TRUE) from the mvtnorm package
  truth = truth + dgamma(x,3,2,log=TRUE)
  truth = truth + dlnorm(x,3,2,log=TRUE)
  truth = truth + plnorm(x,3,2)

  expect_equal(signif(evaluate(), digits = 7), signif(truth, digits = 7))
})
