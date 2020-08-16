test_that("importance_sample works", {

  data("gaussian")

  simulate_prior = function()
  {
    rgamma(1,1,1)
  }

  evaluate_log_likelihood = function(inputs,data)
  {
    sum(dnorm(data,mean=0,sd=inputs,log=TRUE))
  }

  number_of_points = 10

  model = list(simulate_prior = simulate_prior,
               evaluate_log_likelihood = evaluate_log_likelihood,
               data = gaussian)

  set.seed(1)

  output_r = importance_sample(number_of_points,
                               model)

  expect_equal(signif(output_r$log_normalising_constant, digits = 7), signif(-142.5295, digits = 7))
})
