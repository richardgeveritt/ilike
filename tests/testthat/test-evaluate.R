test_that("evaluate_log works for various priors", {

  skip_on_cran()

  Sys.setenv(PKG_CPPFLAGS = paste0("-I", system.file("include", package = "ilike")))

  filename = system.file("extdata", "evaluate_priors.ilike", package = "ilike")

  recipe = compile(filename)

  output = evaluate_log(recipe,
                        parameters = list(x=0.5))

  x = 0.5
  truth = 0
  truth = truth + dunif(x,log=TRUE)
  truth = truth + dunif(x,-2,2,log=TRUE)

  x2 = matrix(0.5,2,1)
  #truth = truth + dunif(x2,matrix(-2,2,1),matrix(2,2,1),log=TRUE)



  truth = truth + dnorm(x,0,1,log=TRUE)


  expect_equal(signif(output, digits = 7), signif(truth, digits = 7))

  # data("gaussian")
  #
  # simulate_prior = function()
  # {
  #   return(list(precision=rgamma(1, 1, 1)))
  # }
  #
  # evaluate_log_likelihood = function(inputs, observed_data)
  # {
  #   return(sum(dnorm(observed_data$data, mean=0, sd=1/sqrt(inputs$precision), log=TRUE)))
  # }
  #
  # number_of_points = 10
  #
  # model = list(simulate_prior = simulate_prior,
  #              evaluate_log_likelihood = evaluate_log_likelihood,
  #              observed_data = list(data=gaussian))
  #
  # algorithm = list(number_of_points=number_of_points)
  #
  # set.seed(1)
  #
  # output_r = importance_sample(model, algorithm)
  #
  # expect_equal(signif(output_r$log_normalising_constant, digits = 7), signif(-141.8712, digits = 7))
})
