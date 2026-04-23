test_is_prior_as_proposal_z <- function()
{
  # IS estimation of normalising constant
  filename = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  recipe = ilike::compile(filename)
  output = ilike::IS(recipe,results_name="is_prior_as_proposal",number_of_importance_points=100000,seed=2)
  return(output)
}

test_is_custom_proposal_z <- function()
{
  # testing both the functions in distributions.h and also the evaluate_log function
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_proposal = system.file("extdata", "test_is_custom_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_proposal))
  output = ilike::IS(recipe,results_name="is_custom",number_of_importance_points=100000,seed=2)
  return(output)
}
