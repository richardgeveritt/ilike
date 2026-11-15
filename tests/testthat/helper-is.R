test_is_prior_as_proposal_z <- function()
{
  # IS estimation of normalising constant
  filename = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  recipe = ilike::compile(filename)
  output = ilike::IS(recipe,results_name="is_prior_as_proposal",number_of_importance_points=1000,seed=2)
  return(output)
}

test_is_custom_proposal_z <- function()
{
  # testing both the functions in distributions.h and also the evaluate_log function
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_proposal = system.file("extdata", "test_is_custom_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_proposal))
  output = ilike::IS(recipe,results_name="is_custom",number_of_importance_points=1000,seed=2)
  return(output)
}

test_is_builtin_z <- function(fixture, results_name, seed = 1)
{
  # IS with a built-in prior and a built-in proposal that do not cancel
  filename = system.file("extdata", fixture, package = "ilike")
  recipe = ilike::compile(filename)
  ilike::IS(recipe, results_name = results_name, number_of_importance_points = 5000, seed = seed)
}
