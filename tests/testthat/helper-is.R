test_is_prior_as_proposal_z <- function()
{
  # IS estimation of normalising constant
  message("[diag] IS prior: before compile")
  filename = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  recipe = ilike::compile(filename)
  message("[diag] IS prior: after compile, before IS()")
  output = ilike::IS(recipe,results_name="is_prior_as_proposal",number_of_importance_points=1000,seed=2)
  message("[diag] IS prior: after IS(), returning")
  return(output)
}

test_is_custom_proposal_z <- function()
{
  # testing both the functions in distributions.h and also the evaluate_log function
  message("[diag] IS custom: before compile")
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_proposal = system.file("extdata", "test_is_custom_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_proposal))
  message("[diag] IS custom: after compile, before IS()")
  output = ilike::IS(recipe,results_name="is_custom",number_of_importance_points=1000,seed=2)
  message("[diag] IS custom: after IS(), returning")
  return(output)
}
