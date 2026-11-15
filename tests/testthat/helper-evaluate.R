evaluate <- function()
{
  # testing both the functions in distributions.h and also the evaluate_log function
  filename = system.file("extdata", "test_evaluate_priors.ilike", package = "ilike")
  recipe = ilike::compile(filename)
  output = ilike::evaluate_log(recipe,parameters = list(x=1))
  return(output)
}

evaluate_cut_model_factors <- function(index)
{
  # evaluate_log on the four-factor cut model, at a point where every factor
  # has a different value, for a chosen set of (1-based) factor numbers
  filename = system.file("extdata", "test_cut_model.ilike", package = "ilike")
  recipe = ilike::compile(filename)
  ilike::evaluate_log(recipe, parameters = list(phi = 0.3, theta = 0.5), index = index)
}
