evaluate <- function()
{
  # testing both the functions in distributions.h and also the evaluate_log function
  filename = system.file("extdata", "test_evaluate_priors.ilike", package = "ilike")
  recipe = ilike::compile(filename)
  output = ilike::evaluate_log(recipe,parameters = list(x=1))
  return(output)
}
