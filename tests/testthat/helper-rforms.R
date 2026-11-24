rforms_files <- function(...)
{
  system.file("extdata", c(...), package = "ilike")
}

# MCMC on the R-expression model with theta written as a non-ASCII name
run_rforms_mcmc <- function(algorithm, results_name, seed = 1)
{
  recipe = ilike::compile(rforms_files("test_rforms_gaussian_posterior.ilike", algorithm))
  initial = list(1)
  names(initial) = "theta"
  ilike::MCMC(recipe, results_name = results_name, initial_values = list(initial), seed = seed)
  ilike::load_mcmc_output(results_name)$Value
}
