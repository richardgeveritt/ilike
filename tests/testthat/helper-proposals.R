run_mcmc_fixture <- function(model, algorithm, results_name, initial_values,
                             model_parameter_list = list(), seed = 2, ...)
{
  # MCMC with the model and algorithm in two files from inst/extdata; further
  # arguments (algorithm_parameter_list, fixed_parameter_list) go to MCMC()
  filenames = system.file("extdata", c(model, algorithm), package = "ilike")
  recipe = ilike::compile(filenames, model_parameter_list = model_parameter_list)
  ilike::MCMC(recipe, results_name = results_name, seed = seed,
              initial_values = list(initial_values),
              model_parameter_list = model_parameter_list, ...)
  ilike::load_mcmc_output(results_name)
}
