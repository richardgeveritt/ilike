#' MCMC
#'
#' @param model A file containing the model, or a pre-compiled model list.
#' @param results_directory The name of the directory to which results will be written.
#' @param number_of_mcmc_iterations The number of importance points.
#' @param number_of_chains (optional) The number of chains.
#' @param initial_values (optional) A list of lists containing the initial values for the chains.
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param seed (optional) The seed for the random number generator.
#' @param parallel_flag (optional) Set to true to perform the importance sampling in parallel, false for serial.
#' @param grain_size (optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Nothing: output can be found in the output_directory.
#' @export
mcmc = function(model,
                results_directory,
                number_of_mcmc_iterations,
                number_of_chains=1,
                initial_values = list(),
                model_parameter_list = list(),
                algorithm_parameter_list = list(),
                seed = NULL,
                parallel_flag = FALSE,
                grain_size = 100000)
{
  if ((is.character(model)) && (length(model) == 1))
    model = parse_ilike_model(model,model_parameter_list)

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  do_mcmc(model,
          model_parameter_list,
          algorithm_parameter_list,
          initial_values,
          number_of_mcmc_iterations,
          number_of_chains,
          parallel_flag,
          grain_size,
          results_directory,
          seed)
}
