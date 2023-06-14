#' MCMC
#'
#' @param model A file containing the model, or a pre-compiled model list.
#' @param results_directory The name of the directory to which results will be written.
#' @param number_of_mcmc_iterations The number of MCMC iterations.
#' @param mcmc_termination_method (optional) The method used to terminate the MCMC runs.
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
                number_of_mcmc_iterations = 1,
                mcmc_termination_method = NULL,
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

  # Sort MCMC termination method.
  if (!is.null(mcmc_termination_method))
  {
    print("MCMC termination method provided: using this method instead of number_of_mcmc_iterations.")

    if ("values" %in% names(mcmc_termination_method))
    {
      mcmc_termination_method[["values"]] = as.list(paste(mcmc_termination_method[["values"]]))
    }
  }
  else
  {
    mcmc_termination_method = list(method="iterations",values=list(as.character(number_of_mcmc_iterations)))
  }

  if (length(initial_values)==0)
  {
    print("No initial values set: initial points drawn from the prior.")
  }
  else
  {
    if (length(initial_values)!=number_of_chains)
    {
      print("Initial values provided: number of chains determined by the number of initial values.")
      number_of_chains = length(initial_values)
    }
  }

  do_mcmc(model,
          model_parameter_list,
          algorithm_parameter_list,
          initial_values,
          mcmc_termination_method,
          number_of_chains,
          parallel_flag,
          grain_size,
          results_directory,
          seed)
}
