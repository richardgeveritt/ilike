#' MCMC
#'
#' @param recipe A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.
#' @param results_name (optional) The name of the directory to which results will be written (default is to not write output to file).
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param number_of_chains (optional) The number of chains.
#' @param initial_values (optional) A list of lists containing the initial values for the chains.
#' @param parallel (optional) Set to true to perform the MCMC chains in parallel, false for serial.
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param fixed_parameter_list (optional) A list containing parameters to condition on.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_recipe_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE, this file is left in the working directory.
#' @param seed (optional) The seed for the random number generator.
#' @param grain_size (optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Nothing: output can be found in the output_directory.
#' @export
MCMC = function(recipe,
                results_name = "",
                results_path = getwd(),
                number_of_chains=1,
                initial_values = list(),
                parallel = FALSE,
                model_parameter_list = list(),
                algorithm_parameter_list = list(),
                fixed_parameter_list = list(),
                external_packages = c(),
                julia_bin_dir="",
                julia_required_libraries=c(),
                verify_cpp_function_types=FALSE,
                keep_temporary_recipe_code=FALSE,
                seed = NULL,
                grain_size = 100000)
{

  if (is.character(recipe))
    recipe = compile(filenames = recipe,
                     model_parameter_list = model_parameter_list,
                     external_packages = external_packages,
                     julia_bin_dir = julia_bin_dir,
                     julia_required_libraries = julia_required_libraries,
                     verify_cpp_function_types = verify_cpp_function_types,
                     keep_temporary_recipe_code = keep_temporary_recipe_code)
  else if (!is.list(recipe))
    stop('"Receipe" argument must be either a compiled ilike recipe, the filename of an ilike file, or a vector of filenames of ilike files.')

  if (!results_name == "")
    results_directory = make_results_directory(results_name,results_path)
  else
    results_directory = ""

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  set.seed(as.numeric(substr(as.character(seed),1,9)))

  # Sort MCMC termination method.
  mcmc_termination_method = get_method(recipe,"mcmc_termination")

  if (is.null(mcmc_termination_method))
  {
    print("No MCMC termination method provided: using default of 10^5 iterations.")
    mcmc_termination_method = list(method="iterations",values=list(as.character(100000)))
  }

  if (length(initial_values)==0)
  {
    stop("Initial values set using a list. Future versions will allow initial points drawn from the prior.")
  }
  else
  {
    if (length(initial_values)!=number_of_chains)
    {
      print("Initial values provided: number of chains determined by the number of initial values.")
      number_of_chains = length(initial_values)
    }
  }

  # MCMC weights method.
  mcmc_weights_method = get_method(recipe,"mcmc_weights")

  if (is.null(mcmc_weights_method))
  {
    mcmc_weights_method = list()
  }

  if (length(fixed_parameter_list)==0)
  {
    do_mcmc(recipe,
            model_parameter_list,
            algorithm_parameter_list,
            initial_values,
            mcmc_termination_method,
            mcmc_weights_method,
            number_of_chains,
            parallel,
            grain_size,
            results_directory,
            seed)
  }
  else
  {
    do_mcmc_with_fixed_params(recipe,
                              model_parameter_list,
                              algorithm_parameter_list,
                              fixed_parameter_list,
                              initial_values,
                              mcmc_termination_method,
                              mcmc_weights_method,
                              number_of_chains,
                              parallel,
                              grain_size,
                              results_directory,
                              seed)
  }
}
