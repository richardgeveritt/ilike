#' Ensemble Kalman filter
#'
#' @param recipe A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.
#' @param results_name The name of the directory to which results will be written.
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param number_of_ensemble_members The number of ensemble members.
#' @param parallel (optional) Set to true to perform the importance sampling in parallel, false for serial.
#' @param enkf_iterations_to_store (optional) The number of iterations of filter output stored in memory as the algorithm is running (cannot be fewer than 2).
#' @param write_to_file_at_each_iteration (optional) Do we write the algorithm output to file at each filtering step (TRUE/FALSE)?
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param fixed_parameter_list (optional) A list containing parameters to condition on.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_recipe_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE,
#' @param seed (optional) The seed for the random number generator.
#' @param grain_size (optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Estimate of the marginal likelihood.
#' @export
ENKF = function(recipe,
                                  results_name,
                                  results_path = getwd(),
                                  number_of_ensemble_members,
                                  parallel = FALSE,
                                  enkf_iterations_to_store = 2,
                                  write_to_file_at_each_iteration = TRUE,
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

  results_directory = make_results_directory(results_name,results_path)

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  set.seed(as.numeric(substr(as.character(seed),1,9)))

  # Sort filter method.
  filter_method = get_method(recipe,"filter")

  if (is.null(filter_method))
  {
    stop("No filter method provided: Kalman filter failed.")
  }

  # # The indices of the factors used as likelihoods in the EnK method.
  # enk_likelihood_index_method = get_method(recipe,"enk_likelihood_index")
  #
  # if (is.null(enk_likelihood_index_method))
  # {
  #   enk_likelihood_index_method = list()
  # }
  #
  # # EnK shifter method.
  # enk_shifter_method = get_method(recipe,"enki_shifter")
  #
  # # Method for shifting the ensemble.
  # if (is.null(enk_shifter_method))
  # {
  #   print("No method set for shifting ensemble in EnK: defaulting to stochastic approach.")
  #   enk_shifter_method = list(method='stochastic')
  # }

  if (length(fixed_parameter_list)==0)
  {
    return(do_ensemble_kalman_filter(recipe,
                                     model_parameter_list,
                                     algorithm_parameter_list,
                                     number_of_ensemble_members,
                                     enkf_iterations_to_store,
                                     write_to_file_at_each_iteration,
                                     parallel,
                                     grain_size,
                                     results_directory,
                                     seed))
  }
  else
  {
    return(do_ensemble_kalman_filter_with_fixed_params(recipe,
                                     model_parameter_list,
                                     algorithm_parameter_list,
                                     fixed_parameter_list,
                                     number_of_ensemble_members,
                                     enkf_iterations_to_store,
                                     write_to_file_at_each_iteration,
                                     parallel,
                                     grain_size,
                                     results_directory,
                                     seed))
  }
}
