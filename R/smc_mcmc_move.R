#' SMC with an MCMC move
#'
#' @param recipe A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.
#' @param results_name The name of the directory to which results will be written.
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param number_of_particles The number of particles.
#' @param parallel (optional) Set to true to perform the importance sampling in parallel, false for serial.
#' @param smc_iterations_to_store (optional) The number of iterations of SMC output stored in memory as the algorithm is running (cannot be fewer than 2).
#' @param write_to_file_at_each_iteration (optional) Do we write the algorithm output to file at each SMC step (TRUE/FALSE)?
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_recipe_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE,
#' @param seed (optional) The seed for the random number generator.
#' @param grain_size (optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Estimate of the marginal likelihood.
#' @export
SMC_with_MCMC = function(recipe,
                         results_name,
                         results_path = getwd(),
                         number_of_particles,
                         parallel = FALSE,
                         smc_iterations_to_store = 2,
                         write_to_file_at_each_iteration = TRUE,
                         model_parameter_list = list(),
                         algorithm_parameter_list = list(),
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

  # Sort MCMC termination method.
  mcmc_termination_method = get_method(recipe,"mcmc_termination")

  if (is.null(mcmc_termination_method))
  {
    print("No MCMC termination method provided: using default of 1 iteration.")
    mcmc_termination_method = list(method="iterations",values=list(as.character(1)))
  }

  # Adaptive resampling method.
  # adaptive_resampling_method = get_method(recipe,"adaptive_resampling")
  #
  # if (is.null(adaptive_resampling_method))
  # {
  #   print("No method set for adaptive resampling: defaulting to resampling whenever the ESS drops below the number of particles.")
  #   adaptive_resampling_method = list(method='ess',values=as.list(paste(number_of_particles)))
  # }

  # Adaptive target method.
  adaptive_target_method = get_method(recipe,"adaptive_target")

  if (is.null(adaptive_target_method))
  {
    print("No method set for adaptive sequence of target distributions in SMC: defaulting to moving directly from one distribution to another.")
    adaptive_target_method = list(method='positive')
  }

  # SMC sequencing.
  smc_sequencer_method = get_smc_sequences(recipe,model_parameter_list)

  if (is.null(smc_sequencer_method))
  {
    print("No SMC sequence set. Defaulting to importance sampling with prior as proposal.")
    smc_sequencer_method = list(types=c("annealing"),variables=c("power"),schedules=c(0,1))
  }

  # SMC termination method.
  smc_termination_method = get_method(recipe,"smc_termination")

  if (is.null(smc_termination_method))
  {
    smc_termination_method = list()
  }

  # MCMC weights method.
  mcmc_weights_method = get_method(recipe,"mcmc_weights")

  if (is.null(mcmc_weights_method))
  {
    mcmc_weights_method = list()
  }

  return(do_smc_mcmc_move(recipe,
                          model_parameter_list,
                          algorithm_parameter_list,
                          number_of_particles,
                          mcmc_termination_method,
                          mcmc_weights_method,
                          smc_sequencer_method,
                          adaptive_target_method,
                          smc_termination_method,
                          smc_iterations_to_store,
                          write_to_file_at_each_iteration,
                          parallel,
                          grain_size,
                          results_directory,
                          seed))
}
