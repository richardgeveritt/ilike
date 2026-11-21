#' SMC with an MCMC move
#'
#' @param recipe A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.
#' @param results_name (optional) The name of the directory to which results will be written (default is to not write output to file).
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param write_factor_output (optional) If TRUE (default), write each point's factor outputs (the iterationK/factorJ_* files used to load nested output, e.g. the states of a particle filter inside particle marginal MH). FALSE writes only the points, weights and evidence.
#' @param number_of_particles The number of particles.
#' @param parallel (optional) If TRUE, the MCMC moves of each chunk of grain_size particles run in parallel threads (the number set by RcppParallel::setThreadOptions). The results are exactly those of the sequential run with the same seed and grain_size. Needs a model written in C++ or with ilike:: built-in functions (not R, Python or Julia), exact or annealed factors, and an MCMC termination with a fixed number of iterations; otherwise the moves are sequential, with a warning. The simulation of the first generation and the choice of each next target run in parallel too when the moves can; resampling and writing output are sequential.
#' @param smc_iterations_to_store (optional) The number of iterations of SMC output stored in memory as the algorithm is running (cannot be fewer than 2).
#' @param write_to_file_at_each_iteration (optional) Do we write the algorithm output to file at each SMC step (TRUE/FALSE)?
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param fixed_parameter_list (optional) A list containing parameters to condition on.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_recipe_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE,
#' @param check_factor_index (optional) If TRUE, after every accepted MCMC move that has a factor index, re-evaluate the exact factors that were left out of the move (neither in its factor index nor its cut factor index) and stop if any changed value, i.e. did not cancel in the acceptance ratio. Slower; for checking a model file. Default FALSE.
#' @param waste_free_chains (optional) If greater than 0 (M), run waste-free SMC (Dau and Chopin, 2022): at every step M particles are resampled and each is run as an MCMC chain of number_of_particles/M states (itself and number_of_particles/M - 1 moves, whatever mcmc_termination says), and all the states form the next generation. Must divide number_of_particles. Needs exact or annealed factors, a move that does not adapt, and mcmc_termination ilike::iterations. Default 0: standard SMC.
#' @param seed (optional) The seed for the random number generator.
#' @param progress (optional) If TRUE, show a progress bar (cli) while the algorithm runs. Default getOption("ilike.progress", interactive()): shown in interactive sessions only.
#' @param grain_size (optional) The number of particles moved together (default 64). Each chunk of particles has its own random-number stream, determined by the seed and the chunk's position, so results depend on the seed and on grain_size when there are more particles than grain_size (and will not depend on whether the run is parallel).
#' @param algorithm (optional) For a recipe in the ilike language: the name of the algorithm section to run (needed when the recipe has several).
#' @return Estimate of the marginal likelihood.
#' @export
SMC_with_MCMC = function(recipe,
                         results_name = "",
                         results_path = getwd(),
                         write_factor_output = TRUE,
                         number_of_particles,
                         parallel = FALSE,
                         smc_iterations_to_store = 2,
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
                         grain_size = 64,
                         check_factor_index = FALSE,
                         waste_free_chains = 0,
                         progress = getOption("ilike.progress", interactive()),
                         algorithm = NULL)
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

  # a recipe in the ilike language: the algorithm to run, and the constants
  # with the file's defaults (dev/ilike-language.md)
  prepared = prepare_recipe(recipe, algorithm, model_parameter_list)
  recipe = prepared$recipe
  model_parameter_list = prepared$model_parameter_list
  apply_algorithm_settings(prepared$settings, names(match.call())[-1], environment(), names(formals(sys.function())))

  # parallel = TRUE: refuse models that call R; pass the number of threads
  recipe = prepare_parallel(recipe, parallel)
  # progress bar for this run (read in C++, apply_output_options)
  recipe$show_progress = isTRUE(progress)

  # Read in C++ (apply_output_options).
  recipe$write_factor_output = write_factor_output
  recipe$waste_free_chains = waste_free_chains

  if (!results_name == "")
    results_directory = make_results_directory(results_name,results_path)
  else
    results_directory = ""

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  set.seed(r_seed(seed))

  # Sort MCMC termination method.
  mcmc_termination_method = get_method(recipe,"mcmc_termination")

  if (is.null(mcmc_termination_method))
  {
    message("No MCMC termination method provided: using default of 1 iteration.")
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
    message("No method set for adaptive sequence of target distributions in SMC: defaulting to moving directly from one distribution to another.")
    adaptive_target_method = list(method='positive')
  }

  # SMC sequencing.
  smc_sequencer_method = get_smc_sequences(recipe,model_parameter_list)

  if (is.null(smc_sequencer_method))
  {
    message("No SMC sequence set. Defaulting to importance sampling with prior as proposal.")
    smc_sequencer_method = default_annealing_sequence()
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

  # Read by make_mcmc in C++.
  recipe$check_factor_index = check_factor_index

  if (length(fixed_parameter_list)==0)
  {
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
  else
  {
    return(do_smc_mcmc_move_with_fixed_params(recipe,
                            model_parameter_list,
                            algorithm_parameter_list,
                            fixed_parameter_list,
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
}
