#' Importance sampler.
#'
#' @param recipe A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.
#' @param results_name (optional) The name of the directory to which results will be written (default is to not write output to file).
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param write_factor_output (optional) If TRUE (default), write each point's factor outputs (the iterationK/factorJ_* files used to load nested output, e.g. the states of a particle filter inside particle marginal MH). FALSE writes only the points, weights and evidence.
#' @param number_of_importance_points The number of importance points.
#' @param parallel (optional) If TRUE, the chunks of grain_size points are simulated and weighted in parallel threads (the number set by RcppParallel::setThreadOptions). The results are exactly those of the sequential run with the same seed and grain_size. Needs a model written in C++ or with ilike:: built-in functions (not R, Python or Julia), a fixed target, and factors that do not estimate their likelihood with an inner algorithm; otherwise the run is sequential, with a warning.
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param fixed_parameter_list (optional) A list containing parameters to condition on.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_recipe_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE,
#' @param seed (optional) The seed for the random number generator.
#' @param progress (optional) If TRUE, show a progress bar (cli) while the algorithm runs. Default getOption("ilike.progress", interactive()): shown in interactive sessions only.
#' @param grain_size (optional) The number of points simulated and weighted together (default 256). Each chunk of points has its own random-number stream, determined by the seed and the chunk's position, so results depend on the seed and on grain_size (and will not depend on whether the run is parallel). Larger chunks use more memory when the likelihood is estimated by an inner algorithm (each point of a chunk keeps its inner particles until the chunk is weighted); a chunk larger than the number of points gives the results of a single generation.
#' @param algorithm (optional) For a recipe in the ilike language: the name of the algorithm section to run (needed when the recipe has several).
#' @return Estimate of the marginal likelihood.
#' @export
IS = function(recipe,
                             results_name = "",
                             results_path = getwd(),
                             write_factor_output = TRUE,
                             number_of_importance_points,
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
                             grain_size = 256,
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

  if (!results_name == "")
    results_directory = make_results_directory(results_name,results_path)
  else
    results_directory = ""

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  set.seed(r_seed(seed))

  if (length(fixed_parameter_list)==0)
  {
    return(do_importance_sampler(recipe,
                                 model_parameter_list,
                                 algorithm_parameter_list,
                                 number_of_importance_points,
                                 parallel,
                                 grain_size,
                                 results_directory,
                                 seed))
  }
  else
  {
    return(do_importance_sampler_with_fixed_params(recipe,
                                 model_parameter_list,
                                 algorithm_parameter_list,
                                 fixed_parameter_list,
                                 number_of_importance_points,
                                 parallel,
                                 grain_size,
                                 results_directory,
                                 seed))
  }

}
