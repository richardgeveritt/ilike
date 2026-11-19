#' Mean-field dynamical system ensemble Kalman inversion (MFDS)
#'
#' The derivative-free method of Huang, Huang, Reich and Stuart (2022,
#' "Efficient derivative-free Bayesian inference for large-scale inverse
#' problems"). The posterior is the fixed point of a filtering problem: at each
#' iteration the ensemble is inflated about its mean (the prediction) and then
#' updated by an ensemble Kalman step with the augmented observation (the data
#' and the prior mean, observed through the simulated data and the parameters),
#' with the noise of both scaled by 1/delta_t. For linear-Gaussian problems the
#' ensemble's mean and covariance converge to the posterior's.
#'
#' The model needs a Gaussian prior on every parameter, given by built-in prior
#' blocks (prior,ilike::norm or prior,ilike::mvnorm), and the likelihood as a
#' simulator: a simulate_data_model block (which includes the observation
#' noise), a data_variable block and optionally a summary_statistics block. An
#' enk_shifter block chooses the update: ilike::stochastic (perturbed
#' observations, the default) or ilike::adjustment (ensemble adjustment). There
#' is no estimate of the evidence; read the ensembles with load_enk_output (one
#' target per iteration).
#'
#' @param recipe A recipe from compile, or the file name(s) of an ilike file.
#' @param results_name The name of the folder for the results (in results_path).
#' @param results_path (optional) The folder in which results_name is made.
#' @param number_of_ensemble_members The number of ensemble members.
#' @param number_of_iterations The number of iterations.
#' @param delta_t (optional) The step, between 0 and 1 (default 0.5, as in the paper).
#' @param iterations_to_store (optional) The number of iterations kept in memory while running (at least 2).
#' @param model_parameter_list (optional) A list of run-time constants for the recipe.
#' @param external_packages (optional) Packages used by R functions in the recipe.
#' @param julia_bin_dir (optional) The folder of the Julia executable.
#' @param julia_required_libraries (optional) Julia packages used by the recipe.
#' @param verify_cpp_function_types (optional) Check the types of C++ functions in the recipe.
#' @param keep_temporary_recipe_code (optional) Keep the generated C++ code.
#' @param seed (optional) The seed of the random number generator.
#' @param algorithm (optional) For a recipe in the ilike language: the name of the algorithm section to run (needed when the recipe has several).
#' @return Invisibly, the folder holding the results.
#' @export
MFDS = function(recipe,
                results_name,
                results_path = getwd(),
                number_of_ensemble_members,
                number_of_iterations,
                delta_t = 0.5,
                iterations_to_store = 2,
                model_parameter_list = list(),
                external_packages = c(),
                julia_bin_dir = "",
                julia_required_libraries = c(),
                verify_cpp_function_types = FALSE,
                keep_temporary_recipe_code = FALSE,
                seed = NULL,
                algorithm = NULL)
{
  # a new run: its start-of-run notes are printed once each (dev/plan.md 3.11)
  reset_informs()

  if (is.character(recipe))
    recipe = compile(filenames = recipe,
                     model_parameter_list = model_parameter_list,
                     external_packages = external_packages,
                     julia_bin_dir = julia_bin_dir,
                     julia_required_libraries = julia_required_libraries,
                     verify_cpp_function_types = verify_cpp_function_types,
                     keep_temporary_recipe_code = keep_temporary_recipe_code)
  else if (!is.list(recipe))
    stop('"recipe" must be a compiled ilike recipe, or the file name(s) of ilike files.')

  # a recipe in the ilike language: the algorithm to run, and the constants
  # with the file's defaults (dev/ilike-language.md)
  prepared = prepare_recipe(recipe, algorithm, model_parameter_list, transform_allowed = TRUE)
  recipe = prepared$recipe
  model_parameter_list = prepared$model_parameter_list
  apply_algorithm_settings(prepared$settings, names(match.call())[-1], environment(), names(formals(sys.function())))

  if (missing(results_name) || (results_name == ""))
    stop("MFDS: results_name is needed (the ensembles are the output).")
  if (!is.numeric(delta_t) || (length(delta_t) != 1) || (delta_t <= 0) || (delta_t >= 1))
    stop("MFDS: delta_t must be a number between 0 and 1.")
  if (!is.numeric(number_of_iterations) || (number_of_iterations < 1) ||
      (number_of_iterations != round(number_of_iterations)))
    stop("MFDS: number_of_iterations must be a positive whole number.")
  if (!is.numeric(number_of_ensemble_members) || (number_of_ensemble_members < 2) ||
      (number_of_ensemble_members != round(number_of_ensemble_members)))
    stop("MFDS: number_of_ensemble_members must be a whole number of at least 2.")

  results_directory = make_results_directory(results_name, results_path)

  if (is.null(seed))
    seed = ilike_rdtsc()
  set.seed(as.numeric(substr(as.character(seed), 1, 9)))

  do_enkmfds(recipe,
             model_parameter_list,
             number_of_ensemble_members,
             number_of_iterations,
             delta_t,
             max(2, iterations_to_store),
             results_directory,
             seed)

  invisible(file.path(results_path, results_name))
}
