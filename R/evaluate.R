#' Evaluate a model specified in a recipe.
#'
#' @param recipe A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.
#' @param results_name (optional) The name of the directory to which results will be written (default is to not write output to file).
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param parameters A list containing the parameters at which to evaluate the model.
#' @param index (optional) The index of the factors to evaluate, where the indexing starts at 1 (default is to evaluate all factors).
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
#' @param algorithm (optional) For a recipe in the ilike language: the name of the algorithm section to run (needed when the recipe has several).
#' @return Estimate of the marginal likelihood.
#' @export
evaluate_log = function(recipe,
                         results_name = "",
                         results_path = getwd(),
                    parameters,
                    index = matrix(0,0,0),
                         model_parameter_list = list(),
                         algorithm_parameter_list = list(),
                    fixed_parameter_list = list(),
                         external_packages = c(),
                         julia_bin_dir="",
                         julia_required_libraries=c(),
                         verify_cpp_function_types=FALSE,
                         keep_temporary_recipe_code=FALSE,
                         seed = NULL,
                         grain_size = 100000,
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
    stop('"Receipe" argument must be either a compiled ilike recipe, the filename of an ilike file, or a vector of filenames of ilike files.')

  # a recipe in the ilike language: the algorithm to run, and the constants
  # with the file's defaults (dev/ilike-language.md)
  prepared = prepare_recipe(recipe, algorithm, model_parameter_list)
  recipe = prepared$recipe
  model_parameter_list = prepared$model_parameter_list
  apply_algorithm_settings(prepared$settings, names(match.call())[-1], environment(), names(formals(sys.function())))

  if (!results_name == "")
    results_directory = make_results_directory(results_name,results_path)
  else
    results_directory = ""

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  set.seed(as.numeric(substr(as.character(seed),1,9)))

  if (length(fixed_parameter_list)==0)
  {
    return(do_evaluate_log(recipe,
                            model_parameter_list,
                            algorithm_parameter_list,
                            parameters,
                            index,
                            results_directory,
                            seed))
  }
  else
  {
    return(do_evaluate_log_with_fixed_params(recipe,
                                 model_parameter_list,
                                 algorithm_parameter_list,
                                 fixed_parameter_list,
                                 parameters,
                                 index,
                                 results_directory,
                                 seed))
  }
}

# The gradient of the log target with respect to one variable at a point (the
# sum of the factors' gradients; index chooses factors, numbered from 1), for
# checking gradients against finite differences (dev/plan.md I4e). Internal.
evaluate_gradient_log = function(recipe, parameters, variable, index = integer(0), model_parameter_list = list(),
                                 algorithm = NULL)
{
  if (is.character(recipe))
    recipe = compile(filenames = recipe, model_parameter_list = model_parameter_list)
  prepared = prepare_recipe(recipe, algorithm, model_parameter_list)
  do_evaluate_gradient_log(prepared$recipe, prepared$model_parameter_list, parameters, variable, as.integer(index))
}

# The gradient of a log density f(parameters, data) with respect to one
# variable by central differences: each element moved by
# h = eps^(1/3) max(1, |x|), the step that balances the truncation error,
# O(h^2), against the rounding error, O(eps / h): the same rule as the C++
# that gradient = autodiff generates (dev/plan.md 3.10). Used by
# check_gradient. Internal.
finite_difference_gradient = function(f, variable, parameters, data = NULL)
{
  x = parameters[[variable]]
  if (is.null(x))
    stop(paste("finite differences: the variable ", variable, " is not in the parameters.", sep = ""), call. = FALSE)
  gradient = x
  storage.mode(gradient) = "double"
  scale = .Machine$double.eps^(1/3)
  moved = parameters
  for (k in seq_along(x))
  {
    h = scale * max(1, abs(x[k]))
    up = x[k] + h
    down = x[k] - h
    moved[[variable]][k] = up
    f_up = f(moved, data)
    moved[[variable]][k] = down
    f_down = f(moved, data)
    moved[[variable]][k] = x[k]
    gradient[k] = (f_up - f_down) / (up - down)
  }
  gradient
}

#' Check gradients against finite differences.
#'
#' Compares the gradient of the log target (or of some of its factors) that a
#' recipe provides, from functions written by the user, from model statements,
#' or from \code{gradient = autodiff}, with central differences of the log
#' target itself, element by element.
#'
#' @param recipe A compiled ilike recipe, or the name of an ilike file.
#' @param parameters A named list: the point at which to check.
#' @param variables (optional) The variables whose gradients are checked (default: all those in \code{parameters}).
#' @param index (optional) The factors to check, numbered from 1 (default: all factors).
#' @param model_parameter_list (optional) Values for the constants of the model.
#' @param algorithm (optional) For a recipe with several algorithm sections: the name of one.
#' @return A data frame with one row per element: \code{variable}, \code{element}, \code{gradient} (the recipe's),
#' \code{finite_difference}, and \code{error}, the absolute difference divided by \code{max(1, |finite_difference|)}.
#' Errors around 1e-6 or smaller are what finite differences can resolve.
#' @export
check_gradient = function(recipe, parameters, variables = names(parameters), index = integer(0),
                          model_parameter_list = list(), algorithm = NULL)
{
  if (is.character(recipe))
    recipe = compile(filenames = recipe, model_parameter_list = model_parameter_list)
  evaluation_index = if (length(index) > 0) as.integer(index) else matrix(0, 0, 0)
  log_target = function(p) evaluate_log(recipe, parameters = p, index = evaluation_index, seed = 1,
                                        model_parameter_list = model_parameter_list, algorithm = algorithm)
  rows = list()
  for (variable in variables)
  {
    gradient = as.vector(evaluate_gradient_log(recipe, parameters, variable, index, model_parameter_list, algorithm))
    finite_difference = as.vector(finite_difference_gradient(function(p, data) log_target(p), variable, parameters))
    rows[[length(rows) + 1]] = data.frame(variable = variable, element = seq_along(finite_difference),
                                          gradient = gradient, finite_difference = finite_difference,
                                          error = abs(gradient - finite_difference) / pmax(1, abs(finite_difference)))
  }
  do.call(rbind, rows)
}
