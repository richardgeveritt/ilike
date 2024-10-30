#' Kalman filter
#'
#' @param model A file containing the model, or a pre-compiled model list.
#' @param results_name The name of the directory to which results will be written.
#' @param results_path (optional) The path in which the results folder will be created (current working directory is the default).
#' @param kf_iterations_to_store (optional) The number of iterations of filter output stored in memory as the algorithm is running (cannot be fewer than 2).
#' @param write_to_file_at_each_iteration (optional) Do we write the algorithm output to file at each filtering step (TRUE/FALSE)?
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_model_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE,
#' @return Estimate of the marginal likelihood.
#' @export
kf = function(model,
                         results_name,
                         results_path = getwd(),
                         kf_iterations_to_store = 2,
                         write_to_file_at_each_iteration = TRUE,
                         model_parameter_list = list(),
                         algorithm_parameter_list = list(),
                         external_packages = c(),
                         julia_bin_dir="",
                         julia_required_libraries=c(),
                         verify_cpp_function_types=FALSE,
                         keep_temporary_model_code=FALSE)
{
  if ((is.character(model)) && (length(model) == 1))
    model = compile(filename = model,
                    model_parameter_list = model_parameter_list,
                    external_packages = external_packages,
                    julia_bin_dir = julia_bin_dir,
                    julia_required_libraries = julia_required_libraries,
                    verify_cpp_function_types = verify_cpp_function_types,
                    keep_temporary_model_code = keep_temporary_model_code)

  results_directory = make_results_directory(results_name,results_path)

  # Sort filter method.
  filter_method = get_method(model,"filter")

  if (is.null(filter_method))
  {
    stop("No filter method provided: Kalman filter failed.")
  }

  #filter_method = list()

  #browser()

  return(do_kalman_filter(model,
                          model_parameter_list,
                          algorithm_parameter_list,
                          kf_iterations_to_store,
                          write_to_file_at_each_iteration,
                          results_directory))
}
