#' Importance sampler.
#'
#' @param model_filename A file containing the model.
#' @param output_directory The name of the directory to which results will be written..
#' @param number_of_importance_points The number of importance points.
#' @param parameter_list A list containing parameters for the model.
#' @param seed The seed for the random number generator.
#' @param parallel_flag Set to true to perform the importance sampling in parallel, false for serial.
#' @param grain_size Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Nothing: output can be found in the output_directory.
#' @export
importance_sample = function(model_filename,
                             results_directory,
                             number_of_importance_points,
                             parameter_list=list(),
                             seed = NULL,
                             parallel_flag = FALSE,
                             grain_size = 1)
{
  model = parse_like_model(filename)

  if (is.null(seed))
  {
    do_importance_sampler(model,
                        parameter_list,
                        number_of_importance_points,
                        parallel_flag,
                        grain_size,
                        results_name)
  }
  else
  {
    do_importance_sampler(model,
                          parameter_list,
                          number_of_importance_points,
                          parallel_flag,
                          grain_size,
                          results_name,
                          seed)
  }
}
