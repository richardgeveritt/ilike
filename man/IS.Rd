% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importance_sampler.R
\name{IS}
\alias{IS}
\title{Importance sampler.}
\usage{
IS(
  recipe,
  results_name = "",
  results_path = getwd(),
  number_of_importance_points,
  parallel = FALSE,
  model_parameter_list = list(),
  algorithm_parameter_list = list(),
  fixed_parameter_list = list(),
  external_packages = c(),
  julia_bin_dir = "",
  julia_required_libraries = c(),
  verify_cpp_function_types = FALSE,
  keep_temporary_recipe_code = FALSE,
  seed = NULL,
  grain_size = 1e+05
)
}
\arguments{
\item{recipe}{A pre-compiled ilike recipe, or an ilike file, or a vector of ilike files.}

\item{results_name}{(optional) The name of the directory to which results will be written (default is to not write output to file).}

\item{results_path}{(optional) The path in which the results folder will be created (current working directory is the default).}

\item{number_of_importance_points}{The number of importance points.}

\item{parallel}{(optional) Set to true to perform the importance sampling in parallel, false for serial.}

\item{model_parameter_list}{(optional) A list containing parameters for the model.}

\item{algorithm_parameter_list}{(optional) A list containing named parameters for the algorithm.}

\item{fixed_parameter_list}{(optional) A list containing parameters to condition on.}

\item{external_packages}{(optional) A vector of names of other R packages the functions rely on.}

\item{julia_bin_dir}{(optional) The directory containing the Julia bin file - only needed if Julia functions are used.}

\item{julia_required_libraries}{(optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.}

\item{verify_cpp_function_types}{(optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.}

\item{keep_temporary_recipe_code}{(optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE,}

\item{seed}{(optional) The seed for the random number generator.}

\item{grain_size}{(optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).}
}
\value{
Estimate of the marginal likelihood.
}
\description{
Importance sampler.
}
