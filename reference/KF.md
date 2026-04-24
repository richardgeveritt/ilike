# Kalman filter

Kalman filter

## Usage

``` r
KF(
  recipe,
  results_name = "",
  results_path = getwd(),
  kf_iterations_to_store = 2,
  write_to_file_at_each_iteration = TRUE,
  model_parameter_list = list(),
  algorithm_parameter_list = list(),
  external_packages = c(),
  julia_bin_dir = "",
  julia_required_libraries = c(),
  verify_cpp_function_types = FALSE,
  keep_temporary_recipe_code = FALSE
)
```

## Arguments

- recipe:

  A pre-compiled ilike recipe, or an ilike file, or a vector of ilike
  files.

- results_name:

  (optional) The name of the directory to which results will be written
  (default is to not write output to file).

- results_path:

  (optional) The path in which the results folder will be created
  (current working directory is the default).

- kf_iterations_to_store:

  (optional) The number of iterations of filter output stored in memory
  as the algorithm is running (cannot be fewer than 2).

- write_to_file_at_each_iteration:

  (optional) Do we write the algorithm output to file at each filtering
  step (TRUE/FALSE)?

- model_parameter_list:

  (optional) A list containing parameters for the model.

- algorithm_parameter_list:

  (optional) A list containing named parameters for the algorithm.

- external_packages:

  (optional) A vector of names of other R packages the functions rely
  on.

- julia_bin_dir:

  (optional) The directory containing the Julia bin file - only needed
  if Julia functions are used.

- julia_required_libraries:

  (optional) Vector of strings, each of which is a Julia packge that
  will be installed and loaded.

- verify_cpp_function_types:

  (optional) If TRUE, check the types of the parameters of user-defined
  .cpp functions. If FALSE (default), types are not checked.

- keep_temporary_recipe_code:

  (optional) If FALSE (default), the .cpp file generated for compilation
  is deleted. If TRUE,

## Value

Estimate of the marginal likelihood.
