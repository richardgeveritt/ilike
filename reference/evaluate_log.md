# Evaluate a model specified in a recipe.

Evaluate a model specified in a recipe.

## Usage

``` r
evaluate_log(
  recipe,
  results_name = "",
  results_path = getwd(),
  parameters,
  index = matrix(0, 0, 0),
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

- parameters:

  A list containing the parameters at which to evaluate the model.

- index:

  (optional) The index of the factors to evaluate, where the indexing
  starts at 1 (default is to evaluate all factors).

- model_parameter_list:

  (optional) A list containing parameters for the model.

- algorithm_parameter_list:

  (optional) A list containing named parameters for the algorithm.

- fixed_parameter_list:

  (optional) A list containing parameters to condition on.

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

- seed:

  (optional) The seed for the random number generator.

- grain_size:

  (optional) Sets a minimum chunk size for parallelisation (see
  https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).

## Value

Estimate of the marginal likelihood.
