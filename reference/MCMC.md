# MCMC

MCMC

## Usage

``` r
MCMC(
  recipe,
  results_name = "",
  results_path = getwd(),
  number_of_chains = 1,
  initial_values = list(),
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

- number_of_chains:

  (optional) The number of chains.

- initial_values:

  (optional) A list of lists containing the initial values for the
  chains.

- parallel:

  (optional) Set to true to perform the MCMC chains in parallel, false
  for serial.

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
  is deleted. If TRUE, this file is left in the working directory.

- seed:

  (optional) The seed for the random number generator.

- grain_size:

  (optional) Sets a minimum chunk size for parallelisation (see
  https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).

## Value

Nothing: output can be found in the output_directory.
