# Parse .ilike file to give a compiled ilike "recipe".

Parse .ilike file to give a compiled ilike "recipe".

## Usage

``` r
compile(
  filenames,
  model_parameter_list = list(),
  R_functions = FALSE,
  external_packages = c(),
  julia_bin_dir = "",
  julia_required_libraries = c(),
  verify_cpp_function_types = FALSE,
  keep_temporary_recipe_code = FALSE,
  nesting_level = 1,
  print_block_ends = FALSE
)
```

## Arguments

- filenames:

  The name (and path) of the .ilike files containing the recipe, stored
  in a vector if there is more than one.

- model_parameter_list:

  (optional) A list containing parameters for the recipe.

- R_functions:

  (optional) If TRUE, returns R functions. If FALSE, returns C++
  functions.

- external_packages:

  (optional) A vector of names of other R packages the functions rely
  on.

- julia_bin_dir:

  (optional) The directory containing the Julia bin file - only needed
  if Julia functions are used.

- julia_required_libraries:

  (optional) Vector of strings, each of which is a Julia package that
  will be installed and loaded.

- verify_cpp_function_types:

  (optional) If TRUE, check the types of the parameters of user-defined
  .cpp functions. If FALSE (default), types are not checked.

- keep_temporary_recipe_code:

  (optional) If FALSE (default), the .cpp file generated for compilation
  is deleted. If TRUE, this file is left in the working directory.

- nesting_level:

  (optional) The level of nesting of the current call to compile. A user
  should always use the default of 1.

- print_block_ends:

  (optional) If TRUE, print the end of each block of code. If FALSE
  (default), do not print the end of each block of code.

## Value

A list containing the recipe.
