# Importance sampler.

Importance sampler.

## Usage

``` r
IS(
  recipe,
  results_name = "",
  results_path = getwd(),
  write_factor_output = TRUE,
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
  grain_size = 256,
  progress = getOption("ilike.progress", interactive()),
  algorithm = NULL
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

- write_factor_output:

  (optional) If TRUE (default), write each point's factor outputs (the
  iterationK/factorJ\_\* files used to load nested output, e.g. the
  states of a particle filter inside particle marginal MH). FALSE writes
  only the points, weights and evidence.

- number_of_importance_points:

  The number of importance points.

- parallel:

  (optional) If TRUE, the chunks of grain_size points are simulated and
  weighted in parallel threads (the number set by
  RcppParallel::setThreadOptions). The results are exactly those of the
  sequential run with the same seed and grain_size. Needs a model
  written in C++ or with ilike:: built-in functions (not R, Python or
  Julia), a fixed target, and factors that do not estimate their
  likelihood with an inner algorithm; otherwise the run is sequential,
  with a warning.

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

- progress:

  (optional) If TRUE, show a progress bar (cli) while the algorithm
  runs. Default getOption("ilike.progress", interactive()): shown in
  interactive sessions only.

- grain_size:

  (optional) The number of points simulated and weighted together
  (default 256). Each chunk of points has its own random-number stream,
  determined by the seed and the chunk's position, so results depend on
  the seed and on grain_size (and will not depend on whether the run is
  parallel). Larger chunks use more memory when the likelihood is
  estimated by an inner algorithm (each point of a chunk keeps its inner
  particles until the chunk is weighted); a chunk larger than the number
  of points gives the results of a single generation.

- algorithm:

  (optional) For a recipe in the ilike language: the name of the
  algorithm section to run (needed when the recipe has several).

## Value

Estimate of the marginal likelihood.
