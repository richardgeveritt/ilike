# SMC with an MCMC move

SMC with an MCMC move

## Usage

``` r
SMC_with_MCMC(
  recipe,
  results_name = "",
  results_path = getwd(),
  write_factor_output = TRUE,
  number_of_particles,
  parallel = FALSE,
  smc_iterations_to_store = 2,
  write_to_file_at_each_iteration = TRUE,
  model_parameter_list = list(),
  algorithm_parameter_list = list(),
  fixed_parameter_list = list(),
  external_packages = c(),
  julia_bin_dir = "",
  julia_required_libraries = c(),
  verify_cpp_function_types = FALSE,
  keep_temporary_recipe_code = FALSE,
  seed = NULL,
  grain_size = 64,
  check_factor_index = FALSE,
  waste_free_chains = 0,
  progress = getOption("ilike.progress", interactive())
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

- number_of_particles:

  The number of particles.

- parallel:

  (optional) If TRUE, the MCMC moves of each chunk of grain_size
  particles run in parallel threads (the number set by
  RcppParallel::setThreadOptions). The results are exactly those of the
  sequential run with the same seed and grain_size. Needs a model
  written in C++ or with ilike:: built-in functions (not R, Python or
  Julia), exact or annealed factors, and an MCMC termination with a
  fixed number of iterations; otherwise the moves are sequential, with a
  warning. The simulation of the first generation and the choice of each
  next target run in parallel too when the moves can; resampling and
  writing output are sequential.

- smc_iterations_to_store:

  (optional) The number of iterations of SMC output stored in memory as
  the algorithm is running (cannot be fewer than 2).

- write_to_file_at_each_iteration:

  (optional) Do we write the algorithm output to file at each SMC step
  (TRUE/FALSE)?

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

  (optional) The number of particles moved together (default 64). Each
  chunk of particles has its own random-number stream, determined by the
  seed and the chunk's position, so results depend on the seed and on
  grain_size when there are more particles than grain_size (and will not
  depend on whether the run is parallel).

- check_factor_index:

  (optional) If TRUE, after every accepted MCMC move that has a factor
  index, re-evaluate the exact factors that were left out of the move
  (neither in its factor index nor its cut factor index) and stop if any
  changed value, i.e. did not cancel in the acceptance ratio. Slower;
  for checking a model file. Default FALSE.

- waste_free_chains:

  (optional) If greater than 0 (M), run waste-free SMC (Dau and Chopin,
  2022): at every step M particles are resampled and each is run as an
  MCMC chain of number_of_particles/M states (itself and
  number_of_particles/M - 1 moves, whatever mcmc_termination says), and
  all the states form the next generation. Must divide
  number_of_particles. Needs exact or annealed factors, a move that does
  not adapt, and mcmc_termination ilike::iterations. Default 0: standard
  SMC.

## Value

Estimate of the marginal likelihood.
