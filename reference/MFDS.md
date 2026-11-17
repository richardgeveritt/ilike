# Mean-field dynamical system ensemble Kalman inversion (MFDS)

The derivative-free method of Huang, Huang, Reich and Stuart (2022,
"Efficient derivative-free Bayesian inference for large-scale inverse
problems"). The posterior is the fixed point of a filtering problem: at
each iteration the ensemble is inflated about its mean (the prediction)
and then updated by an ensemble Kalman step with the augmented
observation (the data and the prior mean, observed through the simulated
data and the parameters), with the noise of both scaled by 1/delta_t.
For linear-Gaussian problems the ensemble's mean and covariance converge
to the posterior's.

## Usage

``` r
MFDS(
  recipe,
  results_name,
  results_path = getwd(),
  number_of_ensemble_members,
  number_of_iterations,
  delta_t = 0.5,
  iterations_to_store = 2,
  model_parameter_list = list(),
  external_packages = c(),
  julia_bin_dir = "",
  julia_required_libraries = c(),
  verify_cpp_function_types = FALSE,
  keep_temporary_recipe_code = FALSE,
  seed = NULL,
  algorithm = NULL
)
```

## Arguments

- recipe:

  A recipe from compile, or the file name(s) of an ilike file.

- results_name:

  The name of the folder for the results (in results_path).

- results_path:

  (optional) The folder in which results_name is made.

- number_of_ensemble_members:

  The number of ensemble members.

- number_of_iterations:

  The number of iterations.

- delta_t:

  (optional) The step, between 0 and 1 (default 0.5, as in the paper).

- iterations_to_store:

  (optional) The number of iterations kept in memory while running (at
  least 2).

- model_parameter_list:

  (optional) A list of run-time constants for the recipe.

- external_packages:

  (optional) Packages used by R functions in the recipe.

- julia_bin_dir:

  (optional) The folder of the Julia executable.

- julia_required_libraries:

  (optional) Julia packages used by the recipe.

- verify_cpp_function_types:

  (optional) Check the types of C++ functions in the recipe.

- keep_temporary_recipe_code:

  (optional) Keep the generated C++ code.

- seed:

  (optional) The seed of the random number generator.

- algorithm:

  (optional) For a recipe in the ilike language: the name of the
  algorithm section to run (needed when the recipe has several).

## Value

Invisibly, the folder holding the results.

## Details

The model needs a Gaussian prior on every parameter, given by built-in
prior blocks (prior,ilike::norm or prior,ilike::mvnorm), and the
likelihood as a simulator: a simulate_data_model block (which includes
the observation noise), a data_variable block and optionally a
summary_statistics block. An enk_shifter block chooses the update:
ilike::stochastic (perturbed observations, the default) or
ilike::adjustment (ensemble adjustment). There is no estimate of the
evidence; read the ensembles with load_enk_output (one target per
iteration).
