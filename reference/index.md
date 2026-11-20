# Package index

## All functions

- [`EnKF()`](https://richardgeveritt.github.io/ilike/reference/EnKF.md)
  : Ensemble Kalman filter

- [`IENKI()`](https://richardgeveritt.github.io/ilike/reference/IENKI.md)
  : Iterative ensemble Kalman inversion

- [`IS()`](https://richardgeveritt.github.io/ilike/reference/IS.md) :
  Importance sampler.

- [`KF()`](https://richardgeveritt.github.io/ilike/reference/KF.md) :
  Kalman filter

- [`MCMC()`](https://richardgeveritt.github.io/ilike/reference/MCMC.md)
  : MCMC

- [`MFDS()`](https://richardgeveritt.github.io/ilike/reference/MFDS.md)
  : Mean-field dynamical system ensemble Kalman inversion (MFDS)

- [`PF()`](https://richardgeveritt.github.io/ilike/reference/PF.md) :
  Particle filter

- [`SMC_with_MCMC()`](https://richardgeveritt.github.io/ilike/reference/SMC_with_MCMC.md)
  : SMC with an MCMC move

- [`compile()`](https://richardgeveritt.github.io/ilike/reference/compile.md)
  : Parse .ilike file to give a compiled ilike "recipe".

- [`constant_velocity_x`](https://richardgeveritt.github.io/ilike/reference/constant_velocity_x.md)
  : 10000 points simulated from a constant velocity linear-Gaussian
  dynamic model (process noise standard deviation = 1).

- [`constant_velocity_y`](https://richardgeveritt.github.io/ilike/reference/constant_velocity_y.md)
  : 10000 noisy observations (Gaussian noise, standard deviation 1) of
  10000 points simulated from a constant velocity linear-Gaussian
  dynamic model (process noise standard deviation = 1).

- [`convert_bugs_model()`](https://richardgeveritt.github.io/ilike/reference/convert_bugs_model.md)
  : Convert a BUGS, JAGS or NIMBLE model to the ilike language

- [`convert_ilike_file()`](https://richardgeveritt.github.io/ilike/reference/convert_ilike_file.md)
  : Convert a recipe from the old .ilike format to the ilike language

- [`convert_stan_model()`](https://richardgeveritt.github.io/ilike/reference/convert_stan_model.md)
  : Convert a Stan program to the ilike language

- [`evaluate_log()`](https://richardgeveritt.github.io/ilike/reference/evaluate_log.md)
  : Evaluate a model specified in a recipe.

- [`gaussian`](https://richardgeveritt.github.io/ilike/reference/gaussian.md)
  : 100 points simulated from a Gaussian distribution.

- [`ilike_pivot_longer()`](https://richardgeveritt.github.io/ilike/reference/ilike_pivot_longer.md)
  :

  Change ilike output from the wide format (one row per point, one
  column per scalar, as given by `wide = TRUE` in the load functions or
  by `ilike_pivot_wider`) to the long format (one row per scalar, with
  columns ParameterName, Dimension and Value), as used by ggsmc.

- [`ilike_pivot_wider()`](https://richardgeveritt.github.io/ilike/reference/ilike_pivot_wider.md)
  : Take output from the functions load_mcmc_output, load_smc_output or
  load_enk_output and convert it to a more standard format (non-tidy
  format with one column per dimension).

- [`load_enk_output()`](https://richardgeveritt.github.io/ilike/reference/load_enk_output.md)
  : Loading ensemble Kalman output into R memory.

- [`load_mcmc_adaptation()`](https://richardgeveritt.github.io/ilike/reference/load_mcmc_adaptation.md)
  : Loading the adaptation of adaptive MCMC moves into R memory.

- [`load_mcmc_diagnostics()`](https://richardgeveritt.github.io/ilike/reference/load_mcmc_diagnostics.md)
  : Load the per-iteration diagnostics of an MCMC run with nuts.

- [`load_mcmc_output()`](https://richardgeveritt.github.io/ilike/reference/load_mcmc_output.md)
  : Loading MCMC output into R memory.

- [`load_smc_output()`](https://richardgeveritt.github.io/ilike/reference/load_smc_output.md)
  : Loading SMC output into R memory.

- [`load_smc_tuning()`](https://richardgeveritt.github.io/ilike/reference/load_smc_tuning.md)
  : Load the tuning values of an SMC run with esjd adaptation.

- [`rdtsc_seed()`](https://richardgeveritt.github.io/ilike/reference/rdtsc_seed.md)
  : Gets a seed for random number generation by using the processor time
  stamp (calling rdtsc in C).

- [`toy_model_rwm_ggmcmc`](https://richardgeveritt.github.io/ilike/reference/toy_model_rwm_ggmcmc.md)
  : 100000 points simulated from a two-dimensional posterior using MCMC,
  in the format required by the ggmcmc package.

- [`toy_model_rwm_ilike`](https://richardgeveritt.github.io/ilike/reference/toy_model_rwm_ilike.md)
  : 100000 points simulated from a two-dimensional posterior using MCMC,
  in the format required by the ilike.output and ggsmc packages.
