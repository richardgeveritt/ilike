# Load the number of particles of SMC^2's filters at each step.

In SMC^2 (an outer `SMC_with_MCMC` whose `sequence = data(...)` takes in
the observations of a factor `pf(inner, ...)`), the setting
`increase = exchange(below = b, multiply = m, max = M)` of `pf` gives
the filters more particles when the acceptance rate of the moves falls
below `b` (the exchange step of Chopin, Jacob and Papaspiliopoulos,
2013). This function reads the record of each step that moved the
particles.

## Usage

``` r
load_smc_inner_particles(results_directory)
```

## Arguments

- results_directory:

  The results directory of the run (`results_name`).

## Value

A data frame with one row per step that moved the particles and filter
with `exchange`: columns step, factor (the factor's number in the
model), acceptance_rate (the fraction of particles the step's moves
changed), particles_before and particles_after (the filters' number of
particles before and after the step's exchange; equal when there was
none).
