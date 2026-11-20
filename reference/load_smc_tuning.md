# Load the tuning values of an SMC run with esjd adaptation.

With `smc_adapt = esjd(...)` on a move, each particle of `SMC_with_MCMC`
moves with its own tuning values, chosen from how far the particles
moved at the previous step. This function reads the summary of those
values written at each step: per tuned move and parameter, the mean,
standard deviation and 10%, 50% and 90% quantiles over the particles.

## Usage

``` r
load_smc_tuning(results_directory)
```

## Arguments

- results_directory:

  The results directory of the run (`results_name`).

## Value

A data frame with one row per move step and tuner: columns step (the
move step: 1 for the first move), tuner (its number, when several moves
or proposals are tuned), variable, and for each tuned parameter (scale,
or step and steps) its mean, sd, q10, q50 and q90; for pretune also
step_max and steps_max (the ranges of the next trial moves); for choose,
kernelK.share (the share of the particles using each proposal).
