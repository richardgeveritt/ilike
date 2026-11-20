# Load the per-iteration diagnostics of an MCMC run with nuts.

A `nuts` move records, at every stored iteration of each chain, the
depth of the tree it built (the trajectory doubled that many times), the
number of leapfrog steps (gradient evaluations), whether the trajectory
diverged (the error in the Hamiltonian passed 1000, a sign that the step
is too large for the curvature of the target there), and the acceptance
statistic (the average acceptance probability over the last doubling,
which the adaptors drive towards their target). Many divergences after
warm-up mean the chain may miss part of the target; a tree depth often
at max_depth means the step is too small or the target poorly scaled.

## Usage

``` r
load_mcmc_diagnostics(results_directory, target = NULL)
```

## Arguments

- results_directory:

  The folder in which the results are stored.

- target:

  (optional) The stored SMC iteration to read (default: the last one,
  which for MCMC holds the chains).

## Value

A data frame with one row per stored iteration: Chain, Iteration
(counting stored iterations, as in load_mcmc_output), tree_depth,
leapfrog_steps, divergent (logical) and accept_stat.
