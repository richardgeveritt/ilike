# Loading the adaptation of adaptive MCMC moves into R memory.

For MCMC moves with an adaptor (an m_adaptor or mh_adaptor block, or an
adapt_m_proposal or adapt_mh_proposal function), each chain records the
adapted parameters of its proposal: at every stored iteration, the
scales and the running estimate of the target's mean (the trace); and at
the end of the chain, also the estimated covariance (its upper triangle)
and, for principal-components moves, the eigenvalues and directions (the
final state). Parameter names are prefixed with "moveK." when the MCMC
consists of several moves.

## Usage

``` r
load_mcmc_adaptation(results_directory, target = NULL)
```

## Arguments

- results_directory:

  The folder in which the results are stored.

- target:

  (optional) The stored SMC iteration to read (default: the last one,
  which for MCMC holds the chains).

## Value

A list with data frames trace (columns Chain, Iteration, Parameter,
Value; Iteration counts stored iterations, as in load_mcmc_output) and
final (columns Chain, Parameter, Value).
