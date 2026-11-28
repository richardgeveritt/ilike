# Posterior expectations from a run with factor estimates that may be negative.

A run whose factor estimates may be negative (normaliser = estimate(...,
debias = roulette(...), block_poisson(...) or taylor(...))) targets the
absolute value of the estimated posterior and records each point's sign.
The expectation of a function of the parameters under the posterior is
then the mean of the function weighted by the signs (Lyne et al. 2015):
sum(s_i f_i)/sum(s_i) for MCMC, and sum(W_i s_i f_i)/sum(W_i s_i) for
importance sampling and SMC, with W_i the points' normalised weights.

## Usage

``` r
signed_mean(
  output,
  parameter,
  dimension = 1,
  f = identity,
  burn_in = 0,
  target = NULL
)
```

## Arguments

- output:

  The output of load_mcmc_output, or of load_smc_output, with its Sign
  column.

- parameter:

  The name of a parameter.

- dimension:

  Its element (default 1).

- f:

  A function of the parameter's values (default: the identity).

- burn_in:

  MCMC: the number of first iterations of each chain to leave out.

- target:

  IS and SMC: the target whose points are used (default: the last).

## Value

A list: the estimate, and the (weighted) mean sign (near 0 makes the
estimate unreliable).
