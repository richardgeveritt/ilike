# Posterior expectations from a signed pseudo-marginal MCMC run.

A chain whose factor estimates may be negative (normaliser =
estimate(..., debias = roulette(...) or block_poisson(...))) targets the
absolute value of the estimated posterior and records each iteration's
sign. The expectation of a function of the parameters under the
posterior is then the mean of the function weighted by the signs (Lyne
et al. 2015): sum(s_i f_i)/sum(s_i).

## Usage

``` r
signed_mean(output, parameter, dimension = 1, f = identity, burn_in = 0)
```

## Arguments

- output:

  The output of load_mcmc_output, with its Sign column.

- parameter:

  The name of a parameter.

- dimension:

  Its element (default 1).

- f:

  A function of the parameter's values (default: the identity).

- burn_in:

  The number of first iterations of each chain to leave out.

## Value

A list: the estimate, and the mean sign (near 0 makes the estimate
unreliable).
