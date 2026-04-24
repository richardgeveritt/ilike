# Take output from the functions load_mcmc_output, load_smc_output or load_enk_output and convert it to a more standard format (non-tidy format with one column per dimension).

Take output from the functions load_mcmc_output, load_smc_output or
load_enk_output and convert it to a more standard format (non-tidy
format with one column per dimension).

## Usage

``` r
ilike_pivot_wider(output, variables = NULL)
```

## Arguments

- output:

  Output from the functions load_mcmc_output, load_smc_output or
  load_enk_output.

- variables:

  (optional) Variables to include in the output (default is all).
