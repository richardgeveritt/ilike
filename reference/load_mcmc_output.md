# Loading MCMC output into R memory.

Loading MCMC output into R memory.

## Usage

``` r
load_mcmc_output(
  results_directory,
  ggmcmc = NULL,
  ilike.output = NULL,
  wide = FALSE
)
```

## Arguments

- results_directory:

  The folder in which the results are stored.

- ggmcmc:

  (optional) Output in tidy format for plotting in ggmcmc package
  (default is FALSE if ilike.output is set to TRUE or not set and TRUE
  otherwise).

- ilike.output:

  (optional) Output can be processed by ilike,output package (default is
  TRUE if ggmcmc is set to FALSE or is not set and FALSE otherwise).
  Choosing TRUE for this argument is incompatible with ggmcmc=TRUE,
  since the two packages require different formatting of the output.

- wide:

  (optional) If TRUE, return one row per point, with one column per
  scalar (named variable_dimension, e.g. theta_1), instead of one row
  per scalar (default FALSE). Much faster for large outputs;
  ilike_pivot_longer converts it to the long format.

## Value

A list containing the MCMC chains.
