# Loading SMC output into R memory.

Loading SMC output into R memory.

## Usage

``` r
load_smc_output(
  results_directory,
  ggmcmc = FALSE,
  ggsmc = TRUE,
  ilike.output = TRUE,
  as.mcmc = FALSE,
  as.enk = FALSE,
  which.targets = NULL,
  directory_prefix = "ilike",
  external_log_weights = c(0),
  external_target_parameters = "",
  nesting_level = 0,
  factor = 1
)
```

## Arguments

- results_directory:

  The folder in which the results are stored.

- ggmcmc:

  (optional) Output in tidy format for plotting in ggmcmc package.

- ggsmc:

  (optional) Output in tidy format for plotting in ggsmc package.

- ilike.output:

  (optional) Output can be processed by ilike,output package (default is
  TRUE). Choosing TRUE for this argument is incompatiable with
  ggmcmc=TRUE, since the two packages require different formatting of
  the output.

- as.mcmc:

  (optional) Output treats particles as different MCMC chains.

- as.enk:

  (optional) Output treats particles as an ensemble.

- which.targets:

  (optional) The indices of the targets to output (defaults to all).

- directory_prefix:

  (optional; for nested output only) The first part of the name of the
  directory within results_directory that contains the results. (default
  is "ilike", giving a directory of results_directory/ilike_smc)

- external_log_weights:

  (optional; for nested output only) The weights of the importance
  points external to the current folder. (default is 1, to be used at
  the top level of nested output)

- external_target_parameters:

  (optional; for nested output only) The parameters of the target
  external to the current folder. (default is "", corresponding to no
  parameters)

- nesting_level:

  (optional; for nested output only) The level of nesting at which to
  extract points. (default is 0, representing the top level of nested
  output)

- factor:

  (optional; for nested output only) The factor from which to extract
  points. (default is 0)

## Value

A list containing the SMC output.
