# Loading ensemble Kalman output into R memory.

Loading ensemble Kalman output into R memory.

## Usage

``` r
load_enk_output(
  results_directory,
  ggsmc = TRUE,
  external_log_weights = c(0),
  external_target_parameters = "",
  nesting_level = 0,
  factor = 1
)
```

## Arguments

- results_directory:

  The folder in which the results are stored.

- ggsmc:

  (optional) Output in tidy format for plotting in ggsmc package.

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

A list containing the ensemble members (called particles).
