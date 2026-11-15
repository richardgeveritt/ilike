# Change ilike output from the wide format (one row per point, one column per scalar, as given by `wide = TRUE` in the load functions or by `ilike_pivot_wider`) to the long format (one row per scalar, with columns ParameterName, Dimension and Value), as used by ggsmc.

Change ilike output from the wide format (one row per point, one column
per scalar, as given by `wide = TRUE` in the load functions or by
`ilike_pivot_wider`) to the long format (one row per scalar, with
columns ParameterName, Dimension and Value), as used by ggsmc.

## Usage

``` r
ilike_pivot_longer(output)
```

## Arguments

- output:

  Output in the wide format.

## Value

The output in the long format: the index columns, then ParameterName,
Dimension and Value; rows ordered point by point, as in the load
functions.
