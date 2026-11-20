# Check gradients against finite differences.

Compares the gradient of the log target (or of some of its factors) that
a recipe provides, from functions written by the user, from model
statements, or from `gradient = autodiff`, with central differences of
the log target itself, element by element.

## Usage

``` r
check_gradient(
  recipe,
  parameters,
  variables = names(parameters),
  index = integer(0),
  model_parameter_list = list(),
  algorithm = NULL
)
```

## Arguments

- recipe:

  A compiled ilike recipe, or the name of an ilike file.

- parameters:

  A named list: the point at which to check.

- variables:

  (optional) The variables whose gradients are checked (default: all
  those in `parameters`).

- index:

  (optional) The factors to check, numbered from 1 (default: all
  factors).

- model_parameter_list:

  (optional) Values for the constants of the model.

- algorithm:

  (optional) For a recipe with several algorithm sections: the name of
  one.

## Value

A data frame with one row per element: `variable`, `element`, `gradient`
(the recipe's), `finite_difference`, and `error`, the absolute
difference divided by `max(1, |finite_difference|)`. Errors around 1e-6
or smaller are what finite differences can resolve.
