# Convert a BUGS, JAGS or NIMBLE model to the ilike language

Translates a model written in the BUGS language (as used by WinBUGS,
OpenBUGS and JAGS) or in NIMBLE's dialect of it into a model in the
ilike language, with a data section holding the data given (which also
tells which variables are data: the others are parameters). Translated:
stochastic relations (~), deterministic relations (\<-, with logit(p)
\<- … and log(x) \<- … on the left), for loops, arithmetic and the
functions exp, log, sqrt, pow, abs, logit, ilogit (and NIMBLE's expit),
sin, cos, tanh, min, max, and inprod and sum over vectors whose lengths
are known from the data; truncation T(a, b) (and BUGS's I(a, b)) for
dnorm. Distributions (with BUGS's parameterisations, or NIMBLE's when
dialect = "nimble"): dnorm, dlnorm, dt, dgamma, dbeta, dunif, dexp,
dchisqr (dchisq), dweib, ddexp, dlogis, dpois, dbern, dbin, dnegbin,
dcat, dmulti, ddirch (ddirich), dmnorm, dmt and dwish (whose precision
matrices must be data: they are inverted into new data, name_inverse).
Missing values in the data, censoring (dinterval), other truncations,
precision matrices that are parameters, and other functions stop with a
message.

An algorithm section is added with a NUTS move on all the parameters
when they are all continuous and unconstrained vectors or scalars.

## Usage

``` r
convert_bugs_model(
  model,
  data = list(),
  dialect = c("jags", "bugs", "nimble"),
  output = NULL,
  overwrite = FALSE
)
```

## Arguments

- model:

  The model: a file name, the code as a character string, or an R
  expression (for example from nimbleCode() or quote()).

- data:

  A list with the data (and constants), by name.

- dialect:

  (optional) "jags" (the default), "bugs" or "nimble": the
  parameterisations of the distributions (JAGS and BUGS share them).

- output:

  (optional) The file to write.

- overwrite:

  (optional) If TRUE, an existing output file is replaced.

## Value

Invisibly, the converted text.
