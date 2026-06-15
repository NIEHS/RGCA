# Mixture Response Calculator Wrapper for Manuscript, with Sampling

A factory method that returns a function that computes the mixture
response. Samples directly from the posterior MCMC chains for the Hill
parameters of slope, sill and EC50. Noise is added to some parameters
according to the random effect variance.

## Usage

``` r
create_mix_calc(idx, par_list, add_RE = TRUE, unit_slopes = FALSE)
```

## Arguments

- idx:

  Specifies which clustering to apply to the parameters

- par_list:

  a data frame with individual chemical dose response parameters

- add_RE:

  boolean to include or exclude random effect variances

- unit_slopes:

  boolean to fix slopes to 1 (but still use slope clustering) Used for
  special case of GCA, where slope = 1

## Value

function to take a concentration vector as input and response as output
