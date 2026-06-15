# Inverse Hill Function Factory

The inverse function is used with a single dose response curve, which is
assumed to be a Hill function. A factory method is used to avoid
repeatedly passing in the curve parameters; once generated, the inverse
calculator only needs a response point y to return an inverse (dose)
required to get that response.

## Usage

``` r
hill_invs_factry(a, b, c, max_R = 1, d = 0)
```

## Arguments

- a:

  sill (max effect)

- b:

  EC50

- c:

  slope

- max_R:

  maximum effect across all chemicals

- d:

  minimum effect

## Value

a function to compute the inverse given a mixture response R

## Examples

``` r
NA
#> [1] NA
```
