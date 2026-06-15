# Hill function

A simple function to compute the Hill function response given parameters
and a concentration

## Usage

``` r
hill_function(a, b, c, conc)
```

## Arguments

- a:

  the maximum effect or sill

- b:

  the EC50

- c:

  the slope

- conc:

  an input concentration or dose

## Value

a real number prediction of the dose response

## Examples

``` r
hill_function(1, 1.5, 2, 3)
#> [1] 0.8
```
