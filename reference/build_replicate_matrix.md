# Builds a design matrix for replicates of dose response data for performing a simultaneous Gibbs update for the sill, sill random effect, and intercept random effect

Builds a design matrix for replicates of dose response data for
performing a simultaneous Gibbs update for the sill, sill random effect,
and intercept random effect

## Usage

``` r
build_replicate_matrix(v_list)
```

## Arguments

- v_list:

  a list with a vector for each replicate. Each vector is a sequence of
  coefficients with length equal to the number of doses (ie samples),
  computed according to the Hill model.

## Value

a matrix built from the vectors of the input to mimic a design matrix
for a linear model
