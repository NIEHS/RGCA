# Maximum likelihood dose response curves

Wrapper function to apply the drc (dose response curve) package to the
input data.

## Usage

``` r
get_mle_curve_fits(y_i, Cx, replicate_sets)
```

## Arguments

- y_i:

  matrix of observed dose responses (columns) for multiple chemicals and
  replicates (rows)

- Cx:

  a matrix of the doses corresponding to the responses

- replicate_sets:

  a list of vectors where each vector contains the indices of replicates
  for a particular chemical

## Value

a matrix of parameters with a row for each chemical and columns
corresponding to ()
