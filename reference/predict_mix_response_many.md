# Predict Responses for a List of Calculators, generic

A similar function to predict_mix_response. This is a convenience
function that takes a list of lists as input and iteratively applies the
elements, which are bootstrapped mixture repsonse calculators, to the
columns of the mixture dose matrix. For example, three different
clustering approaches can be used for the outer list while each inner
list has 50 bootstrapped mixture response estimators. Then the outer
list has length 3 and each inner list has length 50.

## Usage

``` r
predict_mix_response_many(
  n_dose,
  chem_conc_matr,
  bootstrap_calc_list,
  default_entry = 0
)
```

## Arguments

- n_dose:

  number of mixture doses

- chem_conc_matr:

  a matrix where the rows represent the constituent chemicals and the
  columns represent the dose. The column sum is the mixture dose.

- bootstrap_calc_list:

  a list of lists of boostrapped functions, instances produced by the
  factory method mix_function_generator

- default_entry:

  a default entry for the methods, with default=0. Can be set to null NA
  to aid in plotting

## Value

A list of matrices where each matrix has dose by column and the
bootstrapped response prediction as a row.
