# Pull parameter posterior chains from an MCMC object.

Given an MCMC chain fit by RE_MCMC_fit, pull out the chains representing
the parameter posteriors and compute some estimates of the parameters by
applying the summary statistic function to the posterior samples.

## Usage

``` r
pull_parameters(re_chains, summry_stat = stats::median, input_replicates = NA)
```

## Arguments

- re_chains:

  output from RE_MCMC_fit consisting of a list of parameter chains,
  where each chain is a vector of posterior samples

- summry_stat:

  the statistic used to summarize the posterior. Default is median, but
  can be mean or any other similar measure of the center.

- input_replicates:

  The optional list of indices of replicates for each chemical

## Value

a list with posterior samples and summary estimates of the parameters
