# Pull parameters from an MCMC chain

Given an MCMC chain fit by RE_MCMC_fit, estimate the parameter values
using a summary statistic applied to the posterior samples.

## Usage

``` r
pull_summary_parameters(re_chains, summry_stat = stats::median)
```

## Arguments

- re_chains:

  output from RE_MCMC_fit consisting of a list of parameter chains,
  where each chain is a vector of posterior samples

- summry_stat:

  the statistic used to summarize the posterior. Default is median, but
  can be mean or any other similar measure of the center.

## Value

a list of parameter estimates
