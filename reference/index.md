# Package index

## RGCA MCMC functions

Functions related to MCMC parameter estimation

- [`RE_MCMC_fit()`](https://niehs.github.io/RGCA/reference/RE_MCMC_fit.md)
  : Random Effect MCMC fitting
- [`pull_parameters()`](https://niehs.github.io/RGCA/reference/pull_parameters.md)
  : Pull parameter posterior chains from an MCMC object.
- [`pull_summary_parameters()`](https://niehs.github.io/RGCA/reference/pull_summary_parameters.md)
  : Pull parameters from an MCMC chain

## RGCA Plotting functions

Functions to create plots from the manuscript

- [`get_mle_curve_fits()`](https://niehs.github.io/RGCA/reference/get_mle_curve_fits.md)
  : Maximum likelihood dose response curves

## RGCA Mixture calculators

Functions that compute mixture responses

- [`create_mix_calc()`](https://niehs.github.io/RGCA/reference/create_mix_calc.md)
  : Mixture Response Calculator Wrapper for Manuscript, with Sampling
- [`hill_function()`](https://niehs.github.io/RGCA/reference/hill_function.md)
  : Hill function
- [`hill_invs_factry()`](https://niehs.github.io/RGCA/reference/hill_invs_factry.md)
  : Inverse Hill Function Factory
- [`eff_response_opt()`](https://niehs.github.io/RGCA/reference/eff_response_opt.md)
  : Mixture Response Optimization Routine
- [`predict_mix_response_many()`](https://niehs.github.io/RGCA/reference/predict_mix_response_many.md)
  : Predict Responses for a List of Calculators, generic
- [`mix_function_generator()`](https://niehs.github.io/RGCA/reference/mix_function_generator.md)
  : Mixture Response Calculator

## RGCA Simulation functions

Functions related to the simulation study of the manuscript

## RGCA Miscellaneous functions

Other helpful functions

- [`build_replicate_matrix()`](https://niehs.github.io/RGCA/reference/build_replicate_matrix.md)
  : Builds a design matrix for replicates of dose response data for
  performing a simultaneous Gibbs update for the sill, sill random
  effect, and intercept random effect
