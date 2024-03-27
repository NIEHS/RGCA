library(truncnorm)
library(CVXR)
library(coda)
library(ggplot2)
library(scoringRules)
library(tables) # for converting data frames to latex tabulars

#' Hill function
#'
#' A simple function to compute the Hill function response given parameters and
#' a concentration
#'
#' @param a the maximum effect or sill
#' @param b the EC50
#' @param c the slope
#' @param conc an input concentration or dose
#'
#' @return a real number prediction of the dose response
#' @export
#'
#' @examples
hill_function <- function(a, b, c, conc) {
  a / (1 + (b / conc)^c)
}

#' Maximum likelihood dose response curves
#'
#' Wrapper function to apply the drc (dose response curve) package to the input
#' data.
#'
#' @param y_i matrix of observed dose responses (columns) for multiple chemicals
#'   and replicates (rows)
#' @param Cx a matrix of the doses corresponding to the responses
#' @param replicate_sets a list of vectors where each vector contains the
#'   indices of replicates for a particular chemical
#'
#' @return a matrix of parameters with a row for each chemical and columns
#'        corresponding to ()
#' @export
#'
#' @examples
get_mle_curve_fits <- function(y_i, Cx, replicate_sets) {
  # individual fits
  curve_fits <- matrix(0, nrow = length(replicate_sets), ncol = 3)
  for (ri in seq(length(replicate_sets))) {
    drc_yi <- c(y_i[replicate_sets[[ri]], ]) # nolint
    drc_xi <- c(Cx[replicate_sets[[ri]], ]) # nolint
    fit_coeff <- try(drc::drm(drc_yi ~ drc_xi,
                              fct = drc::LL.3())$coefficients,
                     silent = TRUE)
    if (is(fit_coeff, "try-error")) next
    fit_coeff <- fit_coeff[c(2, 3, 1)]
    fit_coeff[3] <- -fit_coeff[3]
    fit_coeff[2] <- fit_coeff[2]
    curve_fits[ri, ] <- fit_coeff
  }
  return(curve_fits)
}



get_MCMC_diagnostics <- function(re_chains, re_chains2, re_chains3) {
  # explicit call to global env
  n_chems <- .GlobalEnv$n_chems
  get_GR_Diagnostic <- function(data_id, input_cols) {
    mcmc1 <- coda::as.mcmc(re_chains[[data_id]][, input_cols])
    mcmc2 <- coda::as.mcmc(re_chains2[[data_id]][, input_cols])
    mcmc3 <- coda::as.mcmc(re_chains3[[data_id]][, input_cols])
    mcmclist <- coda::as.mcmc.list(list(mcmc1, mcmc2, mcmc3))
    coda_test <- coda::gelman.diag(mcmclist)$psrf
    return(coda_test)
  }
  coda_u_RE_sd <- get_GR_Diagnostic("u_RE_sd", seq(ncol(re_chains$u_RE_sd)))
  coda_sigma <- get_GR_Diagnostic("sigma", seq(ncol(re_chains$sigma)))
  coda_slope <- get_GR_Diagnostic("slope_record",
                                  seq(ncol(re_chains$slope_record)))
  coda_sill <- get_GR_Diagnostic("sill_mideffect_record", 1:n_chems)
  # EC50 part, needs log
  chain_sill_mid <- re_chains$sill_mideffect_record
  mcmc1 <- coda::as.mcmc(log(chain_sill_mid[, (1 + n_chems):(2 * n_chems)]))
  mcmc2 <- coda::as.mcmc(log(chain_sill_mid[, (1 + n_chems):(2 * n_chems)]))
  mcmc3 <- coda::as.mcmc(log(chain_sill_mid[, (1 + n_chems):(2 * n_chems)]))
  mcmclist <- coda::as.mcmc.list(list(mcmc1, mcmc2, mcmc3))
  coda_ec50 <- coda::gelman.diag(mcmclist)$psrf
  all_dats <- round(cbind(coda_sill,
                          coda_ec50,
                          coda_slope,
                          coda_sigma,
                          coda_u_RE_sd),
                    digits = 2)
  all_dats <- as.data.frame(all_dats)
  new_names <- array(sapply(c("Sill", "EC50", "Slope", "Sigma", "Variance_U"),
                            FUN = function(sx) {
                              (paste(sx, c("Point est.", "Upper C.I.")))
                            }))
  names(all_dats) <- new_names
  rownames(all_dats) <- chem_map
  return(all_dats)
}





# Some notes on posterior summary stats:
# mode (MAP) is usually the lowest, followed by median and mean due to skewness
# MAP or median makes sense for heavily skewed parameters, but slope was not so
# skewed due to strong prior
modes <- function(dats) {
  if (all(dats == 0)) {
    return(0)
  }
  dens <- density(dats)
  i <- which(diff(sign(diff(dens$y))) < 0) + 1
  max_i <- i[which(dens$y[i] == max(dens$y[i]))]
  dens$x[max_i]
}

#' Pull parameters from an MCMC chain
#'
#' Given an MCMC chain fit by RE_MCMC_fit, estimate the parameter values using a
#' summary statistic applied to the posterior samples.
#'
#' @param re_chains output from RE_MCMC_fit consisting of a list of parameter
#'   chains, where each chain is a vector of posterior samples
#' @param summry_stat the statistic used to summarize the posterior.  Default is
#'   median, but can be mean or any other similar measure of the center.
#'
#' @return a list of parameter estimates
#' @export
#'
#' @examples
pull_summary_parameters <- function(re_chains,
                                    summry_stat = median) {
  burnin <- 5000
  # explicit call to global env
  n_chems <- ncol(re_chains$slope_record)
  re_iter <- nrow(re_chains$slope_record)
  if (re_iter <= burnin) burnin <- 100
  thin_idx <- seq(burnin, re_iter - 100, by = 20)
  chain_sill_ec_param <- re_chains$sill_mideffect_record
  sill_params <- apply(chain_sill_ec_param[thin_idx, 1:n_chems],
                       MARGIN = 2,
                       summry_stat)
  sill_sd <- apply(chain_sill_ec_param[thin_idx, 1:n_chems],
                   MARGIN = 2,
                   sd)
  ec50_params <- apply(chain_sill_ec_param[thin_idx, 1:n_chems + n_chems],
                       MARGIN = 2,
                       summry_stat)
  ec50_sd <- apply(chain_sill_ec_param[thin_idx, 1:n_chems + n_chems],
                   MARGIN = 2,
                   sd)
  u_RE_params <- apply(re_chains$u_RE[thin_idx, ], MARGIN = 2, summry_stat)
  v_RE_params <- apply(re_chains$v_RE[thin_idx, ], MARGIN = 2, summry_stat)
  u_RE_sd_params <- apply(re_chains$u_RE_sd[thin_idx, ],
                          MARGIN = 2,
                          summry_stat)
  v_RE_sd_params <- apply(re_chains$v_RE_sd[thin_idx, ],
                          MARGIN = 2,
                          summry_stat)
  phi_c <- apply(re_chains$slope_record[thin_idx, ],
                 MARGIN = 2,
                 summry_stat)
  return(list(
    "sill_params" = sill_params,
    "sill_sd" = sill_sd,
    "ec50_params" = ec50_params,
    "ec50_stdev" = ec50_sd,
    "u_RE_params" = u_RE_params,
    "v_RE_params" = v_RE_params,
    "u_RE_sd_params" = u_RE_sd_params,
    "v_RE_sd_params" = v_RE_sd_params,
    "slope_params" = phi_c
  ))
}



#' Pull parameter posterior chains from an MCMC object.
#'
#' Given an MCMC chain fit by RE_MCMC_fit, pull out the chains representing the
#' parameter posteriors and compute some estimates of the parameters by applying
#' the summary statistic function to the posterior samples.
#'
#' @param re_chains output from RE_MCMC_fit consisting of a list of parameter
#'   chains, where each chain is a vector of posterior samples
#' @param summry_stat the statistic used to summarize the posterior.  Default is
#'   median, but can be mean or any other similar measure of the center.
#'
#' @return a list with posterior samples and summary estimates of the parameters
#' @export
#'
#' @examples
pull_parameters <- function(re_chains, summry_stat = median,
                            input_replicates = NA) {
  if (is.list(input_replicates)) replicate_sets <- input_replicates
  burnin <- 5000
  # explicit call to global env
  n_chems <- ncol(re_chains$slope_record)
  re_iter <- nrow(re_chains$slope_record)
  if (re_iter <= burnin) burnin <- 100
  thin_idx <- seq(burnin, re_iter - 100, by = 20)
  sill_sample <- re_chains$sill_mideffect_record[thin_idx, 1:n_chems]
  sill_params <- apply(re_chains$sill_mideffect_record[thin_idx, 1:n_chems],
                       MARGIN = 2,
                       summry_stat)
  ec50_sample <- re_chains$sill_mideffect_record[thin_idx, 1:n_chems + n_chems]
  slope_sample <- re_chains$slope_record[thin_idx, 1:n_chems]
  u_RE_center <- sapply(1:n_chems,
                        FUN = function(ridx) {
                          median(re_chains$u_RE[thin_idx,
                                                replicate_sets[[ridx]]])
                        })
  u_RE_sd_est <- apply(re_chains$u_RE_sd[thin_idx, ], MARGIN = 2, summry_stat)
  v_RE_sd_est <- apply(re_chains$v_RE_sd[thin_idx, ], MARGIN = 2, summry_stat)
  # slopes must be esimated for DP clustering
  phi_c <- apply(re_chains$slope_record[thin_idx, ], MARGIN = 2, summry_stat)
  return(list(
    "sill_sample" = sill_sample,
    "sill_params" = sill_params,
    "ec50_sample" = ec50_sample,
    "u_RE_center" = u_RE_center,
    "u_RE_sd" = u_RE_sd_est,
    "v_RE_sd" = v_RE_sd_est,
    "slope_sample" = slope_sample,
    "slope_params" = phi_c
  ))
}

pull_parameters_nimble <- function(nimble_samples, summry_stat = median) {
  sill_sample <- nimble_samples[, 1:n_chems]
  ec50_idx <- grep("theta", colnames(nimble_samples))
  ec50_sample <- nimble_samples[, ec50_idx]
  slope_idx <- grep("beta", colnames(nimble_samples))
  slope_sample <- nimble_samples[, slope_idx]
  u_RE_sd_idx <- grep("sigma_u", colnames(nimble_samples))
  u_RE_sd_est <- nimble_samples[, u_RE_sd_idx]
  v_RE_sd_idx <- grep("sigma_v", colnames(nimble_samples))
  v_RE_sd_est <- nimble_samples[, v_RE_sd_idx]
  u_RE_cent_idx <- setdiff(grep("u", colnames(nimble_samples)), u_RE_sd_idx)
  u_RE_samples <- nimble_samples[, u_RE_cent_idx]
  u_RE_center <- sapply(1:n_chems,
                        FUN = function(ridx) {
                          median(u_RE_samples[, replicate_sets[[ridx]]])
                        })
  # slopes must be estimated for DP clustering
  phi_c <- apply(slope_sample, MARGIN = 2, summry_stat)
  return(list(
    "sill_sample" = sill_sample,
    "ec50_sample" = ec50_sample,
    "u_RE_center" = u_RE_center,
    "u_RE_sd" = u_RE_sd_est,
    "v_RE_sd" = v_RE_sd_est,
    "slope_sample" = slope_sample,
    "slope_params" = phi_c
  ))
}




random_assignment <- function(n_chems, n_clust) {
  # assume either random clustering or use DP
  clust_rand <- sample(1:n_clust, size = n_chems, replace = TRUE)
  unique_assign <- unique(clust_rand)
  unique_id <- sapply(clust_rand,
                      FUN = function(c_entry) which(c_entry == unique_assign))
  return(unique_id)
}

#--------------------- Mix Calulators ------------------- ####
#' These functions are wrappers for the factory method mix_function_generator.
#' The wrappers provide different ways of sampling the parameters.


#' Creates a function that computes mixture response.
#'
#' A factory method that returns a function that computes the mixture response.
#' Uses summary statistics to create a random sample rather than directly
#' sampling from the posterior MCMC chains
#'
#' @param idx Specifies which clustering to apply to the parameters
#' @param par_list Contains all estimated parameters and DP cluster options
#'
#' @return a function that takes a concentration vector as input and yields a
#' predicted response as output
#' ie a "calculator" for the mixture effect given component concentrations
#' @export
#'
create_mix_calc_from_summary <- function(idx, par_list) {
  # bootstrap the cluster slope from top clusters
  cluster_name <- names(par_list$cluster_assign[idx])
  cluster_assign_vec <- as.numeric(strsplit(cluster_name, split = " ")[[1]])
  slope_mean <- par_list$centers
  # For summary, dont need: slope_sd <- par_list$cent_sd
  if (is.matrix(par_list$centers)) {
    slope_mean <- par_list$centers[idx, ]
    # dont need again: slope_sd <- par_list$cent_sd[idx, ]
  }
  # must sample both from same idx since diff idx may have small likelihood
  n_chem_pars <- length(par_list$sill_params)
  tot_sd <- sqrt(par_list$u_RE_sd_params^2 + par_list$sill_sd^2)
  sill_params_boot <- truncnorm::rtruncnorm(
    n = n_chem_pars,
    mean = par_list$sill_params,
    sd = tot_sd
  )
  ec50_params_boot <- truncnorm::rtruncnorm(
    n = n_chem_pars,
    a = 0,
    mean = par_list$ec50_params,
    sd = par_list$ec50_stdev
  )
  #if bootstrap: rnorm(n = n_chem_pars,mean = 0,sd = par_list$v_RE_sd_params)
  max_effect_R <- max(par_list$sill_params, sill_params_boot)
  # test: control has maxR=100
  bootstrap_param_matrix <- as.matrix(cbind(
    "a" = sill_params_boot,
    "b" = ec50_params_boot,
    "c" = slope_mean,
    "max_R" = max_effect_R,
    "d" = 0
  ))
  mixture_calculator <- mix_function_generator(
    bootstrap_param_matrix,
    cluster_assign_vec
  )
  return(mixture_calculator)
}


#' Mixture Response Calculator Wrapper for Manuscript, with Sampling
#'
#' A factory method that returns a function that computes the mixture response.
#' Samples directly from the posterior MCMC chains for the Hill parameters of
#' slope, sill and EC50.  Noise is added to some parameters according to the
#' random effect variance.
#'
#' @param idx Specifies which clustering to apply to the parameters
#' @param par_list a data frame with individual chemical dose response
#'   parameters
#' @param add_RE boolean to include or exclude random effect variances
#' @param unit_slopes boolean to fix slopes to 1 (but still use slope
#'   clustering) Used for special case of GCA, where slope = 1
#'
#' @return function to take a concentration vector as input and response as
#'   output
#' @export
create_mix_calc <- function(idx, par_list, add_RE = TRUE, unit_slopes = FALSE) {
  # bootstrap the cluster slope from top clusters
  cluster_name <- names(par_list$cluster_assign[idx])
  cluster_assign_vec <- as.numeric(strsplit(cluster_name, split = " ")[[1]])
  # sample sill parameter
  n_chem_pars <- ncol(par_list$sill_sample)
  sill_params_samp <- apply(par_list$sill_sample,
                            MARGIN = 2,
                            FUN = function(colx) sample(colx, 1))
  # sample ec50
  ec50_params_samp <- apply(par_list$ec50_sample,
                            MARGIN = 2,
                            FUN = function(colx) sample(colx, 1))
  # sample slope
  slope_params_samp <- apply(par_list$slope_sample,
                             MARGIN = 2,
                             FUN = function(colx) sample(colx, 1))
  if (unit_slopes) slope_params_samp <- slope_params_samp * 0 + 1
  # add iid gauss noise for RE variance
  sill_params_boot <- sill_params_samp
  sill_RE_boot <- rnorm(n_chem_pars,
                        mean = par_list$u_RE_center,
                        sd = par_list$u_RE_sd)
  if (add_RE) sill_params_boot <- sill_params_boot + sill_RE_boot
  # We dont sample intercept RE: assume intercept 0, d=0
  # TODO consider empirical max, or clusterwise max
  max_effect_R <- max(sill_params_boot)
  # generate calculator for
  bootstrap_param_matrix <- as.matrix(cbind(
    "a" = sill_params_boot,
    "b" = ec50_params_samp,
    "c" = slope_params_samp,
    "max_R" = max_effect_R,
    "d" = 0
  ))
  mixture_calculator <- mix_function_generator(
    bootstrap_param_matrix,
    cluster_assign_vec
  )
  return(mixture_calculator)
}

#' Mixture Response Calculator Wrapper for Manuscript, using DP
#'
#' A factory method that returns a function that computes the mixture response.
#' Samples directly from the posterior MCMC chains for the EC50 and sill
#' parameters. Slope is taken from the Dirichlet Process clustering.  Noise is
#' added to some parameters according to the random effect variance.
#'
#'
#' @param idx Specifies which clustering to apply to the parameters
#' @param par_list a data frame with individual chemical dose response
#'   parameters
#' @param add_RE boolean to include or exclude random effect variances
#' @param unit_slopes boolean to fix slopes to 1 (but still use slope
#'   clustering) Used for special case of GCA, where slope = 1
#'
#' @return function to take a concentration vector as input and response as
#'   output
#' @export
create_mix_calc_clustered <- function(idx,
                                      par_list,
                                      add_RE = TRUE,
                                      unit_slopes = FALSE) {
  # bootstrap the cluster slope from top clusters
  cluster_name <- names(par_list$cluster_assign[idx])
  cluster_assign_vec <- as.numeric(strsplit(cluster_name, split = " ")[[1]])
  # sample the cluster value rather than use mean/sd
  slope_mean <- par_list$centers[idx, ]
  # if add noise from DP cluster to centers:
  # sample rnorm(mean = 0, sd = unique(par_list$cent_sd[idx, ])
  # add to mean using slope_mean + slope_noise[cluster_assign_vec]
  if (unit_slopes) slope_mean <- slope_mean * 0 + 1
  # sample values for the sill parameter, can be independent
  n_chem_pars <- ncol(par_list$sill_sample)
  sill_params_samp <- apply(par_list$sill_sample,
                            MARGIN = 2,
                            FUN = function(colx) sample(colx, 1))
  # add iid gauss noise for RE variance
  sill_RE_boot <- rnorm(n_chem_pars,
                        mean = par_list$u_RE_center,
                        sd = par_list$u_RE_sd)
  sill_params_boot <- sill_params_samp
  if (add_RE) sill_params_boot <- sill_params_boot + sill_RE_boot
  # sample ec50
  ec50_params_boot <- apply(par_list$ec50_sample,
                            MARGIN = 2,
                            FUN = function(colx) sample(colx, 1))
  # We dont sample intercept RE: assume intercept 0, d=0
  # TODO consider empirical max, or clusterwise max
  max_effect_R <- max(sill_params_boot)
  # generate calculator for
  bootstrap_param_matrix <- as.matrix(cbind("a" = sill_params_boot,
                                            "b" = ec50_params_boot,
                                            "c" = slope_mean,
                                            "max_R" = max_effect_R,
                                            "d" = 0))
  mixture_calculator <- mix_function_generator(bootstrap_param_matrix,
                                               cluster_assign_vec)
  return(mixture_calculator)
}

#' Mixture Response Calculator Wrapper for Manuscript, matched index
#'
#' A slight variation on the standard create_mix_calc, the parameters are still
#' sampled from the posterior MCMC but they are sampled with a single index so
#' that the set of parameters is feasible. Sampling all the parameters randomly
#' could cause an issue when there is non-identifiability, and an extreme value
#' in one parameter is not offset by a small value in another parameter.
#'
#' @param idx Specifies which clustering to apply to the parameters
#' @param par_list a data frame with individual chemical dose response
#'   parameters
#' @param add_RE boolean to include or exclude random effect variances
#'
#' @return a function to take a concentration vector as input and response as
#'   output
#'
create_mix_calc_sample_row <- function(idx, par_list, add_RE = TRUE) {
  # bootstrap the cluster slope from top clusters
  cluster_name <- names(par_list$cluster_assign[idx])
  cluster_assign_vec <- as.numeric(strsplit(cluster_name, split = " ")[[1]])
  # sample the cluster value rather than use mean/sd
  slope_mean <- par_list$centers[idx, ]
  # if add noise from DP cluster to centers:
  # sample rnorm(mean = 0, sd = unique( par_list$centers[idx, ]))
  # add to mean using slope_mean + slope_noise[cluster_assign_vec]
  # sample values for the sill parameter, can be independent
  n_chem_pars <- ncol(par_list$sill_sample)
  n_samples <- nrow(par_list$sill_sample)
  param_idx <- sample(1:n_samples, n_chem_pars, replace = TRUE)
  sill_params_samp <- rep(0, n_chem_pars)
  ec50_params_boot <- rep(0, n_chem_pars)
  for (pidx in 1:n_chem_pars) {
    sill_params_samp[pidx] <- par_list$sill_sample[param_idx[pidx], pidx]
    ec50_params_boot[pidx] <- par_list$ec50_sample[param_idx[pidx], pidx]
  }
  # add iid gauss noise for RE variance
  sill_RE_boot <- rnorm(n_chem_pars, mean = 0, sd = par_list$u_RE_sd)
  sill_params_boot <- sill_params_samp
  if (add_RE) sill_params_boot <- sill_params_boot + sill_RE_boot
  # We dont sample intercept RE: assume intercept 0, d=0
  # TODO consider empirical max, or clusterwise max
  max_effect_R <- max(sill_params_boot)
  # generate calculator for
  bootstrap_param_matrix <- as.matrix(cbind(
    "a" = sill_params_boot,
    "b" = ec50_params_boot,
    "c" = slope_mean,
    "max_R" = max_effect_R,
    "d" = 0
  ))
  mixture_calculator <- mix_function_generator(
    bootstrap_param_matrix,
    cluster_assign_vec
  )
  return(mixture_calculator)
}



#' Predict Responses for a List of Calculators
#'
#' This is a convenience function that applies the lists of bootstrapped mixture
#' response calculators to the matrix of concentrations that define the doses
#' for the mixture.
#'
#'
#' @param sampled_mix_funs list of boostrapped functions for the main method of
#'   the manuscript, RGCA+DP, instances produced by the factory method
#'   mix_function_generator
#' @param sampled_mix_funs_GCA list of boostrapped functions for the standard
#'   method of GCA
#' @param sampled_mix_funs_IA list of boostrapped functions for the standard
#'   method of IA
#' @param n_dose integer number of doses of the mixture
#' @param chem_conc_matr a matrix where the rows represent the constituent
#'   chemicals and the columns represent the dose.  The column sum is the
#'   mixture dose.
#' @param default_entry  a default entry for the methods, with default=0.  Can
#'   be set to null NA to aid in plotting.
#'
#' @return
#'
predict_mix_response <- function(sampled_mix_funs, sampled_mix_funs_GCA,
                                 sampled_mix_funs_IA, n_dose, chem_conc_matr,
                                 default_entry = 0) {
  n_bootstraps <- length(sampled_mix_funs)
  curve_data <- matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
  curve_data_GCA <- matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
  curve_data_IA <- matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
  for (row_idx in 1:n_dose) {
    # each row of conc_matrix is a mixture dose
    conc_val <- chem_conc_matr[row_idx, ]
    # if there are missing conc values, skip prediction
    if (any(is.na(conc_val))) next
    curve_data[, row_idx] <- sapply(sampled_mix_funs,
                                    FUN = function(x) x(conc_val))
    curve_data_GCA[, row_idx] <- sapply(sampled_mix_funs_GCA,
                                        FUN = function(x) x(conc_val))
    curve_data_IA[, row_idx] <- sapply(sampled_mix_funs_IA,
                                       FUN = function(x) x(conc_val))
  }
  return(list(curve_data, curve_data_GCA, curve_data_IA))
}


#' Predict Responses for a List of Calculators, generic
#'
#' A similar function to predict_mix_response.  This is a convenience function
#' that takes a list of lists as input and iteratively applies the elements,
#' which are bootstrapped mixture repsonse calculators, to the columns of the
#' mixture dose matrix.
#'
#'
#' @param n_dose number of mixture doses
#' @param chem_conc_matr a matrix where the rows represent the constituent
#'   chemicals and the columns represent the dose.  The column sum is the
#'   mixture dose.
#' @param bootstrap_calc_list a list of lists of boostrapped functions,
#'   instances produced by the factory method mix_function_generator
#' @param default_entry a default entry for the methods, with default=0.  Can be
#'   set to null NA to aid in plotting
#'
#' @return
#' @export
predict_mix_response_many <- function(n_dose,
                                      chem_conc_matr,
                                      bootstrap_calc_list,
                                      default_entry = 0) {
  n_bootstraps <- length(bootstrap_calc_list[[1]])
  curve_list <- list()
  for (bootcalc in bootstrap_calc_list) {
    curve_data <- matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
    for (row_idx in 1:n_dose) {
      # each row of conc_matrix is a mixture dose
      conc_val <- chem_conc_matr[row_idx, ]
      # if there are missing conc values, skip prediction
      if (any(is.na(conc_val))) next
      curve_data[, row_idx] <- sapply(bootcalc,
                                      FUN = function(x) x(conc_val))
    }
    curve_list <- c(curve_list, list(curve_data))
  }
  names(curve_list) <- names(bootstrap_calc_list)
  return(curve_list)
}

#' Compute Scores for Mixture Predictions
#'
#' Computes the log likelihood, mean square error, and continuous rank
#' probability score for predictions of the mixture response compared to the
#' observed mixture response.
#'
#' @param mix_df data frame with the observed mixture responses
#' @param mix_idx integer index specifying which observed data should be used
#' @param unplotted_repl integer indices specifying additional observed values
#'   to plot
#' @param curve_data_list list of vectors where each vector is a predicted
#'   curve.  List length matches the number of methods used to make predictions.
#'
#' @return a vector of scores with scores grouped by method.  If there are two
#'   methods and three scores, the output would have the form (M1S1, M1S2, M1S3,
#'   M2S1, M2S2, M2S3).
#' @export
#'
compute_mixpred_scores <- function(mix_df,
                                   mix_idx,
                                   unplotted_repl,
                                   curve_data_list) {
  n_x_values <- ncol(Cx)
  score_matrix <- c()
  # for each set of curve data
  for (curve_data in curve_data_list) {
    llh <- 0
    mse <- 0
    crps <- 0
    # response index is set globally when running read_prepared_Tox21_data()
    resp_idx <- .GlobalEnv$resp_idx
    # for each replicate of data
    for (repl_idx in c(mix_idx, unplotted_repl)) {
      # get the replicate response
      resp_y_repl <- array(unlist(mix_df[repl_idx, resp_idx]))
      for (x_idx in 1:n_x_values) {
        # compute scores for given replicate
        llh <- llh + scoringRules::logs_sample(resp_y_repl[x_idx],
                                               curve_data[, x_idx])
        mse <- mse + (resp_y_repl[x_idx] - median(curve_data[, x_idx]))^2 /
          n_x_values
        crps <- crps + scoringRules::crps_sample(resp_y_repl[x_idx],
                                                 curve_data[, x_idx])
      }
    }
    score_matrix <- cbind(score_matrix, c(llh, mse, crps))
  }
  return(score_matrix)
}


# Not used; compute empirical LLH given bootstrapped curve estimates
get_emp_llh <- function(obs_res, boot_vals) {
  # get optimal bandwidth for KDE
  bw <- bw.SJ(boot_vals)
  # sum equivalent to a kernel density estimate at the point
  p_density <- sum(dnorm(obs_res,
                         mean = boot_vals,
                         sd = bw)) / length(boot_vals)
  return(log(p_density))
}

# Not used; compute CRPS given bootstraped estimates, empirical value
get_kernel_CRPS <- function(obs_res, boot_vals) {
  # get optimal bandwidth for KDE
  bw <- bw.SJ(boot_vals)
  #
  interval <- min(1, bw / 5)
  pdf_xvals <- seq(min(boot_vals) - 5 * bw,
                   max(boot_vals) + 5 * bw,
                   by = interval) # min(1, bw/5 ) )
  # sum equivalent to a kernel density estimate at the point
  pdf_kernel <- sapply(pdf_xvals,
                       FUN = function(x) {
                         sum(dnorm(x, mean = boot_vals, sd = bw)) /
                           length(boot_vals)
                       })
  # approximate CDF integral with summation
  cdf_kernel <- cumsum(pdf_kernel) * interval
  # approximate CRPS integral with summation
  CRPS_val <- sum((cdf_kernel - 1 * (pdf_xvals >= obs_res))^2) * interval
  return(CRPS_val)
}

# not used
get_empr_CRPS <- function(obs_res, boot_vals) {
  cdf_edf <- sapply(sort(boot_vals), FUN = function(x) mean(x > boot_vals))
  CRPS_val <- sum(((cdf_edf - 1 * (sort(boot_vals) >= obs_res))^2)[-1] *
                    diff(sort(boot_vals)))
  return(CRPS_val)
}

# Inverse Hill Functions ####

#' Inverse Hill Function Factory
#'
#' The inverse function is used with a single dose response curve, which is
#' assumed to be a Hill function.  A factory method is used to avoid repeatedly
#' passing in the curve parameters; once generated, the inverse calculator only
#' needs a response point y to return an inverse (dose) required to get that
#' response.
#'
#' @param a sill (max effect)
#' @param b EC50
#' @param c slope
#' @param max_R maximum effect across all chemicals
#' @param d minimum effect
#'
#' @return a function to compute the inverse given a mixture response R
#' @export
#'
#' @examples NA
hill_invs_factry <- function(a, b, c, max_R = 1, d = 0) {
  hilly_inverse <- function(y) {
    # input y: the response to invert
    # force: in case values change before function is used
    force(a)
    force(b)
    force(c)
    if (y == 0) {
      return(0)
    }
    # case 1, y extended to small negative conc, invert slope
    if ((a > 0 && y < 0) || (a < 0 && y > 0)) {
      return(-b / (1 + (-a / y)^(1 / c)))
    }
    # case 2, standard inverse
    if ((a > 0 && y < a) || (a < 0 && y > a)) {
      return(b / (a / y - 1)^(1 / c))
    }
    # case 3, reflected part of the standard inverse
    if ((a > 0 && y < 2 * a) || (a < 0 && y > 2 * a)) {
      return(-2 * b - b / (a / (2 * a - y) - 1)^(1 / c))
    }
    # case 4, reflection of extension, slope inverted
    return(-2 * b + b / (1 + (a / (-2 * a + y))^(1 / c)))
  }
  # the following method change was suggested by Insang Song to reduce
  # cyclomatic complexity
  hilly_inverse_low_cyclo <- function(y) {
    # input y: the response to invert
    # force: in case values change before function is used
    force(a)
    force(b)
    force(c)
    if (a * y == 0) {
      return(0)
    }
    # case 1, y extended to small negative conc, invert slope
    if (a * y < 0) {
      return(-b / (1 + (-a / y)^(1 / c)))
    }
    # case 2, standard inverse
    if (a * (y - a) < 0) {
      return(b / (a / y - 1)^(1 / c))
    }
    # case 3, reflected part of the standard inverse
    if (a * (y - 2 * a) < 0) {
      return(-2 * b - b / (a / (2 * a - y) - 1)^(1 / c))
    }
    # case 4, reflection of extension, slope inverted
    return(-2 * b + b / (1 + (a / (-2 * a + y))^(1 / c)))
  }
  return(hilly_inverse_low_cyclo)
}

#' Mixture Response Optimization Routine
#'
#' An internal factory function that returns a function with one input that can
#' be optimized with respect to that input. The returned function computes the
#' generalized concentration addition formula for the input set of chemicals and
#' input concentrations.
#'
#' @param hill_inverse_list a list of the hill inverse functions with one
#'   function per chemical
#' @param conc_vec a vector of non-negative values with the same length as the
#'   hill_inverse_list
#' @param synergy_const a scaling term for synergy. Difficult to specify so it
#'   is not used and currently set to 0, implying no synergy.
#' @param interval_sign 1 by default and can be set to -1.  Used for
#'   chemicals with a negative sill. Details: For numerical stability around
#'   very small concentrations and responses, the output function exponentiates
#'   the input response, making it impossible to invert a negative response. If
#'   at least one chemical in the mixture has a negative sill, this parameter
#'   should be set to -1 to check if the optimal (predicted) response is
#'   negative.
#'
#' @return a function GCA_over_list that has a response r as a input and a
#'   norm (measure of optimality according to GCA) as output
#' @export
#'
#' @examples
eff_response_opt <- function(hill_inverse_list,
                             conc_vec,
                             synergy_const = 0,
                             interval_sign = 1) {
  GCA_over_list <- function(r) {
    # get inverse for every component of cluster for current effect level r
    invs_list <- lapply(hill_inverse_list,
                        FUN = function(x) {
                          do.call(x, list(interval_sign * exp(r)))
                        })
    invs_vec <- unlist(invs_list)
    # avoid constant function when all 0
    if (all(conc_vec == 0)) {
      return(r)
    }
    # if using synergy: compute the GCA inverse, then normalize
    z_vec <- conc_vec / invs_vec
    z_vec <- z_vec / sum(z_vec)
    # by default synergy_scale = 1 (since synergy_const=0), implying no synergy
    synergy_scale <- exp(synergy_const * abs(prod(z_vec)))
    # the GCA equation balances at the optimal value:
    # [sum_i concentration_i / f_inverse(r) - 1 = 0]
    return(abs(sum(conc_vec / invs_vec) - synergy_scale))
  }
  return(GCA_over_list)
}

#' Mixture Response Calculator
#'
#' Factory method to return a function that predicts the mixture response given
#' the concentrations of the mixture component chemicals. A factory method is
#' used because the uncertainty quantification is based on a bootstrapped
#' collection of feasible parameters, and each set of parameters leads to a
#' different predictor for the mixture effect.  Rather than keeping track of
#' lists of parameters, we immediately convert the sampled parameter set to a
#' function (a "calculator") that takes the mixture concentration as input. The
#' list of calculators can be applied to a
#'
#' @param param_matrix a simple matrix object with rows corresponding to unique
#'   chemicals and columns representing parameters (a, b, c, max_R, d), where a
#'   is the sill b is the EC50 c is the slope value max_R is the maximum sill
#'   across all chemicals (deprecated) d is the minimum response (deprecated)
#' @param clust_assign a vector with integers defining cluster membership. For 3
#'   chemicals, a clust_assign could be c(1,2,1), meaning chemicals 1 and 3 are
#'   clustered together and 2 is by itself.  It is assumed that the cluster
#'   assignments start from 1 and do not skip integers.
#' @param get_counts boolean flag to return a count of the number of solutions
#'   found
#' @param scale_CA boolean flag to apply Concentration Addition assumption of
#'   equal sills
#' @param synergy_const a real number used to adjust the mixture prediction
#'   for synergistic effects
#'
#' @return an instance of the function mix_effect_fun which takes as input the
#'   concentration of the component chemicals and outputs a predict response
#' @export
#'
#' @examples
mix_function_generator <- function(param_matrix,
                                   clust_assign,
                                   get_counts = FALSE,
                                   scale_CA = FALSE,
                                   synergy_const = 0) {
  # mix_effect_fun takes as input the concentration of the component chemicals
  # and outputs a predict response
  mix_effect_fun <- function(true_conc_value) {
    CA_scaling_constant <- 1
    conc_value <- true_conc_value
    # NOT USED: adjust dose for effective dose, ligand competition
    clust_ids <- unique(clust_assign)
    response_per_cluster <- array(0, max(clust_ids))
    sol_per_cluster <- array(0, max(clust_ids))
    # for each group, get concentration addition part
    for (cid in 1:max(clust_ids)) {
      active_chems <- which(clust_assign == cid)
      active_params <- param_matrix[active_chems, ]
      active_conc <- conc_value[active_chems]
      if (all(active_conc == 0)) next
      # if a cluster only has 1 value, make sure to still get a matrix
      if (length(active_chems) == 1) active_params <- matrix(active_params,
                                                             nrow = 1)
      # if we are computing plain concentration addition, rescale sills to 1
      if (scale_CA) {
        CA_scaling_constant <- max(active_params[, 1][which(active_conc > 0)])
        # set all sills to be +-1
        active_params[, 1] <- sign(active_params[, 1])
      }
      # get list of inverse hill functions based on params using a factory
      hill_inverse_list <- apply(active_params,
                                 MARGIN = 1,
                                 function(x) {
                                   do.call(hill_invs_factry, as.list(x))
                                 })
      # assume concent addition and create function to find equivalent dose
      GCA_function <- eff_response_opt(hill_inverse_list,
                                       active_conc,
                                       synergy_const = synergy_const)
      # optim_res is the equivalent dose across all chems, note Log scale
      optim_res <- optimize(GCA_function,
                            interval = c(-20, 10),
                            tol = .Machine$double.eps)
      response_per_cluster[cid] <- exp(optim_res$minimum)
      # for the cases when there is a negative sill, check if the mix response
      # should be negative
      GCA_function_neg <- eff_response_opt(hill_inverse_list,
                                           active_conc,
                                           synergy_const = synergy_const,
                                           interval_sign = -1)
      optim_res_neg <- optimize(GCA_function_neg,
                                interval = c(-20, 10),
                                tol = .Machine$double.eps)
      # compare positive vs negative predicted response objective
      if (optim_res_neg$objective < optim_res$objective) {
        # if negative response provides the optimal value, return negative
        response_per_cluster[cid] <- -exp(optim_res_neg$minimum)
      }
      # count how many solutions were found
      n_sols <- 0
      if (optim_res_neg$objective < 5e-5) n_sols <- n_sols + 1
      if (optim_res$objective < 5e-5) n_sols <- n_sols + 1
      sol_per_cluster[cid] <- n_sols
    }
    # combine CA with IA
    r_max <- 1
    if (ncol(param_matrix) > 3) r_max <- max(param_matrix[, 4])
    total_response <- (1 - prod(1 - response_per_cluster / r_max))
    total_response <- CA_scaling_constant * r_max * total_response
    # synergy with IA not used, could be pnorm(qnorm(1-prod(...))-prod(syn))
    if (get_counts) {
      return(sol_per_cluster[1])
    }
    return(total_response)
  }
  return(mix_effect_fun)
}

#' Correcting Responses assuming Deactivation
#'
#' Given input concentrations, suppose some chemicals "deactivate" other
#' chemicals if they have a higher weight (eg molecular weight, docking score,
#' or affinity).  This function will adjust the true concentrations to reflect
#' how the higher weight chemicals will substitute the lower weight chemicals.
#' If the highest weight chemical has high enough concentration, it can dominate
#' and result in a mixture with only one effective chemical.  We can set rules
#' for which chemicals can substitute others using a matrix A_full.
#'
#' @param conc a vector of postive values representing concentrations (eg
#'   micromol/L)
#' @param R a vector of real-valued weights used in the objective; examples
#'   include docking scores or molecular weights
#' @param A_full a binary matrix.  if A_{ij}= 1, chemical i can substitute
#'   chemical j.  The diagonal must be 1.  If A is the Id matrix, the output
#'   adjusted concentrations will be equal to the input.
#'
#' @return a vector of adjusted concentrations, some of which may be 0
#' @export
#'
#' @examples
#' A_full <- matrix(1, nrow=3, ncol=3)
#' conc <- matrix(c(4,9, 7), nrow=3, ncol=1)
#' binding_affinity <- c(-2, -1, -0.5)
adjust_concentrations <- function(conc, R, A_full) {
  # scale conc for numerical accuracy?
  active_chems <- which(conc > 0)
  n_chem <- length(active_chems)
  if (n_chem == 1) return(conc)
  #normalize concentrations so that the smallest conc has value 1
  scale_factor <- 1 / min(conc[conc > 0])
  conc <- conc * scale_factor
  active_conc <- conc[active_chems]
  active_weight <- R[active_chems]
  # subset the substitution matrix for building constrained optimization
  A <- A_full[active_chems, active_chems]
  K <- CVXR::Variable(rows = n_chem, cols = n_chem)
  j_vec <- matrix(1, nrow = n_chem, ncol = 1)
  Jo <- matrix(1, nrow = n_chem, ncol = n_chem)
  diag(Jo) <- -1
  # total concentration conserved
  constraint_1 <- K %*% j_vec == active_conc
  # allocation bounded by diagonals
  constraint_2 <- (diag(t(K) %*% Jo) <= 0)
  # nonnegativity
  constraint_3 <- K >= 0
  # Ligand-Ligand (L-L) binding respected
  constraint_4 <- K * (1 - A) == 0
  probl <- CVXR::Problem(
    objective = CVXR::Minimize(sum(diag(active_weight) %*% t(K) %*% j_vec)),
    constraints = list(
      constraint_1,
      constraint_2,
      constraint_3,
      constraint_4
    )
  )
  result <- solve(probl)
  if (result$status == "solver_error") {
    warning("Adjusting Conc Failed")
    return(conc / scale_factor)
  }
  Kp <- result$getValue(K)
  new_conc <- rep(0, length(conc))
  new_conc[active_chems] <- base::diag(Kp)
  # machine precision issues, zero out small vals
  new_conc[new_conc < 1e-8] <- 0
  new_conc <- new_conc / scale_factor
  return(new_conc)
}
