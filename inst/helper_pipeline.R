# Imports ####
require(coda)
require(scoringRules)
require(tables)

#  Functions ####
# Not used; compute empirical LLH given bootstrapped curve estimates
get_emp_llh <- function(obs_res, boot_vals) {
  # get optimal bandwidth for KDE
  bw <- stats::bw.SJ(boot_vals)
  # sum equivalent to a kernel density estimate at the point
  p_density <- sum(stats::dnorm(obs_res,
                                mean = boot_vals,
                                sd = bw)) / length(boot_vals)
  return(log(p_density))
}

# Not used; compute CRPS given bootstraped estimates, empirical value
get_kernel_CRPS <- function(obs_res, boot_vals) {
  # get optimal bandwidth for KDE
  bw <- stats::bw.SJ(boot_vals)
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

generate_mix_estimators <- function(responses, doses,
                                    replicate_sets, re_iter = 2.5e4,
                                    clust_iter = 3.5e4, n_top_clust = 20,
                                    n_estimators = 100) {
  # source helper functions and MCMC code
  source("inst/manuscript_plots.R")
  source("inst/dirichlet_MCMC.R")
  y_i <- responses
  Cx <- doses
  # fit random effects model
  re_chains <- RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
  re_par_list <- pull_parameters(re_chains)
  # fit DP clustering model
  cluster_chain <- DP_MCMC_fit((re_par_list$slope_params[1:18]),
                               n_iter = clust_iter)
  # check what top assignments are
  clust_centers_w_prob <- cluster_centers(cluster_chain, n_top = n_top_clust)
  # visualize clusters
  visualize_clusters_blocks(re_par_list$slope_params, clust_centers_w_prob)
  # sample clusterings and pass into the mixture estimation function
  cluster_prob <- clust_centers_w_prob$probs
  centers <- clust_centers_w_prob$centers
  cent_sd <- clust_centers_w_prob$center_sd
  cluster_assign <- clust_centers_w_prob$assign
  DP_par_list <- list(
    "centers" = centers,
    "cluster_assign" = cluster_assign,
    "cent_sd" = cent_sd,
    "cluster_prob" = cluster_prob
  )
  tot_par_list <- c(re_par_list, DP_par_list)
  samp_idx <- sample(1:n_top_clust,
                     size = n_estimators,
                     prob = cluster_prob, replace = TRUE)
  sampled_mix_funs <- sapply(samp_idx,
                             FUN = function(x) {
                               create_mix_calc(x, tot_par_list, add_RE = TRUE)
                             })
  return(list(
    indv_chem_params = re_par_list,
    cluster_results = clust_centers_w_prob,
    sampled_mix_functions = sampled_mix_funs
  ))
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



#' Mixture Response Calculator Wrapper for Manuscript using Dirichlet Process
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
  sill_RE_boot <- stats::rnorm(n_chem_pars,
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



# Tabular stuff for latex ####
format_output_latex <- function() {
  pure_ordering <- sapply(pure_unique_CAS,
                          FUN = function(x) {
                            which(x == CAS_nums[relevant_guide_CASs])
                          })
  chem_map <- sapply(names(pure_ordering),
                     FUN = function(cas_num) {
                       pure_df$Sample.Name[which(pure_df$CAS == cas_num)[1]]
                     })
  # full RE_MCMC specification table
  RE_curve_fits <- as.data.frame(list(
    "Name" = chem_map,
    "Sill" = re_par_summary$sill_params,
    "EC50" = re_par_summary$ec50_params,
    "Sill_RE_sd" = re_par_summary$u_RE_sd_params,
    "Int_RE_sd" = re_par_summary$v_RE_sd_params,
    "Slope" = re_par_summary$slope_params
  ))
  RE_curve_fits <- RE_curve_fits %>%
    mutate_at(vars(Sill, Slope, Sill_RE_sd, Int_RE_sd),
              function(x) (round(x, 2))) %>%
    mutate_at(vars(EC50), function(x) (signif(x, digits = 3)))
  latex.tabular(as.tabular(RE_curve_fits))
  # MLE and RE-MCMC parameter table
  RE_curve_fits <- as.data.frame(list(
    "Name" = chem_map,
    "Sill" = re_par_summary$sill_params,
    "EC50" = re_par_summary$ec50_params,
    "Slope" = re_par_summary$slope_params
  ))
  RE_curve_fits <- RE_curve_fits %>%
    mutate_at(vars(Sill, Slope), function(x) (round(x, 2))) %>%
    mutate_at(vars(EC50), function(x) (signif(x, digits = 3)))
  DRC_curve_fits <- as.data.frame(curve_fits) %>%
    mutate_at(vars(V1, V3), function(x) (round(x, 2))) %>%
    mutate_at(vars(V2), function(x) (signif(x, digits = 3)))
  names(DRC_curve_fits) <- c("Sill (drc)", " EC50 (drc)", "Slope (drc)")
  cbind(RE_curve_fits, DRC_curve_fits)
  tables::latex.tabular(tables::as.tabular(cbind(RE_curve_fits,
                                                 DRC_curve_fits)))
  # Mixture and CAS descriptions
  as.table(score_df$`Mix Desc`, row.names = rownames(score_df))
  tables::latex.tabular(tables::as.tabular(matrix(score_df$`Mix Desc`,
                                                  ncol = 1)))
  tables::latex.tabular(tables::as.tabular(score_df[set_1, c(1, 2)]))
  tables::latex.tabular(tables::as.tabular(score_df[, 2:ncol(score_df)]))
  # simple table for cluster centers, tbc
  intmd_tabular <- tables::as.tabular(as.matrix(clust_centers_w_prob$assign))
  tables::latex.tabular(intmd_tabular)
  names(score_df) <- c(
    "Mix id", names(bootstrap_calc_list),
    names(bootstrap_calc_list),
    names(bootstrap_calc_list)
  )
  score_df
  latex.tabular(as.tabular(score_df[set_1, c(2:8)])) # 1:8, 15:22
}



# MCMC Diag ####
get_MCMC_diagnostics <- function(re_chains, re_chains2, re_chains3) {
  n_chems <- ncol(re_chains$slope_record)
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
  # use chem_map from tox21_prep_data : rownames(all_dats) <- chem_map
  return(all_dats)
}


# Synergy CVX ####

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
#'
#' @examples
#' A_full <- matrix(1, nrow=3, ncol=3)
#' conc <- matrix(c(4,9, 7), nrow=3, ncol=1)
#' binding_affinity <- c(-2, -1, -0.5)
adjust_concentrations <- function(conc, R, A_full) {
  require(CVXR)
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
