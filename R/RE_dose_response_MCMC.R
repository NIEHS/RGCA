################################################################################
# -------  MH/Gibbs sampler for the mixed effect model: ------------------------
# ------------------------------------------------------------------------------
# nolint start
# y = (theta1 + u)/(1 + (theta2/x)^phi_c) + v + epsilon
# epsilon ~ N(0, sigma^2)
# theta1 ~ N(0, 1000)
# theta2 ~ Gamma(1e-3, 1e-3)
# phi_c ~ Gamma(0.3, 0.3)
# u ~ N(0, sigma_u^2)
# v ~ N(0, sigma_v^2)
# sigma^2 ~ InvGamma(0.1, 0.1) # indexed for each chemical, not global!
# sigma_u^2 ~ InvGamma(0.1, 0.1)
# sigma_v^2 ~ InvGamma(0.1, 0.1)
# nolint end
################################################################################


# code to fit a tox mixture model given individual dose response curve
# models clustering through the hill slope parameter


# boolean values deactive chunks of code when sourcing the file
fit_ind_curves <- FALSE
read_in_data <- FALSE

#### Read in data ####
if (read_in_data) source("Desktop/tox_mixtures/code/tox21_prep_data.R")

#### MCMCM ####
# a function to make a list into a vector
mk_vec <- function(arow) unlist(array(arow))
# functions to take dot or matrix products when nulls present:
vprod_NA <- function(a, b) {
  mat_out <- (t(stats::na.omit(as.matrix(a))) %*% stats::na.omit(as.matrix(b)))
  if (dim(a)[2] > 1) {
    return(mat_out)
  }
  return(mat_out[1])
}
qprod_NA <- function(a, b, c) {
  return((t(stats::na.omit(as.matrix(a))) %*%
            diag(stats::na.omit(b)) %*%
            stats::na.omit(as.matrix(c)))[1])
}
oprod_NA <- function(a, b) {
  return(outer(stats::na.omit(a), stats::na.omit(b)))
}

#' Builds a design matrix for replicates of dose response data for performing a
#' simultaneous Gibbs update for the sill, sill random effect, and intercept
#' random effect
#'
#' @param v_list a list with a vector for each replicate. Each vector is a
#'   sequence of coefficients with length equal to the number of doses (ie
#'   samples), computed according to the Hill model.
#'
#' @return a matrix built from the vectors of the input to mimic a design matrix
#'   for a linear model
build_replicate_matrix <- function(v_list) {
  num_repls <- length(v_list)
  num_samples <- length(v_list[[1]])
  # tot dim: reps* samples X 1+2*(reps-1)
  V_mat <- matrix(0,
                  nrow = num_repls * num_samples,
                  ncol = 1 + 2 * (num_repls - 1))
  # sill part, no random effects (identifiability)
  V_mat[1:num_samples, 1] <- v_list[[1]]
  for (ri in 2:num_repls) {
    # careful indexing when constructing matrix
    r_start_idx <- 1 + (ri - 1) * num_samples
    r_end_idx <- r_start_idx + num_samples - 1
    r_idx <- r_start_idx:r_end_idx
    V_mat[r_idx, 1] <- v_list[[ri]] # sill part for current replicate
    V_mat[r_idx, ri] <- v_list[[ri]] # sill random effect part
    V_mat[r_idx, ri + num_repls - 1] <- 1 # intercept random effect coefficient
  }
  return(V_mat)
}

# Fit RE hill model to arbitrary dose response data.
# Input:
# n_c = number of chemicals (with replicates, may be redundant)
# n_d = number of doses
# y_i : the chemical responses as a (n_c, n_d) matrix
# Cx : the dosages for each chemical and response, a (n_c, n_d) matrix
# replicate_sets :  A list of replicate indices for each chemical.


#' Random Effect MCMC fitting
#'
#' Performs a Markov Chain Monte Carlo sampling procedure to create posterior
#' samples for the random effect model specified above.  Uses Gibbs updates when
#' possible and reverts to Metropolis-Hastings with Gaussian random walk
#' proposals where needed.  The EC50 parameter varies over multiple orders of
#' magnitude and involves a log normal proposal distribution.
#'
#'
#' @param y_i a matrix of dose responses for individual chemicals.  Rows are
#'   chemicals where each replicate has a separate row, columns are the dose,
#'   and entries are the response.
#' @param Cx a matrix of the doses given for individual chemicals.  Rows are
#'   chemicals where each replicate has a separate row, columns are the index,
#'   and the entry is the dose.  Should match y_i
#' @param replicate_sets a list of vectors where each vector has the row index
#'   of all replicates of a particular chemical.  The length of the list should
#'   match the number of unique chemicals.
#' @param n_iter the number of iterations, defaults to 10,000
#' @param n_hill_par specifies if the full Hill model with 3 parameters is fit
#'   (default) or if a simplified model with 2 parameters (slope =1) is fit.
#'   Useful for comparing our method to standard GCA, whcih requires slope=1.
#'
#' @return a list with the full sampled chains for the parameters: slope (phi),
#'   sill+ec50 (theta1 and theta2), noise variance (sigma), random effects,
#'   random effect prior variances.
#'
#' @export
RE_MCMC_fit <- function(y_i, Cx, replicate_sets,
                        n_iter = 10000,
                        n_hill_par = 3) {
  # parameter id:
  n_chems <- length(replicate_sets)
  n_replicates <- unlist(lapply(replicate_sets, length))
  n_dose <- ncol(Cx)
  # no clustering here, just do RE
  c_i <- 1:n_chems
  phi_c <- rep(1, n_chems)
  theta_1 <- rep(1, n_chems)
  theta_2 <- rep(1e-6, n_chems)
  u_repl <- sapply(1:n_chems, FUN = function(x) list(rep(0, n_replicates[x])))
  v_repl <- u_repl
  sigma <- rep(1, n_chems) # noise for obs data
  sigma_mh <- .1 # proposal distr for MH
  sigma_u <- rep(10, length(u_repl)) # RE for sill
  sigma_v <- rep(10, length(v_repl)) # RE for offset
  phi_rate <- .3 # parameter for gamma prior on slope
  phi_shape <- .3 # 2nd parameter for gamma prior on slope
  ec50_rate <- 1e-3
  ec50_shape <- 1e-3
  sill_var <- 1000 # approximate non-informative prior for sill ~ N(0, var)
  num_sampling_iters <- 5
  record_assigned_mean <- matrix(nrow = n_iter, ncol = n_chems)
  record_assigned_clust <- matrix(nrow = n_iter, ncol = n_chems)
  record_sigma <- matrix(nrow = n_iter, ncol = n_chems)
  record_tau <- matrix(nrow = n_iter, ncol = 1)
  record_RE_u <- matrix(nrow = n_iter, ncol = sum(n_replicates))
  record_RE_v <- matrix(nrow = n_iter, ncol = sum(n_replicates))
  record_RE_u_sd <- matrix(nrow = n_iter, ncol = length(sigma_u))
  record_RE_v_sd <- matrix(nrow = n_iter, ncol = length(sigma_v))
  record_hill_params <- matrix(nrow = n_iter, ncol = 2 * n_chems)
  # if profiling, begin here: Rprof(line.profiling=TRUE)
  for (iter in 1:n_iter) {
    # update phi
    # If model is 2 parameter, fix phi to 1, do not update
    if (n_hill_par > 2) {
      for (j in 1:n_chems) {
        # consider multiple MH steps
        phi_curr <- phi_c[j]
        for (itern in 1:num_sampling_iters) {
          # MH:  use proposal distr of N(phi_curr, sigma_mh)
          prop_phi_cj <- truncnorm::rtruncnorm(1, a = .1,
                                               mean = phi_curr, sd = sigma_mh)
          # if symmetric proposal, just get alpha from ratio
          total_llh_curr <- stats::dgamma(phi_curr, shape = phi_shape,
                                          rate = phi_rate, log = TRUE) +
            log(truncnorm::dtruncnorm(prop_phi_cj, a = .1,
                                      mean = phi_curr, sd = sigma_mh))
          total_llh_prop <- stats::dgamma(prop_phi_cj, shape = phi_shape,
                                          rate = phi_rate, log = TRUE) +
            log(truncnorm::dtruncnorm(phi_curr, a = .1,
                                      mean = prop_phi_cj, sd = sigma_mh))
          ## MH: allow for negative slopes
          cluster_member_idx <- which(c_i == j)
          for (clm in cluster_member_idx) {
            for (rid in 1:n_replicates[clm]) {
              act_idx <- replicate_sets[[clm]][rid]
              u_RE <- u_repl[[clm]][rid]
              v_RE <- v_repl[[clm]][rid]
              # concatenate all relevant Conc, responses:  [y_1(C), y_2(C)], etc
              active_X <- Cx[act_idx, ]
              active_Y <- y_i[act_idx, ]
              # sample conditional on the y_i assigned to phi
              curr_denom <- 1 + (theta_2[j] / active_X)^phi_curr
              curr_mean <- (theta_1[j] + u_RE) / curr_denom + v_RE
              curr_sd_RE <- sigma[j]
              llh_curr <- sum(stats::dnorm(active_Y, mean = curr_mean,
                                           sd = curr_sd_RE, log = TRUE),
                              na.rm = TRUE)
              total_llh_curr <- total_llh_curr + llh_curr
              prop_denom <- 1 + (theta_2[j] / active_X)^prop_phi_cj
              prop_mean <- (theta_1[j] + u_RE) / (prop_denom) + v_RE
              prop_sd_RE <- sigma[j]
              llh_prop <- sum(stats::dnorm(active_Y, mean = prop_mean,
                                           sd = prop_sd_RE, log = TRUE),
                              na.rm = TRUE)
              total_llh_prop <- total_llh_prop + llh_prop
            }
          }
          accept_thresh <- exp(total_llh_prop - total_llh_curr)
          if (stats::runif(1) < accept_thresh) {
            phi_curr <- prop_phi_cj
          }
        }
        phi_c[j] <- phi_curr
      }
    }
    # update theta_1 and u v RE (va Gibbs): prior is N(0, 1000) and
    # N(0,sig_u^2), likelihood is (y-aV-uV)sigm^-2(y-av-uV)
    for (j in 1:n_chems) {
      act_idx <- replicate_sets[[j]]
      # concatenate all relevantresponses:  [y_1(C), y_2(C)], etc
      active_Y <- c(t(y_i[act_idx, ]))
      # compute the coefficients for the parameters of interest
      V_ri <- lapply(act_idx, FUN = function(x) {
        1 / (1 + (theta_2[j] / Cx[x, ])^phi_c[c_i[j]])
      })
      # build design matrix of coefficients for replicates
      V_mat <- build_replicate_matrix(V_ri)
      prec_RE <- 1 / sigma[j]^2
      lin_terms_post_prec <- vprod_NA(V_mat, V_mat) * prec_RE +
        diag(c(
          1 / sill_var,
          rep(sigma_u[j]^(-2), n_replicates[j] - 1),
          rep(sigma_v[j]^(-2), n_replicates[j] - 1)
        ))
      lin_terms_post_var <- solve(lin_terms_post_prec +
                                    1e-5 * diag(nrow(lin_terms_post_prec)))
      # XY, RE: [from (XX)^-1 XY] for post mean is post_var * [v*prec_RE*y +
      # prior_mean/prior var] XY_term = qprod_NA(theta1_coeff, prec_RE,
      # active_Y)
      XY_term <- vprod_NA(V_mat, active_Y) * prec_RE
      lin_terms_post_mean <- vprod_NA(lin_terms_post_var, XY_term)
      # sample from the post distribution
      lin_term_sample <- t(chol(lin_terms_post_var)) %*%
        stats::rnorm(length(lin_terms_post_mean)) + lin_terms_post_mean
      theta_1[j] <- lin_term_sample[1]
      fill_idx <- 2:n_replicates[j]
      u_repl[[j]][fill_idx] <- lin_term_sample[fill_idx]
      # dont update v_RE yet?
      v_repl[[j]][fill_idx] <- lin_term_sample[fill_idx + n_replicates[j] - 1]
    }
    # update theta_2 (via MH).
    for (j in 1:n_chems) {
      act_idx <- replicate_sets[[j]]
      # concatenate all relevant Conc, repsonses
      active_X <- c(t(Cx[act_idx, ]))
      active_Y <- c(t(y_i[act_idx, ]))
      u_RE_vec <- c(sapply(u_repl[[j]], FUN = function(x) rep(x, n_dose)))
      v_RE_vec <- c(sapply(v_repl[[j]], FUN = function(x) rep(x, n_dose)))
      # for each index, do a few MH steps instead of 1.  as an adaptive type
      # update, use the current value for sd sigmh2 = signif(theta_2[j],1)*5
      sigma_mh_theta_2 <- .1
      # for log scale: lower_bd = -12; upper_bd = -2
      for (itern in 1:num_sampling_iters) {
        log_theta2_prop <- stats::rnorm(1, mean = log(theta_2[j]),
                                        sd = sigma_mh_theta_2)
        theta2_prop <- exp(log_theta2_prop)
        numerator <- theta_1[j] + u_RE_vec
        curr_mean <- numerator /
          (1 + (theta_2[j] / active_X)^phi_c[c_i[j]]) + v_RE_vec
        llh_curr <- sum(stats::dnorm(active_Y, mean = curr_mean,
                                     sd = sigma[j], log = TRUE),
                        na.rm = TRUE) +
          log(stats::dlnorm(theta2_prop,
                            meanlog = log(theta_2[j]),
                            sdlog = sigma_mh_theta_2)) +
          stats::dgamma(theta_2[j], shape = ec50_shape, rate = ec50_rate,
                        log = TRUE)
        prop_mean <- numerator /
          (1 + (theta2_prop / active_X)^phi_c[c_i[j]]) + v_RE_vec
        llh_prop <- sum(stats::dnorm(active_Y, mean = prop_mean,
                                     sd = sigma[j], log = TRUE),
                        na.rm = TRUE) +
          log(stats::dlnorm(theta_2[j],
                            meanlog = log(theta2_prop),
                            sdlog = sigma_mh_theta_2)) +
          stats::dgamma(theta2_prop, shape = ec50_shape, rate = ec50_rate,
                        log = TRUE)
        accept_p <- min(exp(llh_prop - llh_curr), 1)
        if (stats::runif(1) < accept_p) theta_2[j] <- theta2_prop
      }
    }
    # update sigma: for each chem (vs across all data, originally)
    for (j in 1:n_chems) {
      tot_error <- 0
      act_idx <- replicate_sets[[j]]
      active_X <- c(t(Cx[act_idx, ]))
      active_Y <- c(t(y_i[act_idx, ]))
      u_RE_vec <- c(sapply(u_repl[[j]], FUN = function(x) rep(x, n_dose)))
      v_RE_vec <- c(sapply(v_repl[[j]], FUN = function(x) rep(x, n_dose)))
      numerator <- theta_1[j] + u_RE_vec
      denom <- 1 + (theta_2[j] / active_X)^phi_c[c_i[j]]
      sum_sqr_err <- sum((active_Y - numerator / denom - v_RE_vec)^2,
                         na.rm = TRUE)
      tot_error <- tot_error + sum_sqr_err
      sigma_sqr <- 1 / stats::rgamma(1,
                                     shape = .1 + n_dose * length(act_idx) / 2,
                                     rate = .1 + tot_error / 2)
      # if sigma is global, sum all errors else update each term within loop
      sigma[j] <- sqrt(sigma_sqr)
    }
    # Update RE sigma_u, sigma_v with Jeffrey, limit = IG(0,0)
    G_beta_prior_u <- .1
    G_alpha_prior_u <- .1
    G_beta_prior_v <- .1
    G_alpha_prior_v <- .1
    for (j in 1:n_chems) {
      act_idx <- replicate_sets[[j]]
      G_beta_post_u <- G_beta_prior_u + u_repl[[j]] %*% u_repl[[j]] / 2
      G_alpha_post_u <- G_alpha_prior_u + (n_replicates[j] - 1) / 2
      sigma_u[j] <- sqrt(1 / stats::rgamma(1, shape = G_alpha_post_u,
                                           rate = G_beta_post_u))
      G_beta_post_v <- G_beta_prior_v + v_repl[[j]] %*% v_repl[[j]] / 2
      G_alpha_post_v <- G_alpha_prior_v + (n_replicates[j] - 1) / 2
      sigma_v[j] <- sqrt(1 / stats::rgamma(1, shape = G_alpha_post_v,
                                           rate = G_beta_post_v))
    }
    record_RE_u[iter, ] <- mk_vec(u_repl)
    record_RE_v[iter, ] <- mk_vec(v_repl)
    record_RE_u_sd[iter, ] <- sigma_u
    record_RE_v_sd[iter, ] <- sigma_v
    record_sigma[iter, ] <- sigma
    record_assigned_mean[iter, ] <- phi_c[c_i]
    record_assigned_clust[iter, ] <- c_i
    record_hill_params[iter, ] <- c(theta_1, theta_2)
  }
  return(list(
    "slope_record" = record_assigned_mean,
    "sill_mideffect_record" = record_hill_params,
    "sigma" = record_sigma,
    "tau" = record_tau,
    "u_RE" = record_RE_u,
    "v_RE" = record_RE_v,
    "u_RE_sd" = record_RE_u_sd,
    "v_RE_sd" = record_RE_v_sd
  ))
}
#if profiling, end here Rprof(NULL)
# review profiling results summaryRprof("Rprof.out")
