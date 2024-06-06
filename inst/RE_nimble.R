# global variable definition for CRAN: vars appear in Nimble model
globalVariables(c("CM", "a1", "u", "theta1", "x", "beta1", "v"))

################################################################
# -----------  nimble: mixed effect -------------
################################################################
#' Nimble MCMC script to fit a random effect Hill model.
#'
#' Fit a random effect dose response model using Nimble for fast MCMC.  Nimble
#' creates a compiled sampler that can run iterations much faster than the
#' manually implemented version.
#'
#' @param y_i Response data as a matrix where each row is a chemical and each
#'   column is a dose.
#' @param Cx Dose data with dimension equal to y_i.
#' @param replicate_sets A list of arrays where each array indicates the rows of
#'   y_i that correspond to replicates of a single chemical
#' @param n_iter Number of MCMC iterations to run, defaults to 50k
#' @param n_burn Number of initial MCMC iterations to discard, defaults to 10k
#'   or half of n_iter, whichever is lower.
#' @param n_thin Number of MCMC iterations to skip to approximate iid samples,
#'   defaults to 20.  If it is larger than difference in iterations and burn in,
#'   it is reset to 1.
#' @param prior_alpha_sd hyperparameter for the sill, controls prior variance
#' @param prior_theta_shape hyperparameter for the EC50, controls gamma shape
#' @param prior_theta_rate hyperparameter for the EC50, controls gamma rate
#' @param prior_beta_shape hyperparameter for the slope, controls gamma shape
#' @param prior_beta_rate hyperparameter for the slope, controls gamma rate
#' @param prior_sd_noninfo hyperparameter for various noise terms which use a
#'   weakly informative prior, controls spread from 1.  Default is 0.1.
#'
#' @return nimble_samples: a list of MCMC chains with an entry for each
#'   parameter of the model.
run_RE_nimble <- function(y_i, Cx, replicate_sets,
                          n_iter = 5e4, n_burn = 1e4, n_thin = 20,
                          prior_alpha_sd = sqrt(1000),
                          prior_theta_shape = 1e-3,
                          prior_theta_rate = 1e-3,
                          prior_beta_shape = 0.3,
                          prior_beta_rate = 0.3,
                          prior_sd_noninfo = 0.1) {
  if (n_burn >= as.integer(n_iter / 2)) n_burn <- as.integer(n_iter / 2)
  if (n_thin >= n_iter - n_burn) n_thin <- 1
  code <- nimble::nimbleCode({
    for (chm in 1:CM) {
      a1[chm] ~ dnorm(0, sd = prior_alpha_sd)
      theta1[chm] ~ dgamma(shape = prior_theta_shape, rate = prior_theta_rate)
      beta1[chm] ~ dgamma(shape = prior_beta_shape, rate = prior_beta_rate)
      # give each chemical a noise term
      sigma_eps[chm] ~ dinvgamma(shape = prior_sd_noninfo,
                                 rate = prior_sd_noninfo)
      sigma_u[chm] ~ dinvgamma(shape = prior_sd_noninfo,
                               rate = prior_sd_noninfo)
      sigma_v[chm] ~ dinvgamma(shape = prior_sd_noninfo,
                               rate = prior_sd_noninfo)
      for (j in 1:n_reps[chm]) {
        # j is the relative index for replicates
        # replicate_matrix[j, chm]: gets the correct index for the replicate
        u[replicate_matrix[j, chm]] ~
          dnorm(mean = 0, sd = sigma_u[chm] * identifiability_constr[j, chm])
        v[replicate_matrix[j, chm]] ~
          dnorm(mean = 0, sd = sigma_v[chm] * identifiability_constr[j, chm])
        for (i in 1:y_lengths[replicate_matrix[j, chm]]) {
          f[j, i, chm] <- (a1[chm] + u[replicate_matrix[j, chm]]) /
            (1 + (theta1[chm] / x[replicate_matrix[j, chm], i])^beta1[chm]) +
            v[replicate_matrix[j, chm]]
          y[replicate_matrix[j, chm], i] ~ dnorm(f[j, i, chm],
            sd = sigma_eps[chm]
          )
        }
      }
    }
  })

  ## constants, data, and initial values
  nchm <- length(replicate_sets)
  n_reps <- unlist(lapply(replicate_sets, length))
  max_repls <- max(n_reps)
  replicate_matrix <- matrix(unlist(lapply(replicate_sets, function(x) {
    a <- rep(0, max_repls)
    a[1:length(x)] <- x
    return(a)
  })), nrow = max_repls)
  # for random effects, fix first replicate to 0 by setting variance to 0
  identifiability_constr <- replicate_matrix * 0 + 1
  # 0 doesnt work, use small value instead
  identifiability_constr[1, ] <- 1e-6
  y_lengths <- apply(y_i, MARGIN = 1,
                     FUN = function(x) {
                       ifelse(any(is.na(x)), min(which(is.na(x))) - 1, 15)
                     })
  constants <- list(
    y_lengths = y_lengths,
    replicate_matrix = replicate_matrix,
    n_reps = n_reps,
    CM = ncol(replicate_matrix),
    identifiability_constr = identifiability_constr,
    prior_theta_shape = prior_theta_shape,
    prior_theta_rate = prior_theta_rate,
    prior_beta_shape = prior_beta_shape,
    prior_beta_rate = prior_beta_rate,
    prior_sd_noninfo = prior_sd_noninfo,
    prior_alpha_sd = prior_alpha_sd
  )

  data <- list(
    y = y_i,
    x = Cx
  )

  inits <- list(
    beta1 = rep(1, nchm),
    theta1 = rep(1, nchm),
    a1 = rep(1, nchm),
    sigma_eps = rep(1, nchm),
    sigma_u = rep(1, nchm),
    sigma_v = rep(1, nchm),
    u = rep(0, nrow(y_i)),
    v = rep(0, nrow(y_i))
  )
  drcModel <- nimble::nimbleModel(
    code = code, constants = constants, data = data,
    inits = inits, check = FALSE
  )

  ## Ensure we have the nodes needed to simulate new datasets
  dataNodes_drc <- drcModel$getNodeNames(dataOnly = TRUE)
  parentNodes_drc <- drcModel$getParents(dataNodes_drc,
                                         stochOnly = TRUE,
                                         upstream = TRUE)
  ## Ensure we have both data nodes and deterministic intermediates (e.g.,
  ## lifted nodes)
  simNodes_drc <- drcModel$getDependencies(parentNodes_drc, self = FALSE)
  c_drcmodel <- nimble::compileNimble(drcModel)
  mcmc <- nimble::buildMCMC(c_drcmodel, monitors = parentNodes_drc)
  cmcmc <- nimble::compileNimble(mcmc, project = drcModel)
  nimble_samples <- nimble::runMCMC(cmcmc, niter = n_iter,
                                    nburnin = n_burn, thin = n_thin)
  return(nimble_samples)
}


#' Organize output of the Nimble MCMC
#'
#' Fitting a random effect model to data using a Bayesian method like MCMC
#' creates posterior samples.  Our implementation with Nimble creates a large
#' output object that contains posterior samples after burn-in and thinning.
#' This function organizes the samples for convenient downstream usage.
#'
#' @param nimble_samples the object returned by the run_RE_nimble function found
#'   in the inst folder.
#' @param summry_stat The function used to compute the summarizing statistic for
#'   the slope parameter.  Default function is median.
#' @param input_replicates The optional list of indices of replicates for each
#'   chemical.  By default, input_replicates is NA and the method assumes there
#'   exists a global var called replicate_sets
#'
#' @return A named list of arrays where each array is a posterior thinned sample
#'   from the nimble chain.
pull_parameters_nimble <- function(nimble_samples,
                                   summry_stat = stats::median,
                                   input_replicates = NA) {
  if (is.list(input_replicates)) replicate_sets <- input_replicates
  ec50_idx <- grep("theta", colnames(nimble_samples))
  ec50_sample <- nimble_samples[, ec50_idx]
  n_chems <- ncol(ec50_sample)
  sill_sample <- nimble_samples[, 1:n_chems]
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
                          stats::median(u_RE_samples[, replicate_sets[[ridx]]])
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



if (FALSE) {
  # -----------  DP with slopes  -------------
  codeDP <- nimbleCode({
    sigsqrd ~ dinvgamma(shape = 10, rate = .1) # var for each cluster xi[i]
    for (i in 1:nchm) {
      beta1[i] ~ dnorm(nu[i], var = sigsqrd) # slopes=samples from mixture dist.
      nu[i] <- nuTilde[ki[i]] # mean for each cluster xi[i]
    }
    # mixture component parameters drawn from base measures
    for (i in 1:nchm) {
      nuTilde[i] ~ dnorm(0, 1)
    }
    # CRP for clustering studies to mixture components
    ki[1:nchm] ~ dCRP(alpha, size = nchm)
    # hyperparameters
    alpha ~ dgamma(shape = 3 / 2, rate = 1 / 6)
  })
  # if using previous RE model nimble MCMC result
  beta1 <- scale(slopes)[1:18]
  nchm <- length(slopes)
  inits <- list(
    sigsqrd = 1,
    ki = rep(1, nchm),
    alpha = 1,
    nuTilde = rnorm(nchm)
  )
  data <- list(beta1 = beta1)
  constants <- list(nchm = length(beta1))
  dpModel <- nimble::nimbleModel(
    code = codeDP, name = "DPMM", constants = constants,
    data = data, inits = inits
  )
  # slow, not always needed: CdpModel <- compileNimble(dpModel)
  dpMCMC <- nimble::buildMCMC(dpModel,
    monitors = c("nu", "ki", "sigsqrd")
  )
  CdpMCMC <- nimble::compileNimble(dpMCMC, project = dpModel)
  niter <- 5e4
  nburnin <- 1e4
  dpMCMC_samples <- nimble::runMCMC(CdpMCMC, niter = niter, nburnin = nburnin)
  cluster_assign <- dpMCMC_samples[, grep("ki", colnames(dpMCMC_samples))]
  cluster_value <- dpMCMC_samples[, grep("nu", colnames(dpMCMC_samples))]
  for (i in 1:(niter - nburnin)) {
    # relabel for uniqueness
    c_i <- cluster_assign[i, ]
    c_i <- c_i - min(c_i) + 1
    uniq_labs <- unique(c_i)
    relabeled <- sapply(c_i, function(x) which(uniq_labs == x)[1])
    cluster_assign[i, ] <- relabeled
  }
  cluster_chain_nimble <- list(
    "means" = as.matrix(cluster_value),
    "cluster" = as.matrix(cluster_assign)
  )
  cluster_centers(cluster_chain_nimble, n_top = 50)
  # as one line
  samples <- nimble::nimbleMCMC(
    code = codeDP, data = data, inits = inits,
    constants = constants,
    monitors = c("nu", "ki", "sigsqrd", "alpha"),
    thin = 5, niter = 5000, nburnin = 2000, nchains = 1,
    setSeed = TRUE
  )
}
