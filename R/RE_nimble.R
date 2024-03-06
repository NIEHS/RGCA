library(nimble, warn.conflicts = TRUE)
################################################################
# -----------  nimble: mixed effect -------------
################################################################
# -----------  samples + replicates  -------------
#' Nimble MCMC script to fit the manuscript model to the Tox21 data.
#'
#' Fit a random effect dose response model using Nimble for fast MCMC.  Nimble
#' creates a compiled sampler that can run iterations much faster than the
#' manually implemented version.
#'
#' @return nimble_samples: a list of MCMC chains with an entry for each
#'   parameter of the model.
#' @export
run_RE_nimble <- function() {
  code <- nimbleCode({
    for (chm in 1:CM) {
      a1[chm] ~ dnorm(0, sd = sqrt(1000))
      theta1[chm] ~ dgamma(shape = 1e-3, rate = 1e-3)
      beta1[chm] ~ dgamma(shape = .3, rate = .3)
      # give each chemical a noise term
      sigma_eps[chm] ~ dinvgamma(shape = 0.1, rate = 0.1)
      sigma_u[chm] ~ dinvgamma(shape = 0.1, rate = 0.1)
      sigma_v[chm] ~ dinvgamma(shape = 0.1, rate = 0.1)
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
  replicate_matrix <- matrix(unlist(lapply(replicate_sets, function(x) {
    a <- rep(0, 6)
    a[1:length(x)] <- x
    return(a)
  })), nrow = 6)
  n_reps <- unlist(lapply(replicate_sets, length))
  # for random effects, fix first replicate to 0 my setting variance to 0
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
    identifiability_constr = identifiability_constr
  )

  data <- list(
    y = y_i,
    x = Cx
  )

  inits <- list(
    beta1 = rep(1, nchm),
    theta1 = rep(1e-6, nchm),
    a1 = rep(10, nchm),
    sigma_eps = rep(1, nchm),
    sigma_u = rep(1, nchm),
    sigma_v = rep(1, nchm),
    u = rep(0, nrow(y_i)),
    v = rep(0, nrow(y_i))
  )
  drcModel <- nimbleModel(
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
  c_drcmodel <- compileNimble(drcModel)
  mcmc <- buildMCMC(c_drcmodel, monitors = parentNodes_drc)
  cmcmc <- compileNimble(mcmc, project = drcModel)
  nimble_samples <- runMCMC(cmcmc, niter = 50000, nburnin = 10000, thin = 20)
  return(nimble_samples)
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
  dpModel <- nimbleModel(
    code = codeDP, name = "DPMM", constants = constants,
    data = data, inits = inits
  )
  # slow, not always needed: CdpModel <- compileNimble(dpModel)
  dpMCMC <- buildMCMC(dpModel,
    monitors = c("nu", "ki", "sigsqrd")
  )
  CdpMCMC <- compileNimble(dpMCMC, project = dpModel)
  niter <- 5e4
  nburnin <- 1e4
  dpMCMC_samples <- runMCMC(CdpMCMC, niter = niter, nburnin = nburnin)
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
  samples <- nimbleMCMC(
    code = codeDP, data = data, inits = inits,
    constants = constants,
    monitors = c("nu", "ki", "sigsqrd", "alpha"),
    thin = 5, niter = 5000, nburnin = 2000, nchains = 1,
    setSeed = TRUE
  )
}
