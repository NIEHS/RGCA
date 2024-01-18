################################################################################
# -------  Gibbs sampler provided in the Neal 1999 paper for the model: --------
#-------------------------------------------------------------------------------
# y|phi, c ~ N(phi_c, 1)
# Grouping c | p ~ Discrete(p_1,...,p_K)
# phi ~ G_0
# probabilities p ~ Dir(a/K,...,a/K)
################################################################################



# maintenance function to remove proposed slope values that aren't used
# and correct the indices so that there are no skipped integers
remove_unused_phi <- function(c_i, phi_c) {
  # collect current indices
  valid_phi <- 1:length(phi_c)
  valid_ci <- unique(c_i)
  # check which phi is no longer being used
  missing_idx <- setdiff(valid_phi, valid_ci)
  # if missing, drop and shift remaining
  if (length(missing_idx) > 0) {
    phi_c <- phi_c[-missing_idx]
    # decrement all higher indices by 1
    c_i[c_i > missing_idx] <- c_i[c_i > missing_idx] - 1
  }
  return(list(c_i, phi_c))
}



#' Dirichlet Process MCMC
#'
#' Performs clustering on the input data according to a Dirichlet Process
#' Mixture model.  An implementation of Neal 1999.
#'
#' @param y_in a vector of real numbers to be clustered
#' @param n_iter an integer for the number of iterations, defaults to 10,000
#' @param sigma_in a positive integer for the variances of the clusters.  A
#'   large value leads to fewer clusters.
#'
#' @return a list with MCMC samples for the cluster means, cluster assignments,
#'   concentration parameter alpha, and the cluster variance sigma
#' @export
#'
#' @examples
#' y_i <- c(-1.48, -1.4, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
#' DP_MCMC_fit(y_i)
DP_MCMC_fit <- function(y_in, n_iter = 10000, sigma_in = .1) {
  # standardize to avoid changing parameters:
  mu_input <- mean(y_in)
  sd_input <- sqrt(var(y_in))
  y_i <- (y_in - mu_input) / sd_input

  # init: assume 1 cluster, but generate a value of phi_c
  n <- length(y_i)
  c_i <- rep(1, n)
  # initialize group mean at data mean
  phi_c <- c(mean(y_i))
  m <- 3 # num aux
  alpha <- 1 # dirichlet prior concentration
  sigma <- sigma_in # noise variance
  lambda <- 1 # prior var phi_c
  record_assigned_mean <- matrix(nrow = n_iter, ncol = n)
  record_assigned_clust <- matrix(nrow = n_iter, ncol = n)
  record_sigma <- matrix(nrow = n_iter, ncol = 1)
  record_alpha <- matrix(nrow = n_iter, ncol = 1)
  for (iter in 1:n_iter) {
    # update c
    for (j in 1:n) {
      c_i_up <- c_i[j]
      remain_c_i <- c_i[-j]
      remaining_c_id <- sort(unique(remain_c_i))
      k_minus <- length(remaining_c_id)
      # generate aux
      phi_m <- rnorm(m, mean = 0, sd = 1) # vs mean = 1
      # reuse cluster id if current idx is lone member
      last_of_group <- !(c_i_up %in% remaining_c_id)
      if (last_of_group) {
        # phi[c_i]  already in set, but treated as new cluster
        phi_m <- phi_m[-1]
      }
      # sample cluster assignment according to likelihood of obs
      # same regardless of lone cluster
      n_phis <- k_minus + length(phi_m) + last_of_group
      potential_phs <- c(phi_c, phi_m)
      prob_c_vals <- rep(0, n_phis)
      for (k in 1:n_phis) {
        n_ic <- sum(remain_c_i == k)
        if (n_ic == 0) n_ic <- alpha / m
        prob_c_vals[k] <- n_ic / (n - 1 + alpha) *
          dnorm(y_i[j], mean = potential_phs[k], sd = sigma)
      }
      prob_c_vals <- prob_c_vals / sum(prob_c_vals)
      sampling_set <- 1:n_phis
      new_c <- sample(sampling_set, size = 1, prob = prob_c_vals)
      # if chose self, go to next iter
      if (new_c == c_i_up) next
      # otherwise: chose existing or new phi
      c_i[j] <- new_c
      if (new_c > (k_minus + last_of_group)) {
        c_i[j] <- k_minus + last_of_group + 1
        phi_c <- c(phi_c, potential_phs[new_c])
      }
      # might have lost a group, update phi
      up_list <- remove_unused_phi(c_i, phi_c)
      c_i <- up_list[[1]]
      phi_c <- up_list[[2]]
    }

    # update phi
    n_phi <- length(unique(phi_c))
    for (j in 1:n_phi) {
      # sample conditional on the y_i assigned to phi
      y_condition <- y_i[c_i == j]
      post_var <- sigma^2 / (length(y_condition) + lambda) # assume N(0,1) prior
      post_mean <- (sum(y_condition) + phi_c[j] * lambda) / sigma^2 * post_var
      new_phi_val <- rnorm(1, post_mean, sqrt(post_var))
      phi_c[j] <- new_phi_val
    }

    # update sigma
    p <- length(phi_c)
    sigma_sqr <- 1 / rgamma(1,
      shape = 10 + n / 2, # was 10, .1
      rate = 1 + sum((y_i - phi_c[c_i])^2) / 2
    )
    # mean of IG is b/(a-1), var is b/(a+1)
    sigma <- sqrt(sigma_sqr)
    record_sigma[iter] <- sigma

    # update the DP concentration alpha: large alpha = more clusters
    eta <- rbeta(n = 1, alpha + 1, n)
    # prior on alpha ~ G(a= 3/2, b= 1/2)
    alpha_a <- 3 / 2
    alpha_b <- 1 / 6
    mix_prob <- (alpha_a + k - 1) / (alpha_a + k - 1 + n * (alpha_b - log(eta)))
    # the posterior is a mixture; sample given eta
    alpha <- ifelse(runif(1) < mix_prob,
      rgamma(1, shape = alpha_a + p, rate = alpha_b - log(eta)),
      rgamma(1, shape = alpha_a + p - 1, rate = alpha_b - log(eta))
    )

    record_alpha[iter] <- alpha

    # relabel for uniqueness
    uniq_labs <- unique(c_i)
    c_i <- order(uniq_labs)[c_i]
    phi_c <- phi_c[uniq_labs]

    record_assigned_mean[iter, ] <- phi_c[c_i]
    record_assigned_clust[iter, ] <- c_i
  }

  # gibbs type update = no need to check acceptance rate

  # correct for initial standardization
  record_assigned_mean <- record_assigned_mean * sd_input + mu_input
  # record_sigma reflects Z_i rather than y_i, but is not used
  return(list(
    "means" = record_assigned_mean,
    "cluster" = record_assigned_clust,
    "alpha_vals" = record_alpha,
    "sigma_vals" = record_sigma
  ))
}


#' Get the cluster results
#'
#' Provides the top clusters and the corresponding centers given the MCMC chains
#' provided by [DP_MCMC_fit()]
#'
#' @param cluster_chain the list of MCMC samples returned by [DP_MCMC_fit()]
#' @param n_top integer number of clusters to select, starting with the most
#'   frequently occuring
#' @param plot_hist boolean to plot a histogram showing how often a cluster has
#'   k components or centers
#'
#' @return a list with the cluster centers, standard deviations, approximated
#'   probabilities, assignments, and the number of leftover or unused
#'   clusterings
#' @export
#'
#' @examples
#' y_i <- c(-1.48, -1.4, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
#' cluster_chain <- DP_MCMC_fit(y_i)
#' clust_centers_w_prob <- cluster_centers(cluster_chain, n_top = 20)
cluster_centers <- function(cluster_chain, n_top = 5, plot_hist = TRUE) {
  record_assgnd_mean <- cluster_chain$mean
  record_assigned_clust <- cluster_chain$cluster
  # check how many clusters there are per sample, for histogram
  n_clust_over_time <- apply(X = record_assigned_clust, MARGIN = 1, max)
  if (plot_hist) hist(n_clust_over_time)
  # Change cluster ID from vector to a string, for filtering commong clusters
  clust_arrange <- apply(
    X = record_assigned_clust,
    MARGIN = 1,
    FUN = function(x) paste(x, collapse = " ")
  )
  # for each element in cluster, get the mean for the (slope) parameters
  cluster_table <- sort(table(clust_arrange), decreasing = TRUE)
  clust_centers <- matrix(0, nrow = n_top, ncol = ncol(record_assigned_clust))
  clust_sds <- matrix(0, nrow = n_top, ncol = ncol(record_assigned_clust))
  for (cidx in 1:n_top) {
    clust_name <- names(cluster_table[cidx])
    clust_arrange[which(clust_arrange == clust_name)]
    cluster_samples <- record_assgnd_mean[which(clust_arrange == clust_name), ]
    # allow for a top cluster to occur once?
    if (is.null(dim(cluster_samples))) {
      clust_centers[cidx, ] <- cluster_samples
      clust_sds[cidx, ] <- sqrt(var(cluster_samples))
    } else {
      clust_centers[cidx, ] <- colMeans(cluster_samples)
      clust_sds[cidx, ] <- sqrt(diag(var(cluster_samples)))
    }
  }

  print(cluster_table[1:n_top])
  unused <- sum(cluster_table) - sum(cluster_table[1:n_top])
  cluster_prob <- array(cluster_table[1:n_top] / sum(cluster_table[1:n_top]))
  return(list(
    "centers" = clust_centers,
    "center_sd" = clust_sds,
    "probs" = cluster_prob,
    "assign" = cluster_table[1:n_top],
    "unused_samples" = unused
  ))
}



#' Convert output from dirichletprocess package
#'
#' This script converts the fitted model object from the
#' [dirichletprocess::DirichletProcessGaussian()] to match our own
#' [cluster_centers()] output.  This allows us to reuse our visualize functions
#' or use their code for the pipeline.
#'
#' @param dp_fit output from [dirichletprocess::DirichletProcessGaussian()]
#' @param n_top integer number of top clusters to select
#'
#' @return a list with the same attributes as the [cluster_centers()] output
#' @export
#'
#' @examples
convert_DP_pack_obj <- function(dp_fit, n_top = 5) {
  plot(dp_fit$alphaChain[seq(1, clust_iter, by = 10)])

  record_assigned_clust <- dp_fit$labelsChain
  # check how many clusters there are per sample, for histogram
  n_clust_over_time <- unlist(lapply(X = record_assigned_clust,
                                     MARGIN = 1,
                                     max))
  hist(n_clust_over_time)
  # Change cluster ID from vector to a string, for filtering commong clusters
  clust_arrange <- unlist(lapply(record_assigned_clust,
    FUN = function(x) {
      paste(unlist(x), collapse = " ")
    }
  ))
  # for each element in cluster, get the mean for the (slope) parameters
  cluster_table <- sort(table(clust_arrange), decreasing = TRUE)
  clust_centers <- matrix(0, nrow = n_top, ncol = length(dp_fit$data))
  clust_sds <- matrix(0, nrow = n_top, ncol = length(dp_fit$data))
  for (cidx in 1:n_top) {
    clust_name <- names(cluster_table[cidx])
    # check which index matches cluster, get all parameter vectors
    relevant_idx <- which(clust_arrange == clust_name)
    cluster_sample_mat <- rbind(sapply(relevant_idx,
      FUN = function(x) unlist(dp_fit$clusterParametersChain[[x]])
    ))

    # allow for a top cluster to occur once?
    if (is.null(dim(cluster_sample_mat))) {
      clust_centers[cidx, ] <- cluster_sample_mat
      clust_sds[cidx, ] <- sqrt(var(cluster_sample_mat))
    } else {
      clust_centers[cidx, ] <- colMeans(cluster_sample_mat)
      clust_sds[cidx, ] <- sqrt(diag(var(cluster_sample_mat)))
    }
  }

  print(cluster_table[1:n_top])
  unused <- sum(cluster_table) - sum(cluster_table[1:n_top])
  cluster_prob <- array(cluster_table[1:n_top] / sum(cluster_table[1:n_top]))
  return(list(
    "centers" = clust_centers,
    "center_sd" = clust_sds,
    "probs" = cluster_prob,
    "assign" = cluster_table[1:n_top],
    "unused_samples" = unused
  ))
}
