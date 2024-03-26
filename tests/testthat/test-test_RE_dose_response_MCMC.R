test_that("Random Effect MCMC works", {
  # lots of setup:  generate some data
  set.seed(112)
  # sample some coefficients: theta_1, theta_2, phi_c
  n_chems <- 10
  n_replicates <- rep(3, n_chems)
  # dosages observed: generate curves at Cx, add noise
  Cx <- seq(0, 20, by = .8)
  max_dose <- 20
  n_dose <- length(Cx)
  # assume 3-4 clusters
  sill_true <- rnorm(n_chems, mean = 1, sd = .3)
  EC50_true <- rgamma(n_chems, shape = 3, rate = 1)
  phi_c_true_vals <- c(2, 1, 3)
  clust_pattern <- c(2, 5, 3)
  true_clust_assign <- rep(1:length(clust_pattern), clust_pattern)
  slope_true <- rep(phi_c_true_vals, clust_pattern) # slope
  # random effects vary, allow for some larger variances
  u_RE_sd_true <- c(rep(0, 10), rep(rep(c(.1, .3), c(5, 5)), 2))
  # true RE: precompute
  u_RE_true <- rnorm(sum(n_replicates), sd = u_RE_sd_true)
  v_RE_sd_true <- c(rep(0, 10), rep(rep(c(.2, .4), c(5, 5)), 2))
  v_RE_true <- rnorm(sum(n_replicates), sd = u_RE_sd_true)
  # now generate curves for each chem
  y_i <- matrix(nrow = sum(n_replicates), ncol = n_dose + 1)
  for (i in 1:n_chems) {
    for (r in 1:n_replicates[n_chems]) {
      curr_idx <- i + (r - 1) * n_chems
      response_data <- hill_function(
        sill_true[i] + u_RE_true[curr_idx], EC50_true[i],
        slope_true[i], Cx
      ) + v_RE_true[curr_idx] + rnorm(n_dose, sd = .05)
      # save the chem idx for next step
      y_i[curr_idx, ] <- c(i, response_data)
    }
  }
  # visualize generated data
  if(FALSE){
    par(mfrow = c(3, 3))
    for (i in 1:9) {
      plot(y_i[i, -1])
      points(y_i[i + n_chems, -1], pch = 2)
      points(y_i[i + n_chems + 1, -1], pch = 3)
    }
    par(mfrow = c(1, 1))
  }
  # replicates require chem ID
  replicate_sets <- c()
  for (chem in unique(y_i[, 1])) {
    rep_idx <- which(y_i[, 1] == chem)
    replicate_sets <- c(replicate_sets, list(rep_idx))
  }
  # drop the ID col, have replicate set map
  y_i <- y_i[, -1]
  # create matrix to allow for different levels of each chem
  Cx <- do.call(rbind, rep(list(Cx), sum(n_replicates)))
  # now have obs and doses, can run MCMC.  Large list returned, use snapshot
  expect_snapshot(RE_MCMC_fit(y_i,
                              Cx,
                              replicate_sets,
                              n_iter = 100,
                              n_hill_par = 3))
  # for testing, MCMC run is too short to get good estimates.
  # pull_summary_parameters(re_chain_list)
})
