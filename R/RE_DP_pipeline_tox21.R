# for profiling: Rprof(), code, Rprof(NULL) summaryRprof("Rprof.out")

#  Functions: ####
library(ggplot2)
library(reshape2)
library(scoringRules)
library(drc)



sessionInfo()
# investigated checkpoint package but found issues
# rather than rbind, use data.table funcction rbindlist
# MKL rathter than BLAS
# Pipeline start #####

source("R/helper_plots.R")
source("R/helper_calculators.R")
source("R/RE_dose_response_MCMC.R")
source("R/dirichlet_MCMC.R")
source("R/RE_nimble.R")

# read in data
source("R/tox21_prep_data.R")
# y_i = obs responses
# Cx = measured concentrations
# replicate_sets = indices for replicates of each compound

# Run MLE via drc package
run_pipe <- FALSE # avoid long build times
small_run <- TRUE

if (run_pipe) {
  input_df <- read.delim("input/AR-luc.txt")
  # mix guide describes doses of components in each mixture
  mix_guide <- readxl::read_xls("input/AllMixtureComponentsARER.xls")
  read_prepared_Tox21_data(input_df, mix_guide)
  # get dose response parameters from existing pacakge
  curve_fits <- get_mle_curve_fits(y_i, Cx, replicate_sets)
  # MCMC for RE model
  re_iter <- 2.5e4 # 2.5e4
  clust_iter <- 3.5e4 # 3.5e4
  if (small_run) {
    re_iter <- 2.5e3 # 2.5e4
    clust_iter <- 3.5e3 # 3.5e4
  }

  set.seed(102)
  # fit random effects model
  re_chains <- RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
  beepr::beep()
  re_par_list <- pull_parameters(re_chains)
  # if nimble is used, get re_par_list_nimble from pull_parameters_nimble(chain)
  re_par_summary <- pull_summary_parameters(re_chains, summry_stat = median)
  RE_curve_fits <- as.data.frame(list(
    "sill" = re_par_summary$sill_params,
    "ec50" = re_par_summary$ec50_params,
    "slope" = re_par_summary$slope_params
  ))

  # compare RE to MLE
  cbind(RE_curve_fits, curve_fits)
  #Can show the comparison plots using the plot_individual_response command with
  #params (replicate_sets, re_par_summary,curve_fits, RE_curve_fits, Cx, y_i)


  # diagnostic for MCMC ###
  run_diagnostic <- FALSE
  if (run_diagnostic) {
    set.seed(123)
    re_chains2 <- RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
    set.seed(124)
    re_chains3 <- RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
    all_dats <- get_MCMC_diagnostics(re_chains, re_chains2, re_chains3)
    # create a latex table with latex.tabular(as.tabular(all_dats))
  }


  # MCMC - GCA ####
  # fit curves with slope parameter fixed to 1
  set.seed(1025)
  re_chains_2param <- RE_MCMC_fit(y_i,
                                  Cx,
                                  replicate_sets,
                                  n_iter = re_iter,
                                  n_hill_par = 2)
  re_2par_list <- pull_parameters(re_chains_2param)
  re_2par_summary <- pull_summary_parameters(re_chains_2param)

  # fit DP clustering model
  set.seed(131)
  cluster_chain <- DP_MCMC_fit((re_par_list$sill_params[1:18]),
                               n_iter = clust_iter)
  beepr::beep()

  plot(cluster_chain$alpha_vals[seq(1000, clust_iter, by = 20)])

  # check what top assignments are
  n_top <- 20
  clust_centers_w_prob <- cluster_centers(cluster_chain, n_top = n_top)
  # visualize clusters
  visualize_clusters_blocks(re_par_list$sill_params,
                            clust_centers_w_prob,
                            ul = 400,
                            ll = -50)
  # alternatives:  random clustering, kmeans,
  # sample clusterings and pass into the mixture estimation function
  cluster_prob <- clust_centers_w_prob$probs
  centers <- clust_centers_w_prob$centers
  cent_sd <- clust_centers_w_prob$center_sd
  cluster_assign <- clust_centers_w_prob$assign
  DP_par_list <- list(
    "centers" = centers,
    "cluster_assign" = cluster_assign,
    "cent_sd" = cent_sd
  )
  # RGCA+DP on Sill ####
  RGCA_DP_par_list <- c(re_par_list, DP_par_list)
  cluster_prob_sills <- clust_centers_w_prob$probs
  RGCA_DP_sill_par_list <- c(
    re_par_list,
    list("cluster_assign" = cluster_assign)
  )
  # RGCA + random sampling
  # create 20 random clusters: 1 GCA, 1 IA, and then 4x2,5x3,5x4,4x5
  set.seed(1331)
  rand_clust_mat <- c(
    rep(1, 18),
    sample(1:2, 18 * 10, replace = TRUE),
    sample(1:3, 18 * 20, replace = TRUE),
    sample(1:4, 18 * 30, replace = TRUE),
    sample(1:5, 18 * 38, replace = TRUE),
    1:18
  ) |> matrix(ncol = 18, byrow = TRUE)
  rand_clust_assign <- rep(1, nrow(rand_clust_mat))
  names(rand_clust_assign) <-
    apply(rand_clust_mat,
          MARGIN = 1,
          FUN = function(rx) do.call(paste, as.list(rx)))
  randclust_par_list <- list("cluster_assign" = rand_clust_assign)
  RGCA_randclust_par_list <- c(re_par_list, randclust_par_list)
  # Nimble version: c(re_par_list_nimble, randclust_par_list)
  # RGCA with Kmeans
  kmeans_clust_mat <- matrix(0, nrow = 6, ncol = n_chems)
  for (i in 1:6) {
    kmeans_clust_mat[i, ] <- kmeans(RE_curve_fits, i)$cluster
    print(kmeans(RE_curve_fits, i)$betweenss / kmeans(RE_curve_fits, i)$totss)
  }
  kmeans_clust_assign <- rep(1, nrow(kmeans_clust_mat))
  names(kmeans_clust_assign) <-
    apply(kmeans_clust_mat,
          MARGIN = 1,
          FUN = function(rx) {
            do.call(paste, as.list(rx))
          })
  kmeans_clust_list <- list("cluster_assign" = kmeans_clust_assign)
  RGCA_kmeansclust_par_list <- c(re_par_list, kmeans_clust_list)



  # RGCA alone ####
  one_clust_assign <- 1
  names(one_clust_assign) <- do.call(paste, as.list(rep(1, n_chems)))
  RGCA_clust_par_list <- list(
    "centers" = matrix(re_par_summary$slope_params, nrow = 1),
    "cluster_assign" = one_clust_assign,
    "cent_sd" = matrix(rep(0, n_chems), nrow = 1)
  )
  # nimble version: c(re_par_list_nimble, RGCA_clust_par_list)
  RGCA_par_list <- c(re_par_list, RGCA_clust_par_list)




  # RGCA + cluster by agonist ####
  one_clust_assign <- 1
  clust_by_agonist <- rep(1, n_chems)
  clust_by_agonist[AR_agonist_rows] <- 2
  names(one_clust_assign) <- do.call(paste, as.list(clust_by_agonist))
  RGCA_ARER_clust_par_list <- list(
    "centers" = matrix(re_par_summary$slope_params, nrow = 1),
    "cluster_assign" = one_clust_assign,
    "cent_sd" = matrix(rep(0, n_chems), nrow = 1)
  )
  # Nimble version: re_par_list_nimble
  RGCA_ARER_par_list <- c(re_par_list, RGCA_ARER_clust_par_list)
  # Bayesian GCA and CA ####
  GCA_assign <- 1
  names(GCA_assign) <- do.call(paste, as.list(rep(1, n_chems)))
  GCA_clust_list <- list(
    "centers" = matrix(rep(1, n_chems), nrow = 1),
    "cluster_assign" = GCA_assign,
    "cent_sd" = matrix(rep(0, n_chems), nrow = 1)
  )
  GCA_par_list <- c(re_2par_list, GCA_clust_list)
  #  Bayesian IA ####
  # allow different slopes
  IA_assign <- 1
  names(IA_assign) <- do.call(paste, as.list(1:n_chems))
  DP_2parIA_list <- list(
    "centers" = matrix(re_par_summary$slope_params, nrow = 1),
    "cluster_assign" = IA_assign,
    "cent_sd" = matrix(rep(0, n_chems), nrow = 1)
  )
  IA_par_list <- c(re_par_list, DP_2parIA_list)
  beepr::beep(2)

  # Just Curves! ####
  # instead of a list of curves, there should just be 1
  # get parameter means or medians of MAPs

  #  IA single ####
  IA_assign <- 1
  names(IA_assign) <- do.call(paste, as.list(1:n_chems))
  IA_assign_vec <- as.numeric(strsplit(names(IA_assign), split = " ")[[1]])
  param_matrix_IA <- as.matrix(cbind(
    "a" = re_par_summary$sill_params,
    "b" = re_par_summary$ec50_params,
    "c" = re_par_summary$slope_params,
    "max_R" = max(re_par_summary$sill_params),
    "d" = 0
  ))
  IA_calculator <- mix_function_generator(param_matrix_IA, IA_assign_vec)


  # GCA single ####
  GCA_assign <- 1
  names(GCA_assign) <- do.call(paste, as.list(rep(1, n_chems)))
  GCA_assign_vec <- as.numeric(strsplit(names(GCA_assign), split = " ")[[1]])
  param_matrix_GCA <- as.matrix(cbind(
    "a" = re_2par_summary$sill_params,
    "b" = re_2par_summary$ec50_params,
    "c" = re_2par_summary$slope_params,
    "max_R" = max(re_par_summary$sill_params),
    "d" = 0
  ))
  GCA_calculator <- mix_function_generator(param_matrix_GCA, GCA_assign_vec)
  RGCA_single_calculator <- mix_function_generator(param_matrix_IA,
                                                   GCA_assign_vec)
  param_matrix_large_slope <- param_matrix_IA
  param_matrix_large_slope[, 3] <- param_matrix_large_slope[, 3] + .7
  RGCA_big_slope_calculator <- mix_function_generator(param_matrix_large_slope,
                                                      GCA_assign_vec)

  # CA single ####
  param_matrix_CA <- as.matrix(cbind(
    "a" = re_par_summary$sill_params,
    "b" = re_par_summary$ec50_params,
    "c" = re_par_summary$slope_params,
    "max_R" = 1,
    "d" = 0
  ))
  # simple CA calculator:  force all sill params to be equal
  CA_calculator <- mix_function_generator(param_matrix_CA,
                                          GCA_assign_vec,
                                          scale_CA = TRUE)
  # bootstraps ####
  # sample from posterior parameter estimates and clustering assigments
  set.seed(1026)
  n_bootstraps <- 100
  samp_idx <- sample(1:n_top, size = n_bootstraps,
                     prob = cluster_prob, replace = TRUE)
  sampled_mix_funs_RGCA_DP <-
    sapply(samp_idx,
           FUN = function(x) {
             create_mix_calc_clustered(x, RGCA_DP_par_list, add_RE = TRUE)
           })
  random_clustered_RGCA <-
    sapply(1:n_bootstraps,
           FUN = function(x) {
             create_mix_calc(x, RGCA_randclust_par_list, add_RE = TRUE)
           })
  # Nimble:  use RGCA_randclust_par_list_nimble
  kmeans_samp_idx <- sample(1:6, size = n_bootstraps, replace = TRUE)
  sampled_kmeans_clustered_RGCA <-
    sapply(kmeans_samp_idx,
           FUN = function(x) {
             create_mix_calc(x, RGCA_kmeansclust_par_list, add_RE = TRUE)
           })
  samp_idx <- sample(1:n_top, size = n_bootstraps,
                     prob = cluster_prob_sills,
                     replace = TRUE)
  sampled_sill_DPclustered_RGCA <-
    sapply(samp_idx,
           FUN = function(x) {
             create_mix_calc(x, RGCA_DP_sill_par_list, add_RE = TRUE)
           })
  # single cluster assignment versions
  sampled_mix_funs_RGCA <-
    sapply(rep(1, n_bootstraps),
           FUN = function(x) {
             create_mix_calc(x, RGCA_par_list, add_RE = TRUE)
           })
  # use nimble generated parameters: substitute RGCA_par_list_nimble
  sampled_mix_funs_RGCA_ARER <-
    sapply(rep(1, n_bootstraps),
           FUN = function(x) {
             create_mix_calc(x, RGCA_ARER_par_list, add_RE = TRUE)
           })
  sampled_mix_funs_GCA <-
    sapply(rep(1, n_bootstraps),
           FUN = function(x) {
             create_mix_calc(x, GCA_par_list, add_RE = TRUE)
           })
  sampled_mix_funs_IA <-
    sapply(rep(1, n_bootstraps),
           FUN = function(x) {
             create_mix_calc(x, IA_par_list, add_RE = TRUE)
           })
  # select which methods to use for prediction and plotting
  # nolint start
  bootstrap_calc_list <- list(
    "RGCA" = sampled_mix_funs_RGCA,
    # "RGCA_MOA" = sampled_mix_funs_RGCA_ARER,
    "RGCA_kMeans" = sampled_kmeans_clustered_RGCA,
    # "RGCA_DP" = sampled_sill_DPclustered_RGCA,
    "RGCA_mean" = list(RGCA_single_calculator),
    # "Nimble" = sampled_mix_funs_RGCA_nimble,
    # "Nimble_MOA" = sampled_mix_funs_RGCA_ARER_nimble,
    "Random_RGCA" = random_clustered_RGCA,
    # "Random_Nimble" = random_clustered_RGCA_nimble,
    # "RGCA_MOA_alt"= sampled_mix_funs_RGCA_ARER_alt,
    "GCA" = list(GCA_calculator),
    "IA" = list(IA_calculator),
    "CA" = list(CA_calculator)
  )
  # nolint end





  # verify 1-comp mix, sanity and bug check ####
  plot_dummy_mixture(Cx, y_i, re_par_summary, replicate_sets,
    bootstrap_calc_list,
    test_idx = 6
  )



  # Compare mix vs mix constituents
  plot_mix_vs_individuals(mix_idx = 34, ymax = 120)

  # ER/AR mix ####
  # compare to tox21 mixtures 34 109 159
  # good RE?  12, 14, 20  (18:  good for both!)
  # good GCA: 19

  # check 4 ad 12. # 1, 12, 20, 29, 43, 52, 55, 57
  set_1 <- c(1, 5, 10, 12, 20, 25, 29, 30, 31, 32, 43, 45, 50, 52, 55, 57, 62)
  set_4x <- c(1, 12, 20, 29, 43, 52, 55, 57)
  exclude_pure_ER <- c(1, 2, 3, 5, 8, 10, 12, 14, 16, 17, 18, 19, 20, 23, 24,
                       25, 26, 27, 29, 30, 31, 32, 33, 34, 35, 36, 37, 39, 43,
                       44, 45, 46, 47, 49, 50, 51, 52, 55, 56, 57, 61, 62, 66,
                       67, 68, 69)
  subindx <- which(set_1 %in% set_4x)
  set_1minus4x <- setdiff(set_1, set_4x)
  subset_exc <- setdiff(exclude_pure_ER, set_1)
  mix_idx <- 62 # 16, 8, 27!
  binary_mixes <- c(8, 18, 34, 37, 44, 49, 61)
  small_mixes <- c(6, 26, 36, 47, 59, 66, 67, 68)
  # ER subsets
  no_zear_ER <- c(3, 4, 5, 6, 8, 10, 12, 19, 22, 23, 24, 28, 31, 34, 36, 39, 40,
                  41, 43, 44, 45, 46, 47, 48, 53, 55, 56, 61, 63, 65, 67, 69)
  graded_no_zear_ER <- c(3, 6, 8, 22, 44, 47, 53, 63)
  graded_ER <- c(1, 2, 3, 6, 8, 11, 17, 22, 25, 26, 30, 35, 37, 38, 44, 47, 49,
                 51, 53, 54, 57, 58, 62, 63, 64, 66)
  # if saving: pdf(file = "Output/ARluc_5compare_%02d.pdf",onefile=TRUE, 9x5)
  score_matrix <- plot_mixture_response(
    mix_idx, mix_df, mix_conc_df, mix_guide,
    bootstrap_calc_list
  )
  beepr::beep(2)
  # Score results ####
  score_df <- data.frame(score_matrix)
  # drop rows of 0
  score_df <- score_df[apply(score_df,
                             MARGIN = 1,
                             FUN = function(rx) any(rx > 0)), ]
  names(score_df) <- c(
    "Mix id", paste(names(bootstrap_calc_list), c("LLH")),
    paste(names(bootstrap_calc_list), c("MSE")),
    paste(names(bootstrap_calc_list), c("CRPS"))
  )

  ## Violin Plot for Scores ####
  method_names <- c(
    "RGCA",
    "RGCA Sill",
    "RGCA K-Means",
    "RGCA Random",
    "GCA",
    "IA",
    "CA"
  )
  method_levels <- method_names[c(4:1, 5, 7, 6)]
  names(bootstrap_calc_list) <- method_names
  # if saving: pdf(file = "Output/boxplot_slope_invert_set1.pdf", onefile= TRUE,
  # width=8, height = 4)
  plot_scores(score_df,
              bootstrap_calc_list,
              method_levels = method_levels)
  # if saving dev.off()
  mix_descriptions <- sapply(score_df$`Mix id`, FUN = function(x) {
    mix_guide$Description[which(mix_guide$CAS == mix_df$CAS[x])]
  })
  mix_CAS <- sapply(score_df$`Mix id`, FUN = function(x) {
    mix_guide$CAS[which(mix_guide$CAS == mix_df$CAS[x])]
  })
  mix_descriptions <- stringr::str_replace(mix_descriptions, "200%", "4x EC50")
  score_df <- cbind("Mix Desc" = mix_descriptions, score_df)
  score_df <- cbind("Mix CAS" = mix_CAS, score_df)
}
