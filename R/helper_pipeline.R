#  Functions: ####
library(ggplot2)
library(cowplot)
library(reshape2)
library(scoringRules)
library(tables) # for converting data frames to latex tabulars

generate_mix_estimators <- function(responses, doses,
                                    replicate_sets, re_iter = 2.5e4,
                                    clust_iter = 3.5e4, n_top_clust = 20,
                                    n_estimators = 100) {
  # source helper functions and MCMC code
  source("R/helper_plots.R")
  source("R/helper_calculators.R")
  source("R/RE_dose_response_MCMC.R")
  source("R/dirichlet_MCMC.R")
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

format_output_latex <- function() {
  # Tabular stuff for latex ####
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
  latex.tabular(as.tabular(cbind(RE_curve_fits, DRC_curve_fits)))
  # Mixture and CAS descriptions
  as.table(score_df$`Mix Desc`, row.names = rownames(score_df))
  latex.tabular(as.tabular(matrix(score_df$`Mix Desc`, ncol = 1)))
  latex.tabular(as.tabular(score_df[set_1, c(1, 2)]))
  latex.tabular(as.tabular(score_df[, 2:ncol(score_df)]))
  # simple table for cluster centers, tbc
  latex.tabular(as.tabular(as.matrix(clust_centers_w_prob$assign)))
  names(score_df) <- c(
    "Mix id", names(bootstrap_calc_list),
    names(bootstrap_calc_list),
    names(bootstrap_calc_list)
  )
  score_df
  latex.tabular(as.tabular(score_df[set_1, c(2:8)])) # 1:8, 15:22
}
