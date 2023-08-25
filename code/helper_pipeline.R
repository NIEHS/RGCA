#  Functions: ####
library(ggplot2)
library(cowplot)
library(reshape2)
library(scoringRules)
library(tables) # for converting data frames to latex tabulars


generate_mix_estimators = function(responses, doses,replicate_sets, re_iter = 2.5e4, 
                        clust_iter = 3.5e4, n_top_clust = 20, n_estimators = 100){
# source helper functions and MCMC code
  source("Code/helper_plots.R")
  source("Code/helper_calculators.R")
  source("Code/RE_dose_response_MCMC.R")
  source("Code/dirichlet_MCMC.R")

  y_i = responses
  Cx = doses
  # fit random effects model
  re_chains = RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
  re_par_list = pull_parameters(re_chains)
  re_par_summary = pull_summary_parameters(re_chains); #as.data.frame(re_par_summary)
  RE_curve_fits = as.data.frame(list("sill" = re_par_summary$sill_params, 
                                     "ec50" = re_par_summary$ec50_params,
                                     "slope" = re_par_summary$slope_params ))
  
  # fit DP clustering model
  cluster_chain = DP_MCMC_fit((re_par_list$slope_params[1:18]), n_iter = clust_iter)
  # check what top assignments are
  clust_centers_w_prob = cluster_centers(cluster_chain, n_top = n_top_clust)
  # visualize clusters
  #pdf("Output/cluster_vis.pdf", width = 8, height = 5)
  visualize_clusters_blocks(re_par_list$slope_params, clust_centers_w_prob)
  
  #sample clusterings and pass into the mixture estimation function
  cluster_prob = clust_centers_w_prob$probs
  centers = clust_centers_w_prob$centers
  cent_sd = clust_centers_w_prob$center_sd
  cluster_assign = clust_centers_w_prob$assign
  DP_par_list = list("centers" = centers, 
                     "cluster_assign" = cluster_assign,
                     "cent_sd" = cent_sd, 
                     "cluster_prob" = cluster_prob)
  tot_par_list = c(re_par_list, DP_par_list)
  
  
  samp_idx = sample(1:n_top_clust, size = n_estimators, prob =cluster_prob, replace = T)
  sampled_mix_funs = sapply(samp_idx, FUN = function(x) create_mix_calc(x, tot_par_list, add_RE = T))
  
  return(list(indv_chem_params = re_par_list,
              cluster_results = clust_centers_w_prob,
              sampled_mix_functions = sampled_mix_funs))
  
}