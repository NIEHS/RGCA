# for improving speed
#Rprof()
#Rprof(NULL)
#summaryRprof("Rprof.out")

#setwd("/Users/zilberds/Desktop/tox_mixtures/code/")

#  Functions: ####
library(ggplot2)
library(cowplot)
library(reshape2)
library(scoringRules)
library(tables) # for converting data frames to latex tabulars

sessionInfo()
checkpoint::checkpoint()
DT::datatable()
ggplotly: hover-over
rather than rbind, use data.table::rbindlist()
%dopar%, do_parallel(), foreach
MKL rathter than BLAS
# Pipeline start #####
#   fit linear RE model
#   perform dp clustering, save chain
#   bootstrap from chain and plot mixture effect, compare to concentration


source("Code/helper_plots.R")
source("Code/helper_calculators.R")
source("Code/RE_dose_response_MCMC.R")
source("Code/dirichlet_MCMC.R")

# read in data
source("Code/tox21_prep_data.R")
# y_i = obs responses
# Cx = measured concentrations
# replicate_sets = indices for replicates of each compound



# ER-luc problems:  chem #4, idx 14, 42, 47: identifiability issue

# PA-GCA ####
# TODO : check low-concentration bias
# TODO: log transforms for better EC50 estimatipo?
# TODO Add comment about 0 effect chemis:  Necessary!
# TODO Get equivalent data form DP package
# TODO check for pbetter prior for variance: non-informative
# TODO check if there is a burn in?  distribution of top cluster?  No burn in!
# DONE: Null model or estiamte a ratio of parameters that are not mixing well
# DONE  try 100000 iter?  1e6: 400k unique!  not scaling well, max clust was only 2k
# DONE check if slopes are significantly away from 0, Tukey HSD:  means are always sign diff, since sample is large; but 

# Run MLE via drc package
curve_fits = get_mle_curve_fits(y_i, Cx, replicate_sets)



# MCMC for RE model
re_iter = 2.5e4 # 2.5e4
clust_iter =3.5e4  #3.5e4
set.seed(102)
# fit random effects model
re_chains = RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter); beepr::beep()
re_par_list = pull_parameters(re_chains)
re_par_summary = pull_summary_parameters(re_chains); #as.data.frame(re_par_summary)
RE_curve_fits = as.data.frame(list("sill" = re_par_summary$sill_params, 
                                   "ec50" = re_par_summary$ec50_params,
                                   "slope" = re_par_summary$slope_params ))

#compare RE to MLE
cbind(RE_curve_fits, curve_fits)
# plot_individual_response(replicate_sets, re_par_summary,curve_fits, RE_curve_fits, Cx, y_i)


# diagnostic for MCMC ###
run_diagnostic=F
if(run_diagnostic){
  set.seed(123)
  re_chains2 = RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
  set.seed(124)
  re_chains3 = RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter)
  all_dats = get_MCMC_diagnostics(re_chains, re_chains2, re_chains3)
  #latex.tabular(as.tabular(all_dats))
}


# GCA ####
# fit curves with slope parameter fixed to 1
set.seed(1025)
re_chains_2param=RE_MCMC_fit(y_i, Cx, replicate_sets, n_iter = re_iter, n_hill_par = 2)
re_2par_list = pull_parameters(re_chains_2param)



#re_chains_small = re_chains

# fit DP clustering model
# TODO : cluster negative sill as one, remaining with DP?
set.seed(131)
cluster_chain = DP_MCMC_fit((re_par_list$slope_params[1:18]), n_iter = clust_iter)
beepr::beep()

plot(cluster_chain$alpha_vals[seq(1000, clust_iter, by=20)])


# check what top assignments are
n_top = 20
clust_centers_w_prob = cluster_centers(cluster_chain, n_top = n_top)
# visualize clusters
#pdf("Output/cluster_vis.pdf", width = 8, height = 5)
visualize_clusters_blocks(re_par_list$slope_params, clust_centers_w_prob)
#dev.off()

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


# fixed assignments for GCA and IA
GCA_assign = 1
names(GCA_assign) = do.call(paste, as.list(rep(1, n_chems)))
DP_2par_list = list("centers" = matrix(rep(1, n_chems),nrow=1),
                    "cluster_assign" = GCA_assign,
                    "cent_sd" = matrix(rep(0, n_chems), nrow=1), 
                    "cluster_prob" = 1)
GCA_par_list = c(re_2par_list, DP_2par_list)

# for IA, allow different slopes
IA_assign =1 
names(IA_assign) = do.call(paste, as.list( 1:n_chems))
DP_2parIA_list = list("centers" = matrix(re_par_summary$slope_params, nrow = 1),#matrix(rep(1, n_chems),nrow=1),
                      "cluster_assign" = IA_assign,
                      "cent_sd" = matrix(rep(0, n_chems), nrow=1), 
                      "cluster_prob" = 1)
IA_par_list = c(re_par_list, DP_2parIA_list)
beepr::beep(2)

# For RGCA no DP, allow different slopes
GCA_assign = 1
names(GCA_assign) = do.call(paste, as.list(rep(1, n_chems)))
RDP_2par_list = list("centers" = matrix(re_par_summary$slope_params, nrow = 1),
                     "cluster_assign" = GCA_assign,
                     "cent_sd" = matrix(rep(0, n_chems), nrow=1), 
                     "cluster_prob" = 1)
RGCA_par_list = c(re_2par_list, RDP_2par_list)



# bootstraps ####
# sample from posterior parameter estimates and clustering assigments
set.seed(1026)
n_bootstraps = 100
samp_idx = sample(1:n_top, size = n_bootstraps, prob =cluster_prob, replace = T)
sampled_mix_funs = sapply(samp_idx, FUN = function(x) create_mix_calc(x, tot_par_list, add_RE = T))
sampled_mix_funs_GCA = sapply(rep(1, n_bootstraps), FUN = function(x) create_mix_calc(x, GCA_par_list, add_RE = F))
sampled_mix_funs_IA = sapply(rep(1, n_bootstraps), FUN = function(x) create_mix_calc(x, IA_par_list, add_RE = F))
# additional comparisons:  
# use unit slopes, so no RGCA
sampled_mix_funs_GCA_noR = sapply(samp_idx, FUN = function(x) create_mix_calc(x, tot_par_list, add_RE = F, unit_slopes=T))
# only use the top clustering, no DP
sampled_mix_funs_noDP = sapply(rep(1, n_bootstraps), FUN = function(x) create_mix_calc(x, RGCA_par_list, add_RE = F))

bootstrap_calc_list = list("RGD" = sampled_mix_funs, 
                           "G" = sampled_mix_funs_GCA,
                           "I" = sampled_mix_funs_IA, 
                           "DG" = sampled_mix_funs_GCA_noR,
                           "RG" = sampled_mix_funs_noDP)




# verify 1-comp mix ####
plot_dummy_mixture(Cx, y_i, tot_par_list, replicate_sets, 
                   bootstrap_calc_list, test_idx =5)



# Compare mix vs mix constituents
plot_mix_vs_individuals(mix_idx = 34, ymax = 120)

# if using QSAR or docking scores, compute
QSAR=F
if(QSAR){
  get_CAS =function(long_cas){paste(strsplit(long_cas, split = "-")[[1]][1:3],collapse="-")}
  AR_scores = read.csv("/Users/zilberds/Desktop/tox_mixtures/pythonenv/Ligand_Scores_Weights.csv", row.names = 1)
  AR_scores = AR_scores[-which(AR_scores$Name=="Testosterone"),]
  AR_scores$Ligand_id = sapply(AR_scores$Ligand_id, get_CAS)
  dock_scores =AR_scores$Dock_Score
  mol_weight_unordered =AR_scores$Molecular_wt
  
  score_ordering = sapply(pure_unique_CAS, FUN = function(x)
    which(x==AR_scores$Ligand_id))
  mol_weight = mol_weight_unordered[score_ordering]
  cbind(tot_par_list$slope_params, dock_scores[score_ordering], mol_weight)
  
}


# ER/AR mix ####
# compare to tox21 mixtures 34 109 159
# good RE?  12, 14, 20  (18:  good for both!)
# good GCA: 19

#check 4 ad 12. # 1, 12, 20, 29, 43, 52, 55, 57
set_1 = c(1, 5, 10, 12, 20, 25, 29, 30, 31, 32, 43, 45, 50, 52, 55, 57, 62)
set_4x = c( 1, 12, 20, 29, 43, 52, 55, 57)
subindx = which(set_1 %in% set_4x)
mix_idx = 1#16, 8, 27!

binary_mixes = c(8, 34, 37, 61)
small_mixes = c(6, 26, 47 )
#pdf(file = "Output/ARluc_5compare.pdf", onefile= T, width=8, height = 5)
#pdf(file = "Output/ARluc_5compare_%02d.pdf",onefile=F, width=9, height = 5)
#pdf(file = "Output/ARluc_DR_%02d.pdf",onefile=F, width=9, height = 5)
#pdf(file = "Output/ARluc_binary_%02d.pdf",onefile=F, width=9, height = 5)
score_matrix = plot_mixture_response(small_mixes, mix_df, mix_conc_df, mix_guide, 
                                     bootstrap_calc_list)
#dev.off()
beepr::beep(2)



#Score results ####
score_df = data.frame(score_matrix)
# drop rows of 0
score_df = score_df[apply(score_df, MARGIN=1, FUN = function(rx) any(rx>0)),]
names(score_df) = c("Mix id", paste(names(bootstrap_calc_list), c("LLH")),
                    paste(names(bootstrap_calc_list), c("MSE")),
                    paste(names(bootstrap_calc_list), c("CRPS")))

## Violin Plot for Scores ####
#pdf(file = "Output/boxplot_5compare_4xec.pdf", onefile= T, width=8, height = 4)
plot_scores(score_df[subindx,], bootstrap_calc_list)
#dev.off()


mix_descriptions = sapply(score_df$`Mix id`,  FUN = function(x) 
  mix_guide$Description[which(mix_guide$CAS == mix_df$CAS[x])])
mix_descriptions = stringr::str_replace(mix_descriptions, "200%", "4x EC50")
score_df =cbind("Mix Desc" =mix_descriptions , score_df)
#score_df = score_df[,-2]
write.csv(score_df, "Output/ARluc_5comparex.csv")
#matplot(x =Cx_axis_values, t(curve_data), type = "l", log = "x", add=F)
#points(Cx_axis_values, resp_y_values, lwd=3 )
#plot(Cx_axis_values, resp_y_values, lwd=3 ,log="x")



### check which. mixtures have 0 agonists 
no_effect_mix = c()
for(mix_idx in 1:nrow(mix_df)){
  if(mix_df$ReplicateSet[mix_idx]>1) next
  chem_conc_matr = get_conc_matrix(mix_idx)
  if(all(chem_conc_matr[,AR_agonist_rows]==0)) no_effect_mix = c(no_effect_mix, mix_idx)
}

score_df_filt = (score_df[-no_effect_mix,])

# among the cases w/o synergy or antagnoism, best?
best_crps_by_mix = apply(score_df[,8:10],MARGIN=1,  FUN = function(rx) which(rx==min(rx)))
table(best_crps_by_mix)
for(rgc_bet  in as.numeric(names(which(best_crps_by_mix == 1)))){
  print(paste(rgc_bet, mix_guide$Description[ which(mix_guide$CAS == mix_df$CAS[rgc_bet])]))
}




valid_cols = apply(score_df[1:69,2:4], MARGIN=1, FUN = function(rx) all(is.finite(rx)))
apply(score_df[1:69,2:4], MARGIN=2, FUN = function(rx) sum(is.finite(rx)))






# Compare mixtures with and without no-effect chemicals
#pdf(file = "Output/No_effect_EC_EP.pdf", onefile= T, width=10, height = 5)
compare_exclude_include()
#dev.off()


# Tabular stuff for latex ####
pure_ordering = sapply(pure_unique_CAS, FUN = function(x) which(x==CAS_nums[relevant_guide_CASs]))
chem_map = sapply(names(pure_ordering), FUN = function(cas_num) pure_df$Sample.Name[which(pure_df$CAS==cas_num)[1]])

# full RE_MCMC specification table
RE_curve_fits = as.data.frame(list( "Name" =chem_map ,
                                    "Sill" = re_par_summary$sill_params, 
                                    "EC50" = re_par_summary$ec50_params,
                                    "Sill_RE_sd" = re_par_summary$u_RE_sd_params,
                                    "Int_RE_sd" = re_par_summary$v_RE_sd_params,
                                    "Slope" = re_par_summary$slope_params ))
RE_curve_fits = RE_curve_fits %>%
  mutate_at(vars(Sill, Slope,Sill_RE_sd,Int_RE_sd), function(x) (round(x, 2))) %>%
  mutate_at(vars(EC50), function(x) (signif(x, digits=3)))
latex.tabular(as.tabular(RE_curve_fits))






# MLE and RE-MCMC parameter table
RE_curve_fits = as.data.frame(list( "Name" =chem_map ,
                                    "Sill" = re_par_summary$sill_params, 
                                    "EC50" = re_par_summary$ec50_params,
                                    "Slope" = re_par_summary$slope_params ))
RE_curve_fits = RE_curve_fits %>% 
  mutate_at(vars(Sill, Slope), function(x) (round(x, 2))) %>%
  mutate_at(vars(EC50), function(x) (signif(x, digits=3)))

DRC_curve_fits = as.data.frame(curve_fits) %>%
  mutate_at(vars(V1, V3), function(x) (round(x, 2))) %>%
  mutate_at(vars(V2), function(x) (signif(x, digits=3)))

names(DRC_curve_fits) = c("Sill (drc)"," EC50 (drc)", "Slope (drc)")
cbind(RE_curve_fits, DRC_curve_fits)
latex.tabular(as.tabular(cbind(RE_curve_fits, DRC_curve_fits)))



# Mixture and CAS descriptions
as.table(score_df$`Mix Desc`, row.names = rownames(score_df))
latex.tabular(as.tabular(matrix(score_df$`Mix Desc`, ncol=1)))
latex.tabular(as.tabular(score_df[,2:ncol(score_df)]))

#p{0.35\linewidth} | p{0.6\linewidth}

# simple table for cluster centers, tbc
latex.tabular(as.tabular(as.matrix(clust_centers_w_prob$assign)))





