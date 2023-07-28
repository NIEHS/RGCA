#RE_DP_simulation_study
get_emp_llh = function(obs_res, boot_vals){
  bw <- bw.SJ(boot_vals)
  # sum equivalent to a kernel density estimate at the point
  p_density = sum(dnorm(obs_res, mean = boot_vals, sd = bw))/length(boot_vals)
  return(log(p_density))
} 



#### DP-only study ####

#settings
use_random_cluster = F
show_plots = F
n_top=2
mix_size = c(4, 8, 16)
n_big_iters = 25
#record metrics
record_MSE = matrix(0, nrow = n_big_iters*length(mix_size), ncol=4)
record_CRPS = matrix(0, nrow = n_big_iters*length(mix_size), ncol=4)

for(nc in 1:length(mix_size)){
  n_chems = mix_size[nc]
  for( big_loop_iter in 1:n_big_iters){
    #Generate individual curves:
    #random a, random slope, random ec50:
   # n_chems =  4
    a_vec = runif(n_chems, min = 1.5, max = 10)
    slope_vec = runif(n_chems, min = .2, max = 2)
    ec_vec = runif(n_chems, min = .1, max = 20)
    # make one chemical extremely toxic or not, compare Mix vs RE
    # a_vec[1]= 40
    #  ec_vec[1] = .01
    # fit DP clustering model
    cluster_chain = DP_MCMC_fit(slope_vec, n_iter = 5000)
    # allow for true cluster to be unlikely, idx>n_top
    clust_centers_w_prob = cluster_centers(cluster_chain, n_top = n_top, plot_hist = F)
    clust_rand  = as.numeric(strsplit(names(clust_centers_w_prob$assign[sample(1:n_top, 1)]),split = " ")[[1]])
    #assume either random clustering or use DP
    if(use_random_cluster) clust_rand =  sample(1:sample(2:6, 1),size = n_chems, replace = T)
    # generate mixture curve from IA/GCA
    sim_param_matrix = as.matrix(cbind("a" = a_vec, 
                                       "b" = ec_vec, 
                                       "c" = slope_vec,
                                       "max_R" =max(a_vec),
                                       "d" = 0))
    mix_calc = mix_function_generator(sim_param_matrix, 
                                      clust_rand)
    
    
    # generate mix concentrations:  each row of matrix is one mix dose
    #equipotent, equimolar, dominant
    n_samps = 30
    equipot_conc_matrix = matrix(0, nrow = n_samps, ncol = n_chems)
    equimol_conc_matrix = matrix(0, nrow = n_samps, ncol = n_chems)
    
    for(chem_idx in 1:n_chems){
      equipot_conc_matrix[, chem_idx] =ec_vec[chem_idx]/(10^seq(3, -2, length.out = n_samps))
      equimol_conc_matrix[,chem_idx] = 1/(10^seq(4, -2, length.out = n_samps))
    }  
    
    
    
    # compute true curves
    true_response_equipot = sapply(1:n_samps, function(x) mix_calc(equipot_conc_matrix[x,]))
    if(show_plots){
      x_range = range(c(rowSums(equimol_conc_matrix), rowSums(equipot_conc_matrix)))
      plot(rowSums(equipot_conc_matrix), true_response_equipot, log = "x", ylim =c(0, max(a_vec)), xlim = x_range)
      true_response_equimol =sapply(1:n_samps, function(x) mix_calc(equimol_conc_matrix[x,]))
      points(rowSums(equimol_conc_matrix), true_response_equimol, col=2)
    }
    
    
    
    #compare GCA, IA, our method via UQ: CRPS, etc
    
    sim_fixed_params= list("sill_params" = a_vec,  
                           "sill_sd" = a_vec*0, 
                           "ec50_params" = ec_vec,
                           "ec50_stdev" = ec_vec*0,
                           "u_RE_sd_params" = rep(0, n_chems),
                           "v_RE_sd_params" = rep(0, n_chems),
                           "slope_params" = slope_vec)
    
    # DP method: cluster
    #sample clusterings and pass into the mixture estimation function
    cluster_prob = clust_centers_w_prob$probs
    centers = clust_centers_w_prob$centers
    cent_sd = clust_centers_w_prob$center_sd
    cluster_assign = clust_centers_w_prob$assign
    DP_par_list = list("centers" = centers, 
                       "cluster_assign" = cluster_assign,
                       "cent_sd" = cent_sd, 
                       "cluster_prob" = cluster_prob)
    DP_par_list_sim= c(sim_fixed_params, DP_par_list)
    
    
    # GCA/IA 
    
    # fixed assignments for GCA and IA
    GCA_assign = 1
    names(GCA_assign) = do.call(paste, as.list(rep(1, n_chems)))
    DP_2par_list = list("centers" = matrix(rep(1, n_chems),nrow=1),
                        "cluster_assign" = GCA_assign,
                        "cent_sd" = matrix(rep(0, n_chems), nrow=1), 
                        "cluster_prob" = 1)
    GCA_par_list_sim = c(sim_fixed_params, DP_2par_list)
    
    # for IA, assume mean 1?  or all slopes different?
    IA_assign =1 
    names(IA_assign) = do.call(paste, as.list( 1:n_chems))
    DP_2parIA_list = list("centers" = matrix(rep(1, n_chems),nrow=1),#matrix(tot_par_list$slope_params),
                          "cluster_assign" = IA_assign,
                          "cent_sd" = matrix(rep(0, n_chems), nrow=1), 
                          "cluster_prob" = 1)
    IA_par_list_sim = c(sim_fixed_params, DP_2parIA_list)
    
    
    n_bootstraps = 50
    samp_idx = sample(1:n_top, size = n_bootstraps, prob =cluster_prob, replace = T)
    sim_mix_funs = sapply(samp_idx, FUN = function(x) create_mix_calc_from_summary(x, DP_par_list_sim))
    sim_mix_funs_GCA = sapply(rep(1, n_bootstraps), FUN = function(x) create_mix_calc_from_summary(x, GCA_par_list_sim))
    sim_mix_funs_IA = sapply(rep(1, n_bootstraps), FUN = function(x) create_mix_calc_from_summary(x, IA_par_list_sim))
    
    
    ## Plot results for simulated data
    conc_mat = equipot_conc_matrix
    curve_data = matrix(0, ncol = n_samps, nrow = n_bootstraps)
    curve_data_GCA = matrix(0, ncol = n_samps, nrow = n_bootstraps)
    curve_data_IA = matrix(0, ncol = n_samps, nrow = n_bootstraps)
    for(row_idx in 1:n_samps){
      conc_val = conc_mat[row_idx,]#+1e-40
      #Rprof()
      curve_data[,row_idx] = sapply(sim_mix_funs, 
                                    FUN = function(x) x(conc_val))
      #Rprof(NULL)
      #summaryRprof("Rprof.out")
      curve_data_GCA[,row_idx] =  sapply(sim_mix_funs_GCA, 
                                         FUN = function(x) x(conc_val))
      curve_data_IA[,row_idx] =  sapply(sim_mix_funs_IA, 
                                        FUN = function(x) x(conc_val))
    }
    #Cx_axis_values = array(unlist(mix_df[mix_idx,conc_idx_T21_matrix]))
    Cx_axis_values = rowSums(conc_mat)
    
    if(show_plots) make_ggplot(Cx_axis_values, 
                               true_response_equipot,
                               curve_data,
                               curve_data_GCA,
                               curve_data_IA)
    
    # get MSE, Logscore, CRPS
    ## Score 1: LH ####
    n_x_values= ncol(Cx)
    ClGCA_llh = 0; ClGCA_mse = 0; RGCA_crps = 0
    GCA_llh = 0; GCA_mse = 0; GCA_crps = 0
    IA_llh = 0; IA_mse = 0;  IA_crps= 0
    
    for(x_idx in 1:n_x_values){
      # ClGCA_llh = ClGCA_llh+ get_emp_llh(true_response_equipot[x_idx], curve_data[,x_idx])
      # GCA_llh = GCA_llh+ get_emp_llh(true_response_equipot[x_idx], curve_data_GCA[,x_idx])
      # IA_llh = IA_llh+ get_emp_llh(true_response_equipot[x_idx], curve_data_IA[,x_idx])
      RGCA_crps = RGCA_crps + crps_sample(true_response_equipot[x_idx], curve_data[,x_idx])
      GCA_crps = GCA_crps+ crps_sample(true_response_equipot[x_idx], curve_data_GCA[,x_idx])
      IA_crps = IA_crps+ crps_sample(true_response_equipot[x_idx], curve_data_IA[,x_idx])
      
      
      ClGCA_mse=ClGCA_mse + (true_response_equipot[x_idx] - median(curve_data[,x_idx]))^2/n_x_values
      GCA_mse=GCA_mse + (true_response_equipot[x_idx] - median(curve_data_GCA[,x_idx]))^2/n_x_values
      IA_mse=IA_mse + (true_response_equipot[x_idx] - median(curve_data_IA[,x_idx]))^2/n_x_values
    }
    #score_matrix[mix_idx,] = c(mix_idx, ClGCA_llh,GCA_llh,IA_llh, ClGCA_mse, GCA_mse, IA_mse)
    
    curr_idx = big_loop_iter+n_big_iters*(nc-1)
    record_MSE[curr_idx, ]=c(ClGCA_mse, GCA_mse, IA_mse, n_chems)
    record_CRPS[curr_idx, ]=c(RGCA_crps, GCA_crps, IA_crps, n_chems)
    # annotation offset: 10% of max(curve_data_IA)
  }
}
beepr::beep()

pdf(file = "sim_study.pdf", width = 7, height = 3)
sim_study_boxplot(record_MSE, record_CRPS)
dev.off()

# old simple boxplot
# record_CRPS = data.frame(record_CRPS)
# names(record_CRPS) = c("RE+DP", "GCA", "IA")
# pdf(file = "CRPS_2step_4_chems.pdf")
# boxplot(record_CRPS, main = "CRPS, Simulation w 8 Chems")#, ylim = c(0, 20))
# dev.off()

