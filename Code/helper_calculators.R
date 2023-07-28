
#' Hill function
#'
#' A simple function to compute the Hill function response given parameters and
#' a concentration
#'
#' @param a the maximum effect or sill
#' @param b the EC50
#' @param c the slope
#' @param conc an input concentration or dose
#'
#' @return a real number prediction of the dose response
#' @export
#'
#' @examples
hill_function = function(a,b,c,conc){
  a/(1+(b/conc)^c)
}

#' Maximum likelihood dose response curves
#' 
#' Wrapper function to apply the drc (dose response curve) package to the input
#' data.
#'
#' @param y_i matrix of observed dose responses (columns) for multiple chemicals
#'   and replicates (rows)
#' @param Cx a matrix of the doses corresponding to the responses
#' @param replicate_sets a list of vectors where each vector contains the
#'   indices of replicates for a particular chemical
#'
#' @return a matrix of parameters with a row for each chemical and columns corresponding to ()
#' @export
#'
#' @examples
get_mle_curve_fits = function(y_i, Cx, replicate_sets ){
  # group by replicates
  replicate_counts = rep(1:18, unlist(lapply(replicate_sets, length)))
  replicate_vec = unlist(replicate_sets)
  ordering_replicates = order(replicate_vec)
  grouping_vec = replicate_counts[ordering_replicates]
  
  #individual fits
  curve_fits = matrix(0, nrow = length(replicate_sets), ncol=3)
  for(ri in 1:length(replicate_sets)){
    drc_yi = c(y_i[replicate_sets[[ri]],])
    drc_xi = c(Cx[replicate_sets[[ri]],])
    fit_coeff = try(drc::drm(drc_yi~drc_xi, 
                             fct=drc::LL.3())$coefficients, silent = T)
    if(class(fit_coeff) == "try-error") next
    fit_coeff = fit_coeff[c(2,3,1)]
    fit_coeff[3] = -fit_coeff[3]
    fit_coeff[2] = fit_coeff[2]
    curve_fits[ri, ]=fit_coeff
  }
  return(curve_fits)
}



get_MCMC_diagnostics  = function(re_chains, re_chains2, re_chains3) {
  require(coda)
  get_GR_Diagnostic = function(data_id, input_cols){
    mcmc1 = coda::as.mcmc(re_chains[[data_id]][,input_cols])
    mcmc2 = coda::as.mcmc(re_chains2[[data_id]][,input_cols])
    mcmc3 = coda::as.mcmc(re_chains3[[data_id]][,input_cols])
    mcmclist = coda::as.mcmc.list(list(mcmc1, mcmc2, mcmc3))
    coda_test = coda::gelman.diag(mcmclist)$psrf
    return(coda_test)
  }
  
  coda_u_RE = get_GR_Diagnostic("u_RE", which(re_chains$u_RE[1,]!=0))
  coda_u_RE_sd = get_GR_Diagnostic("u_RE_sd", 1:ncols(re_chains$u_RE_sd))
  coda_sigma = get_GR_Diagnostic("sigma", 1:ncols(re_chains$sigma))
  coda_slope = get_GR_Diagnostic("slope_record", 1:ncols(re_chains$slope_record))
  coda_sill = get_GR_Diagnostic("sill_mideffect_record", 1:n_chems)
  
  # EC50 part, needs log
  mcmc1 = coda::as.mcmc(log(re_chains$sill_mideffect_record[, (1+n_chems):(2*n_chems)]))
  mcmc2 = coda::as.mcmc(log(re_chains2$sill_mideffect_record[, (1+n_chems):(2*n_chems)]))
  mcmc3 = coda::as.mcmc(log(re_chains3$sill_mideffect_record[, (1+n_chems):(2*n_chems)]))
  mcmclist = coda::as.mcmc.list(list(mcmc1, mcmc2, mcmc3))
  coda_ec50 = coda::gelman.diag(mcmclist)$psrf
  
  all_dats = round(cbind(coda_sill, coda_ec50,coda_slope,coda_sigma, coda_u_RE_sd ), digits =2)
  all_dats = as.data.frame(all_dats)
  new_names = array(sapply(c("Sill", "EC50", "Slope", "Sigma", "Variance_U"), 
                           FUN = function(sx) (paste(sx, c("Point est.", "Upper C.I.")))))
  names(all_dats)=new_names
  rownames(all_dats) = chem_map
  
  return(all_dats)
}





# Some notes on posterior summary stats:  
# mode (MAP) is usually the lowest value, followed by median and mean due to skewness
# MAP or median makes sense for heavily skewed parameters, but slope was not so skewed due to strong prior
modes <- function(dats){
  if(all(dats == 0)) return(0)
  dens = density(dats)
  i <- which(diff(sign(diff(dens$y))) < 0) + 1
  max_i = i[which(dens$y[i] == max(dens$y[i]))]
  dens$x[max_i]
}

pull_summary_parameters = function(re_chains, summry_stat = median){
  thin_idx = seq(5000, re_iter-100, by=20)
  #sill_params = colMeans(re_chains$sill_mideffect_record[thin_idx, 1:n_chems])
  sill_params = apply(re_chains$sill_mideffect_record[thin_idx, 1:n_chems], MARGIN=2, summry_stat)
  sill_sd = apply(re_chains$sill_mideffect_record[thin_idx, 1:n_chems], MARGIN = 2, sd)
  #insignificant_sill = sapply(1:length(sill_params), FUN = function(idx) abs(sill_params[idx])<1.5*sill_sd[idx] )
  #sill_params[insignificant_sill] = 0
  #ec50_params = colMeans(re_chains$sill_mideffect_record[thin_idx, 1:n_chems+n_chems])
  ec50_params = apply(re_chains$sill_mideffect_record[thin_idx, 1:n_chems+n_chems], MARGIN = 2, summry_stat)
  ec50_sd = apply(re_chains$sill_mideffect_record[thin_idx, 1:n_chems+n_chems], MARGIN=2, sd)
  u_RE_params = apply(re_chains$u_RE[thin_idx,], MARGIN=2, summry_stat)
  v_RE_params = apply(re_chains$v_RE[thin_idx,], MARGIN=2, summry_stat)
  #u_RE_sd_chain= re_chains$u_RE_sd[-which(re_chains$u_RE_sd[,2]>1000),]
  u_RE_sd_params = apply(re_chains$u_RE_sd[thin_idx,], MARGIN=2, summry_stat)
  v_RE_sd_params = apply(re_chains$v_RE_sd[thin_idx,], MARGIN=2, summry_stat)
  phi_c = apply(re_chains$slope_record[thin_idx,], MARGIN=2, summry_stat)
  return(list("sill_params" = sill_params, 
              "sill_sd" = sill_sd, 
              "ec50_params" = ec50_params,
              "ec50_stdev" = ec50_sd,
              "u_RE_params" = u_RE_params,
              "v_RE_params" = v_RE_params,
              "u_RE_sd_params" = u_RE_sd_params,
              "v_RE_sd_params" = v_RE_sd_params,
              "slope_params" = phi_c))
}

pull_parameters = function(re_chains, summry_stat = median){
  thin_idx = seq(5000, re_iter-100, by=20)
  sill_sample = re_chains$sill_mideffect_record[thin_idx, 1:n_chems]
  ec50_sample = re_chains$sill_mideffect_record[thin_idx, 1:n_chems+n_chems]
  u_RE_center = sapply(1:n_chems,FUN = function(ridx) median(re_chains$u_RE[thin_idx,replicate_sets[[ridx]]]))
  u_RE_sd_est = apply(re_chains$u_RE_sd[thin_idx,], MARGIN=2, summry_stat)
  v_RE_sd_est = apply(re_chains$v_RE_sd[thin_idx,], MARGIN=2, summry_stat)
  # slopes must be esimated for DP clustering
  phi_c = apply(re_chains$slope_record[thin_idx,], MARGIN=2, summry_stat)
  return(list("sill_sample" = sill_sample, 
              "ec50_sample" = ec50_sample,
              "u_RE_center" = u_RE_center,
              "u_RE_sd" = u_RE_sd_est,
              "v_RE_sd" = v_RE_sd_est,
              "slope_params" = phi_c))
}




# review chains ####
if(F){
  thin_idx = seq(5000, re_iter-100, by=20)
  paroi = re_chains$sill_mideffect_record[thin_idx, 17+0*n_chems]
  points(paroi, log="y")#, ylim = c(1e-8, 1000))
  
  paroi = re_chains$sill_mideffect_record[thin_idx, 17+0*n_chems]
  paroi_sl = re_chains$slope_record[thin_idx, 17]
  paroi_ec = re_chains$sill_mideffect_record[thin_idx, 17+1*n_chems]
  plot(paroi)
  plot(paroi_ec, log = "y")
  plot(paroi_sl, log="y")
  plot(re_chains$slope_record[thin_idx, 1]/((re_chains$sill_mideffect_record[thin_idx, 1+n_chems]) / 
                                              (re_chains$sill_mideffect_record[thin_idx, 1])), log="y")
  
  plot(log(re_chains$tau[thin_idx]))
  
  
  # compare:  6 vs 1 vs 17 
  #pdf("MCMC_traces.pdf", width = 10, height = 8)
  par(mfrow=c(3,4))
  for( idx in c(1, 6, 17)){
    paroi = re_chains$sill_mideffect_record[thin_idx, idx+0*n_chems]
    paroi_sl = re_chains$slope_record[thin_idx, idx]
    paroi_ec = re_chains$sill_mideffect_record[thin_idx, idx+1*n_chems]
    plot(re_chains$sigma[thin_idx, idx])
    # if(idx == 17){plot(paroi, ylim = c(-.5, 7),main = paste("Sill Trace plot, chemical ", idx), ylab = NA)
    # }else 
      plot(paroi, main = paste("Sill Trace plot, chemical ", idx), ylab = NA)
    plot(paroi_ec, log = "y",  main = paste("EC50 Trace plot (log y), chemical ", idx), ylab = NA)
    plot(paroi_sl, main = paste("Slope Trace plot, chemical ", idx), ylab = NA)
  }
  #dev.off()
  plot(re_chains$sill_mideffect_record[thin_idx, idx+1*n_chems], log= "y")
  plot(re_chains$slope_record[thin_idx, idx], log= "y")
  
  thin_idx = seq(10000, re_iter-100, by=40)
  plot(hill_function(re_chains$sill_mideffect_record[thin_idx, idx+0*n_chems],
                re_chains$sill_mideffect_record[thin_idx, idx+1*n_chems],
                re_chains$slope_record[thin_idx, idx], 1e-10))
  
}



#--------------------- Mix Calulators ------------------- ####
#' These functions are wrappers for the factory method mix_function_generator. 
#' The wrappers provide different ways of sampling the parameters.  


#' Get a function that computes mixture response
#' 
#' A factory method that returns a function that computes the mixture response.
#' Uses summary statistics to create a random sample rather than directly
#' sampling from the posterior MCMC chains
#'
#' @param idx Specifies which clustering to apply to the parameters 
#' @param par_list Contains all estimated parameters and DP cluster options
#'
#' @return a function that takes a concentration vector as input and yields a 
#' predicted response as output 
#' ie a "calculator" for the mixture effect given component concentrations
#' @export
#'
#' @examples
create_mix_calc_from_summary = function(idx, par_list){
  # bootstrap the cluster slope from top clusters
  cluster_name = names(par_list$cluster_assign[idx])
  cluster_assign_vec = as.numeric(strsplit(cluster_name, split = " ")[[1]])
  slope_mean = par_list$centers[idx,]
  slope_sd = par_list$cent_sd[idx,]
  # theta_1, theta_2 bootstraps from mcmc chains?
  #mcmc_bootstrap = sample(1:nrow(sill_chain), 1)
  # must sample both from same idx since diff idx may have small likelihood
  n_chem_pars = length(par_list$sill_params)
  tot_sd = sqrt(par_list$u_RE_sd_params^2+par_list$sill_sd^2)
  sill_params_boot =  truncnorm::rtruncnorm(n = n_chem_pars, a = 0,
                                            mean = par_list$sill_params,
                                            sd = tot_sd)
  ec50_params_boot = truncnorm::rtruncnorm(n = n_chem_pars, a = 0,
                                           mean = par_list$ec50_params,
                                           sd = par_list$ec50_stdev)
  d_boot = rnorm(n = n_chem_pars, 
                 mean=0, 
                 sd = par_list$v_RE_sd_params)
  max_effect_R =  max(par_list$sill_params, sill_params_boot)
  # instead:  max per cluster?
  # n_clust = length(unique(cluster_assign_vec))
  # max_effect_R = rep(0, n_clust)
  # for(i in 1:n_clust) max_effect_R[i] = max(sill_params_boot[which(cluster_assign_vec==i)])
  # test: control has maxR=100
  #max_effect_R = 100
  # generate calculator for 
  bootstrap_param_matrix = as.matrix(cbind("a" = sill_params_boot, 
                                           "b" = ec50_params_boot, 
                                           "c" = slope_mean,
                                           "max_R" =max_effect_R,
                                           "d" = 0
  ))
  mixture_calculator = mix_function_generator(bootstrap_param_matrix, 
                                              cluster_assign_vec)
  return(mixture_calculator)
}


#' Get a function that computes the mixture response.
#' 
#'  A factory method that returns a function that computes the mixture response.
#' Samples directly from the posterior MCMC chains for the EC50 and sill
#' parameters. Slope is taken from the clustering.  Noise is added to some
#' parameters according to the random effect variance.
#'
#'
#' @param idx Specifies which clustering to apply to the parameters 
#' @param par_list a data frame with individual chemical dose response parameters 
#' @param add_RE boolean to include or exclude random effect variances
#' @param unit_slopes boolean to fix slopes to 1 (but still use slope clustering)
#'                    Used for special case of GCA, where slope = 1
#'
#' @return function to take a concentration vector as input and response as output
#' @export
#'
#' @examples
create_mix_calc = function(idx, par_list, add_RE=T, unit_slopes = F){
  # bootstrap the cluster slope from top clusters
  cluster_name = names(par_list$cluster_assign[idx])
  cluster_assign_vec = as.numeric(strsplit(cluster_name, split = " ")[[1]])
  # sample the cluster value rather than use mean/sd
  slope_mean = par_list$centers[idx,]
  slope_sd = par_list$cent_sd[idx,]
  # add noise from DP cluster to centers
  slope_noise = rnorm(length(unique(slope_sd)), mean = 0, sd = unique(slope_sd))
  slope_boot = slope_mean+slope_noise[cluster_assign_vec]
  if(unit_slopes) slope_mean = slope_mean*0+1 
  # sample values for the sill parameter, can be independent
  n_chem_pars = ncol(par_list$sill_sample)
  sill_params_samp = apply(par_list$sill_sample, MARGIN = 2, 
                           FUN =  function(colx) sample(colx, 1)) 
  # add iid gauss noise for RE variance
  sill_RE_boot = rnorm(n_chem_pars, mean = par_list$u_RE_center, sd = par_list$u_RE_sd)
  sill_params_boot= sill_params_samp  
  if(add_RE)sill_params_boot  =sill_params_boot + sill_RE_boot
  # sample ec50
  ec50_params_boot = apply(par_list$ec50_sample, MARGIN = 2, 
                           FUN =  function(colx) sample(colx, 1)) 
  #We dont sample intercept RE: assume intercept 0, d=0
  # TODO consider empirical max, or clusterwise max
  max_effect_R =  max(sill_params_boot)
  # generate calculator for 
  bootstrap_param_matrix = as.matrix(cbind("a" = sill_params_boot, 
                                           "b" = ec50_params_boot, 
                                           "c" = slope_mean,
                                           "max_R" =max_effect_R,
                                           "d" = 0
  ))
  mixture_calculator = mix_function_generator(bootstrap_param_matrix, 
                                              cluster_assign_vec)
  return(mixture_calculator)
}

#' Get a function that computes the mixture response
#'
#' A slight variation on the standard create_mix_calc, the parameters are still
#' sampled from the posterior MCMC but they are sampled with a single index so
#' that the set of parameters is feasible. Sampling all the parameters randomly
#' could cause an issue when there is non-identifiability, and an extreme value
#' in one parameter is not offset by a small value in another parameter.
#'
#' @param idx Specifies which clustering to apply to the parameters
#' @param par_list a data frame with individual chemical dose response
#'   parameters
#' @param add_RE boolean to include or exclude random effect variances
#'
#' @return a function to take a concentration vector as input and response as
#'   output
#' @export
#'
#' @examples
create_mix_calc_sample_row = function(idx, par_list, add_RE=T){
  # bootstrap the cluster slope from top clusters
  cluster_name = names(par_list$cluster_assign[idx])
  cluster_assign_vec = as.numeric(strsplit(cluster_name, split = " ")[[1]])
  # sample the cluster value rather than use mean/sd
  slope_mean = par_list$centers[idx,]
  slope_sd = par_list$cent_sd[idx,]
  # add noise from DP cluster to centers
  slope_noise = rnorm(length(unique(slope_sd)), mean = 0, sd = unique(slope_sd))
  slope_boot = slope_mean+slope_noise[cluster_assign_vec]
  # sample values for the sill parameter, can be independent
  n_chem_pars = ncol(par_list$sill_sample)
  n_samples = nrow(par_list$sill_sample)
  param_idx = sample(1:n_samples, n_chem_pars, replace=T)
  sill_params_samp = rep(0, n_chem_pars)
  ec50_params_boot = rep(0, n_chem_pars)
  for(pidx in 1:n_chem_pars){
    sill_params_samp[pidx] = par_list$sill_sample[param_idx[pidx], pidx]
    ec50_params_boot[pidx] = par_list$ec50_sample[param_idx[pidx], pidx]
  }
  # add iid gauss noise for RE variance
  sill_RE_boot = rnorm(n_chem_pars, mean = 0, sd = par_list$u_RE_sd)
  sill_params_boot= sill_params_samp  
  if(add_RE) sill_params_boot  =sill_params_boot + sill_RE_boot
  #We dont sample intercept RE: assume intercept 0, d=0
  # TODO consider empirical max, or clusterwise max
  max_effect_R =  max(sill_params_boot)
  # generate calculator for 
  bootstrap_param_matrix = as.matrix(cbind("a" = sill_params_boot, 
                                           "b" = ec50_params_boot, 
                                           "c" = slope_mean,
                                           "max_R" =max_effect_R,
                                           "d" = 0
  ))
  mixture_calculator = mix_function_generator(bootstrap_param_matrix, 
                                              cluster_assign_vec)
  return(mixture_calculator)
}



#' Predict the mixture response
#'
#' This is a convenience function that applies the lists of bootstrapped mixture
#' response calculators to the matrix of concentrations that define the doses
#' for the mixture.
#'
#'
#' @param sampled_mix_funs list of boostrapped functions for the main method of
#'   the manuscript, RGCA+DP, instances produced by the factory method
#'   mix_function_generator
#' @param sampled_mix_funs_GCA list of boostrapped functions for the standard
#'   method of GCA
#' @param sampled_mix_funs_IA list of boostrapped functions for the standard
#'   method of IA
#' @param n_dose integer number of doses of the mixture
#' @param chem_conc_matr a matrix where the rows represent the constituent
#'   chemicals and the columns represent the dose.  The column sum is the
#'   mixture dose.
#' @param default_entry  a default entry for the methods, with default=0.  Can
#'   be set to null NA to aid in plotting.
#'
#' @return
#' @export
#'
#' @examples
predict_mix_response = function(sampled_mix_funs, sampled_mix_funs_GCA,  
                                sampled_mix_funs_IA, n_dose, chem_conc_matr,
                                default_entry = 0){
  n_bootstraps = length(sampled_mix_funs)
  curve_data = matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
  curve_data_GCA = matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
  curve_data_IA = matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
  for(row_idx in 1:n_dose){
    # each row of conc_matrix is a mixture dose
    conc_val = chem_conc_matr[row_idx,]#+1e-40
    # if there are missing conc values, skip prediction
    if(any(is.na(conc_val))) next
    
    curve_data[,row_idx] = sapply(sampled_mix_funs, 
                                  FUN = function(x) x(conc_val))
    curve_data_GCA[,row_idx] =  sapply(sampled_mix_funs_GCA, 
                                       FUN = function(x) x(conc_val))
    curve_data_IA[,row_idx] =  sapply(sampled_mix_funs_IA, 
                                      FUN = function(x) x(conc_val))
  }
  return(list(curve_data,curve_data_GCA, curve_data_IA))
}


#' Predict the mixture response
#'
#' A similar function to predict_mix_response.  This is a convenience function
#' that takes a list of lists as input and iteratively applies the elements,
#' which are bootstrapped mixture repsonse calculators, to the columns of the
#' mixture dose matrix.
#'
#'
#' @param n_dose number of mixture doses
#' @param chem_conc_matr a matrix where the rows represent the constituent
#'   chemicals and the columns represent the dose.  The column sum is the
#'   mixture dose.
#' @param bootstrap_calc_list a list of lists of boostrapped functions,
#'   instances produced by the factory method mix_function_generator
#' @param default_entry a default entry for the methods, with default=0.  Can be
#'   set to null NA to aid in plotting
#'
#' @return
#' @export
#'
#' @examples
predict_mix_response_many = function(n_dose, chem_conc_matr, bootstrap_calc_list, default_entry = 0){
  #bootstrapped_calc_list = list(...)
  n_bootstraps = length(bootstrap_calc_list[[1]])
  curve_list = list()
  for (bootcalc in bootstrap_calc_list){
    curve_data = matrix(default_entry, ncol = n_dose, nrow = n_bootstraps)
    for(row_idx in 1:n_dose){
      # each row of conc_matrix is a mixture dose
      conc_val = chem_conc_matr[row_idx,]#+1e-40
      # if there are missing conc values, skip prediction
      if(any(is.na(conc_val))) next
      curve_data[,row_idx] = sapply(bootcalc, 
                                    FUN = function(x) x(conc_val))
    }
    curve_list = c(curve_list, list(curve_data))
  }
  names(curve_list) = names(bootstrap_calc_list)
  return(curve_list)
}

#' Score the predictions
#'
#' Computes the log likelihood, mean square error, and continuous rank
#' probability score for predictions of the mixture response compared to the
#' observed mixture response.
#'
#' @param mix_df data frame with the observed mixture responses
#' @param mix_idx integer index specifying which observed data should be used
#' @param unplotted_repl integer indices specifying additional observed values
#'   to plot
#' @param curve_data_list list of vectors where each vector is a predicted
#'   curve.  List length matches the number of methods used to make predictions.
#'
#' @return a vector of scores with scores grouped by method.  If there are two
#'   methods and three scores, the output would have the form (M1S1, M1S2, M1S3,
#'   M2S1, M2S2, M2S3).
#' @export
#'
#' @examples
compute_mixpred_scores = function(mix_df, mix_idx, unplotted_repl, curve_data_list  ){
  require(scoringRules)
  n_x_values= ncol(Cx)
  score_matrix = c()
  # for each set of curve data
  for(curve_data in curve_data_list){
    llh = 0; mse = 0; crps = 0
    # for each replicate of data
    for(repl_idx in c(mix_idx, unplotted_repl)){
      # get the replicate response
      resp_y_repl =  array(unlist(mix_df[repl_idx,resp_idx])) 
      for(x_idx in 1:n_x_values){
        #compute scores for given replicate
        llh = llh+ scoringRules::logs_sample(resp_y_repl[x_idx], curve_data[,x_idx])
        mse = mse + (resp_y_repl[x_idx] - median(curve_data[,x_idx]))^2/n_x_values
        crps = crps + scoringRules::crps_sample(resp_y_repl[x_idx], curve_data[,x_idx])
      }
    }
    score_matrix=cbind(score_matrix, c(llh, mse, crps))
  }
  #score_matrix_row = c(mix_idx,array(t(score_matrix)))
  return(score_matrix)
}


# Not used; compute empirical LLH given bootstrapped curve estimates
get_emp_llh = function(obs_res, boot_vals){
  # get optimal bandwidth for KDE
  bw <- bw.SJ(boot_vals)
  # sum equivalent to a kernel density estimate at the point
  p_density = sum(dnorm(obs_res, mean = boot_vals, sd = bw))/length(boot_vals)
  return(log(p_density))
} 

# Not used; compute CRPS given bootstraped estimates, empirical value
get_kernel_CRPS = function(obs_res, boot_vals){
  # get optimal bandwidth for KDE
  bw <- bw.SJ(boot_vals)
  #
  interval = min(1, bw/5)
  pdf_xvals = seq(min(boot_vals)-5*bw, max(boot_vals)+5*bw, by =interval)#min(1, bw/5 ) )
  # sum equivalent to a kernel density estimate at the point
  pdf_kernel = sapply(pdf_xvals, 
                      FUN = function(x) sum(dnorm(x, mean = boot_vals, sd = bw))/length(boot_vals))
  # approximate CDF integral with summation
  cdf_kernel = cumsum(pdf_kernel)*interval
  # approximate CRPS integral with summation
  CRPS_val = sum((cdf_kernel-1*(pdf_xvals>=obs_res))^2)*interval
  return(CRPS_val)
} 

# not used
get_empr_CRPS = function(obs_res, boot_vals){
  cdf_edf = sapply(sort(boot_vals), FUN = function(x) mean(x>boot_vals) )
  CRPS_val = sum(((cdf_edf-1*(sort(boot_vals)>=obs_res))^2)[-1]*diff(sort(boot_vals)) )
  return(CRPS_val)
}



# Inverse Hill Functions ####
#' hill_invs_factry is a factory method that returns the inverse of a Hill
#' function given the parameters (ie returns the concentration given an input
#' response).  There are two sub-functions representing two methods of extending
#' the inverse.
#'
#' @param a sill (max effect)
#' @param b EC50 
#' @param c slope
#' @param max_R maximum effect across all chemicals
#' @param d minimum effect
#'
#' @return a function that computes the inverse of a response given the
#'   parameters
#' @export
#'
#' @examples
hill_invs_factry = function(a,b,c, max_R = 1, d=0){
  
  #' Hill inverse computed by reflection for RGCA.
  #'
  #' @param y the response value to invert 
  #'
  #' @return the concentration that corresponds to the input response
  #'
  #' @examples
  hilly_inverse = function(y){
    if(y==0) return(0)
    # case 1, y extended to small negative conc
    if( (a>0 & y<0) |  (a<0 & y>0)){
      return(-b / (1+(-a/y)^c))
      #return(-b / (1+(-b/y)^c)) # why use sill as ec50?
    }
    # case 2, standard inverse
    if( (a>0 & y<a)  | (a<0 & y>a)){
      return(b/(a/y - 1)^(1/c))
    }
    # case 3, reflected part of the standard inverse
    if((a>0 & y<2*a) | (a<0 & y>2*a)){
      return(-2*b-b/(a/(2*a-y)-1)^(1/c))
    }
    # case 4, reflection of extension
    return(-2*b + b/(1+(a/(-2*a + y))^c))
  }

  #' Older version of Hill inverse extension stretches the inverse to cover the
  #' relevant range between d (minimum response) and max_R (maximum response or
  #' sill across all chemicals)
  #'
  #' @param y the response value to invert 
  #'
  #' @return the concentration that corresponds to the input response
  #' 
  #' @examples
  hilly_inverse_stretch =  function(y){
    # GCA case:
    #if(c==1 && y!= d+a) return(b/(a/(y-d)-1))
    
    # negative sill (a<0) case:
    if(a<0){
      if(y>0) return(-b/(1+(-a/y)^c))
      if(y<0 && y>a) return(b/((a/(y-d)-1)^(1/c)))
      # another undefined case?
      if(y<a) return(0)
    }
    
    # general edge or undefined case:
    # regular inverse case: valid for d<y<a
    if(y<a && y>d || y<d && c==1 ) return(b/((a/(y-d)-1)^(1/c)))
    
    #extension for negative y, standard case
    if(y<0) return(-b/(1+(-a/y)^c))
    
    # edge case:  y =a is asymptote
    if(y==a ) return(1e10)
    # undefined case:  y> all responses
    if( y>=max_R  && c!=1)  return(0)
    # partial agonist case, y in [a, max_R]
    return(-2*b-b/((max_R-a)/(max_R-y)-1)^(1/c))
    #partial to match slope = 1 case: no max_R scaling
    #return(-2*b-b/(a/(2*a-y)-1)^(1/c))
  }
  return(hilly_inverse)
}





#' A factory function that returns a function with one input that can be
#' optimized with respect to that input. The returned function computes the
#' generalized concentration addition formula for the input set of chemicals and
#' input concentrations.
#'
#' @param hill_inverse_list a list of the hill inverse functions with one
#'   function per chemical
#' @param conc_vec a vector of non-negative values with the same length as the
#'   hill_inverse_list
#' @param synergy_const a scaling term for synergy. Difficult to specify so it
#'   is not used and currently set to 0, implying no synergy.
#' @param interval_sign 1 by default and can be set to -1.  Used for
#'   chemicals with a negative sill. Details: For numerical stability around
#'   very small concentrations and responses, the output function exponentiates
#'   the input response, making it impossible to invert a negative response. If
#'   at least one chemical in the mixture has a negative sill, this parameter
#'   should be set to -1 to check if the optimal (predicted) response is
#'   negative.
#'
#' @return a function GCA_over_list that has a response r as a input and a
#'   norm (measure of optimality according to GCA) as output
#' @export
#'
#' @examples
eff_response_opt = function(hill_inverse_list, conc_vec, synergy_const = 0, interval_sign= 1){
  
  
  
  
  GCA_over_list = function(r){
    # get inverse for every component of cluster for current effect level r
    invs_list = lapply(hill_inverse_list, 
                       FUN = function(x) do.call(x, list(interval_sign*exp(r))))
    invs_vec = unlist(invs_list)
    # avoid constant function when all 0
    if(all(conc_vec==0)) return(r)
    # some attempts to deal with negative concentrations, not currently needed
    # z_sign = sign(invs_vec)
    # z_vec = conc_vec/abs(invs_vec)
    # z_vec = z_vec/sum(z_vec)*z_sign
    # if using synergy: compute the GCA inverse, then normalize
    z_vec = conc_vec/invs_vec
    z_vec = z_vec/sum(z_vec)
    # by default synergy_scale = 1 (since synergy_const=0), implying no synergy
    synergy_scale = exp(synergy_const*prod(z_vec))
    # the GCA equation should balance at the optimal value:
    # [sum_i concentration_i / f_inverse(r) - 1 = 0]
    return(abs(sum(conc_vec/invs_vec)-synergy_scale))
  }
  
  
  asymmetric_GCA = function(r){
    # get inverse for every component of cluster for current effect level r
    invs_list = lapply(hill_inverse_list, 
                       FUN = function(x) do.call(x, list(exp(r))))
    invs_vec = unlist(invs_list)
    # avoid constant function when all 0
    if(all(conc_vec==0)) return(r)
    # effect undefined for current r value, return a large value
    if(any(invs_vec==0)) return(Inf)
    
    A_scale = diag(1/conc_vec)
    # if a conc is 0, sub 1/conc with 1 
    A_scale[A_scale == Inf] = 1
    A_scale[1,2] = -.1
    A_scale[2,1] = 0#-.1
    A_finv = diag(1/invs_vec)
    A = A_finv %*% A_scale
    synergy_multiplier = exp(A%*%conc_vec)
    
    # combined_gca=  conc_vec %*% A %*% conc_vec
    return(abs(conc_vec/invs_vec*synergy_multiplier -1))
  }
  
  return(GCA_over_list)
}



#' Factory method to return a function that predicts the mixture response given 
#' the concentrations of the mixture component chemicals.
#' A factory method is used because the uncertainty quantification is based on a
#' bootstrapped collection of feasible parameters, and each set of parameters
#' leads to a different predictor for the mixture effect.  Rather than keeping
#' track of lists of parameters, we immediately convert the sampled parameter
#' set to a function (a "calculator") that takes the mixture concentration as
#' input. The list of calculators can be applied to a
#'
#' @param param_matrix a simple matrix object with rows corresponding to 
#' unique chemicals and columns representing parameters (a, b, c, max_R, d), where 
#'          a is the sill
#'          b is the EC50
#'          c is the slope value
#'          max_R is the maximum sill across all chemicals (deprecated)
#'          d is the minimum response (deprecated)
#' @param clust_assign a vector with integers defining cluster membership. 
#' For 3 chemicals, a clust_assign could be c(1,2,1), meaning chemicals 1 and 3 
#' are clustered together and 2 is by itself.  It is assumed that the cluster 
#' assignments start from 1 and do not skip integers.
#'
#' @return an instance of the function mix_effect_fun which takes as input the
#'  concentration of the component chemicals and outputs a predict response
#' @export
#'
#' @examples
mix_function_generator = function(param_matrix, clust_assign){
  
    
  #' mix_effect_fun takes as input the concentration of the component chemicals
  #' and outputs a predict response
  #'
  #' @param true_conc_value a vector with length equal to the number of unique
  #'   chemicals in the mixture and values real non-negative numbers.  Negative
  #'   numbers may be possible but have not been tested
  #'
  #' @return a real number representing the total response of the mixture
  #' @export
  #'
  #' @examples
  mix_effect_fun = function(true_conc_value){
    
    conc_value = true_conc_value
    # NOT USED: adjust dose for effective dose, ligand competition
    clust_ids = unique(clust_assign)
    response_per_cluster = array(0, max(clust_ids))
    r_max = 1
    if(ncol(param_matrix)>3) r_max = max(param_matrix[,4] )
    # for each group, get concentration addition part
    for(cid in 1:max(clust_ids)){
      active_chems = which(clust_assign == cid)
      active_params= param_matrix[active_chems,]
      active_conc = conc_value[active_chems]
      if(all(active_conc==0)) next
      # if a cluster only has 1 value, make sure to still get a matrix
      if(length(active_chems)==1) active_params = matrix(active_params, nrow=1) 
      # get a list of inverse hill functions based on params using a factory method
      hill_inverse_list = apply(active_params, 
                                MARGIN=1, 
                                function(x) do.call(hill_invs_factry,as.list(x)))
      # assume concentration addition and create function to find equivalent dose
      GCA_function = eff_response_opt(hill_inverse_list, active_conc)
      #optim_res is the equivalent dose across all chems, note Log scale
      optim_res = optimize(GCA_function,interval= c(-100,10), 
                           tol = .Machine$double.eps^0.5)
      response_per_cluster[cid] = exp(optim_res$minimum)
      # for the cases when there is a negative sill, check if the mix response
      # should be negative
      GCA_function_neg = eff_response_opt(hill_inverse_list,
                                          active_conc, interval_sign =-1)
      optim_res_neg = optimize(GCA_function_neg,interval= c(-100,10), 
                               tol = .Machine$double.eps^0.5)
      # compare positive vs negative predicted response objective
      if(optim_res_neg$objective<optim_res$objective){
        # if the negative response provides the optimal value, return that
        # prediction
        response_per_cluster[cid] = -exp(optim_res_neg$minimum)
      }
    }
    
    # combine CA with IA
    total_response =r_max*(1-prod(1-response_per_cluster/r_max))
    
    # synergy with IA?  not used
    #z_approx = response_per_cluster/sum(response_per_cluster)
    #total_response_ztrans =qnorm(1-prod(1-response_per_cluster/r_max))-0*prod(z_approx)
    #total_response = r_max*pnorm(total_response_ztrans)
    return(total_response)
  }
  return(mix_effect_fun)
}

# some testing for GCA with antagonists
if(F){
  a1 = .5
  a2 = -1
  theta1 = 1
  theta2 = 1
  test_gca = function(c) {
    c1 = c[1]; c2 = c[2]
    return((c1*a1/theta1 + c2*a2/theta2)/(1+c1/theta1 + c2/theta2))
  }
  
  xy_grid = expand_grid(seq(0,4,by=.1),seq(0,4,by=.1))
  z_val = apply(expand_grid(seq(0,4,by=.1),seq(0,4,by=.1)), MARGIN = 1, test_gca)
  my_tib = cbind(xy_grid, z_val)
  names(my_tib) = c("x", "y","z")
  ggplot(my_tib, aes(x=x, y=y))+geom_tile(aes(fill=z))
}





# 
# assume: more molecules than dock slots
# then: if chem 1 has best score in A matrix, reduce other chem concentrations
# input: 
#c = total concen of each chem
#R = dock_score_vec or molecular weight or some other weighting
#A = ligand-ligand interaction matrix


#' Given input concentrations, suppose some chemicals "deactivate" other
#' chemicals if they have a higher weight (eg molecular weight, docking score,
#' or affinity).  This function will adjust the true concentrations to reflect
#' how the higher weight chemicals will substitute the lower weight chemicals.
#' If the highest weight chemical has high enough concentration, it can dominate
#' and result in a mixture with only one effective chemical.  We can set rules
#' for which chemicals can substitute others using a matrix A_full.
#'
#' @param conc a vector of postive values representing concentrations (eg
#'   micromol/L)
#' @param R a vector of real-valued weights used in the objective; examples
#'   include docking scores or molecular weights
#' @param A_full a binary matrix.  if A_{ij}= 1, chemical i can substitute
#'   chemical j.  The diagonal must be 1.  If A is the Id matrix, the output
#'   adjusted concentrations will be equal to the input.
#'
#' @return a vector of adjusted concentrations, some of which may be 0
#' @export
#'
#' @examples
adjust_concentrations = function(conc, R, A_full){
  require(CVXR)
  # AR_scores = read.csv("Desktop/tox_mixtures/pythonenv/Ligand_Scores_Weights.csv", row.names = 1)
  # AR_scores = AR_scores[complete.cases(AR_scores),]
  # pubchem_AR_scores
  # dock_scores =AR_scores$Dock_Score
  # mol_weight =AR_scores$Molecular_wt
  # 
  # n_chem = nrow(AR_scores)
  # A = diag(n_chem) 
  # A[9,10] = 1; A[10,9] = 1
  # #A= matrix(1, nrow=n_chem, ncol=n_chem)
  # 
  # chem_conc_matr = get_conc_matrix(10)
  # c = matrix(chem_conc_matr[30, 1:n_chem], nrow=n_chem, ncol=1)
  
  # scale conc for numerical accuracy?
  active_chems = which(conc>0)
  n_chem = length(active_chems)
  if(n_chem==1) return (conc)
  
  scale_factor = 1/min(conc[conc>0])
  conc = conc*scale_factor
  active_conc = conc[active_chems]
  active_weight = R[active_chems]
  
  
  A = A_full[active_chems, active_chems]
  K = Variable(rows = n_chem, cols = n_chem)
  j_vec = matrix(1, nrow=n_chem, ncol=1)
  Jo= matrix(1, nrow=n_chem, ncol=n_chem)
  diag(Jo) = -1
  # total concentration conserved
  constraint_1 <-  K%*%j_vec == active_conc
  # allocation bounded by diagonals
  constraint_2 <-  diag(t(K)%*%Jo) <= 0
  # nonnegativity
  constraint_3 <- K >= 0
  # L-L binding respected
  constraint_4 <- K*(1-A)==0
  probl = Problem(objective = Minimize(sum(diag(active_weight)%*%t(K)%*%j_vec)), 
                  constraints = list(constraint_1, 
                                     constraint_2, 
                                     constraint_3, 
                                     constraint_4))
  result = solve(probl)
  if(result$status=="solver_error"){
    warning("Adjusting Conc Failed")
    return(conc/scale_factor)
  }
  Kp = result$getValue(K)
  # Kp
  # diag(Kp)
  # Kp%*%j_vec
  # diag(t(Kp)%*%Jo)
  # sum(t(Kp)%*%R)
  # cbind(conc,new_conc, R)
  new_conc = rep(0, length(conc))
  new_conc[active_chems] = diag(Kp)
  # machine precision issues, zero out small vals
  new_conc[new_conc<1e-8]=0
  new_conc = new_conc/scale_factor
  return(new_conc)
}
# example with tox21 data
# A = diag(n_chem)
# A= matrix(1, nrow=n_chem, ncol=n_chem)
# chem_conc_matr = get_conc_matrix(10)
# c = matrix(chem_conc_matr[30, 1:n_chem], nrow=n_chem, ncol=1)
# adjust_concentrations(c, mol_weight, A)


#--------- Additional plotting stuff -------------- ###########

# Hill inverse testing ####
if(F){
  require(scales)
  a = -.9; b =.5; slope_c = 1.5
  max_R = 2; d=0
  # max_R allows for variable definition of partial agonist
  # rather than assuming partial means < 1 response
  # a: sill parameter
  # b:  EC50 or Disassociation constant
  # c: slope
  # d: minimum effect
  hilly_inverse_test_old =  function(y, c){
    sign_flip = 1
    if(a<0) sign_flip = -1 
    
    # negative a case?
    if(F && a<0){
      # y= sign(a)*a*((sign(a)*b/(x_test) -1)^(-1/slope_c))
      return( sign(a) * b/(1+(-a/y)^c))
    }
    # GCA case:
    if(c==1 && y!= d+a) return(b/(a/(y-d)-1))
    # general edge or undefined case:
    # regular inverse case: valid for d<y<a
    if(y<a && y>d || y<d && c==1 ) return(b/((a/(y-d)-1)^(1/c)))
    # edge case:  y =a is asymptote
    if(y==a ) return(1e10)
    # undefined case:  y> all responses
    if( (y>=max_R  && c!=1))  return(0)
    
    if(y<=0)  return(-b/(1+(-a/y)^c))
    
    
    # partial agonist case, y in [a, max_R]
    return(-2*b-b/((max_R-a)/(max_R-y)-1)^(1/c))
    #partial to match slope = 1 case: no max_R scaling
    #return(-2*b-b/(a/(2*a-y)-1)^(1/c))
  }
  
  hilly_inverse_test = function(y,c){
    # simplifying isnt easy:  y=-.1, a = -1 not equivalent to y=.1, a=1:
    #  sign(a) * y flips around 0, but we flip around a
    # ie y>a if a<0, y<a for a>0 != sign(a) * y <a
    if(y==0) return(0)
    sa =sign(a)
    if(a>0){
      # a>0, y <0 
      if(y<0) return(-b / (1+(-a/y)^c)) # y = (a / (b/x - 1)^(1/c))
      # a>0, y in 0 to a
      if(y<a) return(b/(a/y - 1)^(1/c))
      # a> 0, y in a to 2a
      if(y<2*a) return(-2*b-b/(a/(2*a-y)-1)^(1/c))
      # a>0, y >2a
      return(-2*b + b/(1+(a/(-2*a + y))^c))
    }
    
    if(a<0){
      # y >0 
      if(y>0) return(-b / (1+(-a/y)^c)) # y = (a / (b/x - 1)^(1/c))
      # y in a to 0
      if(y>a) return(b/(a/y - 1)^(1/c))
      #  y<0  in 2a to a
      if(y>2*a) return(-2*b-b/(a/(2*a-y)-1)^(1/c))
      # y <2a
      return(-2*b + b/(1+(a/(-2*a + y))^c))
    }
    
    # 1 
    return(1)
  }
  
  #y_test_old = seq(abs(a)+.02, max_R-.01, by=.001)
  y_test = seq(a+.01, 2*a-.01, by=sign(a)*.001)
  y_test_bey = seq(2*a-.01, 4*a-.01, by=sign(a)*.001)
  y_test_neg = seq(0, -2*a, by=-sign(a)*.001)
  x_test = 1*10^(seq(-8, 8, by=.1)) 
  
  #pdf("RGCA_symmetry_detail.pdf", width = 8 ,height = 5)
  #png("RGCA_symmetry_full_neg.png",width = 8 ,height = 5, units = "in",res = 200)
  x_test = seq(-5,5, by=.01) 
  plot(x_test, a/(1+(b/x_test)^(1)), type = "l", ylim =c(-3,1),#c(-3, 2), #
       lty = 2,col = alpha("black",.95),
       main = "Hill Function Symmetry, Negative Sill", #"GCA Partial Agonist, with Extension", 
       xlab = "Concentration", ylab = "Response")
  
  points(sapply(y_test,function(x) hilly_inverse_test(x,slope_c)),
         y_test,col=5, lwd=3, type = "l")
  points(sapply(y_test_neg,function(x) hilly_inverse_test(x,slope_c)),
         y_test_neg,col=2, lwd=3, type = "l")
  points(sapply(y_test_bey,function(x) hilly_inverse_test(x,slope_c)),
         y_test_bey,col=6, lwd=3, type = "l")
  lines(x_test, a/(1+(b/x_test)^(1)), type = "l", lty = 2)
  lines(x_test, a/(1+(b/x_test)^slope_c),  col =3,  lwd = 3)
  #lines(x_test, -a/(1+(-b/x_test)^slope_c),  col =4,  lwd = 3)
  #lines(x_test, -sign(a)*b*((-sign(a)*a/(x_test) -1)^(-1/slope_c)),  col =4,  lwd = 3)
  #lines(x_test, sign(a)*a*((sign(a)*b/(x_test) -1)^(-1/slope_c)),  col =5,  lwd = 3)
  abline(v = -b, col = "blue", lty=3, lwd = 2)
  abline(h = a, col = "blue", lty=3, lwd = 2)
  text(4, a+sign(a)*.2, paste("Sill: a = ",a))
  text(-b+sign(a)*-.2, -2, paste("-EC50 = ",-b), srt = 90)
  points(0,0, pch=19)
  
  # full symmetry plot
  #text(-b-.2, -.1, paste("-EC50 = ",-b), srt = 90)
  points(-1, 1.8, pch=8)
  points(0,0, pch=19)
  legend("topleft", legend = c("GCA, b=1", 
                               paste("RGCA, b=",slope_c),
                               "RGCA, Local Extension",
                               "RGCA, Reflection",
                               "RGCA, Ext Reflection"), 
         col = c(1, 3,2, 5,6), lty = c(2, 1, 1, 1,1), lwd = c(1,3,3,3,3))
  
  
  # zoomed symmetry plot, with squares
  rect(0, 0, 100, .9, col=alpha("yellow", .5))
  rect(-.5,-100,0,0, col=alpha("orange", .5))
  abline(h =0);abline(v =0)
  legend("bottomright", legend = c("GCA, b=1", 
                                   paste("RGCA, b=",slope_c),
                                   "RGCA, Local Extension"), 
         col = c(1, 3,2, 5,6), lty = c(2, 1, 1, 1,1), lwd = c(1,3,3,3,3))
  
  
  # 
  x_test = seq(-20,10, by=.1) 
  y_1 = seq(0.01, a-.01, by = -.1)
  y_2 = seq(a+.1, 2*a-.01, by = -.1)
  y_3 = seq(2*a+.01,10*a, by = -.1)
  y_4 = seq(-10*a, -.01, by = -.1)
  plot(x_test, a/(1+(b/x_test)^(1)), type = "l", ylim = c(-2,6),
       lty = 2,col = alpha("black",.95), xlim = c(-4,10),
       main ="GCA Partial Agonist, with Extension", 
       xlab = "Concentration", ylab = "Response")
  lines(sapply(y_1,function(x) hilly_inverse_test(x,slope_c)),
        y_1,col=2, lwd=3)
  lines(sapply(y_2,function(x) hilly_inverse_test(x,slope_c)),
        y_2,col=3, lwd=3)
  lines(sapply(y_3,function(x) hilly_inverse_test(x,slope_c)),
        y_3,col=4, lwd=3)
  lines(sapply(y_4,function(x) hilly_inverse_test(x,slope_c)),
        y_4,col=5, lwd=3)
  
}
