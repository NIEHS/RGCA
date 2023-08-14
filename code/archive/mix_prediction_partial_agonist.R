#### Plot Mix ####

# set: theta_1 = a, theta_2=b, phi_c=c
#hilly_reflect_xy =  function(x) b/((a/x)-1)^(1/c)
hill_invs_factry = function(a,b,c){
  hilly_inverse =  function(y){
    if(y<a) return(b/((a/y)-1)^(1/c))
    if(y==a) return(0)
    # regular inverse: good for y in (a, 2a)
    #return(-b/((a/(2*a-y))-1)^(1/c)-2*b)#-2*b
    # scaled inverse: good for y in (a, 1)
    return(-2*b-b/((1-a)/(1-y)-1)^(1/c))
    
  }
  return(hilly_inverse)
}

# given the concentrations, optim to get the total response
eff_response_opt = function(hill_inverse_list, conc_vec){
  GCA_over_list = function(r){
    invs_list = lapply(hill_inverse_list, 
                       FUN = function(x) do.call(x, list(r)))
    invs_vec = unlist(invs_list)
    return(abs(sum(conc_vec/invs_vec) -1))
  }
  
  return(GCA_over_list)
}

mix_function_generator = function(param_matrix, clust_assign){
  mix_effect_fun = function(conc_value){
    clust_ids = unique(clust_assign)
    response_per_cluster = array(0, max(clust_ids))
    # for each group, get concentration addition part
    for(cid in 1:max(clust_ids)){
      active_chems = which(clust_assign == cid)
      active_params= param_matrix[active_chems,]
      active_conc = conc_value[active_chems]
      # get a list of inverse hill functions based on params
      if(length(active_chems)==1) active_params = matrix(active_params, nrow=1) 
      hill_inverse_list = apply(active_params, 
                                MARGIN=1, 
                                function(x) do.call(hill_invs_factry,as.list(x)))
      # use the first element of group as the baseline for total effect
      GCA_function = eff_response_opt(hill_inverse_list, active_conc)
      optim_res = optimize(GCA_function,interval= c(0,1))
      response_per_cluster[cid] = optim_res$minimum
    }
    # combine CA with IA
    total_response = 1-prod(1-response_per_cluster)
    return(total_response)
  }
  return(mix_effect_fun)
}

if(F){
phi_c_est = phi_c_est_vec[sapply(unique(est_clust_assign),function(x) which(est_clust_assign == x)[1])]

est_param_matrix = as.matrix(cbind("a" = theta1_est, 
                                   "b" = theta2_est, 
                                   "c" = phi_c_est[est_clust_assign]))
true_param_matrix = as.matrix(cbind("a" = theta_1_true, 
                                    "b" = theta_2_true, 
                                    "c" = phi_c_true_vals[true_clust_assign]))

est_mix_hill = mix_function_generator(est_param_matrix, est_clust_assign)
true_mix_hill = mix_function_generator(true_param_matrix, true_clust_assign)
# compare to individual curve fits with IA and CA
IA_mix_hill = mix_function_generator(curve_fits, 1:n_chems)
CA_mix_hill = mix_function_generator(curve_fits,rep(1, n_chems) )
true_mix_hill(rep(.1, n_chems))



  #plot some curves as individual chem dosages vary
  ref_point = rep(1e-5, n_chems)
  test_point = rep(.1, n_chems)
  dose_range = seq(1e-5,40, by=.05)
  evaluate_DR_at_point = function(idx, baseline=ref_point, mix_fun = true_mix_hill){
    function(x){baseline[idx]=x; mix_fun(baseline)}
  }
  
  mix_effect_2 = sapply(dose_range, evaluate_DR_at_point(2,test_point))
  mix_effect_10 = sapply(dose_range, evaluate_DR_at_point(10,test_point))
  mix_effect_6 = sapply(dose_range, evaluate_DR_at_point(6,test_point))
  plot(dose_range, mix_effect_2)
  points(dose_range, mix_effect_10, col=2)
  points(dose_range, mix_effect_6, col=3)
  
  # reference effect for individual chem
  mix_effect_3 = sapply(dose_range, evaluate_DR_at_point(3,test_point))
  ref_effect_3 = sapply(dose_range, evaluate_DR_at_point(3,ref_point))
  hill_effect_3 = sapply(dose_range, FUN =function(x){ hill_function(theta_1_true[3],
                                                                     theta_2_true[3],
                                                                     phi_ci_true[3],x)})
  plot(dose_range, mix_effect_3, ylim =c(0,1))
  points(dose_range, hill_effect_3, col=3)
  points(dose_range, ref_effect_3, col=2)
  points(Cx,y_i[3,],pch=1)
  
  # Notice: reference effect messed up for subsequent elements in cluster 
  idx_to_plot = 4
  mix_effect_4 = sapply(dose_range,  evaluate_DR_at_point(idx_to_plot,test_point))
  estmix_effect_4 = sapply(dose_range,  evaluate_DR_at_point(idx_to_plot,test_point, est_mix_hill))
  #ref_effect_4 = sapply(dose_range, evaluate_DR_at_point(idx_to_plot,ref_point))
  #est_ref_effect_4 = sapply(dose_range, evaluate_DR_at_point(idx_to_plot,ref_point, est_mix_hill))
  IA_effect = sapply(dose_range,  evaluate_DR_at_point(idx_to_plot,test_point, IA_mix_hill))
  CA_effect = sapply(dose_range,  evaluate_DR_at_point(idx_to_plot,test_point, CA_mix_hill))
  
  hill_effect_4 = sapply(dose_range, FUN =function(x){ hill_function(theta_1_true[idx_to_plot],
                                                                     theta_2_true[idx_to_plot],
                                                                     phi_ci_true[idx_to_plot],x)})
  
  log_dose_range = log(dose_range)
  plot(log_dose_range, mix_effect_4, type = "l", lwd = 3, ylim = c(.1,.5),
       main = paste("Effect as Chemical varies:", idx_to_plot),
       xlab = "Log Dose", ylab = "Effect")
  #points(dose_range, ref_effect_4, col=2)
  #points(Cx,y_i[idx_to_plot,])
  points(log_dose_range, estmix_effect_4, col=3, type = "l", lwd = 3)
  points(log_dose_range, CA_effect, col=4, type = "l",lty = 2, lwd=2)
  points(log_dose_range, IA_effect, col=2, type = "l", lty = 2, lwd=2)
  
  #confirmed, hill matches ref
  #points(log_dose_range, hill_effect_3, col=4,pch=2)
  #points(Cx,y_i[3,],pch=2)
  legend("topleft", legend = c("True Mix-4, 0.1 base", "DP Mix-4, 0.1 base", "CA, 0.1 base", "IA"),
         col =  c(1,3,4,2), pch = c(NA,NA,NA,NA), lty = c(1,1,2,2), lwd = c(3,3,2,2))
  
  plot(dose_range, mix_effect_6, ylim =c(0,1));  points(dose_range, ref_effect_6, col=3)
  plot(dose_range, mix_effect_10, ylim =c(0,1));  points(dose_range, ref_effect_10, col=3)

  
  #### Simulation Study ####
  hill_true_inverse_list = apply(true_param_matrix, 
                                 MARGIN=1, 
                                 function(x) do.call(hill_invs_factry,as.list(x)))
  hill_est_inverse_list = apply(estim_param_matrix, 
                                MARGIN=1, 
                                function(x) do.call(hill_invs_factry,as.list(x)))
  hill_drc_inverse_list = apply(curve_fits, 
                                MARGIN=1, 
                                function(x) do.call(hill_invs_factry,as.list(x)))
  
  
  
  
  unlist(lapply(hill_drc_inverse_list, FUN = function(x) do.call(x, list(.5))))
  # ideas:
  #simulate constant effects mixtures: EC10, EC20, etc
  
  #recover EC vector given effect percent:
  
  ECx_func= function(EC_level) unlist(lapply(hill_true_inverse_list, 
                                             FUN = function(x) do.call(x, list(EC_level))))
  # create a series of points based on increasing EC:
  EC_levels = seq(.01,1, length.out = 20)
  EC_conc = t(sapply(EC_levels, ECx_func))
  plot(EC_levels, apply(EC_conc, MARGIN=1, true_mix_hill))
  points(EC_levels, apply(EC_conc, MARGIN=1, est_mix_hill),col=2)
  
  
  
  # simulate constant proportions: all at concentration c1, c2, etc (most toxic dominates)
  prop_levels = seq(0.1, 2, by=.1)
  EC_props = t(sapply(prop_levels, FUN=function(x) x*rep(1, n_chems)))
  plot(prop_levels, apply(EC_props, MARGIN=1, true_mix_hill))
  # Larger number of chemicals
  # different SNR cases?  
  
  
  #### Spline inverse testing ####
  x = seq(0, 2, by=.02)
  bet = .5
  my_fun_den = function(z) z^(1/bet)
  plot(x, my_fun_den(x), type = "l");points(x, x/bet+1-(1/bet),col=2,type = "l")
  
  xkn = c(0,.5, 1, 1.5, 2)
  n_spline = length(xkn)-1
  y_knots_den = my_fun_den(xkn)
  knt_slp_den = diff(y_knots_den)/diff(xkn)
  knt_b_den = y_knots_den[1:n_spline]-knt_slp_den*xkn[1:n_spline]
  for(i in 1:n_spline) {
    knot_xs = seq(xkn[i], xkn[i+1],by=.05)
    points(knot_xs, knot_xs*knt_slp_den[i]+knt_b_den[i],col=3)
  }
  
  ### Apply Splines to the numerator
  x = seq(0, 2, by=.02)
  bet = .5; sill_param = .9
  get_numerator_spline_params  = function(bet, sill){
    my_fun = function(z) (sill_param-z)^(1/bet)
    y_knots = my_fun(xkn)
    knt_slp = diff(y_knots)/diff(xkn)
    knt_b = y_knots[1:n_spline]-knt_slp*xkn[1:n_spline]
    return(list(knt_slp, knt_b))
  }
  
  
  spline_param_num = get_numerator_spline_params(bet, sill_param)
  knt_slp=spline_param_num[[1]]
  knt_b = spline_param_num[[2]]
  # checl the taylor approx of numerator:  seems bad
  plot(x, my_fun(x), type = "l");points(x, sill_param-(1/bet)*sill_param^(1/bet-1)*x,col=2,type = "l")
  for(i in 1:n_spline) {
    knot_xs = seq(xkn[i], xkn[i+1],by=.05)
    points(knot_xs, knot_xs*knt_slp[i]+knt_b[i],col=3)
  }
  
  #compare spline inverse for some concentrations?
  theta_mix = c(1,1)
  spline_param_num1 = get_numerator_spline_params(bet, sill = .9)
  knt_slp1=spline_param_num1[[1]]; knt_b1 = spline_param_num1[[2]]
  spline_param_num2 = get_numerator_spline_params(bet, sill = 1.1)
  knt_slp2=spline_param_num2[[1]]; knt_b2 = spline_param_num2[[2]]
  
  R_mix_spline = function(c1, c2, knot_idx){
    
    r_mix_est = (c1*knt_b1[knot_idx]/theta_mix[1]+
                   c2*knt_b2[knot_idx]/theta_mix[1]-
                   knt_b_den[knot_idx])/
      (knt_slp[knot_idx]-
         c1*knt_slp1[knot_idx]-
         c2*knt_slp2[knot_idx])
    return(r_mix_est)
  }
  par(mfrow=c(1,1))
  plot(x, R_mix_spline(.3, x, 3)); for(i in 2:4) points(x, R_mix_spline(.3, x, i),col=i)
  
}