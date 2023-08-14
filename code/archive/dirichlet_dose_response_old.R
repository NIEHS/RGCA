# set: theta_1 = a, theta_2=b, phi_c=c
#hilly_reflect_xy =  function(x) b/((a/x)-1)^(1/c)
hill_invs_factry = function(a,b,c){
  hilly_inverse =  function(y){
    if(y<a) return(b/((a/y)-1)^(1/c))
    if(y==a) return(Inf)
    # regular inverse: good for y in (a, 2a)
    #return(-b/((a/(2*a-y))-1)^(1/c)-2*b)#-2*b
    # scaled inverse: good for y in (a, 1)
    return(-b/(((1-a)/(1+y-2*a))-1)^(1/c)-2*b)#-2*b
  }
  return(hilly_inverse)
}


mix_function_generator = function(theta_1, theta_2, phi_c, clust_assign){
  mix_effect_fun = function(conc_value){
    # use b as EC50?
    param_matrix = as.matrix(cbind("a" = theta_1, 
                                   "b" = theta_2, 
                                   "c" = phi_c[clust_assign]))
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
      EC50_1 = active_params[1,1]/2
      total_effective_dose  = 0
      for( indx in 1:length(active_chems)){
        #compute ratios and inverses, using reflection trick
        EC50 = active_params[indx,1]/2
        TEF = hill_inverse_list[[1]](EC50_1)/hill_inverse_list[[indx]](EC50)
        total_effective_dose = total_effective_dose + TEF*active_conc[indx]
      }
      response =  hill_function(active_params[1,1],
                                active_params[1,2],
                                active_params[1,3],
                                total_effective_dose)
      response_per_cluster[cid] = response
    }
    # combine CA with IA
    total_response = 1-prod(1-response_per_cluster)
    return(total_response)
  }
  return(mix_effect_fun)
}

