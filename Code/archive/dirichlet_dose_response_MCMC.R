# code to fit a tox mixture model given individual dose response curve
# models clustering through the hill slope parameter


#### Generate data ####

hill_function = function(a,b,c,conc) a/(1+(b/conc)^c)

#sample some coefficients: theta_1, theta_2, phi_c
# generate curves at Cx, add noise
set.seed(112)
set.seed(113)

n_chems = 10
# dosages observed
Cx= seq(0,10, by=.4)
n_dose = length(Cx)
#assume 3-4 clusters
theta_1_true = rnorm(n_chems, 1, .3)
theta_2_true = rnorm(n_chems, 3, 1)
phi_c_true = rnorm(4, 1, .5)
clust_pattern= c(2,4,3,1)
true_clust_assign = rep(1:length(clust_pattern), clust_pattern)
phi_ci_true = rep(phi_c_true,  clust_pattern)
# now generate curves for each chem

y_i = matrix(nrow = n_chems, ncol = n_dose )
for(i in 1:n_chems){
  y_i[i,]= hill_function(theta_1_true[i],
                         theta_2_true[i],
                         phi_ci_true[i],Cx) + rnorm(n_dose, sd=.05)
}

par(mfrow=c(3,3))
for(i in 1:9)plot(y_i[i,])






#### MCMCM ####
# apply Gibbs to data provided in the Neal 1999 paper:
#y|phi, c ~ N(phi_c, 1)
# c | p ~ Discrete(p_1,...,p_K)
# phi ~ G_0
# p ~ Dir(a/K,...,a/K)


remove_unused_phi = function(c_i, phi_c){
  # collect current indices
  valid_phi = 1:length(phi_c)
  valid_ci = unique(c_i)
  # check which phi is no longer being used
  missing_idx = setdiff(valid_phi,valid_ci)
  # if missing, drop and shift remaining
  if(length(missing_idx)>0){
    phi_c = phi_c[-missing_idx]
    # decrement all higher indices by 1
    c_i[c_i>missing_idx] = c_i[c_i>missing_idx]-1
  }
  return(list(c_i, phi_c))
}

# could be: computed slopes from Hill?
#y_i =  c(-1.48, -1.4, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
# instead:  curve data, matrix
#y_i = matrix() generated below
n = nrow(y_i)
# init: assume 1 cluster, but generate a value of phi_c
c_i = rep(1,n)
phi_c = c(rnorm(1, mean = 1, sd = .3))
#sig_sqr_c = c(1/rgamma(1,3,1))
m=3 # num aux
alpha=3
sigma = 1
sigma_mh = .1
lambda = 1
n_iter = 20000
num_sampling_iters = 10
record_assigned_mean = matrix(nrow = n_iter, ncol = n)
record_assigned_clust = matrix(nrow = n_iter, ncol = n)
record_sigma = matrix(nrow= n_iter, ncol = 1)
record_hill_params = matrix(nrow = n_iter, ncol=2*n)


for(iter in 1:n_iter){
  #update group assignments
  for(j in 1:n){
    c_i_up  = c_i[j]
    remain_c_i=c_i[-j]
    remaining_c_id = sort(unique(remain_c_i))
    k_minus = length(remaining_c_id)
    
    #generate aux
    phi_m = abs(rnorm(m, mean = 1, sd=.5))
    # reuse cluster id if current idx is lone member
    last_of_group = !(c_i_up %in% remaining_c_id)
    if(last_of_group){
      # phi[c_i]  already in set, but treated as new cluster
      phi_m = phi_m[-1]
    }
    
    #sample c according to likelihood of obs
    n_phis = k_minus+length(phi_m) +last_of_group# same regardless of lone cluster
    potential_phs = c(phi_c, phi_m)
    prob_c_vals = rep(0, n_phis)
    for(k in 1:n_phis){
      n_ic = sum(remain_c_i == k )
      if(n_ic==0) n_ic = alpha/m
      pred_response = theta_1[j]/(1+(theta_2[j]/Cx)^potential_phs[k])
      prob_c_vals[k] = n_ic/(n-1+alpha)*prod(dnorm( y_i[j,], mean = pred_response, sd = sigma))
    }
    prob_c_vals = prob_c_vals/sum(prob_c_vals)
    sampling_set = 1:n_phis
    new_c = sample(sampling_set, size = 1, prob = prob_c_vals)
    # if chose self, go to next iter
    if (new_c == c_i_up) next
    # otherwise: chose existing or new phi
    c_i[j]= new_c
    if(new_c>(k_minus+last_of_group)){
      c_i[j] = k_minus + last_of_group + 1
      phi_c = c(phi_c,  potential_phs[new_c])
    }
    # might have lost a group, update phi
    up_list = remove_unused_phi(c_i, phi_c)
    c_i = up_list[[1]]
    phi_c = up_list[[2]]
  }
  
  #update phi
  n_phi = length(unique(phi_c))
  for(j in 1:n_phi){
    # consider multiple MH steps
    #MH:  use proposal distr of N(phi_curr, sigma_mh)
    prop_phi_cj = abs(rnorm(1,  phi_c[j] , sd = sigma_mh))
    # symmetric proposal, just get alpha from ratio
    total_llh_curr = 0
    total_llh_prop = 0
    cluster_member_idx = which(c_i == j)
    for(clm in cluster_member_idx){
      #sample conditional on the y_i assigned to phi
      curr_mean = theta_1[j]/(1+(theta_2[j]/Cx)^phi_c[j])
      llh_curr = sum(dnorm( y_i[clm,], mean = curr_mean, sd = sigma,log = T))
      total_llh_curr = total_llh_curr + llh_curr
      
      prop_mean = theta_1[j]/(1+(theta_2[j]/Cx)^prop_phi_cj)
      llh_prop = sum(dnorm( y_i[clm,], mean = prop_mean, sd = sigma,log = T))
      total_llh_prop = total_llh_prop + llh_prop  
    }
    total_llh_curr = total_llh_curr+dnorm(phi_c[j], mean=1, sd = .3, log=T)
    total_llh_prop=total_llh_prop + dnorm(prop_phi_cj, mean=1, sd = .3, log=T)
    accept_thresh = exp(total_llh_prop-total_llh_curr)
    if(runif(1)<accept_thresh){
      phi_c[j] = prop_phi_cj
    }
  }
  
  # update theta_1 (va Gibbs)
  for( j in 1:n){
    theta1_coeff = 1/(1+(theta_2[j]/Cx)^phi_c[c_i[j]])
    theta1_post_var = 1/(sigma^(-2)*sum(theta1_coeff*theta1_coeff)+1)
    theta_post_mean = theta1_post_var*(.01+theta1_coeff%*%y_i[j,]*sigma^(-2))
    theta_1[j] = rnorm(1,mean = theta_post_mean,sd = sqrt(theta1_post_var))
  }
  
  # update theta_2 (via MH)
  for( j in 1:n){
    # for each index, do a few MH steps instead of 1 
    for(itern in 1:num_sampling_iters){
      theta2_prop = abs(rnorm(1, mean = theta_2[j], sd = sigma_mh))
      curr_mean = theta_1[j]/(1+(theta_2[j]/Cx)^phi_c[c_i[j]])
      llh_curr = sum(dnorm( y_i[j,], mean = curr_mean, sd = sigma,log = T))+
        dnorm(theta_2[j], mean = 0, sd = 1, log=T)
      
      prop_mean = theta_1[j]/(1+(theta2_prop/Cx)^phi_c[c_i[j]])
      llh_prop = sum(dnorm( y_i[j,], mean = prop_mean, sd = sigma,log = T))+ 
        dnorm(theta2_prop, mean = 0, sd = 1, log=T)
      accept_p =min(exp(llh_prop-llh_curr),1)
      if(runif(1)<accept_p) theta_2[j] = theta2_prop
    }
  }
  
  
  
  
  # update sigma
  # p = length(phi_c)
  # sigma_sqr = 1/rgamma(1, shape = 3+n/2+p/2, 
  #                      scale = 1 + sum((y_i-phi_c[c_i])^2)/2+
  #                        sum(phi_c^2)*lambda/2)
  # sigma= sqrt(sigma_sqr)       
  # record_sigma[iter] =  sigma  
  
  # update sigma
  p = length(phi_c)
  tot_error = 0
  for(j in 1:n){
    sum_sqr_err = sum((y_i[j,]-theta_1[j]/(1+(theta_2[j]/Cx)^phi_c[c_i[j]]))^2)
    tot_error= tot_error + sum_sqr_err
  }
  sigma_sqr = 1/rgamma(1, shape = 5+n_dose*n/2, 
                       rate = .1 + tot_error/2)
  sigma= sqrt(sigma_sqr)       
  record_sigma[iter] =  sigma         
  
  
  
  
  
  # relabel to avoid identical groupings
  uniq_labs = unique(c_i)
  c_i = order(uniq_labs)[c_i]
  phi_c = phi_c[uniq_labs]
  
  record_assigned_mean[iter, ] = phi_c[c_i]
  record_assigned_clust[iter,]= c_i
  record_hill_params[iter,] = c(theta_1, theta_2)
}


#### Plot MCMC ####

plot(record_assigned_mean[seq(1,n_iter,by=10),2])
plot(record_hill_params[seq(1000,n_iter,by=10),12])
plot(record_sigma[seq(1000, n_iter, by=10)])

mu_post = colMeans(record_assigned_mean)
tail(record_assigned_clust)
cbind(mu_post, y_i)

n_clust_over_time = apply(X = record_assigned_clust, MARGIN = 1,max)
hist(n_clust_over_time)

clust_arrange = apply(X = record_assigned_clust, 
                      MARGIN = 1,
                      FUN = function(x)paste(x, collapse=" "))

top_clusters = sort(table(clust_arrange), decreasing = T)[1:5];top_clusters
top_clust = names(top_clusters)[1]
top_clust_idx = which(clust_arrange == top_clust)
phi_c_est = colMeans(record_assigned_mean[top_clust_idx, ])
estimated_thetas = colMeans(record_hill_params[top_clust_idx,])
theta1_est = estimated_thetas[1:10]
theta2_est = estimated_thetas[11:20]
est_clust_assign  = as.numeric(strsplit(top_clust, split = " ")[[1]])

#### Plot Mix ####

# set: theta_1 = a, theta_2=b, phi_c=c
#hilly_reflect_xy =  function(x) b/((a/x)-1)^(1/c)
hill_invs_factry = function(a,b,c){
  hilly_inverse =  function(y){
    if(y<a) return(b/((a/y)-1)^(1/c))
    if(y==a) return(1)
    return(-b/((a/(2*a-y))-1)^(1/c)-2*b)#-2*b
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


est_mix_hill = mix_function_generator(theta1_est, theta2_est, phi_c_est, est_clust_assign)
true_mix_hill = mix_function_generator(theta_1_true, theta_2_true, phi_c_true, true_clust_assign)
true_mix_hill(rep(.1, n_chems))

#plot some curves as individual chem dosages vary
ref_point = rep(0, n_chems)
test_point = rep(.05, n_chems)
dose_range = seq(0,10, by=.05)
evaluate_DR_at_point = function(idx, baseline=ref_point, mix_fun = true_mix_hill)  function(x){baseline[idx]=x; mix_fun(baseline)}

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
points(Cx,y_i[2,],pch=1)

# Notice: reference effect messed up for subsequent elements in cluster 
mix_effect_4 = sapply(dose_range,  evaluate_DR_at_point(4,test_point))
estmix_effect_4 = sapply(dose_range,  evaluate_DR_at_point(4,test_point, est_mix_hill))
ref_effect_4 = sapply(dose_range, evaluate_DR_at_point(4,ref_point))
hill_effect_4 = sapply(dose_range, FUN =function(x){ hill_function(theta_1_true[4],
                                                                   theta_2_true[4],
                                                                   phi_ci_true[4],x)})
plot(dose_range, mix_effect_4, type = "l", lwd = 3)
points(dose_range, estmix_effect_4)
points(dose_range, hill_effect_4, col=3)
points(Cx,y_i[4,])
points(dose_range, ref_effect_4, col=2)
points(dose_range, hill_effect_3, col=4,pch=2)
points(Cx,y_i[3,],pch=2)
legend("topleft", legend = c("True Mix-4, 0.05 base", "Est Mix-4, 0.05 base", "Mix-4, 0 Base", "Indv-4", "Indv-3"),
       col =  c(1,1,2,3,4), pch = c(NA,1,1,1,2), lty = c(1,NA, NA, NA, NA))



plot(dose_range, mix_effect_6, ylim =c(0,1));  points(dose_range, ref_effect_6, col=3)
plot(dose_range, mix_effect_10, ylim =c(0,1));  points(dose_range, ref_effect_10, col=3)

