# code to fit a tox mixture model given individual dose response curve
# models clustering through the hill slope parameter

# 
# # quality checking:  
# -Chain stationarity?  Convergence?  Run 50k+ iterations 
# - allow noise term for each curve
# -zero response chemicals, different agonists: own cluster?
#   top of curve forced to cluster, to avoid GCA?
#   -consider estimating a range of curves instead of the right curve
# abstract
# -alpha via BIC, AIC?  more clusters = more parameters
# #- nonparameteric function inverse for GCA?
# ## Done!
# -inverse function correction
# -replicates:  see v4
# -hyperparameter choice:  added prior. 


#### Generate data ####

hill_function = function(a,b,c,conc) a/(1+(b/conc)^c)

plot(x, hill_function(1,1,1,x)+1.96*.05,type = "l",col=2, log="x")
points(x, hill_function(1,1,1,x))
points(x, hill_function(1,1,1,x)-1.96*.05,type = "l",col=2)
for( i in seq(-.5, .5, by=.05)) points(x, hill_function(1,1,1+i,x), type = "l")
for( i in seq(-.5, .5, by=.05)) points(x, hill_function(1,1+i,1,x), type = "l")


#sample some coefficients: theta_1, theta_2, phi_c
# generate curves at Cx, add noise
set.seed(112)

n_chems = 10
n_replicates = rep(2, n_chems)
# dosages observed
Cx= seq(0,20, by=.8)
max_dose = 20
n_dose = length(Cx)
#assume 3-4 clusters
theta_1_true = rnorm(n_chems, mean=1, sd=.3)#rgamma(n_chems, shape = 3, rate = 3)
theta_2_true = rgamma(n_chems, shape = 3, rate = 1)
phi_c_true_vals = c(2, 1, 3)#rgamma(4, shape = 1.5, rate = .5)
clust_pattern= c(2,5,3)
#theta_1_true[11:12] = 0
true_clust_assign = rep(1:length(clust_pattern), clust_pattern)
phi_ci_true = rep(phi_c_true_vals,  clust_pattern)
# now generate curves for each chem

y_i = matrix(nrow = sum(n_replicates), ncol = n_dose+1 )
for(i in 1:n_chems){
  for(r in 1:n_replicates[n_chems]){
    response_data = hill_function(theta_1_true[i], theta_2_true[i],
                                  phi_ci_true[i],Cx) + rnorm(n_dose, sd=.05)
    y_i[i+(r-1)*n_chems,] = c(i, response_data)
  }
}



par(mfrow=c(3,3))
for(i in c(1,11,2,12,3,13))plot(y_i[i,-1])
par(mfrow=c(1,1))

# replicates require chem ID
replicate_sets = c()
for(chem in unique(y_i[,1])){
  rep_idx = which(y_i[,1]==chem)
  replicate_sets = c(replicate_sets, list(rep_idx))
}

# drop the ID col, have replicate set map
y_i=y_i[,-1]


# create matrix to allow for different levels of each chem
Cx = do.call(rbind, rep(list(Cx), sum(n_replicates)))


#### Fit individual curves ####
curve_fits = t(sapply(1:nrow(y_i), FUN=function(x)  drc::drm(y_i[x,]~Cx[x,], fct=drc::LL.3())$coefficients))
curve_fits = curve_fits[,c(2,3,1)]
curve_fits[,3]=-curve_fits[,3]
colnames(curve_fits) = c('a','b','c')



#### Read in data ####
source("Desktop/tox_mixtures/code/tox21_prep_data.R")

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

init_u_RE = function(repl_matrix, take_RE = T){
  max_curve_values = apply(repl_matrix, 
                           MARGIN = 1,
                           function(x) max(x,na.rm=T))
  offset = max_curve_values - max_curve_values[1]
  if(take_RE) return(offset)
  return(median(max_curve_values))
}

# could be: computed slopes from Hill?
#y_i =  c(-1.48, -1.4, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
# instead:  curve data, matrix
#y_i = matrix() generated above
# replicate_sets: map chemical to replicate idx
#n_chems = n = nrow(y_i)
n = n_chems = length(replicate_sets)
n_replicates = unlist(lapply(replicate_sets, length))
n_dose = ncol(Cx)
# init: assume 1 cluster, but generate a value of phi_c
# c_i = rep(1,n)
# phi_c = c(rgamma(1, shape = 2, rate = 2))
# theta_1 = sapply(1:n_chems, FUN = function(x) init_u_RE(y_i[replicate_sets[[x]],], F))
# theta_2 = rep(1e-7,n)
# u_repl = sapply(1:n_chems, FUN = function(x) init_u_RE(y_i[replicate_sets[[x]],]))
# v_repl = sapply(1:n_chems, FUN = function(x) y_i[replicate_sets[[x]],1])#list(rep(0, n_replicates[x])))
c_i = 1:n_chems
phi_c =  rep(1, n_chems)
theta_1 = rep(1, n_chems)
theta_2 = rep(1, n_chems)
u_repl = sapply(1:n_chems, FUN = function(x) list(rep(0, n_replicates[x])))#
v_repl = u_repl 



#sig_sqr_c = c(1/rgamma(1,3,1))
m=3 # num aux
alpha=1
sigma = 1  # noise for obs data
sigma_mh = .1 # proposal distr for MH
sigma_u = rep(10, length(u_repl)) # RE for sill
sigma_v = rep(10, length(v_repl)) # RE for offset
lambda = 1 # prior variance for phi, not used in this context
phi_rate = 1
phi_shape = 1# 1.5
n_iter = 100000
num_sampling_iters = 5
record_assigned_mean = matrix(nrow = n_iter, ncol = n)
record_assigned_clust = matrix(nrow = n_iter, ncol = n)
record_sigma = matrix(nrow= n_iter, ncol = 1)
record_RE_u = matrix(nrow =n_iter, ncol = sum(n_replicates) )
record_RE_v = matrix(nrow =n_iter, ncol = sum(n_replicates) )
record_RE_u_sd = matrix(nrow = n_iter, ncol = length(sigma_u))
record_RE_v_sd = matrix(nrow = n_iter, ncol = length(sigma_v))
record_hill_params = matrix(nrow = n_iter, ncol=2*n+1)
n_iters_all_0 = 0


#Rprof(line.profiling=T)
for(iter in 1:n_iter){
  # in this version of the DP MCMC, the DP is updated after the parameters
  # as a way to initialize parameters to better values
  #update phi
  n_phi = length(unique(phi_c))
  for(j in 1:n_phi){
    # consider multiple MH steps
    phi_curr = phi_c[j]
    for(itern in 1:num_sampling_iters){
      #MH:  use proposal distr of N(phi_curr, sigma_mh)
      prop_phi_cj = truncnorm::rtruncnorm(1, a=0, mean = phi_curr , sd = sigma_mh)
      # if symmetric proposal, just get alpha from ratio
      total_llh_curr = dgamma(phi_curr, shape = phi_shape, rate =  phi_rate, log=T)+
        log(truncnorm::dtruncnorm(prop_phi_cj, a=0, mean =phi_curr, sd = sigma_mh))
      total_llh_prop = dgamma(prop_phi_cj, shape = phi_shape, rate = phi_rate, log=T)+  
        log(truncnorm::dtruncnorm(phi_curr, a=0, mean =prop_phi_cj, sd = sigma_mh))
      #total_llh_curr = log(dtruncnorm(phi_curr, a=0, mean =1, sd = 1))+
      #total_llh_prop = log(dtruncnorm(prop_phi_cj, a=0, mean =1, sd = 1))+
      cluster_member_idx = which(c_i == j)
      for(clm in cluster_member_idx){
        for(rid in 1:n_replicates[clm]){
          act_idx = replicate_sets[[clm]][rid]
          u_RE = u_repl[[clm]][rid]
          v_RE = v_repl[[clm]][rid]
          # concatenate all relevant Conc, responses:  [y_1(C), y_2(C)], etc
          active_X = Cx[act_idx,] #mk_vec(c(t(Cx[act_idx,])))
          active_Y =  y_i[act_idx,] #mk_vec(c(t(y_i[act_idx,])))
          
          #sample conditional on the y_i assigned to phi
          curr_denom = 1+(theta_2[j]/active_X)^phi_curr
          curr_mean = (theta_1[j]+u_RE)/curr_denom+v_RE
          curr_sd_RE = sigma # sqrt(denom^2*sigma_u^2 + sigma_v^2+sigma^2)
          llh_curr = sum(dnorm( active_Y, mean = curr_mean, sd = curr_sd_RE,log = T), 
                         na.rm = T)
          total_llh_curr = total_llh_curr + llh_curr
          
          prop_denom = 1+(theta_2[j]/active_X)^prop_phi_cj
          prop_mean = (theta_1[j]+u_RE)/(prop_denom)+v_RE
          prop_sd_RE =sigma # sqrt(prop_denom^2*sigma_u^2 + sigma_v^2+sigma^2)
          llh_prop = sum(dnorm(active_Y, mean = prop_mean, sd = prop_sd_RE,log = T), 
                         na.rm = T)
          total_llh_prop = total_llh_prop + llh_prop  
        }
      }
      accept_thresh = exp(total_llh_prop-total_llh_curr)
      if(runif(1)<accept_thresh){
        phi_curr = prop_phi_cj
      }
    }
    phi_c[j] = phi_curr
  }

  # update theta_1, u, v (va Gibbs): prior is N(0, 2) and N(0,sig_u^2), likl is (y-aV-uV)sigm^-2(y-av-uV)
  for( j in 1:n){
    act_idx = replicate_sets[[j]]
    # concatenate all relevantresponses:  [y_1(C), y_2(C)], etc
    active_Y =  c(t(y_i[act_idx,]))
    # build design matrix for replicates
    V_ri = lapply(act_idx,FUN = function(x) 1/(1+(theta_2[j]/Cx[x,])^phi_c[c_i[j]]) )
    V_mat = build_replicate_matrix(V_ri)
    
    #u_RE_vec = c(sapply(u_repl[[j]],FUN = function(x) rep(x,n_dose)))
    #v_RE_vec = c(sapply(v_repl[[j]],FUN = function(x) rep(x,n_dose)))
    lh_prec = 1/ sigma^2
    # Post var, full RE: [V prec_RE V+2^{-1}]^{-1}, V = theta1_coeff
    #theta1_post_var = 1/(na.omit(theta1_coeff)%*% diag(na.omit(prec_RE))%*% na.omit(theta1_coeff)+1/2)
    # Post var, no RE:  
    
    lin_terms_post_var =solve(vprod_NA(V_mat, V_mat)*prec_RE + 
                                diag(c(1/2, 
                                       rep(sigma_u[j]^(-2), n_replicates[j]-1),
                                       rep(sigma_v[j]^(-2), n_replicates[j]-1))))
    # XY, RE: [from (XX)^-1 XY] for post mean is post_var * [v*prec_RE*y + prior_mean/prior var]
    #XY_term = qprod_NA(theta1_coeff, prec_RE, active_Y)
    XY_term = vprod_NA(V_mat, active_Y)*prec_RE
    lin_terms_post_mean = vprod_NA(lin_terms_post_var,XY_term)
    lin_term_sample = t(chol(lin_terms_post_var))%*%rnorm(length(lin_terms_post_mean))+lin_terms_post_mean
    
    theta_1[j] = lin_term_sample[1]
    fill_idx = 2:n_replicates[j]
    u_repl[[j]][fill_idx] = lin_term_sample[fill_idx]
    #dont update v_RE yet?
    v_repl[[j]][fill_idx] = lin_term_sample[fill_idx+ n_replicates[j]-1]
    
  }
  
  # update theta_2 (via MH)
  for( j in 1:n){
    act_idx = replicate_sets[[j]]
    # concatenate all relevant Conc, repsonses
    active_X = c(t(Cx[act_idx,]))
    active_Y = c(t(y_i[act_idx,]))
    u_RE_vec = c(sapply(u_repl[[j]],FUN = function(x) rep(x,n_dose)))
    v_RE_vec = c(sapply(v_repl[[j]],FUN = function(x) rep(x,n_dose)))
    # for each index, do a few MH steps instead of 1 
    for(itern in 1:num_sampling_iters){
      #theta2_prop = abs(rnorm(1, mean = theta_2[j], sd = sigma_mh))
      theta2_prop = truncnorm::rtruncnorm(1,a=0, mean = theta_2[j], sd = sigma_mh)
      numerator = theta_1[j]+u_RE_vec
      curr_mean = numerator/(1+(theta_2[j]/active_X)^phi_c[c_i[j]]) + v_RE_vec
      llh_curr = sum(dnorm(active_Y, mean = curr_mean, sd = sigma,log = T), 
                     na.rm = T)+
        dgamma(theta_2[j], shape= 1.5, rate = .25, log=T)+
        log(truncnorm::dtruncnorm(theta2_prop, a=0, mean = theta_2[j], sd = sigma_mh))
      #dnorm(, mean = 0, sd = 1, log=T)
      
      prop_mean = numerator/(1+(theta2_prop/active_X)^phi_c[c_i[j]]) + v_RE_vec
      llh_prop = sum(dnorm(active_Y, mean = prop_mean, sd = sigma,log = T), 
                     na.rm = T)+ 
        dgamma(theta2_prop, shape = 1.5, rate = .25, log=T)+
        log(truncnorm::dtruncnorm(theta_2[j], a=0, mean = theta2_prop, sd = sigma_mh))
      #dnorm(theta2_prop, mean = 0, sd = 1, log=T)
      accept_p =min(exp(llh_prop-llh_curr),1)
      if(runif(1)<accept_p) theta_2[j] = theta2_prop
    }
  }
  
  
  # update sigma 
  p = length(phi_c)
  tot_error = 0
  for(j in 1:n){
    act_idx = replicate_sets[[j]]
    active_X = c(t(Cx[act_idx,]))
    active_Y = c(t(y_i[act_idx,]))
    u_RE_vec = c(sapply(u_repl[[j]],FUN = function(x) rep(x,n_dose)))
    v_RE_vec = c(sapply(v_repl[[j]],FUN = function(x) rep(x,n_dose)))
    numerator = theta_1[j]+u_RE_vec
    denom = 1+(theta_2[j]/active_X)^phi_c[c_i[j]]
    
    sum_sqr_err = sum((active_Y-numerator/denom-v_RE_vec)^2, na.rm = T)
    tot_error= tot_error + sum_sqr_err
  }
  sigma_sqr = 1/rgamma(1, shape = 5+n_dose*sum(n_replicates)/2, 
                       rate = .1 + tot_error/2)
  sigma= sqrt(sigma_sqr)       
        
  
  
  # Update RE sigma_u, sigma_v
  G_beta_prior_u = .1
  G_alpha_prior_u = 5
  G_beta_prior_v = .1
  G_alpha_prior_v = 5
  for(j in 1:n){
    act_idx = replicate_sets[[j]]
    G_beta_post_u = G_beta_prior_u + u_repl[[j]]%*%u_repl[[j]]/2
    G_alpha_post_u = G_alpha_prior_u + n_replicates[j]/2
    sigma_u[j] = sqrt(1/rgamma(1, shape = G_alpha_post_u, rate = G_beta_post_u))
    
    G_beta_post_v = G_beta_prior_v + v_repl[[j]]%*%v_repl[[j]]/2
    G_alpha_post_v = G_alpha_prior_v + n_replicates[j]/2
    sigma_v[j] = sqrt(1/rgamma(1, shape = G_alpha_post_v, rate = G_beta_post_v))
  }
  
  
  
  
  
  ##### DP update####
  for(j in 1:n){
    c_i_up  = c_i[j]
    remain_c_i=c_i[-j]
    remaining_c_id = sort(unique(remain_c_i))
    k_minus = length(remaining_c_id)
    
    #generate aux
    phi_m = rgamma(m, shape = phi_shape, rate = phi_rate)
    #rtruncnorm(m, a=0, mean=1, sd=2)
    #abs(rnorm(m, mean = 1, sd=.5))
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
    #iterate across cluster options for current chemical
    for(k in 1:n_phis){
      n_ic = sum(remain_c_i == k )
      if(n_ic==0) n_ic = alpha/m
      prob_across_repls = rep(0, n_replicates[j])
      for(rpl in 1:n_replicates[j]){
        # get the index in the data for the current chem and replicate
        rp_idx = replicate_sets[[j]][rpl]
        numerator = theta_1[j]+u_repl[[j]][rpl]
        denom = 1+(theta_2[j]/Cx[rp_idx,])^potential_phs[k]
        # mean effect includes estimated random effect
        pred_response = numerator/denom+v_repl[[j]][rpl]
        sd_RE = sigma # full RE var not needed: sqrt(denom^2*sigma_u^2 + sigma_v^2+sigma^2)
        prob_across_repls[rpl] = sum(dnorm( y_i[rp_idx,], 
                                            mean = pred_response, 
                                            sd = sd_RE,
                                            log = T), 
                                     na.rm = T)
      }
      prob_c_vals[k] = n_ic/(n-1+alpha)*exp(sum(prob_across_repls))
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
  
  
  #update the DP concentration alpha
  eta = rbeta(n=1, alpha+1, n)
  # prior on alpha ~ G(a= 3/2, b= 1/4)
  alpha_a = 3/2
  alpha_b = 1/2
  # k = n_clusts = p 
  mix_prob = (alpha_a+k-1)/(alpha_a+k-1+n*(alpha_b-log(eta)))
  # the posterior is a mixture; sample given eta
  alpha =  ifelse(runif(1)<mix_prob, 
                  rgamma(1,shape = alpha_a+p, rate = alpha_b-log(eta)),
                  rgamma(1, shape = alpha_a+p-1, rate = alpha_b-log(eta)))
  
  
  
  # relabel to avoid identical groupings
  uniq_labs = unique(c_i)
  c_i = order(uniq_labs)[c_i]
  phi_c = phi_c[uniq_labs]
  
  record_RE_u[iter,] = mk_vec(u_repl)
  record_RE_v[iter,] = mk_vec(v_repl)
  record_RE_u_sd[iter,] = sigma_u
  record_RE_v_sd[iter,] = sigma_v
  record_sigma[iter] =  sigma   
  record_assigned_mean[iter, ] = phi_c[c_i]
  record_assigned_clust[iter,]= c_i
  record_hill_params[iter,] = c(alpha, theta_1, theta_2)
}

#Rprof(NULL)
#summaryRprof("Rprof.out")


beepr::beep()


#rearrange u and v to match replicate idx
chem_sorted_repl = mk_vec(u_repl)
repl_idx = unlist(replicate_sets)
u_RE_sort = rep(1,sum(n_replicates))
u_RE_sort[repl_idx] = chem_sorted_repl
#### Plot MCMC ####

plot(record_hill_params[seq(1,n_iter,by=10),2])
plot(record_assigned_mean[seq(1,n_iter,by=10),5])
plot(record_hill_params[seq(1000,n_iter,by=10),3])
plot(record_sigma[seq(1000, n_iter, by=50)])

pori = par(mfrow=c(3,2)); for(i in 1:6) plot(record_RE_u[seq(1000, n_iter, by=20),i])
par(pori)
plot(record_RE_u[seq(1000, n_iter, by=50),14])#sum(n_replicates[1:4])+5])
plot(record_RE_u_sd[seq(1000, n_iter, by=10),4])
plot(record_RE_v[seq(1000, n_iter, by=10),16])
plot(record_RE_v_sd[seq(1000, n_iter, by=10),4])


#compare RE and parameter estimate
#compare RE and parameter estimate
u_all =  colMeans(record_RE_u[seq(1000, n_iter, by=10),])
repl_idx = unlist(replicate_sets)
u_RE_sort = rep(1,sum(n_replicates))
u_RE_sort[repl_idx] = u_all
cbind(u_RE_true, u_RE_sort)

n_clust_over_time = apply(X = record_assigned_clust, MARGIN = 1,max)
hist(n_clust_over_time)



clust_arrange = apply(X = record_assigned_clust, 
                      MARGIN = 1,
                      FUN = function(x)paste(x, collapse=" "))

top_clusters = sort(table(clust_arrange), decreasing = T)[1:5];top_clusters
top_clust = names(top_clusters)[1]
top_clust_idx = which(clust_arrange == top_clust)
phi_c_est_vec = colMeans(record_assigned_mean[top_clust_idx, ])
estimated_thetas = colMeans(record_hill_params[top_clust_idx,])
theta1_est = estimated_thetas[1:length(replicate_sets)+1];# cbind(theta1_est, theta_1_true)
theta2_est = estimated_thetas[1:length(replicate_sets)+1+length(replicate_sets)];# cbind(theta2_est, theta_2_true)
est_clust_assign  = as.numeric(strsplit(top_clust, split = " ")[[1]])

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




