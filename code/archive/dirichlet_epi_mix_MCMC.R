# code to fit a tox mixture model given individual dose response curve
# models clustering through the hill slope parameter


#### Generate data ####

hill_function = function(a,b,c,conc) a/(1+(b/conc)^c)
x = seq(0, 3, by=.02)
plot(x, hill_function(1,1,1,x)+1.96*.05,type = "l",col=2)
points(x, hill_function(1,1,1,x))
points(x, hill_function(1,1,1,x)-1.96*.05,type = "l",col=2)
for( i in seq(-.5, .5, by=.05)) points(x, hill_function(1,1,1+i,x), type = "l")
for( i in seq(-.5, .5, by=.05)) points(x, hill_function(1,1+i,1,x), type = "l")


#sample some coefficients: theta_1, theta_2, phi_c
# generate curves at Cx, add noise
set.seed(112)
set.seed(113)

n_chems = 10
# dosages observed
Cx= seq(0,10, by=.4)
n_dose = length(Cx)
#assume 3-4 clusters
theta_1_true = rgamma(n_chems, shape = 2, rate = 2)
theta_2_true = rgamma(n_chems, shape = 2, rate = 1)
phi_c_true = c(2, 1, 3)#rgamma(4, shape = 1.5, rate = .5)
clust_pattern= c(2,5,3)
true_clust_assign = rep(1:length(clust_pattern), clust_pattern)
phi_ci_true = rep(phi_c_true,  clust_pattern)
# now generate data given true mix components
true_mix_hill = mix_function_generator(theta_1_true, 
                                       theta_2_true, 
                                       phi_c_true, 
                                       true_clust_assign)
#sample uniformly?
n_samples = 100
obs_C = matrix(sample(seq(.01,.4, by=.01), n_chems*n_samples, replace = T), nrow = n_samples)
R_mix = apply(obs_C, MARGIN = 1, true_mix_hill)
plot(R_mix)
cor(obs_C, R_mix)
plot(obs_C[,4], R_mix)
plot(obs_C)



#### MCMCM ####
# apply Gibbs to data provided in the Neal 1999 paper:
#y|phi, c ~ N(phi_c, 1)
# c | p ~ Discrete(p_1,...,p_K)
# phi ~ G_0
# p ~ Dir(a/K,...,a/K)


remove_unused_components = function(c_i, phi_c, theta1_c, theta2_c){
  # collect current indices
  valid_phi = 1:length(phi_c)
  valid_ci = unique(c_i)
  # check which phi is no longer being used
  missing_idx = setdiff(valid_phi,valid_ci)
  # if missing, drop and shift remaining
  if(length(missing_idx)>0){
    phi_c = phi_c[-missing_idx]
    W_c = W_c[-missing_idx,]
    theta1_c = theta1_c[-missing_idx]
    theta2_c = theta2_c[-missing_idx]
    # decrement all higher indices by 1
    c_i[c_i>missing_idx] = c_i[c_i>missing_idx]-1
  }
  return(list(c_i, phi_c, theta1_c, theta2_c))
}

generate_possible_clusters(W_c, W_m,j){
  
}



R_obs = R_mix
x_obs = obs_C
compute_llh_W = function(c_i, W_c, theta1_c, theta2_c, phi_c, rescale = F){
  cluster_ids = unique(c_i)
  response_by_clust = matrix(0, ncol = length(cluster_ids), nrow = nrow(x_obs))
  for(i in cluster_ids){
    active_cluster_idx = which(c_i == i)
    W_vec = rep(0, n_chems)
    W_vec[active_cluster_idx] = W_c[active_cluster_idx]
    #rescale W_vec st at least one term is unscaled
    if(rescale){
      W_vec = W_vec/W_vec[which(abs(W_vec-1) == min((abs(W_vec-1))))]
    }
    eff_doses = W_vec%*%x_obs
    
    clust_resp = sapply(eff_doses, FUN = function(x) hill_function(theta1_c[i], 
                                                                   theta2_c[i], 
                                                                   phi_c[i], x))
    response_by_clust[,i] = clust_resp                    
  }
  pred_resp = apply(response_by_clust, MARGIN = 1, FUN = function(x) 1-prod(1-x))
  llh = sum(dnorm(R_obs, mean = pred_resp, sd = sigma, log=T))
  return(llh)
}



# instead:  mixture data, matrix with Rmix for each concentration vector
#y_i = matrix() generated above
n_samples = nrow(y_i)
# init: assume 1 cluster, but generate a set of values for W
c_i = rep(1,n_chems)
# need to allow for matrix or list structure for W
W_c = rnorm(n_chems)
theta_1 = 1
theta_2 = 1
phi_c = 1
#sig_sqr_c = c(1/rgamma(1,3,1))
m=3 # num aux
alpha=.5
sigma = .1
sigma_mh = .1
lambda = 1 # prior variance for W
phi_rate = 1
phi_shape = 1.5
n_iter = 20000
num_sampling_iters = 5
record_assigned_mean = matrix(nrow = n_iter, ncol = n)
record_assigned_clust = matrix(nrow = n_iter, ncol = n)
record_sigma = matrix(nrow= n_iter, ncol = 1)
record_hill_params = matrix(nrow = n_iter, ncol=2*n)
n_iters_all_0 = 0


for(iter in 1:n_iter){
  #update group assignments
  for(j in 1:n){
    c_i_up  = c_i[j]
    remain_c_i=c_i[-j]
    remaining_c_id = sort(unique(remain_c_i))
    k_minus = length(remaining_c_id)
    
    #clusters are denoted by theta1, theta2, and phi
    #within cluster, have regression params W
    #after shifting W, rescale by W max so there is always 1 term with W=1? 
    #  or rescale by value closest to 1: a constraint
    
    
    
    #generate aux: for W, this is just value of 1?
    # NEED TO GENERATE NEW HILL!  Let TEF = W = 1 in this case
    phi_m = rtruncnorm(m, a=0, mean=1, sd=2)
    theta1_m = rtruncnorm(m, a=0, mean=1, sd=2)
    theta2_m = rtruncnorm(m, a=0, mean=1, sd=2)
    #rgamma(m, shape = phi_shape, rate = phi_rate)
    #abs(rnorm(m, mean = 1, sd=.5))
    # reuse cluster id if current idx is lone member
    last_of_group = !(c_i_up %in% remaining_c_id)
    if(last_of_group){
      # phi[c_i]  already in set, but treated as new cluster
      phi_m = phi_m[-1,]
      theta1_m = theta1_m[-1,]
      theta2_m = theta2_m[-1,]
    }
    
    
    #sample c according to likelihood of obs
    n_phis = k_minus+length(phi_m) +last_of_group# same regardless of lone cluster
    potential_phs = c(phi_c, phi_m)
    potential_theta1 = c(theta1_c, theta1_m)
    potential_theta2 = c(theta2_c, theta2_m)
    
    prob_c_vals = rep(0, n_phis)
    for(k in 1:n_phis){
      # k represents how the jth element is clustered
      temp_c_i = c_i
      temp_c_i[j]=k
      kth_llh = compute_llh_W(temp_c_i, W_c, potential_theta1, potential_theta2, potential_phs)
      
      #  DP posterior weights for llh
      n_ic = sum(remain_c_i == k )
      if(n_ic==0) n_ic = alpha/m
      
      prob_c_vals[k] = n_ic/(n-1+alpha)*kth_llh
    }
    prob_c_vals = prob_c_vals/sum(prob_c_vals)
    sampling_set = 1:n_W_rows
    new_c = sample(sampling_set, size = 1, prob = prob_c_vals)
    # if chose self, go to next iter
    if (new_c == c_i_up) next
    # otherwise: chose existing or new phi
    c_i[j]= new_c
    if(new_c>(k_minus+last_of_group)){
      c_i[j] = k_minus + last_of_group + 1
      phi_c = c(phi_c,  potential_phs[new_c])
      theta1_c = c(theta1_c,  potential_theta1[new_c])
      theta2_c = c(theta2_c,  potential_theta2[new_c])
      
    }
    # might have lost a group, update phi
    up_list = remove_unused_components(c_i, phi_c, theta1_c, theta2_c)
    c_i = up_list[[1]]
    phi_c = up_list[[2]]
    theta1_c = up_list[[3]]
    theta2_c = up_list[[4]]
  }
  
  #update phi
  n_phi = length(unique(phi_c))
  for(j in 1:n_phi){
    # consider multiple MH steps
    phi_curr = phi_c[j]
    for(itern in 1:num_sampling_iters){
      #MH:  use proposal distr of N(phi_curr, sigma_mh)
      prop_phi_cj = rtruncnorm(1, a=0, mean = phi_curr , sd = sigma_mh)
      # symmetric proposal, just get alpha from ratio
      total_llh_curr = dgamma(phi_curr, shape = phi_shape, rate =  phi_rate, log=T)+
        log(dtruncnorm(prop_phi_cj, a=0, mean =phi_curr, sd = sigma_mh))
      total_llh_prop = dgamma(prop_phi_cj, shape = phi_shape, rate = phi_rate, log=T)+  
        log(dtruncnorm(phi_curr, a=0, mean =prop_phi_cj, sd = sigma_mh))
      #total_llh_curr = log(dtruncnorm(phi_curr, a=0, mean =1, sd = 1))+
      #total_llh_prop = log(dtruncnorm(prop_phi_cj, a=0, mean =1, sd = 1))+
      cluster_member_idx = which(c_i == j)
      for(clm in cluster_member_idx){
        #sample conditional on the y_i assigned to phi
        curr_mean = theta_1[j]/(1+(theta_2[j]/Cx)^phi_curr)
        llh_curr = sum(dnorm( y_i[clm,], mean = curr_mean, sd = sigma,log = T))
        total_llh_curr = total_llh_curr + llh_curr
        
        prop_mean = theta_1[j]/(1+(theta_2[j]/Cx)^prop_phi_cj)
        llh_prop = sum(dnorm( y_i[clm,], mean = prop_mean, sd = sigma,log = T))
        total_llh_prop = total_llh_prop + llh_prop  
      }
      accept_thresh = exp(total_llh_prop-total_llh_curr)
      if(runif(1)<accept_thresh){
        phi_curr = prop_phi_cj
      }
    }
    phi_c[j] = phi_curr
  }
  
  # update theta_1 (va Gibbs): prior is N(0, 2)
  for( j in 1:n){
    theta1_coeff = 1/(1+(theta_2[j]/Cx)^phi_c[c_i[j]])
    theta1_post_var = 1/(sigma^(-2)*sum(theta1_coeff*theta1_coeff)+2)
    theta_post_mean = theta1_post_var*(0+theta1_coeff%*%y_i[j,]*sigma^(-2))
    theta_1[j] = rnorm(1,mean = theta_post_mean,sd = sqrt(theta1_post_var))
  }
  
  # update theta_2 (via MH)
  for( j in 1:n){
    # for each index, do a few MH steps instead of 1 
    for(itern in 1:num_sampling_iters){
      #theta2_prop = abs(rnorm(1, mean = theta_2[j], sd = sigma_mh))
      theta2_prop = rtruncnorm(1,a=0, mean = theta_2[j], sd = sigma_mh)
      curr_mean = theta_1[j]/(1+(theta_2[j]/Cx)^phi_c[c_i[j]])
      llh_curr = sum(dnorm( y_i[j,], mean = curr_mean, sd = sigma,log = T))+
        dgamma(theta_2[j], shape= 1.5, rate = .5, log=T)+
        log(dtruncnorm(theta2_prop, a=0, mean = theta_2[j], sd = sigma_mh))
      #dnorm(, mean = 0, sd = 1, log=T)
      
      prop_mean = theta_1[j]/(1+(theta2_prop/Cx)^phi_c[c_i[j]])
      llh_prop = sum(dnorm( y_i[j,], mean = prop_mean, sd = sigma,log = T))+ 
        dgamma(theta2_prop, shape = 1.5, rate = .5, log=T)+
        log(dtruncnorm(theta_2[j], a=0, mean = theta2_prop, sd = sigma_mh))
      #dnorm(theta2_prop, mean = 0, sd = 1, log=T)
      accept_p =min(exp(llh_prop-llh_curr),1)
      if(runif(1)<accept_p) theta_2[j] = theta2_prop
    }
  }
  
  
  
  
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
beepr::beep()
#### Plot MCMC ####

plot(record_assigned_mean[seq(1,n_iter,by=10),5])
plot(record_hill_params[seq(1000,n_iter,by=10),15])
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
