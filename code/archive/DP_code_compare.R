
## TESTING DP compared to alternateiv
cldat=cluster_centers(cluster_chain, n_top = 100)
clust_arrange = apply(X = cluster_chain$cluster, 
                      MARGIN = 1,
                      FUN = function(x)paste(x, collapse=" "))
occured_clust = rep(0, clust_iter)
occured_clust[which(clust_arrange == names(cldat$assign[1]))] = 1
plot(sapply(100:50000, FUN = function(x) mean(occured_clust[(x-99):x]) ), type="l")

# alternative software
scaled_y = scaled(sort(re_par_list$slope_params))
dp_mod_obj = dirichletprocess::DirichletProcessGaussian(scaled_y, 
                                                        g0Priors = c(0,1,20, .5),
                                                        alphaPriors = c(3/2, 1/4))
dp_fit = dirichletprocess::Fit(dp_mod_obj, its = clust_iter)


cluster_table = sort(table(clust_arrange_alt), decreasing = T)
clust_centers = matrix(0, nrow = n_top, ncol = ncol(record_assigned_clust))
clust_sds = matrix(0, nrow = n_top, ncol = ncol(record_assigned_clust))

