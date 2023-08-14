library(nimble, warn.conflicts = T)


################################################################
# -----------  nimble: mixed effect -------------
################################################################

# -----------  samples + replicates  -------------
code <- nimbleCode({
  #CM <- ncol(replicate_matrix)
  for(chm in 1:CM){
    a1[chm] ~ dnorm(0, sd = sqrt(1000))
    theta1[chm] ~ dgamma(shape = 1e-3, rate = 1e-3)
    beta1[chm] ~ dgamma(shape = .3, rate = .3)
    # give each chemical a noise term
    sigma_eps[chm] ~ dinvgamma(shape = 0.1, rate = 0.1)
    sigma_u[chm] ~ dinvgamma(shape = 0.1,rate = 0.1)# dunif(0, 100)
    sigma_v[chm] ~ dinvgamma(shape = 0.1,rate = 0.1)# dunif(0, 100)
    for(j in 1:n_reps[chm]){
      # j is the relative index for replicates
      #replicate_matrix[j, chm]: gets the correct absolute index for the replicate
      u[replicate_matrix[j, chm]] ~ dnorm(mean = 0, sd = sigma_u[chm] *
                                            identifiability_constr[j, chm])
      v[replicate_matrix[j, chm]] ~ dnorm(mean = 0, sd = sigma_v[chm]* 
                                            identifiability_constr[j, chm])
      for (i in 1: y_lengths[replicate_matrix[j, chm]]) {
        f[j, i, chm] <-  (a1[chm]+ u[replicate_matrix[j, chm]])/
          (1+(theta1[chm] / x[replicate_matrix[j, chm],i])^beta1[chm]) +
          v[replicate_matrix[j, chm]]
        y[replicate_matrix[j, chm], i] ~ dnorm(f[j,i, chm], sd = sigma_eps[chm])
      }
    }
  }
})

## constants, data, and initial values
nchm = length(replicate_sets)
replicate_matrix = matrix(unlist(lapply(replicate_sets, function(x) {a = rep(0, 6); a[1:length(x)]=x; return(a)})), nrow = 6)
n_reps = unlist(lapply(replicate_sets, length))
# for random effects, fix first replicate to 0 my setting variance to 0
identifiability_constr = replicate_matrix*0+1
# 0 doesnt work, use small value instead
identifiability_constr[1,] =1e-6
y_lengths = apply(y_i, MARGIN = 1, FUN=function(x) ifelse(any(is.na(x)), min(which(is.na(x)))-1, 15) )
constants <- list(#N = ncol(y_i),
  y_lengths = y_lengths,
  replicate_matrix=replicate_matrix,
  n_reps = n_reps,
  CM = ncol(replicate_matrix),
  identifiability_constr = identifiability_constr)

data <- list(
  y = y_i,
  x =Cx
)

inits <- list(beta1 = rep(1,nchm), 
              theta1 = rep(1e-6, nchm), 
              a1 = rep(10,nchm), 
              sigma_eps = rep(1,nchm),
              sigma_u = rep(1, nchm),
              sigma_v = rep(1, nchm),
              u = rep(0, nrow(y_i)),
              v = rep(0, nrow(y_i)))
drcModel <- nimbleModel(code = code, constants = constants, data = data, 
                        inits = inits, check = FALSE)

## Ensure we have the nodes needed to simulate new datasets
dataNodes_drc <- drcModel$getNodeNames(dataOnly = TRUE)
parentNodes_drc <- drcModel$getParents(dataNodes_drc, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes_drc <- drcModel$getDependencies(parentNodes_drc, self = FALSE)
c_drcmodel  <- compileNimble(drcModel)
mcmc    <- buildMCMC(c_drcmodel, monitors = parentNodes_drc)
cmcmc   <- compileNimble(mcmc, project = drcModel)
samples <- runMCMC(cmcmc, niter = 50000, nburnin = 10000, thin = 20)
t(colMeans(samples))
plot(samples[,3])




# -----------  DP with slopes  -------------

codeDP <- nimbleCode({
  sigsqrd ~ dinvgamma(shape = 10, rate = .1)    # var for each cluster xi[i]
  for(i in 1:nchm) {
    beta1[i] ~ dnorm(nu[i], var = sigsqrd)  # slopes as samples from mixture dist.
    nu[i] <- nuTilde[ki[i]]                 # mean for each cluster xi[i]
  }
  # mixture component parameters drawn from base measures
  for(i in 1:nchm) {
    nuTilde[i] ~ dnorm(0, 1)
  }
  # CRP for clustering studies to mixture components
  ki[1:nchm] ~ dCRP(alpha, size = nchm)
  # hyperparameters
  alpha ~ dgamma(shape = 3/2, rate = 1/6)
})

# if using previous RE model nimble MCMC result
#beta1 = colMeans(samples)[grep("beta",colnames(samples))]
beta1 = scale(slopes)[1:18]
nchm = length(slopes)
inits <- list(sigsqrd = 1, 
              ki = rep(1, nchm),
              alpha = 1, 
              nuTilde = rnorm(nchm))

data <- list(beta1 = beta1)

constants <- list(nchm = length(beta1))

dpModel <- nimbleModel(code = codeDP, name = 'DPMM', constants = constants,
                       data = data, inits = inits)
#CdpModel <- compileNimble(dpModel)

dpMCMC <- buildMCMC(dpModel, 
                    monitors = c("nu", "ki", "sigsqrd"))
#test_time = runMCMC(dpMCMC, niter = 1000)

CdpMCMC <- compileNimble(dpMCMC, project = dpModel)
niter = 5e4; nburnin=1e4
dpMCMC_samples <- runMCMC(CdpMCMC, niter = niter, nburnin = nburnin)


cluster_assign = dpMCMC_samples[,grep("ki", colnames(dpMCMC_samples))]
cluster_value = dpMCMC_samples[,grep("nu", colnames(dpMCMC_samples))]
for( i in 1:(niter-nburnin)){
  # relabel for uniqueness
  c_i = cluster_assign[i,]
  c_i = c_i - min(c_i)+1
  uniq_labs = unique(c_i)
  relabeled = sapply(c_i, function(x) which(uniq_labs == x)[1])
  cluster_assign[i,] = relabeled
  #c_i = uniq_labs[order(lab_map)]
  # phi_c = phi_c[uniq_labs]
  #cluster_value[i, ] = cluster_value[i, c_i]
}


cluster_chain_nimble = list("means" =as.matrix(cluster_value),
                            "cluster"= as.matrix(cluster_assign))
cluster_centers(cluster_chain_nimble, n_top = 50)

# as one line
samples <- nimbleMCMC(code = codeDP, data = data, inits = inits,
                      constants = constants, 
                      monitors = c("nu", "ki", "sigsqrd", "alpha"),
                      thin = 5, niter = 5000, nburnin = 2000, nchains = 1, 
                      setSeed = TRUE)






## nimble test:  pump


pumpCode <- nimbleCode({ 
  # Define relationships between nodes
  for (i in 1:N){
    theta[i] ~ dgamma(alpha,beta)
    lambda[i] <- theta[i]*t[i]
    x[i] ~ dpois(lambda[i])
  }
  # Set priors
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

# Create some contrants, data, and initial values to pass to the model builder
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                         31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))


pumpModel <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                         data = pumpData, inits = pumpInits)

pumpMCMC <- buildMCMC(pumpModel)
test_time = runMCMC(pumpMCMC, niter = 1000)


CpumpModel <- compileNimble(pumpModel)
CpumpMCMC <- compileNimble(pumpMCMC, project = pumpModel)
runMCMC_samples <- runMCMC(CpumpMCMC, nburnin = 1000, niter = 10000)
plot(runMCMC_samples[ , 'alpha'], type = 'l', xlab = 'iteration',  ylab = expression(alpha))





#  ----------- simple example  -------------
code <- nimbleCode({
  ## noninformative priors
  mu ~ dflat()
  sigma ~ dhalfflat()
  ## likelihood
  for(i in 1:n) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
})

data <- list(y = MASS::newcomb)
inits <- list(mu = 0, sigma = 5)
constants <- list(n = length(data$y))

model <- nimbleModel(code = code, data = data, constants = constants, inits = inits)
## Ensure we have the nodes needed to simulate new datasets
dataNodes <- model$getNodeNames(dataOnly = TRUE)
parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- model$getDependencies(parentNodes, self = FALSE)

cmodel  <- compileNimble(model)
mcmc    <- buildMCMC(model, monitors = parentNodes)
cmcmc   <- compileNimble(mcmc, project = model)
samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500)



# -----------  one sample  -------------
code <- nimbleCode({
  a1 ~ dnorm(0, sd = 100)
  theta1 ~ dgamma(shape = 1, rate = 1)
  beta1 ~ dgamma(shape = 1, rate = 1)
  sigma_eps ~ dunif(0, 100)
  for (i in 1:N) {
    f[i] <-  a1/(1+(theta1 / x[i])^beta1)
    y[i] ~ dnorm(f[i], sd = sigma_eps)
  }
})

## constants, data, and initial values
constants <- list(N = length(y_i[1,1:13]))

data <- list(
  y = y_i[1,1:13],
  x =Cx[1,1:13]
)

inits <- list(beta1 = 1, theta1 = 1e-6, a1 = 10, sigma_eps = 1)

## create the model object
drcModel <- nimbleModel(code = code, constants = constants, data = data, 
                        inits = inits, check = FALSE)

## Ensure we have the nodes needed to simulate new datasets
dataNodes_drc <- drcModel$getNodeNames(dataOnly = TRUE)
parentNodes_drc <- drcModel$getParents(dataNodes_drc, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes_drc <- drcModel$getDependencies(parentNodes_drc, self = FALSE)



c_drcmodel  <- compileNimble(drcModel)
mcmc    <- buildMCMC(c_drcmodel, monitors = parentNodes_drc)
cmcmc   <- compileNimble(mcmc, project = drcModel)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 1000)
colMeans(samples)

beepr::beep()



# -----------  replicates  -------------
code <- nimbleCode({
  a1 ~ dnorm(0, sd = 100)
  theta1 ~ dgamma(shape = 1, rate = 1)
  beta1 ~ dgamma(shape = 1, rate = 1)
  sigma_eps ~ dunif(0, 100)
  sigma_u ~ dunif(0, 100)
  sigma_v ~ dunif(0, 100)
  for(j in 1:R){
    u[j] ~ dnorm(0, sd = sigma_u)
    v[j] ~ dnorm(0, sd = sigma_v)
    # each row is a replicate
    for (i in 1:N) {
      f[j, i] <-  (a1+ u[j])/(1+(theta1 / x[j,i])^beta1) + v[j]
      y[j, i] ~ dnorm(f[j,i], sd = sigma_eps)
    }
  }
})

## constants, data, and initial values
curr_ridx = replicate_sets[[1]]
constants <- list(N = ncol(y_i), R =length(curr_ridx))
data <- list(
  y = y_i[curr_ridx,],
  x =Cx[curr_ridx,]
)

inits <- list(beta1 = 1, theta1 = 1e-6, a1 = 10, sigma_eps = 1,
              u = c(0,0,0), v = c(0,0,0))
drcModel <- nimbleModel(code = code, constants = constants, data = data, 
                        inits = inits, check = FALSE)

## Ensure we have the nodes needed to simulate new datasets
dataNodes_drc <- drcModel$getNodeNames(dataOnly = TRUE)
parentNodes_drc <- drcModel$getParents(dataNodes_drc, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes_drc <- drcModel$getDependencies(parentNodes_drc, self = FALSE)
c_drcmodel  <- compileNimble(drcModel)
mcmc    <- buildMCMC(c_drcmodel, monitors = parentNodes_drc)
cmcmc   <- compileNimble(mcmc, project = drcModel)
samples <- runMCMC(cmcmc, niter = 20000, nburnin = 1000, thin = 20)
colMeans(samples)
plot(samples[,4])
beepr::beep()
