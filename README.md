# RGCA
Reflected Generalized Concentration Addition

---
## Summary:  
This code predicts the toxic effect of a mixture of chemicals given the individual toxic effects.  See RGCA_manuscript_pipe.Rmd for a vignette showing how it is applied to a set of Androgen Receptor (AR-luc) dose response data taken from the Tox21 database.

## Background:  
The underlying method is based on Concentration Addition (CA) and Independent Action (IA), two classical methods for predicting the effect of a mixture of chemicals.  CA estimates the effective dose or concentration of a mixture, as if it were a single chemical.  IA adds the response of each chemical, rather than figuring out an effective dose.  These two methods are frequently combined by grouping together chemicals that are assumed to be similar and using CA within the groups and IA across groups.  This two-step or joint strategy is used here.

CA is restricted to normalized responses for 2-parameter functions (eg Hill Model with parameters for the maximum effect and inflection point).  Generalized Concentration Addition removes the restriction of normalized responses, and our Reflected GCA allows for 3-parameter models (adds a slope parameter). 

We provide uncertainty estimation by sampling from the posterior.  The individual toxic response parameters have variance due to noisy data, but the clustering can also introduce uncertainty.  In particular, using a Bayesian procedure such as a Dirichlet Process cluster model can provide a posterior cluster distribution.  Simpler methods such as K-means or random clusters are also available but offer less principled uncertainty estimation.  Our DP sampling method is based on the Polya urn model, which is simple but not very scalable.  The DP clustering model as implemented can be applied to any one parameter.  We have two implementations of the individual dose response estimation with random effects: one written manually and one written with Nimble. 

## Process
0. Prepare data:  dose or concentration by chemical, response by chemical, response by mixture, and mixture definitions
1. Fit individual dose response curves
2. Cluster individual dose responses using K-means, a Dirichlet Process mixture model, or randomly.
3. Sample the parameter and clustering values
3. Predict a mixture given the clustering and component concentrations.
4.  Repeat steps 3 and 4 to create multiple predictions; the uncertainty intervals correspond to 5th and 95 percentiles and the prediction is the pointwise median.
