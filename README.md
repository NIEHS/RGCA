# RGCA-DP
Reflected Generalized Concentration Addition with Dirichlet Process

---
## Summary:  
This code predicts the toxic effect of a mixture of chemicals given the individual toxic effects.  See RGCA_manuscript_pipe.Rmd for a vignette showing how it is applied to a set of Androgen Receptor (AR-luc) dose response data taken from the Tox21 database.

## Background:  
The underlying method is based on Concentration Addition (CA) and Independent Action (IA), two classical methods for predicting the effect of a mixture of chemicals.  CA estimates the effective dose or concentration of a mixture, as if it were a single chemical.  IA essentially adds the response of each chemical, rather than figuring out an effective dose.  These two methods are frequently combined by grouping together chemicals that are assumed to be similar and using CA within the groups and IA across groups, and this strategy is used here.

CA is restricted to normalized responses for 2-parameter functions (eg Hill Model with parameters for the maximum effect and inflection point).  Generalized Concentration Addition removes the restriction of normalized responses, and our Reflected GCA allows for 3-parameter models (adds a slope parameter). 

We provide uncertainty estimation by propagating all variances.  The individual toxic response parameters have variance due to noisy data, but the clustering also introduces uncertainty.  In particular, we cluster chemicals based on the slope parameter fitted to the individual toxic responses using a Dirichlet Process cluster model.  Our sampling method is based on the Polya urn model, which is simple but not very scalable.  The slope parameter is used because it can be suggestive of the mode of action of the chemical.

## Process
0. Prepare data:  dose or concentration by chemical, response by chemical, response by mixture, and mixture definitions
1. Fit individual dose response curves
2. Cluster individual dose responses using slope parameter.
3. Sample the parameter and clustering values
3. Predict a mixture given the clustering and component concentrations.
4.  Repeat steps 3 and 4 to create uncertainty intervals
