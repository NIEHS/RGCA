[![R-CMD-check](https://github.com/kyle-messier/RGCA/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/kyle-messier/RGCA/actions/workflows/check-standard.yaml)
[![lint](https://github.com/kyle-messier/RGCA/actions/workflows/lint.yaml/badge.svg)](https://github.com/kyle-messier/RGCA/actions/workflows/lint.yaml)
[![codecov](https://codecov.io/gh/kyle-messier/RGCA/graph/badge.svg?token=KSNSLY1ULL)](https://codecov.io/gh/kyle-messier/RGCA)
[![cov](https://kyle-messier.github.io/RGCA/badges/coverage.svg)](https://github.com/kyle-messier/RGCA/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

# RGCA
Reflected Generalized Concentration Addition: A geometric, piecewise inverse function for 3+ parameter sigmoidal (e.g. hill) models used in chemical mixture concentration-response modeling

---
## Key Inverse Function:  

$$
\begin{equation} 
f^{-1}(r | \alpha>0, \theta, \beta>0) = 
\begin{cases}
   \frac{- \theta}{1+\left (\frac{-\alpha}{r} \right)^{1/\beta}} & r \in (-\infty, 0)\\
 \theta \left (\frac{\alpha}{r} -1 \right)^{-1/\beta} & r \in[0, \alpha)\\
  -2\theta - \theta\left (\frac{\alpha}{2\alpha - r} -1 \right)^{-1/\beta} & r \in (\alpha, 2\alpha)\\
    -2\theta + \frac{\theta}{1+\left(\frac{\alpha}{r-2\alpha} \right)^{1/\beta}} & r \in (2\alpha, \infty)\\
 \end{cases}
\end{equation}
$$

This inverse provides a wide enough support to satisfy the invertibility requirements of GCA, but with non-unit slopes. The resulting inverse maintains a coarse hyperbolic shape and continuity and is smooth at the transitions.  This procedure is not limited to the Hill function and can be applied to any monotonic dose response function, but the resulting stability may vary.  Note that negative slope parameters for the Hill function are not supported.

## Abstract:  
  Environmental toxicants overwhelmingly occur together as mixtures. The variety of possible chemical interactions makes it difficult to predict the danger of the mixture. In this work, \hl{we propose the novel Reflected Generalized Concentration Addition (RGCA), a piece-wise, geometric technique for sigmoidal dose-responsed inverse functions that extends the use of generalized concentration addition (GCA) for 3+ parameter  models.  Since experimental tests of all relevant mixtures is costly and intractable, we rely only on the individual chemical dose responses. Additionally, RGCA enhances} the classical two-step model for the cumulative effects of mixtures, which assumes a combination of GCA and  independent action (IA).  We explore how various clustering methods can dramatically improve predictions.  We compare our technique to the IA, CA, and GCA models and show in a simulation study that \hl{the two-step approach performs well under a variety of true models. We then apply our method to a challenging data set of individual chemical and mixture responses where the target is an androgen receptor (Tox21 AR-luc). Our results show significantly improved predictions for larger mixtures. Our work complements ongoing efforts to predict environmental exposure to various chemicals and offers a starting point for combining different exposure predictions to quantify a total risk to health.
