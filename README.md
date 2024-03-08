# Zero-Inflated Joint modeling
Perform a Gibbs sampler for hurdle models, which encompass zero-inflated longitudinal models. The package features longitudinal zero-inflated models spanning Poisson, generalized Poisson, and negative binomial distributions. Its main emphasis lies in estimating joint models for zero-inflated longitudinal measurements alongside time-to-event data, encompassing Poisson, generalized Poisson, and negative binomial distributions. The "ZIJM" function facilitates joint modeling, incorporating a proportional hazard sub-model and a piecewise constant baseline hazard, treating the current value as the association between the two sub-models. Meanwhile, the "ZISRE" function engages in joint modeling with a Weibull sub-model by integrating a shared random effects model.

### Installation
To acquire the latest development version of UHM, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/UHM")
```
This will seamlessly fetch and install the most up-to-date version of UHM for your use.

### Reference 
Ganjali, M., Baghfalaki, T. & Balakrishnan, N. (2024). A Unified Bayesian approach for Modeling Zero-Inflated count and continuous outcomes. *Submitted*.
