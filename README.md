# Zero-Inflated Joint modeling
Perform a Gibbs sampler for hurdle models, which encompass zero-inflated longitudinal models. The package features longitudinal zero-inflated models spanning Poisson, generalized Poisson, and negative binomial distributions. Its main emphasis lies in estimating joint models for zero-inflated longitudinal measurements alongside time-to-event data, encompassing Poisson, generalized Poisson, and negative binomial distributions. The "ZIJM" function facilitates joint modeling, incorporating a proportional hazard sub-model and a piecewise constant baseline hazard, treating the current value as the association between the two sub-models. Meanwhile, the "ZISRE" function engages in joint modeling with a Weibull sub-model by integrating a shared random effects model.

### Installation
To acquire the latest development version of ZI, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/ZI")
```
This will seamlessly fetch and install the most up-to-date version of ZI for your use.

### Reference 
Ganjali, M., Baghfalaki, T. & Balakrishnan, N. (2024). Joint Modeling of Zero-Inflated Longitudinal Measurements and Time-to-Event Outcomes with Applications to Dynamic Prediction. *Revised*.
