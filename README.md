# Zero-Inflated Joint modeling
Perform a Gibbs sampler for hurdle models, which encompass zero-inflated longitudinal models. The package features longitudinal zero-inflated models spanning Poisson, generalized Poisson, and negative binomial distributions. Its main emphasis lies in estimating joint models for zero-inflated longitudinal measurements alongside time-to-event data, encompassing Poisson, generalized Poisson, and negative binomial distributions. The "ZIJM" function facilitates joint modeling, incorporating a proportional hazard sub-model and a piecewise constant baseline hazard, treating the current value as the association between the two sub-models. Meanwhile, the "ZISRE" function engages in joint modeling with a Weibull sub-model by integrating a shared random effects model.

### Installation
To acquire the latest development version of ZI, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/ZI")
```
This will seamlessly fetch and install the most up-to-date version of ZI for your use.

### Example Usage

This analysis is presented [here](/Exam1.md)

### Implementation with real data

The code for implementing the applications (Pregnancy Data and AIDS Data) involves fitting a Zero-Inflated Joint Model and computing the risk prediction. This analyses of data are presented [here](https://github.com/tbaghfalaki/ZIJMApp)



### Reference 
Ganjali, M., Baghfalaki, T. & Balakrishnan, N. (2024). Joint Modeling of Zero-Inflated Longitudinal Measurements and Time-to-Event Outcomes with Applications to Dynamic Prediction. SMMR, DOI: 10.1177/09622802241268466. Available at https://journals.sagepub.com/doi/full/10.1177/09622802241268466
