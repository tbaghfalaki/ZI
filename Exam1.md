Getting Started
---------------

```
library(ZI)
```
Loading the data from the package includes both longitudinal data in long format and survival data. It's essential to ensure that the same subject (ID) is present in both datasets.

```
data(long_data_nb)
data(surv_data_nb)
```

Dividing data to 70% training data and 30% validation set:

```
set.seed(2)
  INDTRAIN <- sample(surv_data_nb$id, 0.7 * (dim(surv_data_nb)[1]))
  INDVALID <- surv_data_nb$id[-INDTRAIN]
  dataLong_t <- subset(
    long_data_nb,
    long_data_nb$id %in% INDTRAIN
  )
  dataSurv_t <- subset(
    surv_data_nb,
    surv_data_nb$id %in% INDTRAIN
  )
  names(dataSurv_t)

  dataLong_v <- subset(
    long_data_nb,
    long_data_nb$id %in% INDVALID
  )
  dataSurv_v <- subset(
    surv_data_nb,
    surv_data_nb$id %in% INDVALID
  )
```
We are considering one marker; therefore, we require one fixed effects model for rate and probability, one random effects model for , and another for the survival model.

```
FixedY = Y1 ~ obstime + x1 + x2
RandomY = ~obstime
GroupY = ~ id
FixedZ = ~ obstime + x1
RandomZ = ~obstime
GroupZ = ~id
formSurv = Surv(survtime, death) ~ w1 + w2,
```

This versatile package supports a wide range of distributional assumptions, encompassing Poisson, negative binomial, and generalized Poisson. Further elaboration on the model's specifics can be found in Ganjali et al. (2024).

These joint models are operationalized through two pivotal functions: (1) "ZIJMCV", facilitating joint modeling with a proportional hazard sub-model and a piecewise constant baseline hazard, considering associations based on the current values, and (2) "ZISRE," enabling joint modeling with a Weibull sub-model by incorporating a shared random effects model. At first, we consider the first one:

Finally, we have to use the VS function with the following arguments:

-  FixedY formula for fixed part of longitudinal count model
-  RandomY formula for random part of longitudinal count model
-  GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
-  FixedZ formula for fixed part of longitudinal probability model
-  RandomZ formula for random part of longitudinal probability model
-  GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
-  formSurv formula for survival model
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  obstime the observed time in longitudinal data
-  id the id variable in longitudinal data
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.
-  K Number of nodes and weights for calculating Gaussian quadrature
-  family Family objects provide a convenient way to specify the details of the models. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson joint models, respectively, while the others are utilized for the corresponding models based on their names.

As an example, consider the following command, where this implementation has been performed on training data:

```
Z2 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1 + w2,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", id = "id", n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "NB"
  )
```

A part of its output of the function is as follows:

```
$Estimation
$Estimation$Count_model
                   Est         SD         L_CI       U_CI     Rhat
(Intercept) -0.9393380 0.15128469 -1.224993616 -0.6484760 1.621542
obstime      0.2625674 0.16164251  0.005090672  0.6192264 1.290469
x1          -0.6079522 0.19236562 -1.060006716 -0.2798539 1.072495
x2          -0.6832724 0.08898479 -0.852994882 -0.5144040 1.133840
Dispersion   5.2285879 0.98801470  3.546213139  7.4822129 1.035696

$Estimation$Zero_inflated_model
                   Est        SD       L_CI       U_CI     Rhat
(Intercept) -1.0435121 0.1543119 -1.3578442 -0.7637013 1.025776
obstime     -0.8803696 0.1796026 -1.2730814 -0.5391646 1.056521
x1           1.0631190 0.1720780  0.7092521  1.3984815 1.017474

$Estimation$Survival_model
                    Est         SD       L_CI       U_CI     Rhat
w1            1.0193130 0.09347805  0.8465882  1.2062719 1.019889
w2           -1.1007597 0.16804981 -1.4314663 -0.7755796 1.005280
gamma_lambda -0.5437266 0.08571440 -0.7250236 -0.3866456 1.113134
gamma_pi     -0.4982187 0.09159170 -0.6773262 -0.3245060 1.002459
h1            0.9093232 0.21471229  0.5623703  1.4385323 1.171186
h2            0.8910990 0.17250013  0.5998604  1.2731592 1.126302
h3            0.8472350 0.14652661  0.5958458  1.1625943 1.099862
h4            0.9301828 0.17185870  0.6405757  1.3017021 1.037230
h5            0.8675282 0.39187958  0.3110176  1.7871397 1.001031

$Estimation$D
$Estimation$D$D11
          Intercept     Slope
Intercept 0.9908800 0.5480709
Slope     0.5480709 1.1476204

$Estimation$D$D22
          Intercept     Slope
Intercept 0.8144997 0.2615043
Slope     0.2615043 0.7440443



$DIC
[1] 5073.027

$LPML
[1] -2296.774
```


For the share random effects, the function ZISRE has the following arguments:

-  FixedY formula for fixed part of longitudinal count model
-  RandomY formula for random part of longitudinal count model
-  GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
-  FixedZ formula for fixed part of longitudinal probability model
-  RandomZ formula for random part of longitudinal probability model
-  GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
-  offset the offset or library size for discrete response. If offset=NULL, it is considered without an offset.
-  obstime the observed time in longitudinal data
-  formSurv formula for survival model
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.
-  family Family objects provide a convenient way to specify the details of the models. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson joint models, respectively, while the others are utilized for the corresponding models based on their names.

As an example, consider the following command, where this implementation has been performed on training data:
```
  data(long_data_nb)
  data(surv_data_nb)
  Z1 <- ZISRE(
    FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1, IStructure=TRUE,
    obstime = "obstime", offset = NULL,
    dataLong = long_data_nb, dataSurv = surv_data_nb,
    n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, family = "NB"
  )

```

Dynamic prediction
---------------
##### By considering current value for the association 

To compute dynamic predictions (DP), various functions are available. 

The initial method involves computing DP through a first-order approximation utilizing the *DP_CV* function, which requires the following arguments:
-  object an object inheriting from class ZIJMCV
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.

The second approach entails computing DP via an MCMC approximation, enabling the computation of credible intervals. This is achieved using the *DP_CV_CI* function, which necessitates the following arguments:

-  object an object inheriting from class ZIJMCV
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.

The final method involves generating plots for DP using the *DPplot1* function, which requires the following arguments:

-  object an object inheriting from class ZIJMCV
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  id_new id number for individual who want to plot his/her DP
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  by number: increment of the sequence of DP.
-  Marker_lab the label for the response axis
-  Time_lab the label for the time axis
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.


Example: 
```
 DD <- DP_CV(
    object = Z2, s = 0.5, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )
```

The outputs of this function are as follows:

```
$DP
     id        est
1     5 0.54272784
2    10 0.36723442
3    15 0.21209392
4    16 0.42057850
5    18 0.58191392
6    19 0.93703289
7    21 0.78147087
8    22 0.86166804
.
.
.
144 488 0.79895668
145 492 0.19337471
146 493 0.08209953
147 494 0.42388086
148 495 0.18696540
149 498 0.74811362
150 500 0.66517540

$s
[1] 0.5

$t
[1] 0.5
```
Computing AUC and BS for the predictions
---------------
For this purpose, we use DPCri package <https://github.com/tbaghfalaki/DPCri>.

Computing the criteria using this package is straightforward, as demonstrated by the following commands:

- s the landmark time for prediction
- t the window of prediction for prediction
- Survt the survival time
- CR the indicator for competing risks or censoring
- P the risk predictions
- cause the main cause for prediction


Consider the following command: 

```
library(survival)
library(DPCri)

Criteria(
  s = 0.5, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$death, P = DD$DP$est, cause = 1
)
```
with the following outputs:

```
$Cri
          est         sd
AUC 0.7843750 0.05783069
BS  0.1673699 0.02236938
```

Additionally, for the plot, we have the following visualization for subject 167:

```
DPplot1(Z2,
  s = 1.5, id_new = 167, by = 0.1, mi = 5,
  Marker_lab="Marker", Time_lab="Time",
  n.chains = 1, n.iter = 20, n.burnin = 10,
  dataLong = dataLong_v, dataSurv = dataSurv_v
 )
```


##### By considering shared random effects

To compute dynamic predictions (DP), various functions are available. 

The initial method involves computing DP through a first-order approximation utilizing the *DP_SRE* function, which requires the following arguments:
-  object an object inheriting from class VS
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.

The second approach entails computing DP via an MCMC approximation, enabling the computation of credible intervals. This is achieved using the *DP_SRE_CI* function, which necessitates the following arguments:

-  object an object inheriting from class VS
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.


The final method involves generating plots for DP using the *DPplot2* function, which requires the following arguments:

-  object an object inheriting from class ZISRE
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  id_new id number for individual who want to plot his/her DP
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  by number: increment of the sequence of DP.
-  Marker_lab the label for the response axis
-  Time_lab the label for the time axis
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.


Example: 

```
DD <- DP_SRE(Z1,
    s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )

```
with the following outputs:

```
$DP
     id          est
1     5 0.0108342165
2    10 0.8862672150
3    15 0.0315760312
4    16 0.7846448587
5    18 0.1739335618
6    19 0.0342197666
7    21 0.4936661757
8    22 0.8842215307
.
.
.
143 481 0.1513871921
144 488 0.5682667495
145 492 0.1234556872
146 493 0.3496337934
147 494 0.1883084163
148 495 0.6348319831
149 498 0.6761305704
150 500 0.6103851130

$s
[1] 0.1

$t
[1] 0.5
```



Additionally, for the plot, we have the following visualization for subject 498:

```
DPplot2(Z1,
  s = 0.4, id_new = 498, by = 0.2, mi = 5,
  Marker_lab="Biomarker", Time_lab="Time (week)",
  n.chains = 1, n.iter = 20, n.burnin = 10,
  dataLong = dataLong_v, dataSurv = dataSurv_v
)
```




