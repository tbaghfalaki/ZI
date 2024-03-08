library(survival)
data(long.data)
data(surv.data)

\donttest{
ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = survival::Surv(survtime, death) ~ w1,
  dataLong = long.data, dataSurv = surv.data,
  obstime = "obstime",id="id",
  n.chains = 1,
  n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "Poisson"
)


ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long.data, dataSurv = surv.data,
  obstime = "obstime",id="id",
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "GP"
)



ZIJMCV(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long.data, dataSurv = surv.data,
  obstime = "obstime",id="id",
  n.chains = 1,
  n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "NB"
)
}
