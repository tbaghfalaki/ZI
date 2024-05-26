library(survival)
data(long.data)
data(surv.data)

\donttest{
ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1, IStructure=FALSE,
  dataLong = long.data, dataSurv = surv.data,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "Poisson"
)


ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long.data, dataSurv = surv.data, IStructure=FALSE,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "NB"
)

ZISRE(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = Surv(survtime, death) ~ w1,
  dataLong = long.data, dataSurv = surv.data, IStructure=FALSE,
  n.chains = 2,
  n.iter = 200, n.burnin = 100, n.thin = 1, family = "GP"
)
}
