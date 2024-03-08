data(long.data)
data(surv.data)

ZIJM(
  FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
  formSurv = survival::Surv(survtime, death) ~ w1,
  dataLong = long.data, dataSurv = surv.data,
  obstime = "obstime",
  n.chains = 1,
  n.iter = 200, n.burnin = 100, n.thin = 1, K = 15, family = "NB"
)
