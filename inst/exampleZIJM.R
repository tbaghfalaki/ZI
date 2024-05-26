library(survival)

\donttest{
  data(long_data_nb)
  data(surv_data_nb)

  Z1 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = survival::Surv(survtime, death) ~ w1 + w2,
    dataLong = long_data_nb, dataSurv = surv_data_nb,
    obstime = "obstime", id = "id",
    n.chains = 2,
    n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "Poisson"
  )



  Z2 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1 + w2,
    dataLong = long_data_nb, dataSurv = surv_data_nb,
    obstime = "obstime", id = "id", n.chains = 2,
    n.iter = 100, n.burnin = 50, n.thin = 1, K = 15, family = "NB"
  )




  Z3 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~1, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1 ,
    dataLong = long_data_nb, dataSurv = surv_data_nb,
    obstime = "obstime", id = "id",
    n.chains = 2,
    n.iter = 10, n.burnin = 5, n.thin = 1, K = 15, family = "GP"
  )

}
