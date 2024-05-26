library(survival)


\donttest{

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
}
