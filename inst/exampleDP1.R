library(survival)
rm(list=ls())

\donttest{

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


  Z2 <- ZIJMCV(
    FixedY = Y1 ~ obstime + x1 + x2, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ obstime + x1, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1 + w2,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", id = "id", n.chains = 2,
    n.iter = 20, n.burnin = 10, n.thin = 1, K = 15, family = "NB"
  )


  DD <- DP_CV(
    object = Z2, s = 0.5, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )


  DPplot1(Z2,
          s = 1.1, id_new = 167, by = 0.1, mi = 5,
          Marker_lab="Marker", Time_lab="Time",digits=0,
          n.chains = 1, n.iter = 20, n.burnin = 10,
          dataLong = dataLong_v, dataSurv = dataSurv_v
  )


}

