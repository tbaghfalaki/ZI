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


  Z2 <- ZISRE(
    FixedY = Y1 ~ x1 + x2 + obstime, RandomY = ~obstime, GroupY = ~id,
    FixedZ = ~ x1 + x2 + obstime, RandomZ = ~obstime, GroupZ = ~id,
    formSurv = Surv(survtime, death) ~ w1, IStructure=TRUE,
    dataLong = dataLong_t, dataSurv = dataSurv_t,
    obstime = "obstime", offset = NULL,
    n.chains = 2,
    n.iter = 200, n.burnin = 100, n.thin = 1, family = "NB"
  )

  DD <- DP_SRE(Z2,
               s = 0.1, t = 0.5, offset = NULL, n.chains = 1, n.iter = 2000, n.burnin = 500,
               n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )


  DD <- DP_SRE_CI(Z2,
                  s = 0.1, t = 0.5, offset = NULL, mi=5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
                  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )

  DPplot2(Z2,
          s = 0.4, id_new = 168, by = 0.1, mi = 2,digits=1,
          offset = NULL,
          Marker_lab="Biomarker", Time_lab="Time (week)",
          n.chains = 1, n.iter = 20, n.burnin = 10,
          dataLong = dataLong_v, dataSurv = dataSurv_v
  )

}

