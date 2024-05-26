#' Zero-inflation joint modeling
#'
#' @description
#' Fits zero-inflated hurdle joint modeling with a proportional hazard sub-model and a piecewise constant baseline hazard, considering associations based on the current values. This function accommodates various distributional assumptions, including Gaussian, Gamma, inverse Gaussian, Weibull, exponential, beta, Poisson, negative binomial, logarithmic, Bell, generalized Poisson, and binomial.
#'
#' @details
#' Function using 'JAGS' software to estimate the hurdle joint modeling
#'
#' @param FixedY formula for fixed part of longitudinal count model
#' @param RandomY formula for random part of longitudinal count model
#' @param GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
#' @param FixedZ formula for fixed part of longitudinal probability model
#' @param RandomZ formula for random part of longitudinal probability model
#' @param GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
#' @param formSurv formula for survival model
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param obstime the observed time in longitudinal data
#' @param id the id variable in longitudinal data
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param K Number of nodes and weights for calculating Gaussian quadrature
#' @param family Family objects provide a convenient way to specify the details of the models. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson joint models, respectively, while the others are utilized for the corresponding models based on their names.
#'
#'
#' @importFrom
#'
#' @return
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#' - betaL1 the regression coefficients of rate of the hurdle model
#' - betaL2 the regression coefficients of probability of the hurdle model
#' - betaS the regression coefficients of the survival model
#' - Sigmaa the variance of random effects of the rate
#' - Sigmab the variance of random effects of the probability
#' - gamma_pi the association parameter for linear predictor of the current value of the probability
#' - gamma_lambda the association parameter for linear predictor of the current value of the rate
#' - h the parameters of the piecewise constant baseline hazard
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZIJM.R
#'
#' @md
#' @export

ZIJMCV <- function(FixedY, RandomY, GroupY, FixedZ, RandomZ, GroupZ, formSurv, dataLong, dataSurv,
                   obstime = "obstime", id = "id",
                   n.chains = n.chains,
                   n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, K = 15, family = "Poisson") {
  data_long <- dataLong[unique(c(
    all.vars(GroupY), all.vars(FixedY), all.vars(RandomY),
    all.vars(GroupZ), all.vars(FixedZ), all.vars(RandomZ)
  ))]
  y <- data_long[all.vars(FixedY)][, 1]
  mfX <- stats::model.frame(FixedY, data = data_long)
  X1 <- stats::model.matrix(FixedY, mfX)
  mfU <- stats::model.frame(RandomY, data = data_long)
  Z1 <- stats::model.matrix(RandomY, mfU)
  # id_prim <- as.integer(data_long[id][, 1])

  id <- as.integer(data_long[all.vars(GroupY)][, 1])

  M <- table(id)
  id_prim <- rep(1:length(M), M)


  mfX2 <- stats::model.frame(FixedZ, data = data_long)
  X2 <- stats::model.matrix(FixedZ, mfX2)
  mfU2 <- stats::model.frame(RandomZ, data = data_long)
  Z2 <- stats::model.matrix(RandomZ, mfU2)

  n2 <- length(unique(id_prim))

  n <- length(dataLong[, 1])
  z <- rep(0, n)
  z[y == 0] <- 1
  ## Number of patients and number of longitudinal observations per patient
  n1 <- length(y)
  Nbeta1 <- dim(X1)[2]
  Nbeta2 <- dim(X2)[2]




  Obstime <- obstime
  Xvtime1 <- cbind(id_prim, X1[, colnames(X1) %in% setdiff(colnames(X1), Obstime)])
  Xv1 <- Xvtime1[!duplicated(Xvtime1), -1] ### X of data without time and id replications

  indB1 <- 1:dim(X1)[2]
  indtime1 <- indB1[colnames(X1) %in% Obstime] # index of time


  Xvtime2 <- cbind(id_prim, X2[, colnames(X2) %in% setdiff(colnames(X2), Obstime)])
  Xv2 <- Xvtime2[!duplicated(Xvtime2), -1] ### X of data without time and id replications

  indB2 <- 1:dim(X2)[2]
  indtime2 <- indB2[colnames(X2) %in% Obstime] # index of time

  nindtime1 <- c(1:dim(X1)[2])[-indtime1]
  nindtime2 <- c(1:dim(X2)[2])[-indtime2]

  #####################
  # formSurv=survival::formSurv
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  mfZ <- stats::model.frame(formSurv, data = tmp)
  XS <- stats::model.matrix(formSurv, mfZ)[, -1]
  ########  Gauss-Legendre quadrature (15 points)  ########

  glq <- statmod::gauss.quad(K, kind = "legendre")
  xk <- glq$nodes # Nodes
  wk <- glq$weights # Weights
  K <- length(xk) # K-points
  ################
  peice <- stats::quantile(Time, seq(.2, 0.8, length = 4))

  ###########

  Beta1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.beta(y[i],as[i], bs[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    as[i]<-mu[i]*phi1
    bs[i]<-(1-mu[i])*phi1
      logit(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])






    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
      cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))




    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


 phi1~dgamma(.1,.1)
  phis<-1/phi1

}"



Beta <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i] <- (1-z[i])*(logdensity.beta(y[i],as[i], bs[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    as[i]<-mu[i]*phi1
    bs[i]<-(1-mu[i])*phi1
      logit(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])





    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


 phi1~dgamma(.1,.1)
  phis<-1/phi1

}"




Gamma1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    mu1[i]<-sigma1/mu[i]
      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])




    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


   sigma1~dgamma(.1,.1)
  sigma<-1/sigma1


}"



Gamma <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    mu1[i]<-sigma1/mu[i]
      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])




    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  sigma1~dgamma(.1,.1)
  sigma<-1/sigma1

}"





Weibull1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])




    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


  kappa~dgamma(.1,.1)


}"



Weibull <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])


    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  kappa~dgamma(.1,.1)


}"






Exp1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"


Exp <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])


    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"



IGauss1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


 ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])







    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"


IGauss <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])


    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"

Gaussian1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      mu[i] <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


tau~dgamma(0.01,0.01)
sigma<-1/tau

}"


Gaussian <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])



      mu[i] <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


tau~dgamma(0.01,0.01)
sigma<-1/tau

}"





logar <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-(1-z[i])*(log(-1/log(1-pi[i]))+y[i]*log(pi[i])-log(y[i]))+z[i]*log(lambda[i])+
  (1-z[i])*log(1-lambda[i])


      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"




logar1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


       ll[i]<-(1-z[i])*(log(-1/log(1-pi[i]))+y[i]*log(pi[i])-log(y[i]))+z[i]*log(lambda[i])+
  (1-z[i])*log(1-lambda[i])


      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"

binomial1 <- "model{

m<-max(y)
  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(pi[i])+(m-y[i])*log(1-pi[i])-
                   log(1-pow(1-pi[i],m)))+z[i]*log(lambda[i])+(1-z[i])*log(1-lambda[i])


      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


binomial <- "model{

m<-max(y)
  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(pi[i])+(m-y[i])*log(1-pi[i])-
                   log(1-pow(1-pi[i],m)))+z[i]*log(lambda[i])+(1-z[i])*log(1-lambda[i])

      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"

###########
NB1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-z[i]*log(pi[i]) +
(1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+
y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))


      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Poisson1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


     ll[i]<-z[i]*log(pi[i]) +(1-z[i])*(log(1-pi[i])+logdensity.pois(y[i], lambda[i])-log(1-exp(-lambda[i])))
 #ll[i]<-(1-z[i])*(y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))+z[i]*log(muz[i])+
#(1-z[i])*log(1-muz[i])

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"

GP1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

      phiz~dnorm(0,.001)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


NB <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-z[i]*log(pi[i]) +
(1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+
y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"


Poisson <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

     ll[i]<-z[i]*log(pi[i]) +(1-z[i])*(log(1-pi[i])+logdensity.pois(y[i], lambda[i])-log(1-exp(-lambda[i])))


# ll[i]<-(1-z[i])*(y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))+z[i]*log(muz[i])+
#(1-z[i])*log(1-muz[i])


      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])

    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"

GP <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

      phiz~dnorm(0,.001)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"

Bell1 <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])

    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Bell <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"





Bell1wc <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Bellwc <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k,1])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2]+b[k,2])
    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab[1:Nb2,1:Nb2]<-inverse(Omegab[,])
  Omegab[1:Nb2,1:Nb2]~dwish(V2[,],Nb2)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"




Beta1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.beta(y[i],as[i], bs[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    as[i]<-mu[i]*phi1
    bs[i]<-(1-mu[i])*phi1
      logit(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]






    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
      cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))




        Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])



    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }



  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


 phi1~dgamma(.1,.1)
  phis<-1/phi1

}"



Betat <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i] <- (1-z[i])*(logdensity.beta(y[i],as[i], bs[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    as[i]<-mu[i]*phi1
    bs[i]<-(1-mu[i])*phi1
      logit(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]





    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


 phi1~dgamma(.1,.1)
  phis<-1/phi1

}"




Gamma1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    mu1[i]<-sigma1/mu[i]
      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]




    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

        Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])




    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


   sigma1~dgamma(.1,.1)
  sigma<-1/sigma1


}"



Gammat <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

    mu1[i]<-sigma1/mu[i]
      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]




    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  sigma1~dgamma(.1,.1)
  sigma<-1/sigma1

}"





Weibull1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]




    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

       Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


  kappa~dgamma(.1,.1)


}"



Weibullt <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]


    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  kappa~dgamma(.1,.1)


}"






Exp1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"


Expt <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- -(inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,]))
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]


    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"



IGauss1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


 ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]







    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

        Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])

    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"


IGausst <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]


    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


  lambda~dgamma(.1,.1)
  sigma<-1/lambda

}"

Gaussian1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      mu[i] <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

       Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])

   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])


    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


tau~dgamma(0.01,0.01)
sigma<-1/tau

}"


Gaussiant <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])



      mu[i] <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


  haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }


tau~dgamma(0.01,0.01)
sigma<-1/tau

}"





logart <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-(1-z[i])*(log(-1/log(1-pi[i]))+y[i]*log(pi[i])-log(y[i]))+z[i]*log(lambda[i])+
  (1-z[i])*log(1-lambda[i])


      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


  haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"




logar1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


       ll[i]<-(1-z[i])*(log(-1/log(1-pi[i]))+y[i]*log(pi[i])-log(y[i]))+z[i]*log(lambda[i])+
  (1-z[i])*log(1-lambda[i])


      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])



    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"

binomial1t <- "model{

m<-max(y)
  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(pi[i])+(m-y[i])*log(1-pi[i])-
                   log(1-pow(1-pi[i],m)))+z[i]*log(lambda[i])+(1-z[i])*log(1-lambda[i])


      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])





    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


binomialt <- "model{

m<-max(y)
  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(pi[i])+(m-y[i])*log(1-pi[i])-
                   log(1-pow(1-pi[i],m)))+z[i]*log(lambda[i])+(1-z[i])*log(1-lambda[i])

      logit(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}



  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"

###########
NB1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-z[i]*log(pi[i]) +
(1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+
y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))


      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])



    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Poisson1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


     ll[i]<-z[i]*log(pi[i]) +(1-z[i])*(log(1-pi[i])+logdensity.pois(y[i], lambda[i])-log(1-exp(-lambda[i])))
 #ll[i]<-(1-z[i])*(y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))+z[i]*log(muz[i])+
#(1-z[i])*log(1-muz[i])

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])




    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"

GP1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

      Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])




    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

      phiz~dnorm(0,.001)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


NBt <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-z[i]*log(pi[i]) +
(1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+
y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"


Poissont <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

     ll[i]<-z[i]*log(pi[i]) +(1-z[i])*(log(1-pi[i])+logdensity.pois(y[i], lambda[i])-log(1-exp(-lambda[i])))


# ll[i]<-(1-z[i])*(y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))+z[i]*log(muz[i])+
#(1-z[i])*log(1-muz[i])


      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])




    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"

GPt <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

  ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

      phiz~dnorm(0,.001)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"

Bell1t <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Bellt <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


   haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"





Bell1wct <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))


    Alpha0[k]<- betaS*XS[k]+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])
    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

    betaS~dnorm(0,0.001)


}"


Bellwct <- "model{


  for(i in 1:n){
      cpoinvy[i]<- 1/(0.0001+exp(ll[i]))

    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

    ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])


      log(theta[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+b[id[i]]



    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i]~dnorm(0,Omegab)
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
        cpoinvt[k]<- 1/(0.0001+exp(death[k]*log(haz[k])+logSurv[k]))

     Alpha0[k]<- inprod(betaS[],XS[k,])+gamma_lambda*(inprod(betaL1[nindtime1],Xv1[k,])+a[k,1])+
    gamma_pi*(inprod(betaL2[nindtime2],Xv2[k,])+b[k])

    Alpha1[k]<- gamma_lambda*(betaL1[indtime1]+a[k,2])+gamma_pi*(betaL2[indtime2])


    haz[k]<- ((h[1]*step(s[1]-Time[k]))+
  (h[2]*step(Time[k]-s[1])*step(s[2]-Time[k]))+
  (h[3]*step(Time[k]-s[2])*step(s[3]-Time[k]))+
  (h[4]*step(Time[k]-s[3])*step(s[4]-Time[k]))+
  (h[5]*step(Time[k]-s[4])))*exp(Alpha0[k]+Alpha1[k]*Time[k])
    for(j in 1:K){
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[k,j]<-(xk[j]+1)/2*Time[k]
      wk11[k,j]<- wk[j]*Time[k]/2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k,j]<-  ((h[1]*step(s[1]-xk11[k,j]))+
        (h[2]*step(xk11[k,j]-s[1])*step(s[2]-xk11[k,j]))+
        (h[3]*step(xk11[k,j]-s[2])*step(s[3]-xk11[k,j]))+
        (h[4]*step(xk11[k,j]-s[3])*step(s[4]-xk11[k,j]))+
        (h[5]*step(xk11[k,j]-s[4])))*exp(Alpha0[k]+Alpha1[k]*xk11[k,j])

    }


    logSurv[k]<- -inprod(wk11[k,],chaz[k,])

    #Definition of the survival log-likelihood using zeros trick
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }


  Sigmaa[1:Nb1,1:Nb1]<-inverse(Omegaa[,])
  Omegaa[1:Nb1,1:Nb1]~dwish(V1[,],Nb1)

  Sigmab <-1/Omegab
  Omegab ~dgamma(0.1,0.1)

for(l in 1:J){
    h[l]~dgamma(0.1,0.1)
  }

gamma_lambda~dnorm(0,0.001)
gamma_pi~dnorm(0,0.001)

for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

}"


if (is.matrix(XS) == FALSE) {
  if (family == "Beta") {
    y[y == 1] <- 0.9999

    model.file <- textConnection(Beta1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Beta1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phis", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )



    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }




    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      phi = sim1$sims.list$phis
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }

      names(sim1$mean$phis) <-
        names(sim1$sd$phis) <-
        names(sim1$q2.5$phis) <-
        names(sim1$q97.5$phis) <-
        names(sim1$Rhat$phis) <- "Dispersion"

      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")


      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$phis, sim1$sd$phis, sim1$q2.5$phis, sim1$q97.5$phis, sim1$Rhat$phis)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$phis) <-
        names(sim1$sd$phis) <-
        names(sim1$q2.5$phis) <-
        names(sim1$q97.5$phis) <- "Dispersion"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$phis, sim1$sd$phis, sim1$q2.5$phis, sim1$q97.5$phis)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }




  if (family == "Gamma") {
    model.file <- textConnection(Gamma1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Gamma1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "sigma", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }


    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }
      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Shape"

      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Shape"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  if (family == "Weibull") {
    model.file <- textConnection(Weibull1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Weibull1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1), kappa = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "kappa", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )



    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }



    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      kappa = sim1$sims.list$kappa
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <-
        names(sim1$Rhat$kappa) <- "Shape"

      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa, sim1$Rhat$kappa)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }




      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)

      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <- "Shape"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  if (family == "Exponential") {
    model.file <- textConnection(Exp1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Exp1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1), r = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }




      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }




      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  if (family == "inverse.gaussian") {
    model.file <- textConnection(IGauss1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(IGauss1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1), r = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "sigma", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Shape"

      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }




      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Shape"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  if (family == "Poisson") {
    model.file <- textConnection(Poisson1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Poisson1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }


    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")

    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }




      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }
  #################
  if (family == "Logarithmic") {
    model.file <- textConnection(logar1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(logar1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }


    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")

    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")

      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  if (family == "Binomial") {
    model.file <- textConnection(binomial1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(binomial1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }


    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")

    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }
  #################
  if (family == "Bell") {
    model.file <- textConnection(Bell1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Bell1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }


    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")

    C <- c()
    for (i in 1:n) {
      C[i] <- log(numbers::bell(as.numeric(y[i]))) - lfactorial(y[i])
    }

    if (is.infinite(numbers::bell(max(y))) == TRUE) {
      model.file <- textConnection(Bell1wc)
      if(dim(Z2)[2]==1)(model.file <- textConnection(Bell1wct))

    }

    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, C = C
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, mub1 = rep(0, Nb1),  V1 = diag(1, Nb1), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }
    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  #########################
  if (family == "Gaussian") {
    model.file <- textConnection(Gaussian1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Gaussian1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1), r = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "sigma", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time,
      death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1),
      V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time,
        death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1,indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Variance"

      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Variance"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  if (family == "NB") {
    model.file <- textConnection(NB1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(NB1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1), r = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "r", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1,indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )


    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      r = sim1$sims.list$r
    )


    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }

      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"

      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r, sim1$Rhat$r)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <- "Dispersion"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  if (family == "GP") {
    model.file <- textConnection(GP1)
    if(dim(Z2)[2]==1)(model.file <- textConnection(GP1t))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(1), phiz = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phiz", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }


    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      phi = sim1$sims.list$phiz,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"

      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz, sim1$Rhat$phiz)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- "x"


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }
  ### ????????
} else {
  if (family == "Beta") {
    y[y == 1] <- 0.9999

    model.file <- textConnection(Beta)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Betat))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phis", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }



    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1,
      beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      phi = sim1$sims.list$phis
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$phis) <-
        names(sim1$sd$phis) <-
        names(sim1$q2.5$phis) <-
        names(sim1$q97.5$phis) <-
        names(sim1$Rhat$phis) <- "Dispersion"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$phis, sim1$sd$phis, sim1$q2.5$phis, sim1$q97.5$phis, sim1$Rhat$phis)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$phis) <-
        names(sim1$sd$phis) <-
        names(sim1$q2.5$phis) <-
        names(sim1$q97.5$phis) <- "Dispersion"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$phis, sim1$sd$phis, sim1$q2.5$phis, sim1$q97.5$phis)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }




  if (family == "Gamma") {
    model.file <- textConnection(Gamma)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Gammat))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "sigma", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      sigma = sim1$sims.list$sigma,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Shape"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Shape"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  if (family == "Weibull") {
    model.file <- textConnection(Weibull)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Weibullt))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), kappa = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "kappa", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      kappa = sim1$sims.list$kappa
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <-
        names(sim1$Rhat$kappa) <- "Shape"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa, sim1$Rhat$kappa)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$kappa) <-
        names(sim1$sd$kappa) <-
        names(sim1$q2.5$kappa) <-
        names(sim1$q97.5$kappa) <- "Shape"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$kappa, sim1$sd$kappa, sim1$q2.5$kappa, sim1$q97.5$kappa)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }



  if (family == "Exponential") {
    model.file <- textConnection(Exp)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Expt))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }


    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)




      MM <-
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)


      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)



      MM <-
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  if (family == "inverse.gaussian") {
    model.file <- textConnection(IGauss)
    if(dim(Z2)[2]==1)(model.file <- textConnection(IGausst))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), sigma = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "sigma", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h,
      sigma = sim1$sims.list$sigma
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Shape"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Shape"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  if (family == "Poisson") {
    model.file <- textConnection(Poisson)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Poissont))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
    ### @@@@@@@@@@@@
  }

  if (family == "Logarithmic") {
    model.file <- textConnection(logar)
    if(dim(Z2)[2]==1)(model.file <- textConnection(logart))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
    ### @@@@@@@@@@@@
  }

  if (family == "Binomial") {
    model.file <- textConnection(binomial)
    if(dim(Z2)[2]==1)(model.file <- textConnection(binomialt))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
    ### @@@@@@@@@@@@
  }
  if (family == "Bell") {
    model.file <- textConnection(Bell)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Bellt))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")

    C <- c()
    for (i in 1:n) {
      C[i] <- log(numbers::bell(y[i])) - lfactorial(y[i])
    }
    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K, C = C
    )
    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K,C=C
      )
    }

    if (is.infinite(numbers::bell(max(y))) == TRUE) {
      model.file <- textConnection(Bellwc)
      if(dim(Z2)[2]==1)(model.file <- textConnection(Bellwct))


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
      if(dim(Z2)[2]==1){
        d.jags <- list(
          n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
          indtime1 = indtime1, indtime2 = indtime2,
          X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
          Xv1 = Xv1, Xv2 = Xv2,
          Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
          XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
          s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
        )
      }
    }



    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }


  #########
  if (family == "Gaussian") {
    model.file <- textConnection(Gaussian)
    if(dim(Z2)[2]==1)(model.file <- textConnection(Gaussiant))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), r = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "sigma", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )


    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      sigma = sim1$sims.list$sigma,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <-
        names(sim1$Rhat$sigma) <- "Variance"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma, sim1$Rhat$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$sigma) <-
        names(sim1$sd$sigma) <-
        names(sim1$q2.5$sigma) <-
        names(sim1$q97.5$sigma) <- "Variance"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$sigma, sim1$sd$sigma, sim1$q2.5$sigma, sim1$q97.5$sigma)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  if (family == "NB") {
    model.file <- textConnection(NB)
    if(dim(Z2)[2]==1)(model.file <- textConnection(NBt))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), r = 1,
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "r", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }


    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      r = sim1$sims.list$r,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r, sim1$Rhat$r)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <- "Dispersion"


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r)
      )

      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")




      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }

  if (family == "GP") {
    model.file <- textConnection(GP)
    if(dim(Z2)[2]==1)(model.file <- textConnection(GPt))

    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]
    NbetaS <- dim(XS)[2]

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2),
        betaS = stats::rnorm(NbetaS), phiz = stats::rnorm(1),
        Omegab = diag(Nb2), Omegaa = diag(Nb1), gamma_pi = stats::rnorm(1),
        gamma_lambda = stats::rnorm(1)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "cpoinvy", "cpoinvt", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phiz", "h")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
    )

    if(dim(Z2)[2]==1){
      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death, KF1 = 100000, KF2 = 100000,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1,  z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2, NbetaS = NbetaS,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1,  mub1 = rep(0, Nb1), V1 = diag(1, Nb1),  id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )
    }

    sim1 <- jagsUI::jags(
      data = d.jags,
      parameters.to.save = parameters,
      model.file = model.file,
      n.chains = n.chains,
      parallel = FALSE,
      n.adapt = FALSE,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigmaa = sim1$sims.list$Sigmaa,
      Sigmab = sim1$sims.list$Sigmab,
      gamma_lambda = sim1$sims.list$gamma_lambda,
      gamma_pi = sim1$sims.list$gamma_pi,
      phi = sim1$sims.list$phiz,
      h = sim1$sims.list$h
    )

    if (n.chains > 1) {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <-
        names(sim1$Rhat$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <-
        names(sim1$Rhat$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <-
        names(sim1$Rhat$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <-
        names(sim1$Rhat$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <-
        names(sim1$Rhat$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <-
        names(sim1$Rhat$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        rownames(sim1$Rhat$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$Rhat$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          rownames(sim1$Rhat$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab) <-
          colnames(sim1$Rhat$Sigmab) <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab) <-
          names(sim1$Rhat$Sigmab) <- "Intercept"
      }



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa,Rhat=sim1$Rhat$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab,Rhat=sim1$Rhat$Sigmab)



      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz, sim1$Rhat$phiz)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"

      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda, sim1$Rhat$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi, sim1$Rhat$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h, sim1$Rhat$h)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    } else {
      names(sim1$mean$betaL1) <-
        names(sim1$sd$betaL1) <-
        names(sim1$q2.5$betaL1) <-
        names(sim1$q97.5$betaL1) <- colnames(X1)

      names(sim1$mean$betaL2) <-
        names(sim1$sd$betaL2) <-
        names(sim1$q2.5$betaL2) <-
        names(sim1$q97.5$betaL2) <- colnames(X2)


      names(sim1$mean$betaS) <-
        names(sim1$sd$betaS) <-
        names(sim1$q2.5$betaS) <-
        names(sim1$q97.5$betaS) <- colnames(XS)


      names(sim1$mean$h) <-
        names(sim1$sd$h) <-
        names(sim1$q2.5$h) <-
        names(sim1$q97.5$h) <- c(paste0("h", 1), paste0("h", 2), paste0("h", 3), paste0("h", 4), paste0("h", 5))


      names(sim1$mean$gamma_lambda) <-
        names(sim1$sd$gamma_lambda) <-
        names(sim1$q2.5$gamma_lambda) <-
        names(sim1$q97.5$gamma_lambda) <- "gamma_lambda"


      names(sim1$mean$gamma_pi) <-
        names(sim1$sd$gamma_pi) <-
        names(sim1$q2.5$gamma_pi) <-
        names(sim1$q97.5$gamma_pi) <- "gamma_pi"


      rownames(sim1$mean$Sigmaa) <-
        rownames(sim1$sd$Sigmaa) <-
        rownames(sim1$q2.5$Sigmaa) <-
        rownames(sim1$q97.5$Sigmaa) <-
        colnames(sim1$mean$Sigmaa) <-
        colnames(sim1$sd$Sigmaa) <-
        colnames(sim1$q2.5$Sigmaa) <-
        colnames(sim1$q97.5$Sigmaa) <- c("Intercept", "Slope")


      if(dim(Z2)[2]==2){
        rownames(sim1$mean$Sigmab) <-
          rownames(sim1$sd$Sigmab) <-
          rownames(sim1$q2.5$Sigmab) <-
          rownames(sim1$q97.5$Sigmab) <-
          colnames(sim1$mean$Sigmab) <-
          colnames(sim1$sd$Sigmab) <-
          colnames(sim1$q2.5$Sigmab) <-
          colnames(sim1$q97.5$Sigmab)  <- c("Intercept", "Slope")
      }else{
        names(sim1$mean$Sigmab) <-
          names(sim1$sd$Sigmab) <-
          names(sim1$q2.5$Sigmab) <-
          names(sim1$q97.5$Sigmab)  <- "Intercept"
      }


      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")

      D11 <- list(est=sim1$mean$Sigmaa,sd=sim1$sd$Sigmaa,L=sim1$q2.5$Sigmaa,U=sim1$q97.5$Sigmaa)
      D22 <- list(est=sim1$mean$Sigmab,sd=sim1$sd$Sigmab,L=sim1$q2.5$Sigmab,U=sim1$q97.5$Sigmab)


      MM <- rbind(
        cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz)
      )
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"

      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma_lambda, sim1$sd$gamma_lambda, sim1$q2.5$gamma_lambda, sim1$q97.5$gamma_lambda),
        cbind(sim1$mean$gamma_pi, sim1$sd$gamma_pi, sim1$q2.5$gamma_pi, sim1$q97.5$gamma_pi),
        cbind(sim1$mean$h, sim1$sd$h, sim1$q2.5$h, sim1$q97.5$h)
      )


      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")



      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = list(D11 = D11, D22 = D22))
    }
  }
}


KF1 <- 100000
KF2 <- 100000
DIC <- sim1$DIC - 2 * n1 * KF1 - 2 * n2 * KF2
LPML <- -sum(log(sim1$mean$cpoinvy)) - sum(log(sim1$mean$cpoinvt))




if (family == "Bell") {
  if (is.infinite(numbers::bell(max(y))) == TRUE) {
    DIC <- "It is not possible to compute the DIC because the values of the Bell number are not finite."
    LPML <- "It is not possible to compute the LPML because the values of the Bell number are not finite."
  }
}


list(
  FixedY = FixedY, FixedZ = FixedZ, formSurv = formSurv, GroupY = GroupY, GroupZ = GroupZ,
  RandomY = RandomY, RandomZ = RandomZ, obstime = obstime, id = id, peice = peice,
  family = family, MCMC = MCMC, Estimation = results,
  DIC = DIC, LPML = LPML
)
}
