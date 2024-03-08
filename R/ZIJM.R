#' Zero-inflation joint modeling
#'
#' @description
#' Fits zero-inflated hurdle joint modeling under Poisson, negative binomial, and generalized Poisson distributions.
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
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#' @param K Number of nodes and weights for calculating Gaussian quadrature
#' @param family Family objects provide a convenient way to specify the details of the models. The families include "Poisson", negative binomial (by "NB") and generalized Poisson (by "GP").
#'
#'
#' @importFrom
#'
#' @return
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 25% list of posterior median for each parameter
#' - 50% list of posterior median for each parameter
#' - 75% list of posterior median for each parameter
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
#' - r the over-dispersion parameter of the negative binomial model
#' - phiz the dispersion parameter of the generalized Poisson model
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZIMEM.R
#'
#' @md
#' @export

ZIJM <- function(FixedY, RandomY, GroupY, FixedZ, RandomZ, GroupZ, formSurv, dataLong, dataSurv,
                 obstime = "obstime",
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
  id_prim <- as.integer(data_long[all.vars(GroupY)][, 1])

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
  # delta <- nnet::class.ind(arules::discretize(Time, method = "fixed", c(0, peice, max(Time))))




  ###########
  NB1 <- "model{
  KF<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF

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
    phi2[k]<-1000000-death[k]*log(haz[k])-logSurv[k]
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
  KF<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF

     ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
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
    phi2[k]<-1000000-death[k]*log(haz[k])-logSurv[k]
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
  KF<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF

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
    phi2[k]<-1000000-death[k]*log(haz[k])-logSurv[k]
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
  KF<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF

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
    phi2[k]<-1000000-death[k]*log(haz[k])-logSurv[k]
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
  KF<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF

     ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])

    a[i,1:Nb1]~dmnorm(mub1[],Omegaa[,])
    b[i,1:Nb2]~dmnorm(mub2[],Omegab[,])
  }

  # Survival and censoring times
  # Hazard function
  for(k in 1:n2){
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
    phi2[k]<-1000000-death[k]*log(haz[k])-logSurv[k]
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
  KF<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF

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
    phi2[k]<-1000000-death[k]*log(haz[k])-logSurv[k]
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

  if (is.matrix(XS) == FALSE) {
    if (family == "Poisson") {
      model.file <- textConnection(Poisson1)

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

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )

      sim1 <- R2jags::jags(
        data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
      )
    }

    if (family == "NB") {
      model.file <- textConnection(NB1)

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

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "r", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )

      sim1 <- R2jags::jags(
        data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
      )
    }

    if (family == "GP") {
      model.file <- textConnection(GP1)
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

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phiz", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )

      sim1 <- R2jags::jags(
        data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
      )
    }
  } else {
    if (family == "Poisson") {
      model.file <- textConnection(Poisson)

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

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )

      sim1 <- R2jags::jags(
        data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
      )
    }

    if (family == "NB") {
      model.file <- textConnection(NB)

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

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "r", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )

      sim1 <- R2jags::jags(
        data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
      )
    }

    if (family == "GP") {
      model.file <- textConnection(GP)
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

      parameters <- c("betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phiz", "h")


      d.jags <- list(
        n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = Time, death = death,
        indtime1 = indtime1, indtime2 = indtime2,
        X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
        Xv1 = Xv1, Xv2 = Xv2,
        Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), V1 = diag(1, Nb1), V2 = diag(1, Nb2), id = id_prim,
        XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
        s = peice, J = length(peice) + 1, xk = xk, wk = wk, K = K
      )

      sim1 <- R2jags::jags(
        data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
      )
    }
  }






  ###############
  sim1
}
