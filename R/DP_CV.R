#'  Dynamic prediction
#'
#' @description
#' Dynamic prediction for ZIJMCV
#'
#'
#' @details
#' Estimate DP for joint modeling based on ZIJMCV
#'
#' @param object an object inheriting from class ZIJMCV
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param s the landmark time for prediction
#' @param t the window of prediction for prediction
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
#'
#'
#' @importFrom stats quantile rnorm model.frame model.matrix
#'
#' @return
#' - mu.vect list of posterior mean for each parameter
#' - sd.vect list of standard error for each parameter
#' - 2.5% list of posterior mode for each parameter
#' - 97.5% list of posterior median for each parameter
#' - Rhat Gelman and Rubin diagnostic for all parameter
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#' @example inst/exampleDP1.R
#'
#' @md
#' @export

DP_CV <- function(object, s = s, t = t, n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                  n.thin = max(1, floor((n.iter - n.burnin) / 1000)), dataLong, dataSurv) {
  Dt <- t
  KK <- 1000000

  FixedY <- object$FixedY
  FixedZ <- object$FixedZ
  RandomY <- object$RandomY
  RandomZ <- object$RandomZ

  GroupY <- object$GroupY
  GroupZ <- object$GroupZ
  id <- object$id
  formSurv <- object$formSurv
  obstime <- object$obstime
  nmark <- object$nmark
  mu1 <- object$mu1
  peice <- object$peice
  family <- object$family
  #######

  ########### univariate_jm_random_effect_estimation



  Beta1b <- "model{


  for(i in 1:n){
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

  #Omegaa[1:Nb1,1:Nb1] <- inverse(Sigmaa[1:Nb1,1:Nb1])
 # Omegab[1:Nb2,1:Nb2] <- inverse(Sigmab[1:Nb2,1:Nb2])

phi1<-1/phis

}"



  Betab <- "model{


  for(i in 1:n){
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


#Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
 # Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])

phi1<-1/phis


}"




  Gamma1b <- "model{


  for(i in 1:n){
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




#Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
 # Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])




}"



  Gammab <- "model{


  for(i in 1:n){
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





 # Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
 # Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])



}"





  Weibull1b <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1


 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])




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
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}




#Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
#  Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])




}"



  Weibullb <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

 ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])


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
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


#Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
 # Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])



}"






  Exp1b <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])



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
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}


#Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
 # Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])



}"


  Expb <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+KF1

ll[i] <- (1-z[i])*(logdensity.exp(y[i],mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])



      log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(a[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],1:Nb2],Z2[i,])


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
    phi2[k]<-KF2-death[k]*log(haz[k])-logSurv[k]
    zeros2[k]~dpois(phi2[k])
}






#Omegaa[,] <-inverse(Sigmaa[1:Nb1,1:Nb1])
# Omegab[,] <-inverse(Sigmab[1:Nb2,1:Nb2])



}"



  IGauss1b <- "model{


  for(i in 1:n){
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




}"


  IGaussb <- "model{


  for(i in 1:n){
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



}"

  Gaussian1b <- "model{


  for(i in 1:n){
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





}"


  Gaussianb <- "model{


  for(i in 1:n){
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







}"





  logarb <- "model{


  for(i in 1:n){
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





}"



  logar1b <- "model{


  for(i in 1:n){
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




}"



  binomial1b <- "model{

m<-max(y)
  for(i in 1:n){
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






}"


  binomialb <- "model{

m<-max(y)
  for(i in 1:n){
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







}"

  ###########
  NB1b <- "model{


  for(i in 1:n){
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






}"


  Poisson1b <- "model{


  for(i in 1:n){
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





}"

  GP1b <- "model{


  for(i in 1:n){
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





}"


  NBb <- "model{


  for(i in 1:n){
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




}"


  Poissonb <- "model{


  for(i in 1:n){
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





}"

  GPb <- "model{


  for(i in 1:n){
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




}"

  Bell1b <- "model{


  for(i in 1:n){
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




}"


  Bellb <- "model{


  for(i in 1:n){
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




}"





  Bell1wcb <- "model{


  for(i in 1:n){
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




}"


  Bellwcb <- "model{


  for(i in 1:n){
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

}"





Beta1tb <- "model{


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


 phi1 <-1/phis

}"



Betatb <- "model{


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


  phi1 <-1/phis


}"




Gamma1tb <- "model{


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




}"



Gammatb <- "model{


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





}"





Weibull1tb <- "model{


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



}"



Weibulltb <- "model{


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






}"






Exp1tb <- "model{


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



}"


Exptb <- "model{


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




}"



IGauss1tb <- "model{


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



}"


IGausstb <- "model{


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



}"

Gaussian1tb <- "model{


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



}"


Gaussiantb <- "model{


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





}"





logartb <- "model{


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


}"




logar1tb <- "model{


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





}"

binomial1tb <- "model{

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


}"


binomialtb <- "model{

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



}"

###########
NB1tb <- "model{


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




}"


Poisson1tb <- "model{


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




}"

GP1tb <- "model{


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


}"


NBtb <- "model{


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



}"


Poissontb <- "model{


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



}"

GPtb <- "model{


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



}"

Bell1tb <- "model{


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



}"


Belltb <- "model{


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



}"





Bell1wctb <- "model{


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


}"


Bellwctb <- "model{


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



}"
  ############################
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  mfZ <- stats::model.frame(formSurv, data = tmp)
  XS <- stats::model.matrix(formSurv, mfZ)[, -1]
  ########  Gauss-Legendre quadrature (15 points)  ########
  K <- 15
  glq <- statmod::gauss.quad(K, kind = "legendre")
  xk <- glq$nodes # Nodes
  wk <- glq$weights # Weights
  K <- length(xk) # K-points
  ################
  time_new=dataLong[obstime]
  data_Long_s <- dataLong[time_new <= s, ]
  # data_Long_s <- dataLong[dataLong$obstime <= s, ]
  data_long <- data_Long_s[unique(c(
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


  if (family == "Beta") {
    y[y == 1] <- 0.9999


    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]



    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")
    # "betaL1", "betaL2", "betaS", "Sigmaa", "Sigmab", "gamma_pi", "gamma_lambda", "phis", "h")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
    Sigmab <- mean(object$MCMC$Sigmab)
    }else{
    Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    phis <- mean(object$MCMC$phi)
    h <- apply(object$MCMC$h, 2, mean)

    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Beta1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(Beta1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(Betab)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(Betatb)
      Omegab = 1/Sigmab
      }
    }


    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda, phis = phis,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )



    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }

  if (family == "Gamma") {
    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    sigma <- mean(object$MCMC$sigma)
    h <- apply(object$MCMC$h, 2, mean)



    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Gamma1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(Gamma1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(Gammab)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(Gammatb)
      Omegab = 1/Sigmab
      }
    }






    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      sigma1 = 1 / sigma,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )




    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }

  ####
  if (family == "Weibull") {
    Nbeta1 <- dim(X1)[2]
    Nbeta2 <- dim(X2)[2]

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    kappa <- mean(object$MCMC$kappa)
    h <- apply(object$MCMC$h, 2, mean)



    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Weibull1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(Weibull1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(Weibullb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(Weibulltb)
      Omegab = 1/Sigmab
      }
    }




    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      kappa = kappa,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )


    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }

  if (family == "Exponential") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)




    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Exp1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(Exp1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(Expb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(Exptb)
      Omegab = 1/Sigmab
      }
    }



    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )


    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }





  if (family == "inverse.gaussian") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)
    sigma <- mean(object$MCMC$sigma)



    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(IGauss1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(IGauss1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(IGaussb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(IGausstb)
      Omegab = 1/Sigmab
      }
    }



    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h, lambda = 1 / sigma,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )

    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }


  if (family == "Poisson") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)




    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Poisson1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(Poisson1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(Poissonb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(Poissontb)
      Omegab = 1/Sigmab
      }
    }


    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )

    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }
  #################
  if (family == "Logarithmic") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)


    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(logar1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){
      model.file <- textConnection(logar1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(logarb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){
      model.file <- textConnection(logartb)
      Omegab = 1/Sigmab
      }
    }


    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )


    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }


  if (family == "Binomial") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)



    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(binomial1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(binomial1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(binomialb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(binomialtb)
      Omegab = 1/Sigmab
      }
    }

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )


    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }
  #################
  if (family == "Bell") {
    C <- c()
    for (i in 1:length(y)) {
      C[i] <- log(numbers::bell(y[i])) - lfactorial(y[i])
    }


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)






    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Bell1b)
      betaS <- mean(object$MCMC$beta3)
      if (is.infinite(numbers::bell(max(y))) == TRUE) {
        model.file <- textConnection(Bell1wcb)
      }
      if(dim(Z2)[2]==1){
        model.file <- textConnection(Bell1tb)
        Omegab = 1/Sigmab
        if (is.infinite(numbers::bell(max(y))) == TRUE) {
          model.file <- textConnection(Bell1wctb)
        }
      }

    } else {
      model.file <- textConnection(Bellb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if (is.infinite(numbers::bell(max(y))) == TRUE) {
        model.file <- textConnection(Bellwcb)
      }
      if(dim(Z2)[2]==1){
        model.file <- textConnection(Belltb)
        Omegab = 1/Sigmab
        if (is.infinite(numbers::bell(max(y))) == TRUE) {
          model.file <- textConnection(Bellwctb)
        }
      }
    }

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K, C = C
    )

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
      DIC = FALSE
    )


    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }

  #########################
  if (family == "Gaussian") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }
    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)
    sigma <- mean(object$MCMC$sigma)




    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(Gaussian1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(Gaussian1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(Gaussianb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(Gaussiantb)
      Omegab = 1/Sigmab
      }
    }


    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h, tau = 1 / sigma,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2),
      death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2),
      id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )


    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }

  if (family == "NB") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }
    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)
    r <- mean(object$MCMC$r)



    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(NB1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(NB1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(NBb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(NBtb)
      Omegab = 1/Sigmab
      }
    }


    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h, r = r,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )

    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }

  if (family == "GP") {
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]

    i.jags <- function() {
      list(
        a = matrix(0, nTime, Nb1), b = matrix(0, nTime, Nb1)
      )
    }

    parameters <- c("a", "b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    Sigmaa <- apply(object$MCMC$Sigmaa, c(2, 3), mean)

    if(dim(Z2)[2]==1){
      Sigmab <- mean(object$MCMC$Sigmab)
    }else{
      Sigmab <- apply(object$MCMC$Sigmab, c(2, 3), mean)
    }

    gamma_pi <- mean(object$MCMC$gamma_pi)
    gamma_lambda <- mean(object$MCMC$gamma_lambda)
    h <- apply(object$MCMC$h, 2, mean)
    phiz <- mean(object$MCMC$phi)



    Omegab = solve(Sigmab)
    if (is.matrix(XS) == FALSE) {
      model.file <- textConnection(GP1b)
      betaS <- mean(object$MCMC$beta3)
      if(dim(Z2)[2]==1){model.file <- textConnection(GP1tb)
      Omegab = 1/Sigmab
      }
    } else {
      model.file <- textConnection(GPb)
      betaS <- apply(object$MCMC$beta3, 2, mean)
      if(dim(Z2)[2]==1){model.file <- textConnection(GPtb)
      Omegab = 1/Sigmab
      }
    }

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omegaa = solve(Sigmaa), Omegab = Omegab, gamma_pi = gamma_pi, gamma_lambda = gamma_lambda,
      h = h, phiz = phiz,
      n = n1, zeros = rep(0, n1), n2 = n2, zeros2 = rep(0, n2), y = y, Time = rep(s, n2), death = death, KF1 = 100000, KF2 = 100000,
      indtime1 = indtime1, indtime2 = indtime2,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Xv1 = Xv1, Xv2 = Xv2,
      Nb1 = Nb1, Nb2 = Nb2, mub1 = rep(0, Nb1), mub2 = rep(0, Nb2), id = id_prim,
      XS = XS, nindtime1 = nindtime1, nindtime2 = nindtime2,
      s = peice, xk = xk, wk = wk, K = K
    )

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
      DIC = FALSE
    )
    a_hat <- sim1$mean$a
    b_hat <- sim1$mean$b

    a_sim <- sim1$sims.list$a
    b_sim <- sim1$sims.list$b
  }
  a <- a_hat
  b <- b_hat
  #######################
  ################################
  ############################################
  step <- function(x) {
    z <- 0
    if (x >= 0) (z <- 1)
    z
  }
  inprod <- function(a, b) {
    z <- a %*% b
    z
  }
  Alpha0 <- Alpha1 <- Surv_d <- c()
  chaz <- matrix(0, n2, K)
  for (k in 1:n2) {
    if (is.matrix(XS) == FALSE) {
      if(dim(Z2)[2]==1){
        Alpha0[k] <- betaS * XS[k] + gamma_lambda * (inprod(betaL1[nindtime1], Xv1[k, ]) + a[k, 1]) +
          gamma_pi * (inprod(betaL2[nindtime2], Xv2[k, ]) + b[k])
      }else{
        Alpha0[k] <- betaS * XS[k] + gamma_lambda * (inprod(betaL1[nindtime1], Xv1[k, ]) + a[k, 1]) +
          gamma_pi * (inprod(betaL2[nindtime2], Xv2[k, ]) + b[k, 1])
      }
    } else {
      if(dim(Z2)[2]==1){
        Alpha0[k] <- inprod(betaS[], XS[k, ]) + gamma_lambda * (inprod(betaL1[nindtime1], Xv1[k, ]) + a[k, 1]) +
          gamma_pi * (inprod(betaL2[nindtime2], Xv2[k, ]) + b[k])
      }else{
        Alpha0[k] <- inprod(betaS[], XS[k, ]) + gamma_lambda * (inprod(betaL1[nindtime1], Xv1[k, ]) + a[k, 1]) +
          gamma_pi * (inprod(betaL2[nindtime2], Xv2[k, ]) + b[k, 1])
      }
    }

    if(dim(Z2)[2]==1){
      Alpha1[k] <- gamma_lambda * (betaL1[indtime1] + a[k, 2]) + gamma_pi * (betaL2[indtime2])
    } else{
      Alpha1[k] <- gamma_lambda * (betaL1[indtime1] + a[k, 2]) + gamma_pi * (betaL2[indtime2] + b[k, 2])
    }




    xk11 <- wk11 <- c()

    for (j in 1:K) {
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[j] <- (xk[j] + 1) / 2 * s
      wk11[j] <- wk[j] * s / 2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k, j] <- ((h[1] * step(peice[1] - xk11[j])) +
        (h[2] * step(xk11[j] - peice[1]) * step(peice[2] - xk11[j])) +
        (h[3] * step(xk11[j] - peice[2]) * step(peice[3] - xk11[j])) +
        (h[4] * step(xk11[j] - peice[3]) * step(peice[4] - xk11[j])) +
        (h[5] * step(xk11[j] - peice[4]))) * exp(Alpha0[k] + Alpha1[k] * xk11[j])
    }


    Surv_d[k] <- exp(-inprod(wk11, chaz[k, ]))
  }
  ########
  Surv_n <- c()
  chaz <- matrix(0, n2, K)
  for (k in 1:n2) {
    xk11 <- wk11 <- c()

    for (j in 1:K) {
      # Scaling Gauss-Kronrod/Legendre quadrature
      xk11[j] <- (xk[j] + 1) / 2 * (s + Dt)
      wk11[j] <- wk[j] * (s + Dt) / 2
      #  Hazard function at Gauss-Kronrod/Legendre nodes
      chaz[k, j] <- ((h[1] * step(peice[1] - xk11[j])) +
        (h[2] * step(xk11[j] - peice[1]) * step(peice[2] - xk11[j])) +
        (h[3] * step(xk11[j] - peice[2]) * step(peice[3] - xk11[j])) +
        (h[4] * step(xk11[j] - peice[3]) * step(peice[4] - xk11[j])) +
        (h[5] * step(xk11[j] - peice[4]))) * exp(Alpha0[k] + Alpha1[k] * xk11[j])
    }

    Surv_n[k] <- exp(-inprod(wk11, chaz[k, ]))
  }
  Surv_d[Surv_d == 0] <- 0.000001
  DP <- 1 - Surv_n / Surv_d
  #####################
  DP_last <- cbind(unique(id), DP)
  colnames(DP_last) <- c("id", "est")
  DP_last <- data.frame(DP_last)

  list(DP = DP_last, s = s, t = Dt)
}
