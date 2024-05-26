#'  Dynamic prediction
#'
#' @description
#' Dynamic prediction for ZISRE
#'
#'
#' @details
#' Estimate DP for joint modeling based on VS
#'
#' @param object an object inheriting from class VS
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param s the landmark time for prediction
#' @param t the window of prediction for prediction
#' @param offset the offset or library size for discrete response. If offset=NULL, it is considered without an offset.
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
#' @example inst/exampleDP2.R
#'
#' @md
#' @export

DP_SRE <- function(object, s = s, t = t, offset = NULL, n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                   n.thin = max(1, floor((n.iter - n.burnin) / 1000)), dataLong, dataSurv) {
  Dt <- t
  KK <- 1000000

  FixedY <- object$FixedY
  FixedZ <- object$FixedZ
  RandomY <- object$RandomY
  RandomZ <- object$RandomZ
  IStructure <-  object$IStructure


  GroupY <- object$GroupY
  GroupZ <- object$GroupZ
  id <- object$id
  formSurv <- object$formSurv
  obstime <- object$obstime
  peice <- object$peice
  family <- object$family

  #offset=dataLong[offset]


  ###########
  NBb <- "model{

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))

      log(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])
  }
    #####
   for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }


}"

  ###########
  Poissonb <- "model{


    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))


     log(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }


  }"

  ###########
  Bellb <- "model{


    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

       ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])+C[i]-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(theta[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }



  }"


  Bellwcb <- "model{


    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

       ll[i]<-(1-z[i])*(y[i]*log(theta[i])+1-exp(theta[i])-(log(1-exp(1-exp(theta[i])))))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

     log(theta[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }



  }"
  ###########
  logarb <- "model{


    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i]<-(1-z[i])*(log(-1/log(1-pi[i]))+y[i]*log(pi[i])-log(y[i]))+z[i]*log(lambda[i])+
  (1-z[i])*log(1-lambda[i])

     logit(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])




    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }




  }"

  ###########
  binomialb <- "model{

  m<-max(y)

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i]<-(1-z[i])*(loggam(m+1)-loggam(y[i]+1)-loggam(m-y[i]+1)+y[i]*log(pi[i])+(m-y[i])*log(1-pi[i])-
                   log(1-pow(1-pi[i],m)))+z[i]*log(lambda[i])+(1-z[i])*log(1-lambda[i])


     logit(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }



  }"
  ###########
  GPb <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- offset[i]+inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])

    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }



  }"

  Expb <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i] <- (1-z[i])*(logdensity.exp(y[i],lambda[i]))+z[i]*log(pi[i])+(1-z[i])*log(1-pi[i])


     log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }


  }"


  ######
  Gammab <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.gamma(y[i],sigma1, mu1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

 mu1[i]<-sigma1/mu[i]
     log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }

  }"
  #####

  Betab <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.beta(y[i],a1[i], b1[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])

a1[i]<-mu[i]*phi11
    b1[i]<-(1-mu[i])*phi11

     logit(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }



  }"
  #####
  Weibullb <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

ll[i] <- (1-z[i])*(logdensity.weib(y[i],kappa, mu[i]))+z[i]*log(muz[i])+(1-z[i])*log(1-muz[i])


     log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }



  }"

  Gaussianb <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i] <- (1-z[i])*(logdensity.norm(y[i], mu[i],tau))+z[i]*log(muz[i])+
(1-z[i])*log(1-muz[i])

    mu[i] <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }


  }"


  IGaussb <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

 ll[i]<- (1-z[i])*(0.5*(log(lambda)-log(2*3.14)-3*log(0.000000001+y[i]))-
 0.5*lambda*pow((y[i]-mu[i]),2)/(pow(mu[i],2)*(0.000000001+y[i]))+log(1-muz[i]))+z[i]*log(muz[i])

     log(mu[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(muz[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }


  }"

  ############################


  time_new=dataLong[obstime]
  data_Long_s <- dataLong[time_new <= s, ]
  data_long <- data_Long_s[unique(c(
    all.vars(GroupY), all.vars(FixedY), all.vars(RandomY),
    all.vars(GroupZ), all.vars(FixedZ), all.vars(RandomZ)
  ))]
  y <- data_long[all.vars(FixedY)][, 1]
  mfX <- stats::model.frame(FixedY, data = data_long)
  X1 <- stats::model.matrix(FixedY, mfX)
  mfU <- stats::model.frame(RandomY, data = data_long)
  Z1 <- stats::model.matrix(RandomY, mfU)
  # id_prim <- as.integer(data_Long_s[id][, 1])
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


  if (length(offset) == 0) {
    offset <- rep(0, n1)
  } else {
    offset <- dataLong[offset][, 1]
  }

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


  ### $$$$$$$$$$$
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)
  nTime <- length(Time) # number of subject having Time
  # design matrice
  mfZ <- stats::model.frame(formSurv, data = tmp)
  XS <- stats::model.matrix(formSurv, mfZ)



  surt.cen <- rep(0, n2)
  for (i in 1:n2) {
    if (death[i] == 1) (surt.cen[i] <- 10000)
    if (death[i] == 0) ((surt.cen[i] <- Time[i]) & (Time[i] <- NA))
  }
  ########################################################################################
  if (family == "Exponential") {
    model.file <- textConnection(Expb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
    Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
if(dim(Z2)[2]>1){
  Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
}
if(dim(Z2)[2]==1){

  Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

}

}



    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p,
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }
  ##
  #####
  if (family == "Beta") {
    model.file <- textConnection(Betab)
    y[y == 1] <- 0.99999

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    phi1 <- mean(object$MCMC$phi)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, phi11 = 1 / phi1,
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }
  ###
  if (family == "Weibull") {
    model.file <- textConnection(Weibullb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    kappa <- mean(object$MCMC$kappa)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, kappa = kappa,
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }
  ##
  if (family == "inverse.gaussian") {
    model.file <- textConnection(IGaussb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    sigma <- mean(object$MCMC$sigma)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, lambda = 1 / sigma,
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }
  ###
  if (family == "Gamma") {
    model.file <- textConnection(Gammab)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))
      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    sigma <- mean(object$MCMC$sigma)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, sigma1 = 1 / sigma,
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }
  ####
  if (family == "Gaussian") {
    model.file <- textConnection(Gaussianb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    sigma <- mean(object$MCMC$sigma)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, tau = 1 / sigma,
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }

  if (family == "Poisson") {
    model.file <- textConnection(Poissonb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)


    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p,
      n = n1, zeros = rep(0, n1), n2 = n2,  y = y, offset = offset, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }
  ####################### $$$$$$$$$$$$$
  if (family == "Bell") {
    C <- c()
    for (i in 1:length(y)) {
      C[i] <- log(numbers::bell(y[i])) - lfactorial(y[i])
    }

    model.file <- textConnection(Bellb)
    if (is.infinite(numbers::bell(max(y))) == TRUE) {
      model.file <- textConnection(Bellwcb)
    }


    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))
}
      }

    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p,
      n = n1, zeros = rep(0, n1), n2 = n2,  y = y, offset = offset, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      C = C,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }


  ####################### $$$$$$$$$$$$$
  if (family == "Logarithmic") {
    model.file <- textConnection(logarb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }}

    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p,
      n = n1, zeros = rep(0, n1), n2 = n2,  y = y, offset = offset, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }

  if (family == "Binomial") {
    model.file <- textConnection(binomialb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }}

    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p,
      n = n1, zeros = rep(0, n1), n2 = n2,  y = y, offset = offset, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
    )
    ###############################################
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }

  if (family == "NB") {
    model.file <- textConnection(NBb)

    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)

    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))

      }
}
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    r <- mean(object$MCMC$r)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, r = r,
      n = n1, zeros = rep(0, n1), n2 = n2,  y = y, offset = offset, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
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

    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }

  ############

  if (family == "GP") {
    model.file <- textConnection(GPb)
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    i.jags <- function() {
      list(
        b = matrix(0, n2, Nb1 + Nb2)
      )
    }

    parameters <- c("b")

    betaL1 <- apply(object$MCMC$beta1, 2, mean)
    betaL2 <- apply(object$MCMC$beta2, 2, mean)
    betaS <- apply(object$MCMC$beta3, 2, mean)
    if(IStructure==FALSE){
      Sigma <- apply(object$MCMC$Sigma, c(2, 3), mean)
    }else{
      if(dim(Z2)[2]>1){
        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),apply(object$MCMC$D2, c(2, 3), mean))
      }
      if(dim(Z2)[2]==1){

        Sigma <- Matrix::bdiag(apply(object$MCMC$D1, c(2, 3), mean),mean(object$MCMC$D2))
}
      }
    gamma <- apply(object$MCMC$gamma, 2, mean)
    p <- mean(object$MCMC$p)
    phiz <- mean(object$MCMC$r)

    d.jags <- list(
      betaL1 = betaL1, betaL2 = betaL2, betaS = betaS, Omega = solve(Sigma), gamma = gamma,
      p = p, phiz = phiz,
      n = n1, zeros = rep(0, n1), n2 = n2,  y = y, offset = offset, Time = rep(s, n2), surt.cen = surt.cen, K = 100000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), id = id_prim,
      XS = XS
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
    b_hat <- sim1$mean$b
    b_sim <- sim1$sims.list$b
  }

  ###############################

  inprod <- function(a, b) {
    z <- a %*% b
    z
  }
  ###############################
  b <- b_hat
  mut <- c()
  for (k in 1:n2) {
    mut[k] <- inprod(betaS[], XS[k, ]) + inprod(gamma[], b[k, ])
  }
  r <- mean(object$MCMC$p)


  Num <- 1 - pweibull((s + Dt), r, exp(-mut / r))
  Den <- 1 - pweibull(s, r, exp(-mut / r))

  Den[Den == 0] <- 0.0001


  DP <- 1 - Num / Den


  DP_last <- cbind(unique(id), DP)
  colnames(DP_last) <- c("id", "est")
  DP_last <- data.frame(DP_last)

  list(DP = DP_last, s = s, t = Dt)
}
