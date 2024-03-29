#' Zero-inflation joint modeling
#'
#' @description
#' Fits zero-inflated hurdle joint models under Poisson, negative binomial, and generalized Poisson distributions by considering shared random effects model.
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
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#' @param n.thin integer specifying the thinning of the chains; default is 1.
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
#' - Sigma the variance of random effects
#' - gamma the association parameters
#' - p the scale parameter of Weibull model
#' - r the over-dispersion parameter of the negative binomial model
#' - phiz the dispersion parameter of the generalized Poisson model
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZISRE.R
#'
#' @md
#' @export

ZISRE <- function(FixedY, RandomY, GroupY, FixedZ, RandomZ, GroupZ, formSurv, dataLong, dataSurv,
                  n.chains = n.chains,
                  n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, family = "Poisson") {
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






  #####################
  #formSurv=survival::formSurv
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
  ###########
  NB <- "model{


  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])
  }
    #####
   for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }

 Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))

   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

}"

  ###########
  Poisson <- "model{


    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))


     log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])


    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)
    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"

  ###########
  GP <- "model{

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-y[i]*log(1+phiz*lambda[i])+(y[i]-1)*log(1+phiz*y[i])- loggam(y[i]+1)-
        lambda[i]*(1+phiz*y[i])/(1+phiz*lambda[i])-log(1-exp(-lambda[i]/(1+phiz*lambda[i]))))

      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])

    }

for(k in 1:n2){
    Time[k] ~ dweib(p,mut[k])
    is.censored[k]~dinterval(Time[k],surt.cen[k])
    log(mut[k])<-inprod(betaS[],XS[k,])+inprod(gamma[],b[k,])
    b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])
   }
    for(l in 1:NbetaS){
  betaS[l]~dnorm(0,0.001)
  }

  p~ dgamma(0.1,0.1)

    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }
      phiz~dnorm(0,.001)

  Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))

   for(k in 1:(Nb1+Nb2)){
gamma[k]~dnorm(0,0.001)
}

  }"


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
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "Sigma", "gamma", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 1000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
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
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p
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



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
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



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
    ###############################################
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
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "Sigma", "gamma", "r", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 1000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
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
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
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
        names(sim1$Rhat$betaS) <- colnames(XS)



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")


      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r, sim1$Rhat$r)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
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



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$r) <-
        names(sim1$sd$r) <-
        names(sim1$q2.5$r) <-
        names(sim1$q97.5$r) <-
        names(sim1$Rhat$r) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p),
        cbind(sim1$mean$r, sim1$sd$r, sim1$q2.5$r, sim1$q97.5$r)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }

    #########
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
        betaS = stats::rnorm(NbetaS), phiz = 1,
        Omega = diag(Nb1 + Nb2), gamma = stats::rnorm(Nb1 + Nb2)
      )
    }

    parameters <- c("betaL1", "betaL2", "betaS", "Sigma", "gamma", "phiz", "p")


    d.jags <- list(
      n = n1, zeros = rep(0, n1), n2 = n2, y = y, Time = Time, surt.cen = surt.cen, K = 1000,
      X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      NbetaS = NbetaS,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, Nb1 + Nb2), V = diag(1, Nb1 + Nb2), id = id_prim,
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
      DIC = TRUE
    )

    MCMC <- list(
      beta1 = sim1$sims.list$betaL1, beta2 = sim1$sims.list$betaL2,
      beta3 = sim1$sims.list$betaS,
      Sigma = sim1$sims.list$Sigma,
      gamma = sim1$sims.list$gamma,
      p = sim1$sims.list$p,
      r = sim1$sims.list$phiz
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



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <-
        names(sim1$Rhat$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <-
        names(sim1$Rhat$p) <- "Scale"


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        rownames(sim1$Rhat$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <-
        colnames(sim1$Rhat$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2, sim1$Rhat$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1, sim1$Rhat$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")


      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS, sim1$Rhat$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma, sim1$Rhat$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p, sim1$Rhat$p),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz, sim1$Rhat$phiz)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI", "Rhat")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
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



      names(sim1$mean$gamma) <-
        names(sim1$sd$gamma) <-
        names(sim1$q2.5$gamma) <-
        names(sim1$q97.5$gamma) <- c(colnames(Z1), colnames(Z2))


      names(sim1$mean$p) <-
        names(sim1$sd$p) <-
        names(sim1$q2.5$p) <-
        names(sim1$q97.5$p) <- "Scale"


      names(sim1$mean$phiz) <-
        names(sim1$sd$phiz) <-
        names(sim1$q2.5$phiz) <-
        names(sim1$q97.5$phiz) <-
        names(sim1$Rhat$phiz) <- "Dispersion"


      rownames(sim1$mean$Sigma) <-
        rownames(sim1$sd$Sigma) <-
        rownames(sim1$q2.5$Sigma) <-
        rownames(sim1$q97.5$Sigma) <-
        colnames(sim1$mean$Sigma) <-
        colnames(sim1$sd$Sigma) <-
        colnames(sim1$q2.5$Sigma) <-
        colnames(sim1$q97.5$Sigma) <- c(colnames(Z1), colnames(Z2))



      ZPM <- cbind(sim1$mean$betaL2, sim1$sd$betaL2, sim1$q2.5$betaL2, sim1$q97.5$betaL2)
      colnames(ZPM) <- c("Est", "SD", "L_CI", "U_CI")



      MM <- cbind(sim1$mean$betaL1, sim1$sd$betaL1, sim1$q2.5$betaL1, sim1$q97.5$betaL1)
      colnames(MM) <- c("Est", "SD", "L_CI", "U_CI")



      SM <- rbind(
        cbind(sim1$mean$betaS, sim1$sd$betaS, sim1$q2.5$betaS, sim1$q97.5$betaS),
        cbind(sim1$mean$gamma, sim1$sd$gamma, sim1$q2.5$gamma, sim1$q97.5$gamma),
        cbind(sim1$mean$p, sim1$sd$p, sim1$q2.5$p, sim1$q97.5$p),
        cbind(sim1$mean$phiz, sim1$sd$phiz, sim1$q2.5$phiz, sim1$q97.5$phiz)
      )

      colnames(SM) <- c("Est", "SD", "L_CI", "U_CI")

      D <- sim1$mean$Sigma

      results <- list(Count_model = MM, Zero_inflated_model = ZPM, Survival_model = SM, D = D)
    }
  }

  ############

  K <- 1000
  DIC <- sim1$DIC - 2 * n * K

  list(
    FixedY = FixedY, FixedZ = FixedZ, formSurv = formSurv,
    RandomY = RandomY, RandomZ = RandomZ,
    family = family, MCMC = MCMC, Estimation = results,
    DIC = DIC
  )
}
