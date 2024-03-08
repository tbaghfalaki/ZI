#' Zero-inflation mixed-effects models
#'
#' @description
#' Fits zero-inflated hurdle mixed-effects models under Poisson, negative binomial, and generalized Poisson distributions.
#'
#' @details
#' Function using 'JAGS' software to estimate the hurdle regression mixed-effects models
#'
#' @param FixedY formula for fixed part of longitudinal count model
#' @param RandomY formula for random part of longitudinal count model
#' @param GroupY formula specifying the cluster variable for Y (e.g. = ~ subject)
#' @param FixedZ formula for fixed part of longitudinal probability model
#' @param RandomZ formula for random part of longitudinal probability model
#' @param GroupZ formula specifying the cluster variable for Z (e.g. = ~ subject)
#' @param data dataset of observed variables.
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
#' - Sigma the variance of random effects
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

ZIMEM <- function(FixedY, RandomY, GroupY, FixedZ, RandomZ, GroupZ, data, n.chains = n.chains,
                  n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, family = "Poisson") {
  data_long <- data[unique(c(
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

  n <- length(data[, 1])
  z <- rep(0, n)
  z[y == 0] <- 1
  ## Number of patients and number of longitudinal observations per patient
  n1 <- length(y)
  Nbeta <- dim(X1)[2]
  ###########
  NB <- "model{
  K<-1000

  for(i in 1:n){
    zeros[i]~dpois(phi[i])
    phi[i]<-  - ll[i]+K

    ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+loggam(r+y[i])-loggam(r)-loggam(y[i]+1)+r*log(r/(r+lambda[i]))+y[i]*log(lambda[i]/(lambda[i]+r))-log(1-pow(r/(r+lambda[i]),r)))



      log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])



  }
       for(k in 1:n2){

        b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])

       }

  r~ dgamma(0.1,0.1)

  for(l in 1:Nbeta1){
    betaL1[l]~dnorm(0,0.001)
  }

  for(l in 1:Nbeta2){
    betaL2[l]~dnorm(0,0.001)
  }

 Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


}"

  ###########
  Poisson <- "model{
    K<-1000

    for(i in 1:n){
      zeros[i]~dpois(phi[i])
      phi[i]<-  - ll[i]+K

      ll[i]<-z[i]*log(pi[i]) +
 (1-z[i])*(log(1-pi[i])+y[i]*log(lambda[i])-lambda[i] - loggam(y[i]+1)-log(1-exp(-lambda[i])))


     log(lambda[i]) <- inprod(betaL1[],X1[i,])+inprod(b[id[i],1:Nb1],Z1[i,])
      logit(pi[i]) <-  inprod(betaL2[],X2[i,])+inprod(b[id[i],(Nb1+1):(Nb1+Nb2)],Z2[i,])




    }

for(k in 1:n2){

        b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])

       }



    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }


    Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))

  }"

  ###########
  GP <- "model{
    K<-1000

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

        b[k,1:(Nb1+Nb2)]~dmnorm(mub[],Omega[,])

       }


    for(l in 1:Nbeta1){
      betaL1[l]~dnorm(0,0.001)
    }

    for(l in 1:Nbeta2){
      betaL2[l]~dnorm(0,0.001)
    }
      phiz~dnorm(0,.001)

  Sigma[1:(Nb1+Nb2),1:(Nb1+Nb2)]<-inverse(Omega[,])
  Omega[1:(Nb1+Nb2),1:(Nb1+Nb2)]~dwish(V[,],(Nb1+Nb2))


  }"
  Nbeta1 <- dim(X1)[2]
  Nbeta2 <- dim(X2)[2]


  if (family == "Poisson") {
    model.file <- textConnection(Poisson)
    i.jags <- function() {
      list(betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2))
    }

    parameters <- c("betaL1", "betaL2", "Sigma")
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    d.jags <- list(
      n = n1, n2 = n2, zeros = rep(0, n1), y = y, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, (Nb1 + Nb2)), V = diag(1, (Nb1 + Nb2)), id = id_prim
    )

    sim1 <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
      n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
    )
  }
  ###############
  if (family == "NB") {
    model.file <- textConnection(NB)
    i.jags <- function() {
      list(betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2))
    }

    parameters <- c("betaL1", "betaL2", "Sigma", "r")
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    d.jags <- list(
      n = n1, n2 = n2, zeros = rep(0, n1), y = y, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1, Nbeta2 = Nbeta2,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, (Nb1 + Nb2)), V = diag(1, (Nb1 + Nb2)), id = id_prim
    )

    sim1 <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
      n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
    )
  }
  ###############
  if (family == "GP") {
    model.file <- textConnection(GP)
    i.jags <- function() {
      list(betaL1 = stats::rnorm(Nbeta1), betaL2 = stats::rnorm(Nbeta2))
    }

    parameters <- c("betaL1", "betaL2", "Sigma", "phiz")
    Nb1 <- dim(Z1)[2]
    Nb2 <- dim(Z2)[2]
    d.jags <- list(
      n = n1, n2 = n2, zeros = rep(0, n1), y = y, X1 = X1, X2 = X2, Z1 = Z1, Z2 = Z2, z = z, Nbeta1 = Nbeta1,
      Nbeta2 = Nbeta2,
      Nb1 = Nb1, Nb2 = Nb2, mub = rep(0, (Nb1 + Nb2)), V = diag(1, (Nb1 + Nb2)), id = id_prim
    )

    sim1 <- R2jags::jags(
      data = d.jags, inits = i.jags, parameters, n.chains = n.chains,
      n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = model.file
    )
  }
  ###############
  sim1
}
