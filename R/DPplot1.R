#'  Dynamic prediction plot for one individual
#'
#' @description
#' Dynamic prediction plot for ZIJMCV
#'
#'
#' @details
#' Estimate DP for joint modeling based on ZIJMCV
#'
#' @param object an object inheriting from class ZIJMCV
#' @param dataLong data set of observed longitudinal variables.
#' @param dataSurv data set of observed survival variables.
#' @param s the landmark time for prediction
#' @param id_new id number for individual who want to plot his/her DP
#' @param mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
#' @param by number: increment of the sequence of DP.
#' @param Marker_lab the label for the response axis
#' @param Time_lab the label for the time axis
#' @param digits integer indicating the number of decimal places for Time axis
#' @param n.chains the number of parallel chains for the model; default is 1.
#' @param n.iter integer specifying the total number of iterations; default is 1000.
#' @param n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
#'
#'
#' @return
#' - A plot for DP
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}
#'
#' @example inst/exampleDP1.R
#'
#' @md
#' @export
#'
DPplot1 <- function(object, s = s, id_new = id_new, mi = mi,
                    Marker_lab="Marker", Time_lab="Time",digits=digits,
                    by = 0.1, n.chains = n.chains, n.iter = n.iter, n.burnin = floor(n.iter / 2),
                    dataLong, dataSurv) {

  FixedY <- object$FixedY
  FixedZ <- object$FixedZ
  RandomY <- object$RandomY
  RandomZ <- object$RandomZ
  GroupY <- object$GroupY
  GroupZ <- object$GroupZ
  obstime=object$obstime
  id=object$id
  time_new <- dataLong[obstime]

  data_Long_s <- dataLong[time_new <= s, ]
  data_long <- data_Long_s[unique(c(
    all.vars(GroupY), all.vars(FixedY), all.vars(RandomY),
    all.vars(GroupZ), all.vars(FixedZ), all.vars(RandomZ)
  ))]
  y <- data_long[all.vars(FixedY)][, 1]
  id_prime <- as.integer(data_long[all.vars(GroupY)][, 1])
  time_new <- as.numeric(dataLong[obstime][, 1])

  Data_new <- data_long[id_prime == id_new, ]

  y_new <- Data_new[all.vars(FixedY)][, 1]
  time_y <- as.numeric(Data_new[obstime][, 1])

  #################### survival time #########
  formSurv <- object$formSurv
  tmp <- dataSurv[all.vars(formSurv)]
  Time <- tmp[all.vars(formSurv)][, 1] # matrix of observed time such as Time=min(Tevent,Tcens)
  death <- tmp[all.vars(formSurv)][, 2] # vector of event indicator (delta)

  Time_new <- data.frame(unique(id_prime), Time)
  colnames(Time_new) <- c("id", "surv")
  surt <- Time_new[Time_new$id == id_new, ][2]


  Dt <- seq(s, max(time_new), by = by)

  est_M <- matrix(0, length(Dt), 4)

  for (kk in 1:length(Dt)) {
    DD <- DP_CV_CI(
      object = object, s = s, t = Dt[kk]-s, mi = mi, n.chains = 1, n.iter = n.iter, n.burnin = n.burnin,
      n.thin = 1, dataLong = dataLong, dataSurv = dataSurv
    )

    est_M[kk, ] <- as.numeric(DD$DP[DD$DP$id == id_new, ])
  }

  Tab1 <- cbind(Dt, 1 - est_M[, -1])
  colnames(Tab1) <- c("obstime", "mean", "uzi", "l")

  ytab1 <- cbind(time_y, y_new)
  colnames(ytab1) <- c("obstime", "y")


  TT <- merge(Tab1, ytab1, by = "obstime", all = TRUE)

  #xlab <- seq(from = min(time_new), to = max(time_new), length = 5)
  xlab <- seq(from = floor(min(time_y)), to = ceiling(max(time_new)), length = 5)

  par(mar = c(5, 5, 2, 5))

  plot(TT[, 1], TT[, 5],
       type = "p", pch = 20, xaxt = "n", col = "black",
       xlim=c(min(time_y),max(time_new)),
       ylab = Marker_lab, xlab = Time_lab,
       main = ""
  )
  axis(side = 1, round(xlab, digits =digits), labels = TRUE)

  abline(v = surt, col = "red", lty = 3)
  abline(v = s, col = "blue", lty = 3)

  par(new = T)
  plot(TT[, 1:2], type = "l", pch = 16, col = "purple3", axes=F, xlab = NA, ylab = NA, cex = 1, ylim = c(0, 1)) # axes=F,
  lines(TT[, c(1, 3)], type = "l", col = "purple3", lty = 2, pch = 16, xlab = NA, ylab = NA, cex = 1)
  lines(TT[, c(1, 4)], type = "l", col = "purple3", lty = 2, pch = 16, xlab = NA, ylab = NA, cex = 1)

  axis(side = 4)
  mtext(side = 4, line = 3, "Probability of Event")
}
