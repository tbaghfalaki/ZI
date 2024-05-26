#' Simulated survival data
#'
#' Simulated survival data were generated in which "survtime" follows a survival time with a censoring indicator called "death" and a continuous covariate, "w1".
#'
#'
#' @name surv_data_b
#' @format A data frame which contains id, survtime, death, w1.
#' \describe{
#'   \item{id}{patients identifier}
#'   \item{survtime}{survival time (the response variable)}
#'   \item{death}{censoring indicator, 1=observed, 0=censored}
#'   \item{w1}{a continuous covariate}
#'   \item{w2}{a binary covariate}
#' }
#' @seealso \code{\link{UHJM}}
"surv_data_b"


