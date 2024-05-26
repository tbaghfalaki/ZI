#' Simulated longitudinal data
#'
#' Simulated data were generated, where "Y1" follows a hurdle exponential distribution, and "obstime" represents the observed time with two covariates, x1 and x2.
#'
#'
#' @name long_data_e
#' @format A data frame which contains id, obstime, x1, x2, and Y1.
#' \describe{
#'   \item{id}{patients identifier}
#'   \item{obstime}{observed times for longitudinal measurements}
#'   \item{Y1}{the longitudinal response}
#'   \item{x1}{a continuous covariate}
#'   \item{x2}{a Binary covariate}
#' }
#' @seealso \code{\link{UHJM}}
"long_data_e"
