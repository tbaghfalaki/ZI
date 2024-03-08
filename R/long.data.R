#' Simulated longitudinal data
#'
#' Simulated data were generated in which "Y1" follows a Poisson distribution, "obstime" represents the observed time with two covariates x1 and x2.
#'
#'
#' @name long.data
#' @format A data frame which contains id, obstime, x1, x2, and Y1.
#' \describe{
#'   \item{id}{patients identifier}
#'   \item{obstime}{observed times for longitudinal measerements}
#'   \item{Y1}{the longitudinal response}
#'   \item{x1}{a continoues covariate}
#'   \item{x2}{a Binary covariate}
#' }
#' @seealso \code{\link{ZI}}
"long.data"
