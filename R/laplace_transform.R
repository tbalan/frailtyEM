#' Laplace transform
#'
#' @param x A vector of positive values where to calculate the Laplace transform
#' @param .distribution An \code{emfrail_distribution} object. See \code{?emfrail_distribution}.
#'
#' @return Returns a vector of the same length as \code{x} with the Laplace transform.
#'
#' @details This is a simple function which calculates the Laplace transform for the gamma, positive stable and PVF distribution.
#' It is intended to be used to calculate marginal quantities from an \code{emfrail} object.
#'
#' @export
#'
#' @examples
laplace_transform <- function(x, .distribution) {
  # if(missing(.distribution) & missing())
  if(!inherits(.distribution, "emfrail_distribution"))
    stop(".distribution argument misspecified; see ?emfrail_distribution()")

  getpars <- dist_to_pars(.distribution$dist, log(.distribution$frailtypar), .distribution$pvfm)

  if(getpars$dist == 0L) {
    L <- with(getpars, (bbeta / (bbeta + x))^alpha)
  }

  if(getpars$dist == 1L) {
    L <- with(getpars, exp(-1 * x^bbeta))
  }

  if(getpars$dist == 2L) {
    L <- with(getpars, exp(-alpha * sign(.distribution$pvfm) * (1 - (bbeta / (bbeta + x))^.distribution$pvfm )))
  }

  L

}
