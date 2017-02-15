#' Laplace transform calculation
#'
#' @param x A vector of positive values where to calculate the Laplace transform
#' @param .distribution An \code{emfrail_distribution} object. See \code{?emfrail_distribution}.
#'
#' @return A vector of the same length as \code{x} with the Laplace transform of \code{x}
#'
#' @details This is a simple function which calculates the Laplace transform for the gamma, positive stable and PVF distribution.
#' It is intended to be used to calculate marginal quantities from an \code{emfrail} object.
#' Note that the \code{left_truncation} argument is ignored here;
#' the marginal survival or hazard are given for the Laplace transform of a baseline subject entered at time 0.
#' #'
#' @examples
#' dist <- emfrail_distribution()
#' laplace_transform(rexp(10), dist)
#'
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
