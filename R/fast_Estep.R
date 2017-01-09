#' Fast fitting of the E step
#'
#' This function calculates the E step in a quicker way, without taking all the derivatives of the Laplace transform
#' For a data set with \code{K} clusters,
#' @param c Vector of length \code{K} of cumulative hazards, i.e. total accumulated hazards within a cluster
#' @param delta Vector of integers of length \code{K} of the number of events for each cluster
#' @param alpha,bbeta Parameters of the frailty distribution
#' @param pvfm Parameter for the PVF distribution, only matters in that case
#' @param dist One of 0 (for gamma), 1 (for stable) or 2 (for PVF)
#'
#' @return A \code{K x 3} matrix where the first column and the second column are the numerators
#' and the denominators of the frailty fraction (without the Laplace transform) and the
#' last column is the log(denominator) + log-Laplace transform, i.e. the log-likelihood contribution
#'
#'
#' @export
fast_Estep <- function(c, c_lt = 0, delta, alpha, bbeta, pvfm, dist) {

  if(dist != 0) stop("no fast option available here")
  res <- matrix(0, length(delta), 3)

  if(dist==0) {
    bbeta <- bbeta + c_lt
    res[,3] <- alpha * log(bbeta) - (alpha + delta)*log(bbeta + c) + lgamma(alpha + delta) - lgamma(alpha)
    res[,1] <- (alpha + delta)
    res[,2] <- (bbeta + c)
  }

  res
}
