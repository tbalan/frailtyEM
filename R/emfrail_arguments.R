
#' Distribution parameters for emfrail
#'
#' @param dist One of 'gamma', 'stable' or 'pvf'.
#' @param frailtypar A starting value for the 'outer' maximization with respect to the frailty parameter \eqn{\theta}. Must be positive.
#' @param pvfm Only relevant if \code{dist = 'pvf'} is used. It determines which PVF distribution should be used. Must be  larger than -1 and not equal to 0.
#'
#' @return A list of arguments suitable for \code{emfrail()}!
#' @export
#'
#' @note The \code{frailtypar} argument must be positive.
#' @seealso \code{\link{emfrail}, \link{emfrail_control}}
#' @examples
#' emfrail_distribution()
#' emfrail_distribution(dist = 'pvf', frailtypar = 1.5, pvfm = 0.5)
emfrail_distribution <- function(dist = "gamma", frailtypar = 2, pvfm = -1/2) {

    if (!(dist %in% c("gamma", "stable", "pvf")))
        stop("frailty distribution must be one of gamma, stable, pvf")
    if (length(frailtypar) != 1)
        stop("specify exactly 1 parameter (theta) for the frailty")
    if (frailtypar <= 0)
        stop("frailty parameter (theta) must be positive")
    if (dist == "pvf" & (pvfm < -1 | pvfm == 0))
        stop("pvfm must be >-1 and not equal to 0")


    res <- list(dist = dist, frailtypar = frailtypar, pvfm = pvfm)
    attr(res, "class") <- c("emfrail_distribution")
    return(res)
}




#' Control parameters for emfrail
#'
#' @param eps Convergence criterion for the inner loops (the EM algorithm) for fixed frailty parameters
#' @param maxit Maximum number of iterations in the inner loops (the EM algorithm)
#' @param opt_fit Logical. Whether the outer optimization should be carried out. If \code{FALSE}, then the frailty parameter is treated as fixed.
#' @param verbose Logical. Whether to print out information about what is going on during maximization.
#' @param fast_fit Logical. Whether to go the faster route in the E step and calculate directly (when possible)
#' @param zerotol Lower limit for 1/frailtypar (variance in case of gamma / pvf). Below this value, the variance is taken to be 0 and the
#' EM is not actually performed. The log-likelihood returned by em_fit is that for the Cox model.
#' @param opt_control A list with arguments to be sent to the maximizer. Currently ignored, will be useful in the future.
#'
#' @return A list of arguments suitable for \code{emfrail()}!
#' @export
#' @seealso \code{\link{emfrail}, \link{emfrail_distributon}}
#' @examples
#' emfrail_control()
#' emfrail_control(eps = 10e-7)
emfrail_control <- function(eps = 0.001, maxit = Inf, opt_fit = TRUE, verbose = FALSE, fast_fit = TRUE,
                            zerotol = 1e-4,
                            opt_control = list(method = "bobyqa",
    itnmax = NULL, control = list())) {
    # calculate SE as well

    # Here some checks

    res <- list(eps = eps, maxit = maxit, opt_fit = opt_fit, zerotol = zerotol, verbose = verbose, fast_fit = fast_fit, opt_control = opt_control)
    attr(res, "class") <- c("emfrail_control")
    res
}
