#' Control parameters for emfrail
#'
#' @param eps Convergence criterion for the inner loops (the EM algorithm) for fixed frailty parameters
#' @param maxit Maximum number of iterations in the inner loops (the EM algorithm)
#' @param opt_fit Logical. Whether the outer optimization should be carried out.
#' If \code{FALSE}, then the frailty parameter is treated as fixed and the \code{emfrail} function returns only log-likelihood. See details.
#' @param verbose Logical. Whether to print out information about what is going on during maximization.
#' @param fast_fit Logical. Whether to try to calculate the E step directly, when possible. See details.
#' @param se_fit Logical. Whether to calculate the variance / covariance matrix.
#' @param se_adj Logical. Whether to calculate the adjusted variance / covariance matrix (needs \code{se_fit == TRUE})
#' @param ca_test Logical. Should the Commenges-Andersen test be calculated?
#' @param only_ca_test Logical. Should ONLY the Commenges-Andersen test be calculated?
#' @param lik_ci Logical. Should likelihood-based confidence interval be calculated for the frailty parameter?
#' @param opt_control A list with arguments for the likelihood-based confidence interval and the maximizer.
#' The list should contain two length 2 vectors \code{interval} and \code{interval_stable} that are used in calculating likelihood-based
#' confidence intervals. These are the edges, on the scale of \eqn{\theta}, of the parameter space where to look for the
#' ends of these confidence intervals.
#'
#' @return An object of the type \code{emfrail_control}.
#' @export
#'
#' @details The \code{fast_fit} option make a difference when the distribution is gamma (with or without left truncation) or
#' inverse Gaussian, i.e. pvf with m = -1/2 (without left truncation). For all the other scenarios, the fast_fit option will
#' automatically be changed to FALSE. When the number of events in a cluster / individual is not very small, the cases for which
#' fast fitting is available will show an improvement in performance.
#'
#' @seealso \code{\link{emfrail}}, \code{\link{emfrail_distribution}}, \code{\link{emfrail_pll}}
#' @examples
#' emfrail_control()
#' emfrail_control(eps = 10e-7)
#'

emfrail_control <- function(eps = 0.0001, maxit = Inf, opt_fit = TRUE, verbose = FALSE, fast_fit = TRUE,
                            se_fit = TRUE, se_adj = TRUE, ca_test = TRUE, only_ca_test = FALSE,
                            lik_ci = TRUE,
                            opt_control = list(interval = c(-7, 20),
                                               interval_stable = c(0, 20))) {
    # calculate SE as well

    # Here some checks

  if(isTRUE(only_ca_test) & !isTRUE(ca_test)) stop("control: if only_ca_test is TRUE then ca_test must be TRUE as well")
  if(is.null(opt_control$interval)) stop("opt_control must be a list which contains a named element interval")
  if(length(opt_control$interval) != 2) stop("interval must be of length 2")
  if(opt_control$interval[1] < -7 | opt_control$interval[2] > 20) warning("extreme values for interval, there might be some numerical trouble")

    res <- list(eps = eps,
                maxit = maxit,
                opt_fit = opt_fit,
                verbose = verbose,
                fast_fit = fast_fit,
                se_fit = se_fit,
                se_adj = se_adj,
                ca_test = ca_test,
                only_ca_test = only_ca_test,
                lik_ci = lik_ci,
                opt_control = opt_control)
    attr(res, "class") <- c("emfrail_control")
    res
}





#' Distribution parameters for emfrail
#'
#' @param dist One of 'gamma', 'stable' or 'pvf'.
#' @param theta A starting value for the 'outer' maximization with respect to the frailty parameter \eqn{\theta}. Must be >0.
#' @param pvfm Only relevant if \code{dist = 'pvf'} is used. It determines which PVF distribution should be used. Must be  larger than -1 and not equal to 0.
#' @param left_truncation Logical. Whether the data set represents left truncated survival times.
#'
#' @return An object of the type \code{emfrail_distribution}, which is mostly used to denote the
#' supported frailty distributions in a consistent way.
#' @export
#'
#' @details The \code{theta} argument must be positive. In the case of gamma or PVF, this is the inverse of
#'  the frailty variance, i.e. the larger the \code{theta} is,
#'  the closer the model is to a Cox model. For the positive stable distribution, the \eqn{\gamma} parameter of the Laplace trnasform is
#'  \eqn{\theta / (1 + \theta)}, with the \eqn{alpha} parameter fixed to 1.
#'
#' @seealso \code{\link{emfrail}, \link{emfrail_control}}
#' @examples
#' emfrail_distribution()
#' # Compound Poisson distribution:
#' emfrail_distribution(dist = 'pvf', theta = 1.5, pvfm = 0.5)
#' # Inverse Gaussian distribution:
#' emfrail_distribution(dist = 'pvf')
emfrail_distribution <- function(dist = "gamma", theta = 2, pvfm = -1/2, left_truncation = FALSE) {

  if (!(dist %in% c("gamma", "stable", "pvf")))
    stop("frailty distribution must be one of gamma, stable, pvf")
  if (length(theta) != 1)
    stop("specify exactly 1 parameter (theta>0) for the frailty")
  if (theta <= 0)
    stop("frailty parameter (theta) must be positive")
  if (dist == "pvf" & (pvfm < -1 | pvfm == 0))
    stop("pvfm must be >-1 and not equal to 0")

  if(!is.logical(left_truncation)) stop("left_truncation must be TRUE or FALSE")
  res <- list(dist = dist, theta = theta, pvfm = pvfm, left_truncation = left_truncation)
  attr(res, "class") <- c("emfrail_distribution")
  return(res)
}



