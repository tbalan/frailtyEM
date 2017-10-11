#' Control parameters for emfrail
#'
#' @param opt_fit Logical. Whether the outer optimization should be carried out.
#' If \code{FALSE}, then the frailty parameter is treated as fixed and the \code{emfrail} function returns only log-likelihood. See details.
#' @param se Logical. Whether to calculate the variance / covariance matrix.
#' @param se_adj Logical. Whether to calculate the adjusted variance / covariance matrix (needs \code{se == TRUE})
#' @param ca_test Logical. Should the Commenges-Andersen test be calculated?
#' @param only_ca_test Logical. Should ONLY the Commenges-Andersen test be calculated?
#' @param lik_ci Logical. Should likelihood-based confidence interval be calculated for the frailty parameter?
#' @param lik_ci_intervals This list should contain two length 2 vectors \code{interval} and \code{interval_stable} that are used in calculating likelihood-based
#' confidence intervals. These are the edges, on the scale of \eqn{\theta}, of the parameter space where to look for the
#' ends of these confidence intervals.
#' @param nlm_control A list of named arguments to be sent to \code{nlm} for the outer optimization.
#' @param inner_control A list of parameters for the inner optimization. See details.
#'
#' @return An object of the type \code{emfrail_control}.
#' @export
#'
#' @details
#' The \code{nlm_control} argument should not overalp with \code{hessian}, \code{f} or \code{p}.
#'
#' The \code{inner_control} argument should be a list with the following items:
#' \itemize{
#' \item{\code{eps}}{ A criterion for convergence of the EM algorithm (difference between two consecutive values of the log-likelihood)}
#' \item{\code{maxit}}{ The maximum number of iterations between the E step and the M step}
#' \item{\code{fast_fit}}{ Logical, whether the closed form formulas should be used for the E step when available}
#' \item{\code{verbose}}{ Logical, whether details of the optimization should be printed}
#' \item{\code{lower_tol}}{ A "lower" bound for \eqn{\theta}; after this treshold, the algorithm returns the limiting log-likelihood of the no-frailty model. For example,
#' a value of 20 means that the maximum likelihood for \eqn{\theta} will be \eqn{\exp(20)}. For a frailty variance, this is approx \eqn{2 \times 10^{-9}}}
#' \item{\code{lik_tol}}{ For values higher than this, the algorithm returns a warning when the log-likelihood decreases between EM steps. Technically, this should not happen, but
#' if the parameter \eqn{\theta} is somewhere really far from the maximum, numerical problems might lead in very small likelihood decreases.
#' }}
#' The \code{fast_fit} option make a difference when the distribution is gamma (with or without left truncation) or
#' inverse Gaussian, i.e. pvf with m = -1/2 (without left truncation). For all the other scenarios, the fast_fit option will
#' automatically be changed to FALSE. When the number of events in a cluster / individual is not very small, the cases for which
#' fast fitting is available will show an improvement in performance.
#'
#' The starting value of the outer optimization may be set in the \code{distribution} argument.
#'
#' @seealso \code{\link{emfrail}}, \code{\link{emfrail_dist}}, \code{\link{emfrail_pll}}
#' @examples
#' emfrail_control()
#' emfrail_control(inner_control = list(eps = 1e-7))
#'

emfrail_control <- function(opt_fit = TRUE,
                            se = TRUE,
                            se_adj = TRUE,
                            ca_test = TRUE,
                            only_ca_test = FALSE,
                            lik_ci = TRUE,
                            lik_ci_intervals = list(interval = c(-3, 20),
                                                    interval_stable = c(0, 20)),
                            nlm_control = list(stepmax = 1),
                            inner_control = list(eps = 0.0001,
                                                 maxit = Inf,
                                                 fast_fit = TRUE,
                                                 verbose = FALSE,
                                                 lower_tol = 20,
                                                 lik_tol = 1)
) {
  # calculate SE as well

  # Here some checks

  if(isTRUE(only_ca_test) & !isTRUE(ca_test))
    stop("control: if only_ca_test is TRUE then ca_test must be TRUE as well")

  if(isTRUE(lik_ci)) {
    if(is.null(lik_ci_intervals))
      stop("opt_control must be a list which contains a named element interval")
    if(length(lik_ci_intervals$interval) != 2)
      stop("interval must be of length 2")
    if(lik_ci_intervals$interval[1] < -7 | lik_ci_intervals$interval[2] > 20)
      warning("extreme values for interval, there might be some numerical trouble")
  }

  # make sure the defaults of these function are the same as those from the input!
  inner_c <- function(eps = 0.0001,
                      maxit = Inf,
                      fast_fit = TRUE,
                      verbose = FALSE,
                      lower_tol = 20,
                      lik_tol = 1) {
    list(eps = eps,
         maxit = maxit,
         fast_fit = fast_fit,
         verbose = verbose,
         lower_tol = lower_tol,
         lik_tol = lik_tol)
  }

  inner_control <- do.call(inner_c, inner_control)

  res <- list(opt_fit = opt_fit,
              se = se,
              se_adj = se_adj,
              ca_test = ca_test,
              only_ca_test = only_ca_test,
              lik_ci = lik_ci,
              lik_ci_intervals = lik_ci_intervals,
              nlm_control = nlm_control,
              inner_control = inner_control)
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
#' @return An object of the type \code{emfrail_dist}, which is mostly used to denote the
#' supported frailty distributions in a consistent way.
#' @export
#'
#' @details The \code{theta} argument must be positive. In the case of gamma or PVF, this is the inverse of
#'  the frailty variance, i.e. the larger the \code{theta} is,
#'  the closer the model is to a Cox model. When \code{dist = "pvf"} and \code{pvfm = -0.5}, the inverse Gaussian
#'  distribution is obtained. For the positive stable distribution, the \eqn{\gamma} parameter of the Laplace transform is
#'  \eqn{\theta / (1 + \theta)}, with the \eqn{alpha} parameter fixed to 1.
#'
#' @seealso \code{\link{emfrail}, \link{emfrail_control}}
#' @examples
#' emfrail_dist()
#' # Compound Poisson distribution:
#' emfrail_dist(dist = 'pvf', theta = 1.5, pvfm = 0.5)
#' # Inverse Gaussian distribution:
#' emfrail_dist(dist = 'pvf')
emfrail_dist <- function(dist = "gamma", theta = 2, pvfm = -1/2, left_truncation = FALSE) {

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
  attr(res, "class") <- c("emfrail_dist")
  return(res)
}



