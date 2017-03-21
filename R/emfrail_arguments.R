#' Control parameters for emfrail
#'
#' @param eps Convergence criterion for the inner loops (the EM algorithm) for fixed frailty parameters
#' @param maxit Maximum number of iterations in the inner loops (the EM algorithm)
#' @param opt_fit Logical. Whether the outer optimization should be carried out.
#' If \code{FALSE}, then the frailty parameter is treated as fixed and the \code{emfrail} function returns only log-likelihood. See details.
#' @param verbose Logical. Whether to print out information about what is going on during maximization.
#' @param fast_fit Logical. Whether to try to calculate the E step directly, when possible. See details.
#' @param zerotol Lower limit for 1/frailtypar (variance in case of gamma / pvf). Below this value, the variance is taken to be 0 and the
#' EM is not actually performed. The log-likelihood returned to the maximizer is that for the Cox model.
#' @param se_fit Logical. Whether to calculate the variance / covariance matrix.
#' @param se_adj Logical. Whether to calculate the adjusted variance / covariance matrix (needs \code{se_fit == TRUE})
#' @param ca_test Logical. Should the Commenges-Andersen test be calculated?
#' @param only_ca_test Logical. Should ONLY the Commenges-Andersen test be calculated?
#' @param opt_control A list with arguments to be sent to the maximizer.
#'
#' @return An object of the type \code{emfrail_control}.
#' @export
#'
#' @details The \code{fast_fit} option make a difference when the distribution is gamma (with or without left truncation) or
#' inverse Gaussian, i.e. pvf with m = -1/2 (without left truncation). For all the other scenarios, the fast_fit option will
#' automatically be changed to FALSE. When the number of events in a cluster / individual is not very small, the cases for which
#' fast fitting is available will show an improvement in performance.
#'
#' The \code{zerotol} option defaults to \code{1e-04}, which in practical terms means, for example, that for
#' the gamma / pvf distribution, a frailty variance below \code{1e-04} can not be detected.
#'
#' @seealso \code{\link{emfrail}}, \code{\link{emfrail_distribution}}, \code{\link{emfrail_pll}}
#' @examples
#' emfrail_control()
#' emfrail_control(eps = 10e-7)
#'
#' \dontrun{
#' # A data set with very small heterogeneity
#' set.seed(10)
#' x <- sample(c(0,1/2), 500, TRUE)
#' tstart <- rep(0, 500)
#' tstop <- rexp(1 * exp(x))
#' status <- rep(1, 500)
#' id <- rep(1:100, each = 5)
#' dat <- data.frame(id, tstart, tstop, status, x)
#'
#' # What coxph does:
#' library(survival)
#' m_cph <- coxph(Surv(tstart, tstop, status) ~ x + frailty(id), dat, ties = "breslow")
#' m_cph$history
#' # For the frailty variance, the program tries: 0, 1, 0.5, 0.005, 0.00005, etc. Stops at 5e-7.
#'
#' m_ft <- emfrail(dat, Surv(tstart, tstop, status) ~ x + cluster(id))
#' m_ft
#'
#' # The algorithm gives as frailty parameter 10587.88,
#' # which means frailty variance 1/10587.88 = 9.44e-05
#' # That is because by default, zerotol = 1e-04,
#' # which is the point where the algorithm decides that the frailty is 0.
#'
#' # If you want the exact value of the estimate,
#' # increase the precision so the point stays somewhat interior to the parameter space:
#' m_ft_08 <- emfrail(dat, Surv(tstart, tstop, status) ~ x + cluster(id),
#'                    .control = emfrail_control(zerotol = 1e-8))
#' # This gives a more precise estimate, 5845410 so frailty variance 1/5845410 = 1.71e-07
#' }
emfrail_control <- function(eps = 0.0001, maxit = Inf, opt_fit = TRUE, verbose = FALSE, fast_fit = TRUE,
                            zerotol = 1e-4, se_fit = TRUE, se_adj = TRUE, ca_test = TRUE, only_ca_test = FALSE,
                            opt_control = list(method = "bobyqa",
    itnmax = NULL, control = list())) {
    # calculate SE as well

    # Here some checks

  if(isTRUE(only_ca_test) & !isTRUE(ca_test)) stop("control: if only_ca_test is TRUE then ca_test must be TRUE as well")
    res <- list(eps = eps,
                maxit = maxit,
                opt_fit = opt_fit,
                zerotol = zerotol,
                verbose = verbose,
                fast_fit = fast_fit,
                se_fit = se_fit,
                se_adj = se_adj,
                ca_test = ca_test,
                only_ca_test = only_ca_test,
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
#'  \eqn{1 - \theta / (1 + \theta)}, with the \eqn{alpha} parameter fixed to 1.
#'
#' @seealso \code{\link{emfrail}, \link{emfrail_control}}
#' @examples
#' emfrail_distribution()
#' # Compound Poisson distribution:
#' emfrail_distribution(dist = 'pvf', theta = 1.5, pvfm = 0.5)
#' # Inverse Gaussian distribution:
#' emfrail_distribution(dist = 'pvf')
emfrail_distribution <- function(dist = "gamma", theta, pvfm = -1/2, left_truncation = FALSE) {

  if (!(dist %in% c("gamma", "stable", "pvf")))
    stop("frailty distribution must be one of gamma, stable, pvf")
  if(missing(theta)) {
    if(dist == "stable") theta <- 0.5 else theta <- 2
  }
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



