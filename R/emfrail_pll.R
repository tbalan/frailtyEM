#' Profile log-likelihood calculation
#'
#' @param .data Same as in \code{emfrail}
#' @param .formula Same as in \code{emfrail}
#' @param .distribution Same as in \code{emfrail}
#' @param .values A vector of values on where to calculate the profile likelihood. See details.
#'
#' @return The profile log-likelihood at the specific value of the frailty parameter
#' @export
#'
#' @details It is of interest sometimes to see the profile log-likelihood. The scale is that of \code{frailtypar} as defined in \code{emfrail_distribution()}.
#' For the gamma and pvf frailty, that is the inverse of the frailty variance.
#' @examples
#'
#' fr_var <- c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
#' profloglik <- emfrail_pll(rats,
#'                           Surv(rep(0, nrow(rats)), time, status) ~ rx + sex + cluster(litter),
#'                           .values = 1/fr_var)
#' plot(fr_var, profloglik, xlab = "frailty variance", ylab = "profile log-likelihood")
#'
#' # check with coxph:
#' profloglik_cph<- sapply(fr_var, function(th)
#'   coxph(data =  rats, formula = Surv(time, status) ~ rx + sex + frailty(litter, theta = th),
#'                  method = "breslow")$history[[1]][[3]])
#'
#' lines(fr_var, profloglik_cph, col = 2)
#'
emfrail_pll <- function(.data, .formula,
                     .distribution = emfrail_distribution(),
                     .values) {
  sapply(.values, function(fp) {
    -emfrail(.data = .data,
            .formula = .formula,
            .distribution =  emfrail_distribution(dist = .distribution$dist,
                                                  frailtypar = fp,
                                                  pvfm = .distribution$pvfm,
                                                  left_truncation = .distribution$left_truncation),
            .control = emfrail_control(opt_fit = FALSE))
  })

}
