#' Profile log-likelihood
#'
#' @param .data Same as in \code{emfrail}
#' @param .formula Same as in \code{emfrail}
#' @param .distribution Same as in \code{emfrail}
#' @param .values A vector of values on where to calculate the profile likelihood. These values are the values of "theta"
#'
#' @return The maximum log-likelihood at the specific value of theta
#' @export
#'
#' @examples
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
