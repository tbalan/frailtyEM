#' Predicted hazard and survival curves from an \code{emfrail} object
#'
#' @param fit An \code{emfrail} fit object
#' @param lp A vector of linear predictor values at which to calculate the curves.
#' @param quantity The quantity to be calculated for the values of \code{lp}
#' @param type The type of the quantity (conditional  / marginal)
#' @param conf_int The type of the confidence interval (adjusted / regular)
#'
#' @return A data frame with the column \code{time} and several other columns according to the input.
#' By default, for each \code{lp} it will give the following columns: \code{cumhaz}, \code{survival},
#' \code{cumhaz_m}, \code{survival_m} for the cumulative hazard and survival, conditional and marginal.
#'
#'
#' @export
#'
#' @include laplace_transform.R
#'
#' @examples
#' #' m1 <- emfrail(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#' .distribution = emfrail_distribution(dist = "gamma"))
#' predict.emfrail(m1)
#'
predict.emfrail <- function(fit,
                            lp = c(0),
                            quantity = c("cumhaz", "survival"),
                            type = c("conditional", "marginal"),
                            conf_int = c("regular", "adjusted")) {

  ncoef <- length(fit$res$coef)

  varH <- fit$vcov[(ncoef + 1): nrow(fit$vcov), (ncoef+ 1): nrow(fit$vcov)]
  varH_adj <- fit$vcov_adj[(ncoef + 1): nrow(fit$vcov), (ncoef+ 1): nrow(fit$vcov)]

  varH_time <- numeric(nrow(varH))
  for(i in 1:nrow(varH)) {
    varH_time[i] = sum(varH[1:i, 1:i])
  }

  varH_adj_time <- numeric(nrow(varH))
  for(i in 1:nrow(varH)) {
    varH_adj_time[i] = sum(varH_adj[1:i, 1:i])
  }

  # The "core" is to calculate the cumulative hazard and the confidence band for it
  time <- fit$res$haz$time
  cumhaz <- fit$res$haz$cumhaz

  se_chz <- sqrt(varH_time)
  se_chz_adj <- sqrt(varH_adj_time)
  # lower_chz <- pmax(0, cumhaz - 1.96*se_chz)
  # upper_chz <- cumhaz + 1.96*se_chz
  #
  # lower_chz_adj <- pmax(0, cumhaz - 1.96*se_chz_adj)
  # upper_chz_adj <- cumhaz + 1.96*se_chz_adj

  # Now calculate that for different LPs

  ret <- do.call(rbind, lapply(lp, function(x) data.frame(time = time,
                                              cumhaz = cumhaz * exp(x),
                                              se_chz = exp(x) * se_chz,
                                              se_chz_adj = exp(x) * se_chz_adj,
                                              lp = as.factor(x))))


  chz_to_surv <- function(x) exp(-x)
  surv_to_chz <- function(x) -1 * log(x)

  scenarios <- expand.grid(quantity, type, conf_int)

  survival <- function(chz) exp(-chz)
  survival_m <- function(chz) laplace_transform(chz, fit$est_dist)
  cumhaz_m <- function(chz) -1 * log(laplace_transform(chz, fit$est_dist))

  # cumhaz_m(cumhaz)
  # first, add confidence bands for the cumulative hazard
  bounds <- "cumhaz"

  if(length(conf_int) > 0) {
    if("regular" %in% conf_int) {
      bounds <- c(bounds, "cumhaz_l", "cumhaz_r")
      ret$cumhaz_l <- pmax(0, cumhaz - 1.96*se_chz)
      ret$cumhaz_r <- cumhaz + 1.96*se_chz
    }

    if("adjusted" %in% conf_int) {
      bounds <- c(bounds, "cumhaz_l_a", "cumhaz_r_a")
        ret$cumhaz_l_a <- pmax(0, cumhaz - 1.96*se_chz_adj)
        ret$cumhaz_r_a <- cumhaz + 1.96*se_chz_adj
    }
  }


  #survival_m(1.5)
  #fit$est_dist$dist == "gamma"
  # bounds is the names of columns that we want to transform
  if("survival" %in% quantity & "conditional" %in% type) {
    ret <- cbind(ret, as.data.frame(lapply(ret[bounds], survival), col.names = sub("cumhaz", "survival", bounds)))
  }

  if("survival" %in% quantity & "marginal" %in% type) {
    ret <- cbind(ret, as.data.frame(lapply(ret[bounds], survival_m), col.names = sub("cumhaz", "survival_m", bounds)))
  }

  if("cumhaz" %in% quantity & "marginal" %in% type) {
    ret <- cbind(ret, as.data.frame(lapply(ret[bounds], cumhaz_m), col.names = sub("cumhaz", "cumhaz_m", bounds)))
  }

  if(!("cumhaz" %in% quantity & "conditional" %in% type)) {
    cols_to_keep <- which(!(names(ret) %in% c(bounds, "se_chz", "se_chz_adj")))
    ret <- ret[,cols_to_keep]
  }

  ret
}
