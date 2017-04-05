#' Predicted hazard and survival curves from an \code{emfrail} object
#'
#' @importFrom stats .getXlevels delete.response drop.terms as.formula
#' @param object An \code{emfrail} fit object
#' @param lp A vector of linear predictor values at which to calculate the curves. Default is 0 (baseline).
#' @param newdata A data frame with the same variable names as those that appear in the \code{emfrail} formula, used to calculate the \code{lp} (optional).
#' @param quantity Can be \code{"cumhaz"} and/or \code{"survival"}. The quantity to be calculated for the values of \code{lp}.
#' @param type Can be \code{"conditional"} and/or \code{"marginal"}. The type of the quantity to be calculated.
#' @param conf_int Can be \code{"regular"} and/or \code{"adjusted"}. The type of confidence interval to be calculated.
#' @param ... Ignored
#'
#' @return A data frame with the column \code{time} and several other columns according to the input.
#' By default, for each \code{lp} it will give the following columns: \code{cumhaz}, \code{survival},
#' \code{cumhaz_m}, \code{survival_m} for the cumulative hazard and survival, conditional and marginal.
#'
#' @details There are two ways of specifying for which individuals to calculate the cumulative hazard or survival curve.
#' One way is directly, for certain values of the linear predictor, via \code{lp},
#' and the second way is via \code{newdata}, a \code{data.frame} where the column names are the same as in the original data.
#' If \code{newdata} is specified, then the \code{lp} argument is ignored.
#'
#' The names of the columns in the returned object are as follows: \code{time} represents the unique event time points
#' from the data set, \code{lp} is the value of the linear predictor (as specified in the input or as calculated from the lines of \code{newdata}). If
#' \code{newdata} is specified, columns repeating each line of \code{newdata} are also added.
#'  The two "quantities" that can be returned are
#' named \code{cumhaz} and \code{survival}. If we denote each quantity with \code{q}, then the columns with the marginal estimates
#' are named \code{q_m}. The confidence intervals contain the name of the quantity (conditional or marginal) followed by \code{_l} or \code{_r} for
#' the lower and upper bound. The bounds calculated with the adjusted standard errors have the name of the regular bounds followed by
#' \code{_a}. For example, the adjusted lower bound for the marginal survival is in the column named \code{survival_m_l_a}.
#'
#' The \code{emfrail} only gives the Breslow estimates of the  baseline hazard \eqn{\lambda_0(t)} at the
#' event time points, conditional on the frailty. Let \eqn{\lambda(t)} be the baseline hazard for a linear predictor of interest.
#' The estimated conditional cumulative hazard is then
#' \eqn{\Lambda(t) = \sum_{s= 0}^t \lambda(s)}. The variance of \eqn{\Lambda(t)} can be calculated from the (maybe adjusted)
#' variance-covariance matrix.
#'
#' The conditional survival is obtained by the usual expression \eqn{S(t) = \exp(-\Lambda(t))}. The marginal survival
#' is given by
#' \deqn{\bar S(t) = E \left[\exp(-\Lambda(t)) \right] = \mathcal{L}(\Lambda(t)),}
#' i.e. the Laplace transform of the frailty distribution calculated in \eqn{\Lambda(t)}.
#'
#' The marginal hazard is obtained as \deqn{\bar \Lambda(t) = - \log \bar S(t).}
#'
#' The only standard errors that are available from \code{emfrail} are those for \eqn{\lambda_0(t)}. From this,
#' standard errors of \eqn{\log \Lambda(t)} may be calculated. On this scale, the symmetric confidence intervals are built, and then
#' moved to the desired scale.
#' The following two issues arise: (1) the linear predictor is taken as fixed,
#' i.e. the variability in the estimation of the regression coefficient is not taken into account and (2) this construction of the confidence intervals
#' is a bit dodgy, in the sense that this is not what happens in the survival package. Several options will be added in future versions.
#'
#'
#'
#' @export
#'
#'
#' @examples
#' kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
#' m1 <- emfrail(.data =  kidney,
#'               .formula = Surv(time, status) ~  sex + age  + cluster(id))
#'
#' pred <- predict(m1)
#'
#' names(pred)
#'
#' # Plot baseline cumulative hazard: note that is for someone aged 0!
#' plot_pred(m1)
#'
#' # More realistic:
#' plot_pred(m1, newdata = data.frame(sex = "female", age = mean(kidney$age)))
#'
#' # Plot survival
#' plot_pred(m1,
#'           newdata = data.frame(sex = "female", age = mean(kidney$age)),
#'           quantity = "survival", conf_int = "none")
#'\dontrun{
#' # Plot cumulative hazard with confidence intervals, ggplot2
#' library(ggplot2)
#' ggplot(pred, aes(x = time, y = cumhaz)) +
#'   geom_step() +
#'   geom_ribbon(aes(ymin = cumhaz_l, ymax = cumhaz_r), alpha = 0.2) +
#'   geom_ribbon(aes(ymin = cumhaz_l_a, ymax = cumhaz_r_a), alpha = 0.2) +
#'   ggtitle("Baseline cumulative hazard with confidence intervals")
#'
#' # For two individuals: with sex 1 and sex 0
#' pred2 <- predict(m1, newdata = data.frame(sex = c("female", "male"), age = c(44, 44)))
#' # Plot the conditional & survival of two individuals
#' ggplot(pred2, aes(x = time, y = survival, group = sex)) +
#'   geom_step(aes(colour = sex)) + ggtitle("Conditional survival")
#'
#' ggplot(pred2, aes(x = time, y = survival_m, group = sex)) +
#'   geom_step(aes(colour = sex)) + ggtitle("Marginal survival")
#'
#' # Plot the conditional and the marginal survival in the same place
#' library(dplyr)
#' library(tidyr)
#' pred2 %>%
#'   gather(key = variable, value = survival, survival, survival_m) %>%
#'   mutate(variable = ifelse(variable == "survival", "conditional", "marginal")) %>%
#'   ggplot(aes(x = time, y = survival, colour = sex, linetype = variable)) +
#'   geom_step() + ggtitle("Survival by sex")
#'
#' # The hazard ratio
#' hr_conditional <- pred2$cumhaz[pred2$sex == "female"] / pred2$cumhaz[pred2$sex == "male"]
#' hr_marginal <- pred2$cumhaz_m[pred2$sex == "female"] / pred2$cumhaz_m[pred2$sex == "male"]
#' time <- pred2$time[pred2$sex == "male"]
#'
#' plot(time, hr_marginal, type = "s", col = 2, main = "Hazard ratio female vs male")
#' lines(time, hr_conditional, type = "s")
#' legend(c("conditional", "marginal"), x = "bottomleft", col = c(1,2), lty = 1)
#' # The marginal hazard ratio in the case of gamma frailty shrinks towards 1
#' # With positive stable, this plot would be two parallel lines
#'
#' # Or easier, in this way:
#' plot_hr(m1, newdata = data.frame(sex = c("female", "male"), age = c(44, 44)))
#' }
#' @seealso \code{\link{plot_pred}}, \code{\link{plot_hr}}
predict.emfrail <- function(object,
                            lp = c(0),
                            newdata = NULL,
                            quantity = c("cumhaz", "survival"),
                            type = c("conditional", "marginal"),
                            conf_int = c("regular", "adjusted"),
                            ...) {


  if(!is.null(newdata)) {
    if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame")
    #message("newdata specified, ignoring lp")

    mdata <- attr(object, "metadata")
    mf <- model.frame(mdata[[1]], data = newdata, xlev = mdata[[2]])
    mm <- try(model.matrix(mdata[[1]], mf)[,-1])
    if(inherits(mm, "try-error")) stop("newdata probably misspecified")
    lp <- as.numeric(mm %*% object$inner_m$coef)
    lp_all <- cbind(newdata, lp)
  } else
    lp_all <- data.frame(lp = lp)



  fit <- object
  est_dist <- fit$.distribution
  est_dist$frailtypar <- exp(fit$outer_m$minimum)


  ncoef <- length(fit$inner_m$coef)
  varH <-   fit$inner_m$Vcov[(ncoef + 1): nrow(fit$inner_m$Vcov), (ncoef+ 1): nrow(fit$inner_m$Vcov)]
  varH_adj <- fit$vcov_adj[(ncoef + 1): nrow(fit$vcov_adj), (ncoef+ 1): nrow(fit$vcov_adj)]


  loghaz <- log(fit$inner_m$haz)



  time <- fit$inner_m$tev
  cumhaz <- cumsum(fit$inner_m$haz)

  xs <- lapply(seq_along(fit$inner_m$haz), function(x) text1 <- paste0("x", x))
  for(i in 2:length(xs)) {
    xs[[i]] = paste0(xs[[i-1]], " + ", xs[[i]])
  }
  forms <- lapply(xs, function(x) as.formula(paste0("~log(", x, ")")))

  # These are the SE of log cumulative hazard
  se_logH <- msm::deltamethod(g = forms,
                   mean = fit$inner_m$haz,
                   cov = varH,
                   ses = TRUE)

  se_logH_adj <- msm::deltamethod(g = forms,
                                   mean = fit$inner_m$haz,
                                   cov = varH_adj,
                                   ses = TRUE)

  # these are the variances at every time point
  # varH_time <- numeric(nrow(varH))
  # for(i in 1:nrow(varH)) {
  #   varH_time[i] = sum(varH[1:i, 1:i])
  # }
  #
  # varH_adj_time <- numeric(nrow(varH))
  # for(i in 1:nrow(varH)) {
  #   varH_adj_time[i] = sum(varH_adj[1:i, 1:i])
  # }

  # dodgy stuff: get the

  # The "core" is to calculate the cumulative hazard and the confidence band for it



  # se_chz <- sqrt(varH_time)
  # se_chz_adj <- sqrt(varH_adj_time)

  # lower_chz <- pmax(0, cumhaz - 1.96*se_chz)
  # upper_chz <- cumhaz + 1.96*se_chz
  #
  # lower_chz_adj <- pmax(0, cumhaz - 1.96*se_chz_adj)
  # upper_chz_adj <- cumhaz + 1.96*se_chz_adj

  # Now calculate that for different LPs
  #lp_all <- data.frame(lp_all, row.names = NULL)

  mintime <- max(0, c(min(time)-1))

  #row.names(lp_all) <- NULL

  ret <- do.call(rbind,
                 lapply(split(lp_all, 1:nrow(lp_all)), function(x) cbind(time = c(mintime,time),
                                                                         cumhaz = c(0, cumhaz*exp(x$lp)),
                                                                         se_logchz = c(0, se_logH),
                                                                         se_logchz_adj = c(0, se_logH_adj),
                                                                         x,
                                                                         row.names = NULL))
  )
  ret



  # ret <- do.call(rbind,
  #                lapply(split(lp_all, 1:nrow(lp_all)), function(x) cbind(time = c(mintime,time),
  #                                             cumhaz = c(0, cumhaz * exp(x$lp)),
  #                                             se_chz = c(0, exp(x$lp) * se_chz),
  #                                             se_chz_adj = c(0, exp(x$lp) * se_chz_adj),
  #                                             x,
  #                                             row.names = NULL))
  #                )


  chz_to_surv <- function(x) exp(-x)
  surv_to_chz <- function(x) -1 * log(x)

  scenarios <- expand.grid(quantity, type, conf_int)

  survival <- function(chz) exp(-chz)
  survival_m <- function(chz) laplace_transform(chz, est_dist)
  cumhaz_m <- function(chz) -1 * log(laplace_transform(chz, est_dist))

  # cumhaz_m(cumhaz)
  # first, add confidence bands for the cumulative hazard
  bounds <- "cumhaz"

  if(length(conf_int) > 0) {
    if("regular" %in% conf_int) {
      bounds <- c(bounds, "cumhaz_l", "cumhaz_r")
      ret$cumhaz_l <- pmax(0, exp(log(ret$cumhaz) - 1.96*ret$se_logchz))
      ret$cumhaz_r <- exp(log(ret$cumhaz) + 1.96*ret$se_logchz)
    }

    if("adjusted" %in% conf_int) {
      bounds <- c(bounds, "cumhaz_l_a", "cumhaz_r_a")
        ret$cumhaz_l_a <- pmax(0, exp(log(ret$cumhaz) - 1.96*ret$se_logchz_adj))
        ret$cumhaz_r_a <- exp(log(ret$cumhaz) + 1.96*ret$se_logchz_adj)
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
