#' Summary for \code{emfrail} objects
#'
#' @param object An object of class \code{emfrail}
#' @param ... Ignored
#'
#' @return An object of class \code{emfrail_summary},
#' with some more human-readable results from an \code{emfrail} object.
#'
#' @details
#' Regardless of
#' the fitted model, one can expect the following
#' fields in this object: \code{est_dist} (an object of class \code{emfrail_distribution}) with the estimated
#' distribution, \code{loglik} (a named vector with the log-likelihoods of the no-frailty model, the frailty model,
#' the likelihood ratio test statistic and the p-value of the one-sided likelihood ratio test), \code{theta} (a named vector with
#' the estimated value of the parameter \eqn{\theta}, the standard error, and the limits of a 95% cofindence interval) and \code{frail}, which
#' is a data frame with the following columns: \code{id} (cluster identifier), \code{z} (empirical Bayes frailty estimates), and optional
#' \code{lower_q} and \code{upper_q} as the 2.5% and 97.5% quantiles of the posterior distribution of the frailties (only for gamma distribution).
#'
#' For the the PVF or gamma distributions, the field \code{fr_var} contains a transformation of \code{theta} to correspond to the
#' frailty variance.
#' The fields \code{pvf_pars} and \code{stable_pars} are for quantities that are calculated only when the distribution is PVF or stable.
#' If the model contains covariates, the field \code{coefmat} contains the corresponding details.
#'
#' @export
#'
#' @seealso \code{\link{predict.emfrail}, \link{plot_emfrail}}
#'
#' @examples
#' data("bladder")
#' coxph(Surv(start, stop, status) ~ treatment + frailty(id), data = bladder1, method = "breslow")
#' mod_gamma <- emfrail(bladder1, Surv(start, stop, status) ~ treatment + cluster(id))
#' summary(mod_gamma)
#'
#' # plot the Empirical Bayes estimates of the frailty
#' # easy way:
#' hist_frail(mod_gamma)
#'
#' # a fancy graph:
#' sum_mod <- summary(mod_gamma)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Create a plot just with the points
#' pl1 <- sum_mod$frail %>%
#'   arrange(z) %>%
#'   mutate(x = 1:n()) %>%
#'   ggplot(aes(x = x, y = z)) +
#'   geom_point()
#'
#' # If the quantiles of the posterior distribution are
#' # known, then error bars can be added:
#' if(!is.null(sum_mod$frail$lower_q))
#'   pl1 <- pl1 + geom_errorbar(aes(ymin = lower_q, ymax = upper_q), alpha = 0.5)
#'
#' pl1
#'
#' # The plot can be made interactive!
#' # ggplot2 gives a warning about the "id" aesthetic, just ignore it
#' pl2 <- sum_mod$frail %>%
#'   arrange(z) %>%
#'   mutate(x = 1:n()) %>%
#'   ggplot(aes(x = x, y = z)) +
#'   geom_point(aes(id = id))
#'
#' if(!is.null(sum_mod$z$lower_q))
#'   pl2 <- pl2 + geom_errorbar(aes(ymin = lower_q, ymax = upper_q, id = id), alpha = 0.5)
#'
#' library(plotly)
#' ggplotly(pl2)
#'
#' # Proportional hazards test
#' off_z <- log(sum_mod$frail$z)[match(bladder1$id, sum_mod$frail$id)]
#'
#' zph1 <- cox.zph(coxph(Surv(start, stop, status) ~ treatment + cluster(id), data = bladder1))
#'
#' # no sign of non-proportionality
#' zph2 <- cox.zph(coxph(Surv(start, stop, status) ~ treatment + offset(off_z), data = bladder1))
#'
#' zph2
#' # the p-values are even larger; the frailty "corrects" for proportionality.

summary.emfrail <- function(object, ...) {

  # Calculate the following: estimated distribution of the frailty at time 0
  fit <- object
  est_dist <- fit$.distribution
  est_dist$frailtypar <- exp(fit$outer_m$p1)

  loglik <- c(L0 = fit$loglik_null,
              L = -fit$outer_m$value[1],
              LRT = 2 * (-fit$outer_m$value[1] - fit$loglik_null),
              pval = pchisq(2 * (-fit$outer_m$value[1] - fit$loglik_null), df = 1, lower.tail = FALSE)/2)

  # this is the frailtypar actually.
  theta <- with(fit$outer_m, exp(p1))
  se_theta <- with(fit$outer_m,
                   msm::deltamethod(~exp(x1),
                                    mean = p1,
                                    cov = 1/attr(fit$outer_m, "details")[[3]]))


  # CI is symmetric on log(theta)
  ci_theta_low <- exp(with(fit, outer_m$p1 - 1.96 * sqrt(1/attr(outer_m, "details")[[3]])))
  ci_theta_high <- exp(with(fit, outer_m$p1 + 1.96 * sqrt(1/attr(outer_m, "details")[[3]])))

  if(theta > 9000) {
    ci_theta_low <- theta
    ci_theta_high <- Inf
  }
  # for gamma and pvf theta is 1/variance
  # for stable the L.T. is exp(- c^(1 - theta / (theta + 1)))

  fr_var <- se_fr_var <- ci_frvar_low <- ci_frvar_high <- NULL
  tau <- tau_stab <- se_tau <- attenuation <- ci_tau_high <- ci_tau_low <- e_log_y <- NULL
  se_tau_stab <- ci_tau_high_stab <- ci_tau_low_stab <- NULL
  coefmat <- NULL
  mass_at_0 <- NULL

  if(est_dist$dist != "stable") {
    fr_var <- 1/ theta
    se_fr_var <- with(fit$outer_m,
                      msm::deltamethod(~1/exp(x1),
                                       mean = p1,
                                       cov = 1/attr(fit$outer_m, "details")[[3]]))
    ci_frvar_low <- 1/ci_theta_high
    ci_frvar_high <- 1/ci_theta_low
  } else {

    # Kendall's tau
    tau_stab <- theta / (theta + 1)

    se_tau_stab <- with(fit$outer_m,
                   msm::deltamethod(~exp(x1) / (exp(x1) + 1),
                                    mean = p1,
                                    cov = 1/attr(fit$outer_m, "details")[[3]]))

    attenuation <- 1 - theta  # this is the gamma of the distribution as well

    ci_tau_low_stab <- ci_theta_high / (1 + ci_theta_high)
    ci_tau_high_stab <- ci_theta_low / (1 + ci_theta_low)

    e_log_y <-  (-1) * (1 / (1-theta) - 1) * digamma(1)

  }



  # The coefficients
  coef <- fit$inner_m$coef
  se_coef <- with(fit$inner_m, sqrt(diag(Vcov)[seq_along(coef)]))
  se_coef_adj <- with(fit, sqrt(diag(vcov_adj)[seq_along(inner_m$coef)]))


  # The frailty estimates
  z <- with(fit$inner_m, estep[,1] / estep[,2])
  names(z) <- rownames(fit$inner_m$nev_id)
  # dist_pars <- dist_to_pars(est_dist$dist, log(est_dist$frailtypar), est_dist$pvfm)


  # For the gamma we also know the distribution
  lower_q <- upper_q <- NULL

  if (est_dist$dist == "gamma") {
    # the EB frailties have mean z and variance...

    # Kendall's tau
    tau <- 1 / (1 + 2 * theta)

    se_tau <- with(fit$outer_m,
                   msm::deltamethod(~1 / (1 + 2 * exp(x1)),
                                    mean = p1,
                                    cov = 1/attr(fit$outer_m, "details")[[3]]))
    ci_tau_low <- 1 / (1 + 2 * ci_theta_low)
    ci_tau_high <- 1 / (1 + 2 * ci_theta_high)


    shape <- est_dist$frailtypar + fit$inner_m$nev_id
    rate <- est_dist$frailtypar + fit$inner_m$Cvec
    var_z <- shape / rate^2
    lower_q <- stats::qgamma(0.025,
                      shape = shape,
                      rate = rate)
    upper_q <- stats::qgamma(0.975,
                      shape = shape,
                      rate = rate)
  }
  z <- data.frame(id = rownames(fit$inner_m$nev_id),
                  z = z)

  if(!is.null(lower_q)) {
    z$lower_q <- as.numeric(lower_q)
    z$upper_q <- as.numeric(upper_q)
  }

  if(length(fit$inner_m$coef) > 0) {
    coefmat <- list(
      coef = fit$inner_m$coef,
      "exp(coef)" = exp(fit$inner_m$coef),
      "se(coef)" = sqrt(diag(fit$inner_m$Vcov)[seq_along(fit$inner_m$coef)]),
      "adjusted se" = sqrt(diag(fit$vcov_adj)[seq_along(fit$inner_m$coef)] ))

    coefmat$z <- coefmat$coef / coefmat$`se(coef)`
    coefmat$p <-  1 - pchisq(coefmat$z^2, df = 1)
  }



  dist_pars <- with(est_dist, dist_to_pars(dist, log(frailtypar), pvfm))

  if(est_dist$dist == "pvf" & est_dist$pvfm > 0) {
    mass_at_0 = with(dist_pars, exp(-alpha))
  }


  ret <- list(est_dist = est_dist,
       loglik = loglik,
       ca_test = object$ca_test,
       theta = c(theta = theta,
                 se_theta = se_theta,
                 ci_theta_low = ci_theta_low,
                 ci_theta_high = ci_theta_high),
       fr_var = c(fr_var = fr_var,
                  se_fr_var = se_fr_var,
                  ci_frvar_low = ci_frvar_low,
                  ci_frvar_high = ci_frvar_high),
       gamma_pars = c(tau = tau,
                      se_tau = se_tau,
                      ci_tau_high = ci_tau_high,
                      ci_tau_low = ci_tau_low),
       pvf_pars = c(mass_at_0 = mass_at_0),
       stable_pars = c(tau = tau_stab,
                       se_tau = se_tau_stab,
                       ci_tau_high = ci_tau_high_stab,
                       ci_tau_low = ci_tau_low_stab,
                       attenuation = attenuation,
                       e_log_y = e_log_y),
       coefmat = do.call(cbind,coefmat),
       frail = z
       )

  #attr(ret, "class") <- "emfrail_summary"
  class(ret) <- "emfrail_summary"
  ret
}



# ppl1 <- data.frame(id = fit$res$z$id,
#            z = fit$res$z$z,
#            shape = fit$est_dist$frailtypar + fit$res$z$nev,
#            rate = fit$est_dist$frailtypar + fit$res$z$Lambda) %>%
#   arrange(z) %>%
#   ggplot(aes(x = seq_along(id), y = z)) + geom_point() +
#   geom_errorbar(aes(ymin = qgamma(0.025, shape = shape, rate = rate), ymax = qgamma(0.975, shape = shape, rate = rate),
#                     id = id, gamma_shape = shape, gamma_rate = rate))
#
#
# ggplotly(ppl1)
