#' Fitting shared frailty models with the EM algorithm
#'
#' @importFrom survival Surv
#' @importFrom stats approx coef model.frame model.matrix pchisq printCoefmat nlm uniroot cor
#' @importFrom magrittr "%>%"
#' @importFrom Rcpp evalCpp
#' @useDynLib frailtyEM, .registration=TRUE
#' @include em_fit.R
#' @include emfrail_aux.R
#'
#' @param .data A data frame in which the formula argument can be evaluated
#' @param .formula A formula that contains on the left hand side an object of the type \code{Surv}
#' and on the right hand side a \code{+cluster(id)} statement. Optionally, also a \code{+terminal()} statement
#' may be added, and then a score test for association between the event process and the result in the specified
#' column is performed. See details.
#' @param .distribution An object as created by \code{\link{emfrail_distribution}}
#' @param .control An object as created by \code{\link{emfrail_control}}
#' @return An object of the class \code{emfrail}, that is in fact a list which contains (1) the object returned by the
#' "outer maximization", (2) the object with all the estimates returned by the "inner maximization",
#' (3) the log-likelihood of the Cox model without frailty, (4) the variance-covariance matrix adjusted for the uncertainty in the
#' outer maximization, (5) the results of the Commenges-Andersen test for heterogeneity and (6,7,8) are copies of the original input arguments: .formula, .distribution and .control.
#' Two attributes are also present, \code{class} for determining the object type and \code{metadata} which contains some information that is used to
#' process the input for \code{predict.emfrail()}.
#' @export
#'
#' @details The \code{emfrail} function fits shared frailty models for processes which have intensity
#' \deqn{\lambda(t) = z \lambda_0(t) \exp(\beta' \mathbf{x})}
#' with a non-parametric (Breslow) baseline intensity \eqn{\lambda_0(t)}.
#' The distribution of \eqn{z} is usually described by one parameter \eqn{\theta}.
#' The family of supported distributions can be one of gamma, positive stable or PVF (power-variance-family).
#' For details, see the vignette and \code{\link{emfrail_distribution}} .
#'
#' The algorithm is described in detail in the vignette. The short version is as follows:
#' The objective is to maximize a marginal likelihood with respect to 3 sets of parameters, \eqn{\theta > 0}
#' (the parameter of the frailty distribution), \eqn{\beta} (a vector of regression coefficients) and \eqn{\lambda_0}
#' (a vector of "jumps" of the cumulative hazards at each event time point).
#' The essence of the algorithm relies on the observation that, if \eqn{\theta} would be fixed, then the maximization
#' with respect to \eqn{\beta, \theta} can be done with an EM algorithm. This procedure has been described, for the
#' gamma frailty, in Nielsen (1992).
#' The "inner" problem is, for a fixed of \eqn{\theta}, to calculate
#'\deqn{\widehat{L}(\theta) = \mathrm{max}_{\beta, \lambda_0} L(\beta, \lambda_0 | \theta).}
#' which is done via an EM algorithm. The "outer" problem is to calculate
#' \deqn{\widehat{L} = \mathrm{max}_{\theta} \widehat{L}(\theta).}
#'
#' The "inner" problem is solved relying on \code{agreg.fit} from the \code{survival} package in the M step
#' and with an E step which is done in \code{R} when there are closed form solutions and in \code{C++} when the
#' solutions are calculated with a recursive algorithm. The information matrix (given fixed \eqn{\theta}) is calculated using Louis' formula.
#' The "outer" problem is solved relying on the maximizers in the \code{nlm} function, which also provides an estimate
#' of the Hessian matrix. The remaining elements of the information matrix are approximated numerically.
#'
#' Several options are available to the user. In the \code{.distribution} argument, the frailty distribution and the
#' starting value for \eqn{\theta} can be specified, as well as whether the data consists of left truncated
#' survival times or not (the default). In the \code{.control} argument, several parameters control the maximization.
#' The convergence criterion for the EM can be specified (or a maximum number of iterations). The "outer" procedure can be
#' avoided altogether and then \code{emfrail} calculates only \eqn{\widehat{L}(\theta)} at the given starting value.
#'
#' The Commenges-Andersen test is a score test for heterogeneity which test the null hypothesis that all
#' frailties are equal to 1. This only uses the information from the initial proportional hazards model
#' without frailty.
#'
#' The score test for dependent censoring test is detailed in Balan et al (2016). If significant, it may indicate
#' that dependent censoring is present. A common scenario is when recurrent events are stopped by a terminal event.
#'
#' @note Some possible problems may appear when the maximum likelihood estimate lies on the border of the parameter space.
#' Usually, this will happen when the "outer" parameter MLE is infinity (i.e. variance 0 in case of gamma and PVF).
#' For small enough values of \eqn{1/\theta} the log-likelihood
#' of the Cox model is returned to avoid such problems. This option can be tweaked in \code{emfrail_control()}.
#'
#' @seealso \code{\link{plot_emfrail}} and \code{\link{ggplot_emfrail}} for plot functions directly available, \code{\link{emfrail_pll}} for calculating \eqn{\widehat{L}(\theta)} at specific values of \eqn{\theta},
#' \code{\link{summary.emfrail}} for transforming the \code{emfrail} object into a more human-readable format and for
#' visualizing the frailty (empirical Bayes) estimates,
#' \code{\link{predict.emfrail}} for calculating and visalizing conditional and marginal survival and cumulative
#' hazard curves.
#'
#' @examples
#' dat <- survival::rats
#'
#'
#' m1 <- emfrail(.data =  dat,
#'               .formula = Surv(time, status) ~  rx + sex + cluster(litter))
#' m1
#' summary(m1)
#' \dontrun{
#' # for the Inverse Gaussian distribution
#' m2 <- emfrail(.data =  dat,
#'               .formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'               .distribution = emfrail_distribution(dist = "pvf"))
#' m2
#'
#' # for the PVF distribution with m = 0.75
#' m3 <- emfrail(.data =  dat,
#'               .formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'               .distribution = emfrail_distribution(dist = "pvf", pvfm = 0.75))
#' m3
#'
#' # for the positive stable distribution
#' m4 <- emfrail(.data =  dat,
#'               .formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'               .distribution = emfrail_distribution(dist = "stable"))
#' m4
#'}
#' # Compare marginal log-likelihoods
#' \dontrun{
#' models <- list(m1, m2, m3, m4)
#'
#' logliks <- lapply(models,
#'                   function(x) -x$outer_m$value)
#'
#' names(logliks) <- lapply(models,
#'                          function(x) with(x$.distribution,
#'                                           ifelse(dist == "pvf",
#'                                                  paste(dist, "/", pvfm),
#'                                                  dist))
#' )
#'
#' logliks
#' }
#' # Draw the profile log-likelihood
#' \dontrun{
#' fr_var <- seq(from = 0.01, to = 1.4, length.out = 20)
#'
#' # For gamma the variance is 1/theta (see parametrizations)
#' pll_gamma <- emfrail_pll(.data =  dat,
#'                          .formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'                          .values = 1/fr_var )
#'  plot(fr_var, pll_gamma,
#'      type = "l",
#'      xlab = "Frailty variance",
#'      ylab = "Profile log-likelihood")
#' }
#' \dontrun{
#' # The same can be done with coxph, where variance is refered to as "theta"
#' pll_cph <- sapply(fr_var, function(fr)
#'   coxph(data =  dat, formula = Surv(time, status) ~ rx + sex + frailty(litter, theta = fr),
#'         method = "breslow")$history[[1]][[3]])
#'
#' # Same for inverse gaussian
#' pll_if <- emfrail_pll(.data =  dat,
#'                       .formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'                       .distribution = emfrail_distribution(dist = "pvf"),
#'                       .values = 1/fr_var )
#'
#' # Same for pvf with a psoitive pvfm parameter
#' pll_pvf <- emfrail_pll(.data =  dat,
#'                        .formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'                        .distribution = emfrail_distribution(dist = "pvf", pvfm = 1.5),
#'                        .values = 1/fr_var )
#'
#' miny <- min(c(pll_gamma, pll_cph, pll_if, pll_pvf))
#' maxy <- max(c(pll_gamma, pll_cph, pll_if, pll_pvf))
#'
#' plot(fr_var, pll_gamma,
#'      type = "l",
#'      xlab = "Frailty variance",
#'      ylab = "Profile log-likelihood",
#'      ylim = c(miny, maxy))
#' points(fr_var, pll_cph, col = 2)
#' lines(fr_var, pll_if, col = 3)
#' lines(fr_var, pll_pvf, col = 4)
#'
#' legend(legend = c("gamma (emfrail)", "gamma (coxph)", "inverse gaussian", "pvf, m=1.5"),
#'        col = 1:4,
#'        lty = 1,
#'        x = 0,
#'        y = (maxy + miny)/2)
#' }
#' # Recurrent events
#' mod_rec <- emfrail(bladder1, Surv(start, stop, status) ~ treatment + cluster(id))
#' # The warnings appear from the Surv object, they also appear in coxph.
#'
#' summary(mod_rec)
#'
#' # Create a histogram of the estimated frailties
#'
#' hist_frail(mod_rec)
#'
#' # or, with ggplot:
#' \dontrun{
#' library(ggplot2)
#' sum_mod_rec <- summary(mod_rec)
#'
#' ggplot(sum_mod_rec$frail, aes(x = z)) +
#'   geom_histogram() +
#'   ggtitle("Estimated frailties")
#'
#' # Plot the frailty estimates with quantiles of the estimated distribution
#' ggplot(sum_mod_rec$frail, aes(x = id, y = z)) +
#'   geom_point() +
#'   ggtitle("Estimated frailties") +
#'   geom_errorbar(aes(ymin = lower_q, ymax = upper_q))
#'
#' # We can do the same plot but with the ordered frailties:
#' ord <- order(sum_mod_rec$frail$z)
#' # we need new x coordinates for that:
#' ordering <- 1:length(ord)
#'
#' ggplot(sum_mod_rec$frail[ord,], aes(x = ordering, y = z)) +
#'   geom_point() +
#'   ggtitle("Estimated frailties") +
#'   geom_errorbar(aes(ymin = lower_q, ymax = upper_q))
#'
#' # How do we know which id is which one now?
#' # We can make an interactive plot with ggplotly
#' # To add text to elements we add id in aes()
#'
#' library(plotly)
#' ggplotly(
#'   ggplot(sum_mod_rec$frail[ord,], aes(x = ordering, y = z)) +
#'     geom_point(aes(id = id)) +
#'     ggtitle("Estimated frailties") +
#'     geom_errorbar(aes(ymin = lower_q, ymax = upper_q, id = id))
#' )
#'}
#'
#' # Plot marginal and conditional curves
#' # For recurrent events, the survival is not very meaningful
#'
#' plot_pred(mod_rec, quantity = "cumhaz")
#' #The strong frailty "drags down" the intensity
#'
#'
#'
#' # Left truncation
#'
#' # N.B. This code takes a longer time to run
#'
#' \dontrun{
#' # We simulate some data with truncation times
#' set.seed(1)
#' x <- sample(c(0,1), 5 * 300, TRUE)
#' u <- rep(rgamma(300,1,1), each = 5)
#' stime <- rexp(5*300, rate = u * exp(0.5 * x))
#' status <- ifelse(stime > 5, 0, 1)
#' stime[status] <- 5
#'
#' ltime <- runif(5 * 300, min = 0, max = 2)
#' d <- data.frame(id = rep(1:300, each = 5),
#'                 x = x,
#'                 stime = stime,
#'                 u = u,
#'                 ltime = ltime,
#'                 status = 1)
#' d_left <- d[d$stime > d$ltime,]
#'
#' # The worst thing that can be done is to
#' # Ignore the left truncation:
#' mod_1 <- emfrail(d_left,
#'                  Surv(stime, status)~ x + cluster(id))
#'
#' # The so-and-so thing is to consider the delayed entry time,
#' # But do not "update" the frailty distribution accordingly
#' mod_2 <- emfrail(d_left,
#'                  Surv(ltime, stime, status)~ x + cluster(id))
#'
#' # This is identical with
#' mod_cox <- coxph(Surv(ltime, stime, status)~ x + frailty(id), data = d_left)
#'
#'
#' # The correct thing is to update the frailty.
#' mod_3 <- emfrail(d_left,
#'                  Surv(ltime, stime, status)~ x + cluster(id),
#'                  .distribution = emfrail_distribution(left_truncation = TRUE))
#'
#' summary(mod_1)
#' summary(mod_2)
#' summary(mod_3)
#' }


emfrail <- function(.data,
                    .formula,
                    .distribution = emfrail_distribution(),
                    .control = emfrail_control()) {

  if(!inherits(.distribution, "emfrail_distribution"))
    stop(".distribution argument misspecified; see ?emfrail_distribution()")

  if(!inherits(.control, "emfrail_control"))
    stop(".control argument misspecified; see ?emfrail_control()")

  if(isTRUE(.control$inner_control$fast_fit)) {
    if(!(.distribution$dist %in% c("gamma", "pvf"))) {
      #message("fast_fit option only available for gamma and pvf with m=-1/2 distributions")
      .control$inner_control$fast_fit <- FALSE
    }

    # version 0.5.6, the IG fast fit gets super sensitive at small frailty variance...
    if(.distribution$dist == "pvf")
      .control$inner_control$fast_fit <- FALSE

  }


  Call <- match.call()


  if(missing(.formula)  | missing(.data)) stop("Missing arguments")

  cluster <- function(x) x
  terminal <- function(x) x

  mf <- model.frame(.formula, .data)


  # Identify the cluster and the ID column
  pos_cluster <- grep("cluster", names(mf))
  if(length(pos_cluster) != 1) stop("misspecified or non-specified cluster")
  id <- mf[[pos_cluster]]

  pos_terminal <- grep("terminal", names(mf))
  if(length(pos_terminal) > 1) stop("misspecified terminal()")

  Y <- mf[[1]]
  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {
    # warning("Y not in (tstart, tstop) format; taking tstart = 0")
    # Y <- cbind(rep(0, nrow(Y)), Y)
    # attr(Y, "dimnames") <- list(NULL, c("start", "stop", "status"))
    # attr(Y, "type") <- "counting"
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }
  if(attr(Y, "type") != "counting") stop("use Surv(tstart, tstop, status)")


  # get the model matrix
  X1 <- model.matrix(.formula, .data)
  # this is necessary because when factors have more levels, pos_cluster doesn't correspond any more
  pos_cluster_X1 <- grep("cluster", colnames(X1))
  pos_terminal_X1 <- grep("terminal", colnames(X1))
  X <- X1[,-c(1, pos_cluster_X1, pos_terminal_X1), drop=FALSE]
  # note: X has no attributes, in coxph it does.

  # some stuff for creating the C vector, is used all along.
  # mcox also works with empty matrices, but also with NULL as x.
  mcox <- survival::agreg.fit(x = X, y = Y, strata = NULL, offset = NULL, init = NULL,
                              control = survival::coxph.control(),
                              weights = NULL, method = "breslow", rownames = NULL)

  # the "baseline" case // this will stay constant

  if(length(X) == 0) {
    newrisk <- 1
    exp_g_x <- matrix(rep(1, length(mcox$linear.predictors)), nrow = 1)
    g <- 0
    g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))

  } else {
    x2 <- matrix(rep(0, ncol(X)), nrow = 1, dimnames = list(123, dimnames(X)[[2]]))
    x2 <- (scale(x2, center = mcox$means, scale = FALSE))
    newrisk <- exp(c(x2 %*% mcox$coefficients) + 0)
    exp_g_x <- exp(mcox$coefficients %*% t(X))
    g <- mcox$coefficients
    g_x <- t(mcox$coefficients %*% t(X))

  }

  explp <- exp(mcox$linear.predictors) # these are with centered covariates

  nev_id <- rowsum(Y[,3], id) # nevent per id or am I going crazy


  # for the baseline hazard how the fuck is that gonna happen?
  # Idea: nrisk has the sum of elp who leave later at every tstop
  # esum has the sum of elp who enter at every tstart
  # indx groups which esum is right after each nrisk;
  # the difference between the two is the sum of elp really at risk at that time point.


  nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
  esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))

  # the stuff that won't change
  death <- (Y[, ncol(Y)] == 1)
  nevent <- as.vector(rowsum(1 * death, Y[, ncol(Y) - 1])) # per time point
  time <- sort(unique(Y[,2])) # unique tstops

  # this gives the next entry time for each unique tstop (not only event)
  etime <- c(0, sort(unique(Y[, 1])),  max(Y[, 1]) + min(diff(time)))
  indx <- findInterval(time, etime, left.open = TRUE) # left.open  = TRUE is very important

  # this gives for every tstart (line variable) after which event time did it come
  # indx2 <- findInterval(Y[,1], time, left.open = FALSE, rightmost.closed = TRUE)
  indx2 <- findInterval(Y[,1], time)

  time_to_stop <- match(Y[,2], time)
  order_id <- findInterval(id, unique(id))

  atrisk <- list(death = death, nevent = nevent, nev_id = nev_id,
                 order_id = order_id, time = time, indx = indx, indx2 = indx2,
                 time_to_stop = time_to_stop)

  nrisk <- nrisk - c(esum, 0,0)[indx]

  haz <- nevent/nrisk * newrisk


  basehaz_line <- haz[atrisk$time_to_stop]

  cumhaz <- cumsum(haz)

  cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
  cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
  cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  * explp / newrisk

  Cvec <- rowsum(cumhaz_line, atrisk$order_id)

  ca_test <- NULL

  if(isTRUE(.control$ca_test)) ca_test <- ca_test_fit(mcox, X, atrisk, exp_g_x, cumhaz)
  if(isTRUE(.control$only_ca_test)) return(ca_test)


  if(isTRUE(.distribution$left_truncation)) {
    #indx2 <- findInterval(Y[,1], time, left.open = TRUE)
    cumhaz_tstop <- cumsum(haz)
    cumhaz_tstart <- c(0, cumhaz_tstop)[indx2 + 1]

    Cvec_lt <- rowsum(cumhaz_tstart, order_id)
    # Cvec_lt <- tapply(X = cumhaz_tstart,
    #                   INDEX = id,
    #                   FUN = sum)
  } else Cvec_lt <- 0 * Cvec


  # a fit just for the log-likelihood;
  if(!isTRUE(.control$opt_fit)) {
    return(em_fit(logfrailtypar = log(.distribution$theta),
           dist = .distribution$dist, pvfm = .distribution$pvfm,
           Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
           mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
           Cvec = Cvec, lt = .distribution$left_truncation,
           Cvec_lt = Cvec_lt, se = FALSE,
           inner_control = .control$inner_control))
  }



  # With the stable distribution, a problem pops up for small values, i.e. very large association (tau large)
  # So there is another interval...
  if(.distribution$dist == "stable") {
    .control$lik_ci_intervals$interval <- .control$lik_ci_intervals$interval_stable
  }

  # add a bit to the interval so that it gets to the Cox likelihood, if it is at that end of the parameter space

  # Maybe try nlm as well. Looks alright!

  # outer_m <- optimize(f = em_fit,
  #                     interval = .control$opt_control$interval + c(0, 0.1),
  #                     dist = .distribution$dist, pvfm = .distribution$pvfm,
  #          Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
  #          mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
  #          Cvec = Cvec, lt = .distribution$left_truncation,
  #          Cvec_lt = Cvec_lt,
  #          .control = .control)

  # Hessian
  # hess <- numDeriv::hessian(func = em_fit,
  #         x = outer_m$minimum,
  #         dist = .distribution$dist, pvfm = .distribution$pvfm,
  #         Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
  #         mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
  #         Cvec = Cvec, lt = .distribution$left_truncation,
  #         Cvec_lt = Cvec_lt,
  #         .control = .control)

  outer_m <- nlm(f = em_fit,
                 p = 2, hessian = TRUE,
                 dist = .distribution$dist, pvfm = .distribution$pvfm,
                 Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                 mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                 Cvec = Cvec, lt = .distribution$left_truncation,
                 Cvec_lt = Cvec_lt, se = FALSE,
                 inner_control = .control$inner_control)

  # do.call(nlm, c(list(f = em_fit, p = 2, hessian = TRUE, dist = .distribution$dist, pvfm = .distribution$pvfm,
  #                Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
  #                mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
  #                Cvec = Cvec, lt = .distribution$left_truncation,
  #                Cvec_lt = Cvec_lt,
  #                .control = .control), .control$opt_control))

    # likelihood-based confidence intervals
  theta_low <- theta_high <- NULL
  if(isTRUE(.control$lik_ci)) {

  lower_llik <- em_fit(.control$lik_ci_intervals$interval[1],
                       dist = .distribution$dist,
                       pvfm = .distribution$pvfm,
                       Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                       mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                       Cvec = Cvec, lt = .distribution$left_truncation,
                       Cvec_lt = Cvec_lt, se = FALSE,
                       inner_control = .control$inner_control)

  upper_llik <- em_fit(.control$lik_ci_intervals$interval[2],
         dist = .distribution$dist,
         pvfm = .distribution$pvfm,
         Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
         mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
         Cvec = Cvec, lt = .distribution$left_truncation,
         Cvec_lt = Cvec_lt, se = FALSE,
         inner_control = .control$inner_control)

  theta_low <- uniroot(function(x, ...) outer_m$minimum - em_fit(x, ...) + 1.92,
                       interval = c(.control$lik_ci_intervals$interval[1], outer_m$estimate),
                       f.lower = outer_m$minimum - lower_llik + 1.92, f.upper = 1.92,
                       tol = .Machine$double.eps^0.1,
                       dist = .distribution$dist,
                       pvfm = .distribution$pvfm,
                       Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                       mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                       Cvec = Cvec, lt = .distribution$left_truncation,
                       Cvec_lt = Cvec_lt, se = FALSE,
                       inner_control = .control$inner_control,
                       maxiter = 100)$root


  # this says that if I can't get a significant difference on the right side
  # then screw this it's infinity
  if(upper_llik  - outer_m$minimum < 1.92) theta_high <- Inf else
    theta_high <- uniroot(function(x, ...) outer_m$minimum - em_fit(x, ...) + 1.92,
                          interval = c(outer_m$estimate, .control$lik_ci_intervals$interval[2]),
                          f.lower = 1.92, f.upper = outer_m$minimum - upper_llik + 1.92,
                          extendInt = c("downX"),
                          dist = .distribution$dist,
                          pvfm = .distribution$pvfm,
                          Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                          mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                          Cvec = Cvec, lt = .distribution$left_truncation,
                          Cvec_lt = Cvec_lt, se  = FALSE,
                          inner_control = .control$inner_control)$root
  }

  # outer_m$minimum
  # em_fit(0.45,
  #        dist = .distribution$dist,
  #        pvfm = .distribution$pvfm,
  #        Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
  #        mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
  #        Cvec = Cvec, lt = .distribution$left_truncation,
  #        Cvec_lt = Cvec_lt,
  #        .control = .control)

  #message("Calculating final fit with information matrix...")

  if(isTRUE(.control$se))  {
    inner_m <- em_fit(logfrailtypar = outer_m$estimate,
                      dist = .distribution$dist, pvfm = .distribution$pvfm,
                      Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                      mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                      Cvec = Cvec, lt = .distribution$left_truncation,
                      Cvec_lt = Cvec_lt, se = TRUE,
                      inner_control = .control$inner_control,
                      return_loglik = FALSE)
  } else
    inner_m <- em_fit(logfrailtypar = outer_m$estimate,
                      dist = .distribution$dist, pvfm = .distribution$pvfm,
                      Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                      mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                      Cvec = Cvec, lt = .distribution$left_truncation,
                      Cvec_lt = Cvec_lt, se = FALSE,
                      inner_control = .control$inner_control,
                      return_loglik = FALSE)


  # adjusted standard errors
  if(isTRUE(.control$se) & isTRUE(.control$se_adj)) {

    # absolute value should be redundant. but sometimes the "hessian" might be 0.
    # in that case it might appear negative; this happened only on Linux...
    # h <- as.numeric(sqrt(abs(1/(attr(outer_m, "details")[[3]])))/2)
    h<- as.numeric(sqrt(abs(1/outer_m$hessian))/2)
    lfp_minus <- max(outer_m$estimate - h , outer_m$estimate - 5)
    lfp_plus <- min(outer_m$estimate + h , outer_m$estimate + 5)

    message("Calculating adjustment for information matrix...")


    final_fit_minus <- em_fit(logfrailtypar = lfp_minus,
                              dist = .distribution$dist, pvfm = .distribution$pvfm,
                              Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                              mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                              Cvec = Cvec, lt = .distribution$left_truncation,
                              Cvec_lt = Cvec_lt, se = TRUE,
                              inner_control = .control$inner_control,
                              return_loglik = FALSE)

    final_fit_plus <- em_fit(logfrailtypar = lfp_plus,
                             dist = .distribution$dist, pvfm = .distribution$pvfm,
                             Y = Y, Xmat = X, atrisk = atrisk, basehaz_line = basehaz_line,
                             mcox = list(coefficients = g, loglik = mcox$loglik),  # a "fake" cox model
                             Cvec = Cvec, lt = .distribution$left_truncation,
                             Cvec_lt = Cvec_lt, se = TRUE,
                             inner_control = .control$inner_control, return_loglik = FALSE)


    # instructional: this should be more or less equal to the
    # -(final_fit_plus$loglik + final_fit_minus$loglik - 2 * inner_m$loglik)/h^2

    # se_logtheta^2 / (2 * (final_fit$loglik -final_fit_plus$loglik ))

    deta_dtheta <- (c(final_fit_plus$coef, final_fit_plus$haz) -
                      c(final_fit_minus$coef, final_fit_minus$haz)) / (2*h)

    #adj_se <- sqrt(diag(deta_dtheta %*% (1/(attr(opt_object, "details")[[3]])) %*% t(deta_dtheta)))

    # vcov_adj = inner_m$Vcov + deta_dtheta %*% (1/(attr(outer_m, "details")[[3]])) %*% t(deta_dtheta)
    vcov_adj = inner_m$Vcov + deta_dtheta %*% (1/outer_m$hessian) %*% t(deta_dtheta)

  } else vcov_adj = matrix(NA, nrow(inner_m$Vcov), nrow(inner_m$Vcov))



  if(length(pos_terminal_X1) > 0 & .distribution$dist == "gamma") {
    Y[,3] <- X1[,pos_terminal_X1]

    Mres <- survival::agreg.fit(x = X, y = Y, strata = NULL, offset = NULL, init = NULL,
                        control = survival::coxph.control(),
                        weights = NULL, method = "breslow", rownames = NULL)$residuals
    Mres_id <- rowsum(Mres, atrisk$order_id)

    theta <- exp(outer_m$minimum)

    fr <- with(inner_m, estep[,1] / estep[,2])

    numerator <- theta + inner_m$nev_id
    denominator <- numerator / fr

    lfr <- digamma(numerator) - log(denominator)

    lfr2 <- (digamma(numerator))^2 + trigamma(numerator) - (log(denominator))^2 - 2 * log(denominator) * lfr

    # score test 1 I think
    r <- cor(lfr, Mres_id)
    tr <- r* sqrt((length(fr) - 2) / (1 - r^2))
    p.cor <- pchisq(tr^2, df = 1, lower.tail = F)

    cens_test = c(tstat = tr, pval = p.cor)
  } else cens_test = NULL

  res <- list(outer_m = list(objective = outer_m$minimum,
                             minimum = outer_m$estimate,
                             hess = outer_m$hessian,
                             ltheta_low = theta_low,
                             ltheta_high = theta_high),
              inner_m = inner_m,
              loglik_null = mcox$loglik[length(mcox$loglik)],
              # mcox = mcox,
              vcov_adj = vcov_adj,
              ca_test = ca_test,
              cens_test = cens_test,
              .formula = .formula,
              .distribution = .distribution,
              .control = .control
              )


  # these are things that make the predict work and other methods
  terms_2 <- delete.response(attr(mf, "terms"))
  pos_cluster_2 <- grep("cluster", attr(terms_2, "term.labels"))
  terms <- drop.terms(terms_2, pos_cluster_2)
  myxlev <- .getXlevels(terms, mf)
  attr(res, "metadata") <- list(terms, myxlev)
  attr(res, "call") <-  Call
  attr(res, "class") <- "emfrail"


  res


}
