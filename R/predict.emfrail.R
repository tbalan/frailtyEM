#' Predicted hazard and survival curves from an \code{emfrail} object
#'
#' @importFrom stats .getXlevels delete.response drop.terms as.formula
#' @param object An \code{emfrail} fit object
#' @param lp A vector of linear predictor values at which to calculate the curves. Default is 0 (baseline).
#' @param newdata A data frame with the same variable names as those that appear in the \code{emfrail} formula, used to calculate the \code{lp} (optional).
#' @param quantity Can be \code{"cumhaz"} and/or \code{"survival"}. The quantity to be calculated for the values of \code{lp}.
#' @param type Can be \code{"conditional"} and/or \code{"marginal"}. The type of the quantity to be calculated.
#' @param conf_int Can be \code{"regular"} and/or \code{"adjusted"}. The type of confidence interval to be calculated.
#' @param individual Logical. Are the observations in \code{newdata} from the same individual? See details.
#' @param ... Ignored
#'
#' @return If \code{lp} has length 1, \code{newdata} has 1 row or \code{individual == TRUE}, then the return object
#' is a \code{data.frame}; otherwise, it is a list of \code{data.frame}s, where each element is the prediction
#' for the corresponding element of \code{lp} or row of \code{newdata}.
#' By default, for each \code{lp} it will give the following columns: \code{cumhaz}, \code{survival},
#' \code{cumhaz_m}, \code{survival_m} for the cumulative hazard and survival, conditional and marginal.
#'
#' @details
#' This method calculates predicted cumulative hazard and survival curves for fixed covariate values.
#' There are two ways of inputing this values. The first is via \code{lp},
#' where the function expects a vector of values, of the linear predictor or via \code{newdata},
#' where the function expects a \code{data.frame} with columns of the covariates matching those in the
#' \code{emfrail} fit. If \code{newdata} is specified, the at-risk indicator can be specified by including
#' two columns named \code{tstart} and \code{tstop} specifying the time at risk.
#'
#' If \code{individual == FALSE}, each value of \code{lp} or each row of \code{newdata} will result in a
#' data frame with the required predictions. If \code{individual == TRUE}, then the rows of \code{newdata}
#' are assumed to come from the same individual and the function will return a single data frame.
#'  In this case, the specification of the \code{tstart} and \code{tstop} columns is mandatory.
#'
#' The names of the columns in the returned data frames are as follows: \code{time} represents the unique event time points
#' from the data set, \code{lp} is the value of the linear predictor (as specified in the input or as calculated from the lines of \code{newdata}).
#'
#' The two "quantities" that can be returned are
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
#' The linear predictor is taken as fixed,
#' i.e. the variability in the estimation of the regression coefficient is not taken into account.
#'
#'
#'
#' @export
#'
#'
#' @examples
#' kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
#' m1 <- emfrail(formula = Surv(time, status) ~  sex + age  + cluster(id),
#'               data =  kidney,)
#'
#' pred <- predict(m1, lp = 0)
#'
#' names(pred)
#'
#' # Plot baseline cumulative hazard: note that is for someone aged 0!
#' plot(m1, type = "pred", lp = 0)
#'
#' # More realistic:
#' plot(m1, type = "pred", newdata = data.frame(sex = "female", age = mean(kidney$age)))
#'
#' # Plot survival
#' plot(m1,
#' type = "pred",
#'           newdata = data.frame(sex = "female", age = mean(kidney$age)),
#'           quantity = "survival", conf_int = "none")
#' \dontrun{
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
#' plot(m1, type = "hr", newdata = data.frame(sex = c("female", "male"), age = c(44, 44)))
#' }
#' @seealso \code{\link{plot.emfrail}}, \code{\link{autoplot.emfrail}}
predict.emfrail <- function(object,
                            lp = NULL,
                            newdata = NULL,
                            quantity = c("cumhaz", "survival"),
                            type = c("conditional", "marginal"),
                            conf_int = c("regular", "adjusted"),
                            individual = FALSE,
                            ...) {


  if(is.null(newdata) & is.null(lp)) stop("lp or newdata must be specified")

  if(!is.null(lp) & !is.null(newdata)) stop("specify either lp or newdata")

  if(!is.null(lp)) {
    if(isTRUE(individual)) stop("if lp is specified, individual must be FALSE")
    tstart <- 0
    tstop <- Inf
  }

  if(!is.null(newdata)) {
    if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame")

    mdata <- attr(object, "metadata")
    mf <- model.frame(mdata[[1]], data = newdata, xlev = mdata[[2]])
    mm <- try(model.matrix(mdata[[1]], mf)[,-1])
    if(inherits(mm, "try-error")) stop("newdata probably misspecified")
    lp <- as.numeric(mm %*% object$coef)

    if(!("tstart" %in% names(newdata))) {
      tstart = 0
      tstop = Inf
    } else
        if(!("tstop" %in% names(newdata)))
          stop("if tstart is specified then also tstop must be specified") else {
            if(any(newdata$tstop <= newdata$tstart)) stop("tstop must be larger than tstart")

            # check if there is an overlap
            if(individual == TRUE) {
              if(any(diff(as.numeric(as.matrix(newdata[c("tstart", "tstop")]))) < 0))
                stop("(tstart,tstop) must be ordered and without overlap")
            }
            tstart <- newdata$tstart
            tstop <- newdata$tstop
          }

    # at this point I should have a vector of lp
    # if individual == TRUE then I should also have tstart tstop for each lp
  }

    list_haz <- mapply(FUN = function(lp, haz, tstart, tstop, tev) {
      cbind(tev[tstart <= tev & tev < tstop], lp, haz[tstart <= tev & tev < tstop] * exp(lp))
    }, lp, list(object$hazard), tstart, tstop, list(object$tev), SIMPLIFY = FALSE)

    if(isTRUE(individual)) {
      res <- as.data.frame(do.call(rbind, list_haz))
      names(res) <- c("time", "lp", "cumhaz")
      res$cumhaz <- cumsum(res$cumhaz) # a bit dodgy
      attr(res, "bounds") <- "cumhaz"
      res <- list(res)
    } else
      res <- lapply(list_haz, function(x) {res <- as.data.frame(x)
                    names(res) <- c("time", "lp", "cumhaz")
                    res$cumhaz <- cumsum(res$cumhaz)
                    attr(res, "bounds") <- "cumhaz"
              res})

    # now res a list of data frames


  ncoef <- length(object$coef)


  # Here we start putting a bunch of standard errors and stuff inside
  # for each lp we will get a data frame with a "bounds" attribute.
  # this attributie is meant to keep track of which columns we have in the data frame
  if(("regular" %in% conf_int) | ("adjusted" %in% conf_int)) {
    varH <-   object$var[(ncoef + 1):nrow(object$var), (ncoef + 1):nrow(object$var)]
    varH_adj <- object$var_adj[(ncoef + 1):nrow(object$var_adj), (ncoef + 1):nrow(object$var_adj)]

    res <- lapply(res, function(x) {
      times_res <- match(x$time, object$tev)

      loghaz <- log(object$hazard[times_res])

      # Now here I build a bunch of formulas - that is to get the confidence intervals with the delta mehtod
      xs <- lapply(seq_along(times_res), function(x) text1 <- paste0("x", x))
      for(i in 2:length(xs)) {
        xs[[i]] = paste0(xs[[i-1]], " + ", xs[[i]])
      }
      forms <- lapply(xs, function(x) as.formula(paste0("~log(", x, ")")))

      # attr(x, "bounds") <- "cumhaz"

      if("regular" %in% conf_int) {

        x$se_logH <- msm::deltamethod(g = forms,
                                      mean = object$hazard[times_res],
                                      cov = varH[times_res, times_res],
                                      ses = TRUE)

        attr(x, "bounds") <- c(attr(x, "bounds"), "cumhaz_l", "cumhaz_r")
        x$cumhaz_l <- pmax(0, exp(log(x$cumhaz) - 1.96*x$se_logH))
        x$cumhaz_r <- exp(log(x$cumhaz) + 1.96*x$se_logH)
      }

      if("adjusted" %in% conf_int) {
        x$se_logH_adj <- msm::deltamethod(g = forms,
                                          mean = object$hazard[times_res],
                                          cov = varH_adj[times_res, times_res],
                                          ses = TRUE)

        attr(x, "bounds") <- c(attr(x, "bounds"), "cumhaz_l_a", "cumhaz_r_a")
        x$cumhaz_l_a <- pmax(0, exp(log(x$cumhaz) - 1.96*x$se_logH_adj))
        x$cumhaz_r_a <- exp(log(x$cumhaz) + 1.96*x$se_logH_adj)
      }



      x
    })
  }


  est_dist <- object$distribution
  est_dist$frailtypar <- exp(object$logtheta)

  # the way to convert the cumulative hazard to all sort of quantities
  chz_to_surv <- function(x) exp(-x)
  surv_to_chz <- function(x) -1 * log(x)
  survival <- function(chz) exp(-chz)
  survival_m <- function(chz) laplace_transform(chz, est_dist)
  cumhaz_m <- function(chz) -1 * log(laplace_transform(chz, est_dist))

  scenarios <- expand.grid(quantity, type, conf_int)

  res <- lapply(res, function(x) {

    # bounds gets kicked out by cbind

    bounds <- attr(x, "bounds")
    if("survival" %in% quantity & "conditional" %in% type) {
      x <- cbind(x,
                 as.data.frame(lapply(x[bounds], survival),
                               col.names = sub("cumhaz", "survival", bounds)))
    }

    if("survival" %in% quantity & "marginal" %in% type) {
      x <- cbind(x,
                 as.data.frame(lapply(x[bounds], survival_m),
                               col.names = sub("cumhaz", "survival_m", bounds)))
    }

    if("cumhaz" %in% quantity & "marginal" %in% type) {
      x <- cbind(x,
                 as.data.frame(lapply(x[bounds], cumhaz_m),
                               col.names = sub("cumhaz", "cumhaz_m", bounds)))
    }

    if(!("cumhaz" %in% quantity & "conditional" %in% type)) {
      cols_to_keep <- which(!(names(x) %in% c(bounds, "se_logH", "se_logH_adj")))
      x <- x[,cols_to_keep]
    }
    x

  })


  if(length(res) == 1) res <- res[[1]]
  res
}
