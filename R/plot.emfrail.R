#' Plots for emfrail objects
#'
#' @importFrom graphics abline hist legend lines plot
#' @param x \code{emfrail} object, typically result of \code{emfrail()}
#' @param type One (or more) of \code{hist} for a histogram of the estimated frailty values,
#'  \code{hr} for a plot of the conditional and marginal hazard ratio between two cases and
#'  \code{pred} for the predicted conditional and marginal cumulative hazard or survival for one case
#' @param newdata A \code{data.frame} with values of the covariates. For \code{type == "hr"} the hazard ratio
#' between the first two rows of \code{newdata} is calculated. For \code{type == "pred"} the prediction
#' for the first row of \code{newdata} is calculated.
#' @param lp A numeric vector of values of the linear predictor, each corresponding to a case. For \code{type == "hr"} the hazard ratio
#' between the first two values of \code{lp} is calculated. For \code{type == "pred"} the prediction
#' for the first value of \code{lp} is calculated.
#' @param quantity For \code{type == "pred"} the predicted quantity; see \code{\link{predict.emfrail}}
#' @param marg_cond For \code{type == "pred"} the type of predicted quantity; see \code{\link{predict.emfrail}}
#' @param conf_int For \code{type == "pred"} the type of confidence intervals; see \code{\link{predict.emfrail}}
#' @param ... Further arguments to be passed to the plot function
#'
#' @return Nothing
#' @export
#'
#' @seealso \code{\link{predict.emfrail}}, \code{\link{summary.emfrail}}, \code{\link{autoplot.emfrail}}.
#'
#' @examples
#' mod_rec <- emfrail(bladder1, Surv(start, stop, status) ~ treatment + number + cluster(id))
#' summary(mod_rec)
#'
#' # cumulative hazard
#' # Note: this individual has number = 0, which does not exist in the data
#' plot(mod_rec, type = "pred", lp = 0)
#'
#' # survival, although not very meaningful with recurrent events
#' plot(mod_rec, type = "pred", lp = 0, quantity = "survival")
#'
#'
#' # For an individual with number == 2
#' plot(mod_rec, type = "pred",
#' newdata = data.frame(treatment = "placebo", number = 2))
#'
#' # hazard ratio between an individual with 0 and with 2 recurrences at baseline
#' # the marginal hazard ratio is "pulled" towards 1:
#'
#' plot(mod_rec, type = "hr",
#' newdata = data.frame(treatment = c("placebo", "placebo"), number = c(0, 2)))
#'
#' # hazard ratio with the stable distribution:
#' \dontrun{
#' mod_rec_stab <- emfrail(bladder1,
#'                         Surv(start, stop, status) ~ treatment + number + cluster(id),
#'                         .distribution = emfrail_distribution(dist = "stable"))
#'
#' plot(mod_rec_stab, type = "hr",
#' newdata = data.frame(treatment = c("placebo", "placebo"), number = c(0, 2)))
#'
#' # histogram of frailty estimates
#' plot(mod_rec_stab, type = "hist")
#' }
plot.emfrail <- function(x, type = c("hist", "hr", "pred"), newdata = NULL, lp = NULL, quantity = "cumhaz", marg_cond = "both", conf_int = "adjusted", ...) {

  if("hist" %in% type) {
    sobj <- summary.emfrail(x)
    hist(sobj$frail$z, xlab = "Frailties", main = "Fraily estimates", ... )
  }

  if("hr" %in% type) {

    if(is.null(newdata))
      if(missing(lp)) stop("lp or newdata should be specified") else
        if(length(lp) < 2) stop("lp must be at least of length 2 or newdata must have 2 rows with type = hr") else {
          if(length(lp) > 2) warning("just the first 2 lp values are used for type = hr")
          if(lp[1] == lp[2]) stop("values of lp should be different")
        }


    if(!is.null(newdata))
      if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame") else
        if(nrow(newdata) < 2) stop("newdata must have 2 rows with type = hr") else {
          if(nrow(newdata) > 2) warning("just the first 2 newdata rows are used for type = hr")
          if(isTRUE(all.equal(newdata[1,], newdata[2,]))) stop("the rows of newdata must be different")
        }


    p1 <- predict.emfrail(x,
                          lp = lp[1],
                          newdata = newdata[1,],
                          quantity = "cumhaz",
                          conf_int = NULL)

    p2 <- predict.emfrail(x,
                          lp = lp[2],
                          newdata = newdata[2,],
                          quantity = "cumhaz",
                          conf_int = NULL)

    hr_cond <- p1$cumhaz / p2$cumhaz
    hr_mar <- p1$cumhaz_m / p2$cumhaz_m

    hr_cond[1] <- hr_cond[2]  # that's cause it's 0/0
    hr_mar[1] <- hr_mar[2]
    ylim <- c(min(c(hr_cond, hr_mar, 0.95)), max(c(hr_cond, hr_mar, 1.05)))

    plot(p1$time, hr_cond, type = "s", ylim = ylim,
         xlab = "time",
         ylab = "hazard ratio", ...)
    lines(p1$time, hr_mar, type = "s", col = 2)
    abline(a = 1, b = 0, lty = 2, col = "gray")

    if(hr_mar[length(hr_mar)] > 1) pos <- "bottomright" else
      pos <- "topright"
    legend(x = pos, legend = c("conditional", "marginal"), col = 1:2, lty = 1)
  }

  if("pred" %in% type) {

    if(is.null(newdata) & is.null(lp)) stop("newdata or lp have to be specified for type = pred")

    if(!is.null(newdata))
      if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame") else
        if(nrow(newdata) < 1) stop("newdata must have eat least 1 row for type = pred") else
          if(nrow(newdata) > 1) warning("just the first row is used for type = pred")

    if(is.null(newdata) & !is.null(lp))
      if(length(lp) < 1) stop("lp should have at least length 1") else
        if(length(lp) > 1) warning("just the first value of lp is used for type = pred")
    if(!(quantity %in% c("cumhaz", "survival")) | length(quantity) != 1)
      stop("quantity must be cumhaz or survival for type = pred")

    if(!(marg_cond %in% c("conditional", "marginal", "both")) | length(marg_cond) != 1)
      stop("marg_cond must be conditional, marginal, or both for type = pred")

    if(!(conf_int %in% c("adjusted", "regular", "none")) | length(conf_int) != 1)
      stop("conf_int must be adjusted, regular or none for type = pred")

    if(marg_cond == "both")
      type_pred <- c("conditional", "marginal") else
        type_pred <- marg_cond

      p1 <- predict.emfrail(x,
                            lp = lp[1],
                            newdata = newdata[1,],
                            quantity = quantity,
                            type = type_pred,
                            conf_int = conf_int)

      # determine columns to be plotted by name,
      # possibilities: quantity is just one of survival or cumhaz
      # they may have an _m or not at the end,

      if(marg_cond == "conditional") {
        if(quantity == "cumhaz") {
          # if(conf_int == "adjusted") ylim <- c(0, max(p1$cumhaz_r_a)) else
          #   if(conf_int == "regular") ylim <- c(0, max(p1$cumhaz_r)) else
          #     ylim <- c(0, max(p1$cumhaz_m))

          with(p1, plot(time, cumhaz,
                        type = "s",
                        #main = "Cumulative hazard",
                        ylab = "cumhaz",
                        xlab = "time",
                        #ylim = ylim,
                        ...))
          if(conf_int == "adjusted") {
            with(p1, lines(time, cumhaz_l_a, lty = 2))
            with(p1, lines(time, cumhaz_r_a, lty = 2))
          }
          if(conf_int == "regular") {
            with(p1, lines(time, cumhaz_l, lty = 2))
            with(p1, lines(time, cumhaz_r, lty = 2))
          }
        }

        if(quantity == "survival") {
          with(p1, plot(time, survival,
                        type = "s",
                        #main = "Survival",
                        ylab = "S(t)",
                        xlab = "time",
                        ylim = c(0,1),
                        ...))
          if(conf_int == "adjusted") {
            with(p1, lines(time, survival_l_a, lty = 2))
            with(p1, lines(time, survival_r_a, lty = 2))
          }
          if(conf_int == "regular") {
            with(p1, lines(time, survival_l, lty = 2))
            with(p1, lines(time, survival_r, lty = 2))
          }
        }
      }

      if(marg_cond == "marginal") {
        if(quantity == "cumhaz") {
          # if(conf_int == "adjusted") ylim <- c(0, max(p1$cumhaz_m_r)) else
          #   if(conf_int == "regular") ylim <- c(0, max(p1$cumhaz_m_r_a)) else
          #     ylim <- c(0, max(p1$cumhaz_m))

          with(p1, plot(time, cumhaz_m,
                        type = "s",
                        #main = "Cumulative hazard",
                        ylab = "cumhaz",
                        xlab = "time",
                        # ylim = ylim,
                        ...))
          if(conf_int == "adjusted") {
            with(p1, lines(time, cumhaz_m_l_a, lty = 2))
            with(p1, lines(time, cumhaz_m_r_a, lty = 2))
          }
          if(conf_int == "regular") {
            with(p1, lines(time, cumhaz_m_l, lty = 2))
            with(p1, lines(time, cumhaz_m_r, lty = 2))
          }
        }

        if(quantity == "survival") {
          with(p1, plot(time, survival_m,
                        type = "s",
                        #main = "Survival",
                        ylab = "S(t)",
                        xlab = "time",
                        ylim = c(0,1),
                        ...))
          if(conf_int == "adjusted") {
            with(p1, lines(time, survival_m_l_a, lty = 2))
            with(p1, lines(time, survival_m_r_a, lty = 2))
          }
          if(conf_int == "regular") {
            with(p1, lines(time, survival_m_l, lty = 2))
            with(p1, lines(time, survival_m_r, lty = 2))
          }
        }


      }

      if(marg_cond == "both") {
        if(quantity == "cumhaz") {
          # if(conf_int == "adjusted") ylim <- c(0, max(p1$cumhaz_r_a)) else
          #   if(conf_int == "regular") ylim <- c(0, max(p1$cumhaz_r)) else
          #     ylim <- c(0, max(p1$cumhaz))

          with(p1, plot(time, cumhaz,
                        type = "s",
                        #main = "Cumulative hazard",
                        ylab = "cumhaz",
                        xlab = "time",
                        # ylim = ylim,
                        col = 1,
                        ...))
          with(p1, lines(time, cumhaz_m,
                         type = "s",
                         col = 2,
                         ...))

          if(conf_int == "adjusted") {
            with(p1, lines(time, cumhaz_m_l_a, lty = 2, col = 2))
            with(p1, lines(time, cumhaz_m_r_a, lty = 2, col = 2))
            with(p1, lines(time, cumhaz_l_a, lty = 2))
            with(p1, lines(time, cumhaz_r_a, lty = 2))
          }

          if(conf_int == "regular") {
            with(p1, lines(time, cumhaz_m_l, lty = 2, col = 2))
            with(p1, lines(time, cumhaz_m_r, lty = 2, col = 2))
            with(p1, lines(time, cumhaz_l, lty = 2))
            with(p1, lines(time, cumhaz_r, lty = 2))
          }


          legend(x = "topleft", legend = c("conditional", "marginal"), col = 1:2, lty = 1)
        }

        if(quantity == "survival") {
          with(p1, plot(time, survival_m,
                        type = "s",
                        #main = "Survival",
                        ylab = "S(t)",
                        xlab = "time",
                        col = 2,
                        ylim = c(0,1),
                        ...))
          with(p1, lines(time, survival,
                         type = "s",
                         col = 1,
                         ...))

          if(conf_int == "adjusted") {
            with(p1, lines(time, survival_m_l_a, lty = 2, col = 2))
            with(p1, lines(time, survival_m_r_a, lty = 2, col = 2))
            with(p1, lines(time, survival_l_a, lty = 2))
            with(p1, lines(time, survival_r_a, lty = 2))
          }
          if(conf_int == "regular") {
            with(p1, lines(time, survival_m_l, lty = 2, col = 2))
            with(p1, lines(time, survival_m_r, lty = 2, col = 2))
            with(p1, lines(time, survival_l, lty = 2))
            with(p1, lines(time, survival_r, lty = 2))
          }

          legend(x = "topright", legend = c("conditional", "marginal"), col = 1:2, lty = 1)
        }


      }

  }

}
