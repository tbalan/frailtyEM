#' @name plot_emfrail
#' @rdname plot_emfrail
#'
#' @title Plot functions for \code{emfrail} objects
#'
#' @param object An \code{emfrail} object
#' @param lp The value(s) of the linear predictor. For \code{plot_pred} this should have length 1 and for \code{plot_hr} length 2.
#' @param newdata A \code{data.frame} from each each line corresponds to a value of the linear predictor (optional). For \code{plot_pred} this should have 1 row and for \code{plot_hr} 2 rows.
#' @param quantity One of "cumhaz" (for cumulative hazard) or "survival"
#' @param type One of "conditional", "marginal" or "both"
#' @param conf_int One of "adjusted", "regular", or "none
#' @param ... Parameters passed on to plot functions
#'
#' @details These functions exist mostly for conveince. They are in fact simple wrappers that use \code{predict.emfrail} or \code{summary.emfrail} on
#' \code{object}, extract some quantities of interest, and plot them.
#'
#' @seealso \code{\link{predict.emfrail}}, \code{\link{summary.emfrail}}, \code{\link{ggplot_emfrail}}.
#' @return Nothing.
#' @importFrom graphics abline hist legend lines plot
#' @examples
#' mod_rec <- emfrail(bladder1, Surv(start, stop, status) ~ treatment + number + cluster(id))
#' summary(mod_rec)
#'
#' # cumulative hazard
#' # Note: this individual has number = 0, which does not exist in the data
#' plot_pred(mod_rec)
#'
#' # survival, although not very meaningful with recurrent events
#' plot_pred(mod_rec, quantity = "survival")
#'
#'
#' # For an individual with number == 2
#' plot_pred(mod_rec, newdata = data.frame(treatment = "placebo", number = 2))
#'
#' # hazard ratio between an individual with 0 and with 2 recurrences at baseline
#' # the marginal hazard ratio is "pulled" towards 1:
#'
#' plot_hr(mod_rec, newdata = data.frame(treatment = "placebo", number = c(0, 2)))
#'
#' # hazard ratio with the stable distribution:
#' \dontrun{
#' mod_rec_stab <- emfrail(bladder1,
#'                         Surv(start, stop, status) ~ treatment + number + cluster(id),
#'                         .distribution = emfrail_distribution(dist = "stable"))
#'
#' plot_hr(mod_rec_stab, newdata = data.frame(treatment = "placebo", number = c(0, 2)))
#'
#' # histogram of frailty estimates
#' hist_frail(mod_rec_stab)
#' }
NULL

#' \code{plot_pred} plots predicted cumulative hazard or survival curves, marginal and / or conditional, with or without confidence intervals.
#' @rdname plot_emfrail
#' @name plot_pred
#' @export
#'
plot_pred <- function(object, lp = 0, newdata = NULL, quantity = "cumhaz", type = "both", conf_int = "adjusted", ...) {

  if(is.null(newdata)) {
    if(length(lp) != 1) stop("lp must be a number")
  }

  if(!is.null(newdata))
    if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame") else
      if(nrow(newdata) != 1) stop("newdata must have exactly 1 row")

  if(!(quantity %in% c("cumhaz", "survival")) | length(quantity) != 1)
    stop("quantity must be cumhaz or survival ")

  if(!(type %in% c("conditional", "marginal", "both")) | length(type) != 1)
    stop("type must be conditional, marginal, or both ")

  if(!(conf_int %in% c("adjusted", "regular", "none")) | length(conf_int) != 1)
    stop("conf_int must be adjusted, regular or none")

  if(type == "both")
    type_pred <- c("conditional", "marginal") else
      type_pred <- type

    p1 <- predict.emfrail(object,
                          lp = lp,
                          newdata = newdata,
                          quantity = quantity,
                          type = type_pred,
                          conf_int = conf_int)

    # determine columns to be plotted by name,
    # possibilities: quantity is just one of survival or cumhaz
    # they may have an _m or not at the end,
    if(lp == 0) baseline <- TRUE else baseline <- FALSE

    if(type == "conditional") {
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

    if(type == "marginal") {
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

    if(type == "both") {
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

#' \code{plot_hr} plots the estimated marginal and conditional hazard ratio between two units with different linear predictor values.
#' @rdname plot_emfrail
#' @name plot_hr
#' @export
#'
plot_hr <- function(object, lp, newdata = NULL, ...) {

  if(is.null(newdata))
    if(missing(lp)) stop("lp or newdata should be specified") else
      if(length(lp) != 2) stop("lp must be of length 2, or newdata must have 2 rows") else
        if(lp[1] == lp[2]) stop("values of lp should be different")

  if(!is.null(newdata))
    if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame") else
      if(nrow(newdata) != 2) stop("newdata must have exactly 2 rows") else
        if(isTRUE(all.equal(newdata[1,], newdata[2,]))) stop("the rows of newdata must be different")

  p1 <- predict.emfrail(object,
                        lp = lp[1],
                        newdata = newdata[1,],
                        quantity = "cumhaz",
                        conf_int = NULL)

  p2 <- predict.emfrail(object,
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


#' \code{hist_frail} plots a histogram of the estimated frailties.
#' @rdname plot_emfrail
#' @name hist_frail
#'
#' @export
#'
#' @examples
#'
hist_frail <- function(object, ...) {
  sobj <- summary.emfrail(object)
  hist(sobj$frail$z, xlab = "Frailties", main = "Fraily estimates", ... )
}


