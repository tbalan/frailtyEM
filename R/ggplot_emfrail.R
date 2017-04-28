#' @name ggplot_emfrail
#' @rdname ggplot_emfrail
#'
#' @title Plot functions for \code{emfrail} objects using \code{ggplot2}
#'
#' @param object An \code{emfrail} object
#' @param lp The value(s) of the linear predictor. For \code{ggplot_pred} this should have length 1 and for \code{ggplot_hr} length 2.
#' @param newdata A \code{data.frame} from each each line corresponds to a value of the linear predictor (optional). For \code{ggplot_pred} this should have 1 row and for \code{plot_hr} 2 rows.
#' @param quantity One of "cumhaz" (for cumulative hazard) or "survival"
#' @param type One of "conditional", "marginal" or "both"
#' @param conf_int One of "adjusted", "regular", or "none
#' @param ... Parameters passed on to plot functions
#'
#' @details These functions exist mostly for conveince. They are in fact simple wrappers that use \code{predict.emfrail} or \code{summary.emfrail} on
#' \code{object}, extract some quantities of interest, and plot them. In \code{ggplot_frail},
#' if the object was fitted with a gamma distribution,
#' then quantiles of the posterior distribution of the random effects are also plotted.
#'
#' @importFrom ggplot2 ggplot geom_step geom_path geom_point geom_histogram geom_abline geom_errorbar aes_string ylim xlab ylab scale_colour_manual aes_string geom_hline
#' @seealso \code{\link{predict.emfrail}}, \code{\link{summary.emfrail}}, \code{\link{plot_emfrail}}.
#' @return An object of class \code{ggplot}
#' @examples
#' mod_rec <- emfrail(bladder1, Surv(start, stop, status) ~ treatment + number + cluster(id))
#' summary(mod_rec)
#'
#' # cumulative hazard
#' # Note: this individual has number = 0, which does not exist in the data
#' ggplot_pred(mod_rec)
#'
#' # survival, although not very meaningful with recurrent events
#' \dontrun{
#' ggplot_pred(mod_rec, quantity = "survival")
#' }
#'
#' # For an individual with number == 2
#' ggplot_pred(mod_rec, newdata = data.frame(treatment = "placebo", number = 2))
#'
#' # hazard ratio between an individual with 0 and with 2 recurrences at baseline
#' # the marginal hazard ratio is "pulled" towards 1:
#'
#' ggplot_hr(mod_rec, newdata = data.frame(treatment = "placebo", number = c(0, 2)))
#'
#' \dontrun{
#' # hazard ratio with the stable distribution:
#' mod_rec_stab <- emfrail(bladder1,
#'                         Surv(start, stop, status) ~ treatment + number + cluster(id),
#'                         .distribution = emfrail_distribution(dist = "stable"))
#'
#' ggplot_hr(mod_rec_stab, newdata = data.frame(treatment = "placebo", number = c(0, 2)))
#'
#' # histogram of frailty estimates
#' gghist_frail(mod_rec_stab)
#'
#' # plot of the frailty estimates
#' ggplot_frail(mod_rec_stab)
#' }
NULL

#' \code{ggplot_pred} plots predicted cumulative hazard or survival curves, marginal and / or conditional, with or without confidence intervals.
#' @rdname ggplot_emfrail
#' @name ggplot_pred
#' @export
#'
ggplot_pred <- function(object,
                        lp = 0,
                        newdata = NULL,
                        quantity = "cumhaz", type = "both", conf_int = "adjusted", ...) {

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

        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "cumhaz")) +
          geom_step()

        if(conf_int == "adjusted") {
          plot1 <-
            plot1 +
            geom_path(aes_string(y = "cumhaz_l_a"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <-
            plot1 +
            geom_path(aes_string(y = "cumhaz_l"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r"), lty = 2)
        }
      }

      if(quantity == "survival") {
        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "survival")) +
          geom_step() +
          ylim(c(0,1))

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_l_a"), lty = 2) +
            geom_path(aes_string(y = "survival_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_l"), lty = 2) +
            geom_path(aes_string(y = "survival_r"), lty = 2)
        }
      }
    }

    if(type == "marginal") {
      if(quantity == "cumhaz") {

        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "cumhaz_m")) +
          geom_step()

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l_a"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r"), lty = 2)
        }
      }

      if(quantity == "survival") {
        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "survival_m")) +
          geom_step() +
          ylim(c(0,1)) +
          ylab("S(t)")

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l_a"), lty = 2) +
            geom_path(aes_string(y = "survival_m_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l"), lty = 2) +
            geom_path(aes_string(y = "survival_m_r"), lty = 2)
        }
      }


    }

    if(type == "both") {
      if(quantity == "cumhaz") {

        plot1 <- p1 %>%
          ggplot(aes_string(x = "time")) +
          geom_step(aes_string(y = "cumhaz", col = shQuote("1"))) +
          geom_step(aes_string(y = "cumhaz_m", col = shQuote("2"))) +
          scale_colour_manual(name = "type",
                              values=c("black", "red"),
                              labels = c("conditional", "marginal"))


        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_l_a", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r_a", col = shQuote("1")), lty = 2)

          # with(p1, lines(time, cumhaz_m_l_a, lty = 2, col = 2))
          # with(p1, lines(time, cumhaz_m_r_a, lty = 2, col = 2))
          # with(p1, lines(time, cumhaz_l_a, lty = 2))
          # with(p1, lines(time, cumhaz_r_a, lty = 2))
        }

        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_l", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r", col = shQuote("1")), lty = 2)

        }

      }

      if(quantity == "survival") {
        plot1 <- p1 %>%
          ggplot(aes_string(x = "time")) +
          geom_step(aes_string(y = "survival", col = shQuote("1"))) +
          geom_step(aes_string(y = "survival_m", col = shQuote("2"))) +
          ylim(c(0,1)) +
          ylab("S(t)") +
          scale_colour_manual(name = "type",
                              values=c("black", "red"),
                              labels = c("conditional", "marginal"))
        # with(p1, plot(time, survival_m,
        #               type = "s",
        #               #main = "Survival",
        #               ylab = "S(t)",
        #               xlab = "time",
        #               col = 2,
        #               ylim = c(0,1),
        #               ...))
        # with(p1, lines(time, survival,
        #                type = "s",
        #                col = 1,
        #                ...))
        #
        if(conf_int == "adjusted") {
          plot1 +
            geom_path(aes_string(y = "survival_m_l_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_m_r_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_l_a", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "survival_r_a", col = shQuote("1")), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 +
            geom_path(aes_string(y = "survival_m_l", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_m_r", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_l", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "survival_r", col = shQuote("1")), lty = 2)
        }

        # legend(x = "topright", legend = c("conditional", "marginal"), col = 1:2, lty = 1)
      }


    }
    plot1


}

#' \code{ggplot_hr} plots the estimated marginal and conditional hazard ratio between two units with different linear predictor values.
#' @rdname ggplot_emfrail
#' @name ggplot_hr
#' @export
#'
ggplot_hr <- function(object, lp, newdata = NULL, ...) {

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

  plot1 <- p1 %>%
    ggplot(aes_string(x = "time")) +
    geom_hline(aes_string("yintercept" = 1), col = "gray") +
    geom_step(aes_string(y = "hr_cond", col = shQuote("1"))) +
    geom_step(aes_string(y = "hr_mar", col = shQuote("2"))) +
    scale_colour_manual(name = "type",
                        values=c("black", "red"),
                        labels = c("conditional", "marginal")) +
    ylab("hazard ratio")

  plot1


}


#' \code{gghist_frail} plots a histogram of the estimated frailties.
#' @rdname ggplot_emfrail
#' @name gghist_frail
#'
#' @export
#'
#' @examples
#'
gghist_frail <- function(object, ...) {
  sobj <- summary.emfrail(object)
  plot1 <- sobj$frail %>%
    ggplot(aes_string(x = "z")) +
    geom_histogram() +
    xlab("frailties")

  plot1
}



#' \code{ggplot_frail} plots the ordered estimates of the frailties.
#' @rdname ggplot_emfrail
#' @name ggplot_frail
#'
#' @export
#'
#' @examples
#'
ggplot_frail <- function(object, ...) {
  sobj <- summary.emfrail(object)

  plot1 <- sobj$frail[order(sobj$frail$z),] %>%
    ggplot(aes_string(x = 1:length(sobj$frail$z), y = "z")) +
    geom_point(aes_string(id = "id"))

  if(sobj$est_dist$dist == "gamma")
    plot1 <- plot1 +
      geom_errorbar(aes_string(ymin = "lower_q", ymax = "upper_q", id = "id"))

  plot1
}


