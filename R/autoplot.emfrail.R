#' @importFrom ggplot2 autoplot
#' @export
#' @rdname autoplot.emfrail
autoplot <- autoplot

#' Plots for emfrail objects using \code{ggplot2}
#' @importFrom ggplot2 autoplot ggplot geom_step geom_path geom_point geom_histogram geom_abline geom_errorbar aes_string ylim xlab ylab scale_colour_manual aes_string geom_hline
#' @param object \code{emfrail} object, typically result of \code{emfrail()}
#' @param type One (or more) of \code{hist} for a histogram of the estimated frailty values,
#'  \code{hr} for a plot of the conditional and marginal hazard ratio between two cases,
#'  \code{pred} for the predicted conditional and marginal cumulative hazard or survival for one case,
#'  \code{frail} for a plot of the ordered frailty estimates with confidence intervals, where available.
#' @param newdata A \code{data.frame} with values of the covariates. For \code{type == "hr"} the hazard ratio
#' between the first two rows of \code{newdata} is calculated. For \code{type == "pred"} the prediction
#' for the first row of \code{newdata} is calculated.
#' @param lp A numeric vector of values of the linear predictor, each corresponding to a case. For \code{type == "hr"} the hazard ratio
#' between the first two values of \code{lp} is calculated. For \code{type == "pred"} the prediction
#' for the first value of \code{lp} is calculated.
#' @param quantity For \code{type == "pred"} the predicted quantity; see \code{\link{predict.emfrail}}
#' @param marg_cond For \code{type == "pred"} the type of predicted quantity; see \code{\link{predict.emfrail}}
#' @param conf_int For \code{type == "pred"} the type of confidence intervals; see \code{\link{predict.emfrail}}
#' @param ... Further arguments to be passed on to `ggplot`
#'
#' @return A list of \code{ggplot2} objects corresponding to the required plots, or one \code{ggplot2} if only one plot is selected
#' @export autoplot.emfrail
#' @export
#' @seealso \code{\link{predict.emfrail}}, \code{\link{summary.emfrail}}, \code{\link{plot.emfrail}}.
#'
#' @examples
#' mod_rec <- emfrail(Surv(start, stop, status) ~ treatment + number + cluster(id), bladder1)
#' summary(mod_rec)
#' # cumulative hazard
#' # Note: this individual has number = 0, which does not exist in the data
#' autoplot(mod_rec, type = "pred", lp = 0)
#'
#' # survival, although not very meaningful with recurrent events
#' autoplot(mod_rec, type = "pred", lp = 0, quantity = "survival")
#'
#'
#' # For an individual with number == 2
#' autoplot(mod_rec, type = "pred", newdata = data.frame(treatment = "placebo", number = 2))
#'
#' # hazard ratio between an individual with 0 and with 2 recurrences at baseline
#' # the marginal hazard ratio is "pulled" towards 1:
#'
#' autoplot(mod_rec, type = "hr",
#' newdata = data.frame(treatment = c("placebo", "placebo"), number = c(0, 2)))
#'
#' # hazard ratio with the stable distribution:
#' \dontrun{
#' mod_rec_stab <- emfrail(Surv(start, stop, status) ~ treatment + number + cluster(id),
#'                         bladder1,
#'                         distribution = emfrail_distribution(dist = "stable"))
#'
#' autoplot(mod_rec_stab, type = "hr",
#' newdata = data.frame(treatment = c("placebo", "placebo"), number = c(0, 2)))
#'
#' # histogram of frailty estimates
#' autoplot(mod_rec_stab, type = "hist")
#' }

autoplot.emfrail <- function(object,
                             type = c("hist", "hr", "pred", "frail"),
                             newdata = NULL, lp = NULL,
                             quantity = "cumhaz",
                             marg_cond = "both",
                             conf_int = "adjusted", ...) {

  res <- vector("list", length(type))
  i <- 1

  if("hist" %in% type) {
    sobj <- summary.emfrail(object)
    plot1 <- sobj$frail %>%
      ggplot(aes_string(x = "z")) +
      geom_histogram() +
      xlab("frailties")
    res[[i]] <- plot1
    i <- i+1
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

    plot2 <- p1 %>%
      ggplot(aes_string(x = "time")) +
      geom_hline(aes_string("yintercept" = 1), col = "gray") +
      geom_step(aes_string(y = "hr_cond", col = shQuote("1"))) +
      geom_step(aes_string(y = "hr_mar", col = shQuote("2"))) +
      scale_colour_manual(name = "type",
                          values=c("black", "red"),
                          labels = c("conditional", "marginal")) +
      ylab("hazard ratio")

    res[[i]] <- plot2
    i <- i+1
  }

  if("pred"  %in% type) {


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

    p1 <- predict.emfrail(object,
                          lp = lp,
                          newdata = newdata,
                          quantity = quantity,
                          type = type_pred,
                          conf_int = conf_int)

    # determine columns to be plotted by name,
    # possibilities: quantity is just one of survival or cumhaz
    # they may have an _m or not at the end,
    # if(lp == 0) baseline <- TRUE else baseline <- FALSE

    if(marg_cond == "conditional") {
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

    if(marg_cond == "marginal") {
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

    if(marg_cond == "both") {
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
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_m_r_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_l_a", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "survival_r_a", col = shQuote("1")), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_m_r", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_l", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "survival_r", col = shQuote("1")), lty = 2)
        }

        # legend(x = "topright", legend = c("conditional", "marginal"), col = 1:2, lty = 1)
      }


    }

    res[[i]] <- plot1
    i <- i+1

  }

  if("frail" %in% type) {
    sobj <- summary.emfrail(object)

    plot1 <- sobj$frail[order(sobj$frail$z),] %>%
      ggplot(aes_string(x = 1:length(sobj$frail$z), y = "z")) +
      geom_point(aes_string(id = "id"))

    if(sobj$est_dist$dist == "gamma")
      plot1 <- plot1 +
      geom_errorbar(aes_string(ymin = "lower_q", ymax = "upper_q", id = "id"))

    res[[i]] <- plot1
  }

  if(length(res) == 1) res <- res[[1]]

  res
}
