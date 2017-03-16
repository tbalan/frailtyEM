#' ca_test_fit
#'
#' @param mcox An object as returned by  \code{agreg.fit}
#' @param X The \code{X} matrix of covariates as obtained from \code{model.matrix} without the cluster or intercept
#' @param atrisk A list with some fields from which the at-risk indicator can be deduced, as calculated in \code{emfrail()}
#' @param exp_g_x A vector of exponentiated linear predictor for each row of the data
#' @param cumhaz An estimate of the cumulative hazard at each time point in the data.
#'
#' @return A list with the test statistic, variance, and p-value
#'
#' @details This is my implementation of Commenges & Andersen (1995) test for heterogeneity, with a few adjustments to make it work with recurrent events data (?).
#' It could be made much faster and more efficient, but since it's not such an essential part, I'll let someone else do it.
#' @keywords internal
#'
if(getRversion() >= "2.15.1")  utils::globalVariables(".")
ca_test_fit <- function(mcox, X, atrisk, exp_g_x, cumhaz) {

  # numerator

  S0_t <-
    lapply(seq_along(atrisk$time), function(x)
      exp_g_x[atrisk$indx2 < x & x <= atrisk$time_to_stop ]) %>%
    lapply(sum) %>%
    do.call(c, .)
  #
  #
  # S0 <- (nrisk/newrisk)[atrisk$time_to_stop]
  #
  #
  #

  # the basic ingredients
  # pij_t is a list of length (tev) where at each event time point we have the sum of
  # elp of people at risk at that time point
  # all divided by S0_t
  pij_t <-
    lapply(seq_along(atrisk$time), function(x)
      exp_g_x * as.numeric(atrisk$indx2 < x & x <= atrisk$time_to_stop )) %>%
    lapply(as.numeric) %>%
    mapply(function(a,b) a/ b, ., S0_t, SIMPLIFY = FALSE)

  pi_t <- lapply(pij_t, function(x) rowsum(x, atrisk$order_id)) %>%
    lapply(as.numeric)

  T_stat <- sum(rowsum(mcox$residuals , atrisk$order_id)  ^2) -
    sum(atrisk$death) +
  sum(pi_t %>%
    lapply(function(x) sum(x^2)) %>%
    do.call(c, .) *
    atrisk$nevent )


  # this is the martingale path of each line

  mt_ij <- mapply(function(pos_left, pos_right) {

    if(pos_left == 0) {
      cumhaz[pos_right:length(cumhaz)] <- cumhaz[pos_right]
    return(-cumhaz)}


    cumhaz[(pos_left+1):length(cumhaz)] <- cumhaz[(pos_left+1):length(cumhaz)] - cumhaz[pos_left]
    cumhaz[1:(pos_left+1)] <- 0
    cumhaz[pos_right:length(cumhaz)] <- cumhaz[pos_right]

    cumhaz
    # stop it after pos_right

  }, atrisk$indx2, atrisk$time_to_stop, SIMPLIFY = FALSE) %>%
    mapply(function(a,b) a*b, ., exp_g_x, SIMPLIFY = FALSE) %>%
    mapply(function(a, pos, delta) {
      a[pos:length(a)] <- delta + a[pos]
      a
    }, ., atrisk$time_to_stop, atrisk$death, SIMPLIFY  = FALSE)



  # mt_ij <- lapply(atrisk$time_to_stop, function(pos) {
  #   cumhaz[pos:length(cumhaz)] <- cumhaz[pos]
  #   -cumhaz
  # }) %>%
  #   mapply(function(a,b) a*b, ., exp_g_x, SIMPLIFY = FALSE) %>%
  #   mapply(function(a, pos, delta) {
  #     a[pos:length(a)] <- delta + a[pos]
  #     a
  #   }, ., atrisk$time_to_stop, atrisk$death, SIMPLIFY  = FALSE)

  # mt_ij replaces the first time with 0 and kicks out the last one...

  # mt_ij %>%
  #   lapply(function(x) cbind(atrisk$time, atrisk$nevent,  x))

  m0t_ij <- #mt_ij
    lapply(mt_ij, function(x) c(0, x[-length(x)]))

  mij_t <- m0t_ij %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    as.list()

  mi_t <- lapply(mij_t, function(x) rowsum(x, atrisk$order_id)) %>%
    lapply(as.numeric)


  mp_t <- mapply(function(a,b) a * b, mi_t, pi_t, SIMPLIFY = FALSE) %>%
    lapply(sum)

  pi_t

  pp_t <- pi_t %>%
    lapply(function(x) x^2)  %>%
    lapply(sum)


  qi_t <- mapply(function(a,b,c,d) 2 * (a - b - c + d), mi_t, mp_t, pi_t, pp_t, SIMPLIFY = FALSE)

  mi_t

  # Main shit in V
  V1 <- qi_t %>%
    lapply(function(x) x^2) %>%
    mapply(function(a,b,c) a * b * c, ., pi_t, atrisk$nevent, SIMPLIFY = FALSE) %>%
    do.call(sum, .)

  # second part in V
  theta2i_t <- pij_t %>%
    lapply(function(x) x * X) %>%
    lapply(as.data.frame) %>%
    lapply(function(x) rowsum.data.frame(x, atrisk$order_id))


  theta_h <- atrisk$nevent %>%
    mapply(function(a,b,c) a * b * c, ., qi_t, theta2i_t, SIMPLIFY = FALSE) %>%
    do.call(rbind, .) %>%
    apply(., 2, sum)

  V2 <- t(theta_h) %*% mcox$var %*% theta_h

  V <- V1 + V2

  c(tstat = T_stat, var = V, pval = pchisq(T_stat^2 / V, 1, lower.tail = FALSE))
}


#
#
# ca_test <- function(.data, .formula) {
#
#   Call <- match.call()
#
#
#   if(missing(.formula)  | missing(.data)) stop("Missing arguments")
#   # if(missing(.distribution)) {
#   #   .distribution <-  emfrail_distribution()
#   #   print("default distribution")
#   # }
#   # if(missing(.control)) .control <- emfrail_control()
#
#   cluster <- function(x) x
#   mf <- model.frame(.formula, .data)
#
#   # Identify the cluster and the ID column
#   pos_cluster <- grep("cluster", names(mf))
#   if(length(pos_cluster) != 1) stop("misspecified or non-specified cluster")
#   id <- mf[[pos_cluster]]
#
#
#   Y <- mf[[1]]
#   if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
#   if(ncol(Y) != 3) {
#     # warning("Y not in (tstart, tstop) format; taking tstart = 0")
#     # Y <- cbind(rep(0, nrow(Y)), Y)
#     # attr(Y, "dimnames") <- list(NULL, c("start", "stop", "status"))
#     # attr(Y, "type") <- "counting"
#     Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
#   }
#   if(attr(Y, "type") != "counting") stop("use Surv(tstart, tstop, status)")
#
#
#   # get the model matrix
#   X1 <- model.matrix(.formula, .data)
#   # this is necessary because when factors have more levels, pos_cluster doesn't correspond any more
#   pos_cluster_X1 <- grep("cluster", colnames(X1))
#   X <- X1[,-c(1, pos_cluster_X1), drop=FALSE]
#   # note: X has no attributes, in coxph it does.
#
#
#   # some stuff for creating the C vector, is used all along.
#   # mcox also works with empty matrices, but also with NULL as x.
#   mcox <- survival::agreg.fit(x = X, y = Y, strata = NULL, offset = NULL, init = NULL,
#                               control = survival::coxph.control(),
#                               weights = NULL, method = "breslow", rownames = NULL)
#
#   # the "baseline" case // this will stay constant
#
#   if(length(X) == 0) {
#     newrisk <- 1
#     exp_g_x <- matrix(rep(1, length(mcox$linear.predictors)), nrow = 1)
#     g <- 0
#     g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))
#
#   } else {
#     x2 <- matrix(rep(0, ncol(X)), nrow = 1, dimnames = list(123, dimnames(X)[[2]]))
#     x2 <- (scale(x2, center = mcox$means, scale = FALSE))
#     newrisk <- exp(c(x2 %*% mcox$coefficients) + 0)
#     exp_g_x <- exp(mcox$coefficients %*% t(X))
#     g <- mcox$coefficients
#     g_x <- t(mcox$coefficients %*% t(X))
#
#   }
#
#   explp <- exp(mcox$linear.predictors) # these are with centered covariates
#
#
#   nev_id <- rowsum(Y[,3], id)
#
#   nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
#   esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
#
#   # the stuff that won't change
#   death <- (Y[, ncol(Y)] == 1)
#   nevent <- as.vector(rowsum(1 * death, Y[, ncol(Y) - 1])) # per time point
#   time <- sort(unique(Y[,2])) # unique tstops
#
#   # this gives the next entry time for each unique tstop (not only event)
#   etime <- c(0, sort(unique(Y[, 1])),  max(Y[, 1]) + min(diff(time)))
#   indx <- findInterval(time, etime, left.open = TRUE) # left.open  = TRUE is very important
#
#   # this gives for every tstart (line variable) after which event time did it come
#   # indx2 <- findInterval(Y[,1], time, left.open = FALSE, rightmost.closed = TRUE)
#   indx2 <- findInterval(Y[,1], time)
#
#   time_to_stop <- match(Y[,2], time)
#   order_id <- findInterval(id, unique(id))
#
#   atrisk <- list(death = death, nevent = nevent, nev_id = nev_id,
#                  order_id = order_id, time = time, indx = indx, indx2 = indx2,
#                  time_to_stop = time_to_stop)
#
#   nrisk <- nrisk - c(esum, 0,0)[indx]
#
#   # esum %>% length
#   # indx %>% length
#   #
#   # esum
#   # indx
#   # # last indx is 548
#   # esum
#
#   haz <- nevent/nrisk * newrisk
#
#
#   basehaz_line <- haz[atrisk$time_to_stop]
#
#   cumhaz <- cumsum(haz)
# }
