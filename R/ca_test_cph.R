#' Commenges-Andersen test for heterogeneity
#'
#' @param object A \code{coxph} object with a \code{cluster()} statement in the right-hand side of the formula.
#'
#' @return A named vector containing the test statistic, variance, and p-value
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom survival basehaz
#' @examples
#' mcox1 <- coxph(Surv(start, stop, status==1) ~ treatment + cluster(id),
#' bladder1, model = TRUE, x = TRUE)
#' ca_test(mcox1)

ca_test <- function(object) {

  # Check input
  if(!inherits(object, "coxph"))
    warning("input should be a coxph object")
  if(is.null(object$model))
    stop("object should be created with model=TRUE")
  if(length(grep("cluster", names(object$model)))==0)
    stop("formula should have a +cluster() statement")
  if(is.null(object$x))
    stop("object should be created with x=TRUE")


  # the real linear predictors
  lp <- object$linear.predictors +
    as.numeric(object$means %*% object$coefficients)
  elp <- exp(lp)

  Y <- object$y
  if(ncol(Y) == 2) Y <- cbind(0, Y)
  Y <- as.data.frame(cbind(Y[,1], Y[,2], Y[,3], elp = elp))

  if(length(grep("strata", names(object$model)))==0)
    Y$strata <- rep(1, nrow(Y)) else
      Y$strata <- object$model[,grep("strata", names(object$model))]

  Ysp <- split(Y, Y$strata)
  time <- sort(unique(Y[,2]))
  timesp <- lapply(Ysp, function(x) sort(unique(x[,2])))


  # this gives for every tstart (line variable) after which event time did it come
  # indx2 <- mapply(function(x,y) findInterval(x[,1], y),Ysp, timesp, SIMPLIFY = FALSE)
  indx2 <- lapply(Ysp, function(x) findInterval(x[,1], time))

  # This gives up t0 which time point does each line last
  # time_to_stop <- mapply(function(x,y) match(x[,2], y),
  #                        Ysp, timesp, SIMPLIFY = FALSE)
  time_to_stop <- lapply(Ysp, function(x) match(x[,2], time))


  # Time for S0_ht = for each strata, sum of Yij(t) elp(ij) within that strata
  S0_ht <- mapply(function(b,c,d) {
    lapply(seq_along(time), function(x) {
      b$elp[c < x & x <= d] %>%
        Reduce("+", .)
    })
  }, Ysp, indx2, time_to_stop, SIMPLIFY = FALSE) %>%
    lapply(function(z) do.call(c, z))


  # Time for pij(s) = Y_{ij}(s) exp(beta' x_{ij}) / S0(s)

  pij_ht_rowh <- mapply(function(b,c,d) {
    lapply(seq_along(time), function(x) {
      b$elp * as.numeric(c < x & x <= d)
    })
  }, Ysp, indx2, time_to_stop, SIMPLIFY = FALSE) %>%
    mapply(function(a,b) {
      mapply(function(c,d) c/d, a, b, SIMPLIFY = FALSE)
    }, ., S0_ht, SIMPLIFY = FALSE)

  # need pij_ht
  # pij (strata)(all_times)(alllines_within_strata)
  id <- object$model[,grep("cluster", names(object$model))]

  order_id <- split(findInterval(id, unique(id)), Y$strata)
  # order_id (strata)(alllines_within_strata)

  pij_hrowh_t <- pij_ht_rowh %>%
    lapply(function(x) do.call(cbind, x)) %>%
    lapply(function(x) split(x, 1:nrow(x))) %>%
    lapply(function(x) do.call(rbind, x))

  # This is each for each strata, a matrix with each row a cluster that exists in that strata
  # and each column one time point
  pi_ht <- mapply(function(a,b) {
    rowsum(a, b, reorder = FALSE)
  }, pij_hrowh_t, order_id, SIMPLIFY = FALSE)


  # Now to calculate the T statistic
  M <- split(object$residuals, Y$strata)
  #(strata)
  Mi_strata <- mapply(function(a,b) {
    rowsum(a, b)
  }, M, order_id, SIMPLIFY = FALSE)

  Mi_1 <- Mi_strata %>%
    lapply(as.data.frame) %>%
    lapply(tibble::rownames_to_column) %>%
    do.call(rbind, .)

  Mi <- as.numeric(rowsum(Mi_1$V1, Mi_1$rowname))
    # with(rowsum(V1, rowname)) %>%
    # as.numeric()  # already by cluster


  death <- (Y[, 3] == 1)
  deathsp <- split(death, Y$strata)

  # number of events per time point within each strata
  nevent <- lapply(Ysp, function(x) {
    lapply(time, function(tp) x[,2] == tp & x[,3] ==1)
  }) %>%
    lapply(function(x) {
      do.call(c, lapply(x, sum))
    })

  thirdterm <- pi_ht %>%
    lapply(function(x) x^2) %>%
    lapply(function(x) apply(x, 2, sum)) %>%
    mapply(function(a,b) a * b, ., nevent, SIMPLIFY = FALSE) %>%
    lapply(sum) %>%
    do.call("+", .)

  T_stat <- sum(Mi^2) - sum(death) + thirdterm

  # Now the variance of T_stat

  # without strata: Hi(S) = 2 {Mi(s-) - sum_l=1^n Ml(s-) pl(s) - pi(s) + sum_l=1^n pl^2(s)}
  # with strata:
  # Hih(s) = 2{Mi(s-) - sum_l=1^n Ml(s-) plh(s) - pih(s) + sum_l=1^n plh^2(s)}

  # Martingale path for each line

  # Calculate baseline cumulative hazard

  BH <- basehaz(object, centered = FALSE)
  if(is.null(BH$strata)) BH$strata <- rep(1, nrow(BH))
  BHsp <- split(BH, BH$strata)



  alltp <- sort(unique(BH$time))

  # pos_leftsp <- mapply(function(a,b) findInterval(a[,1], b$time), Ysp, BHsp, SIMPLIFY = FALSE)
  # pos_rightsp <- mapply(function(a,b) match(a[,2], b$time), Ysp, BHsp, SIMPLIFY = FALSE)
  pos_leftsp <- lapply(Ysp, function(x) findInterval(x[,1], alltp))
  pos_rightsp <- lapply(Ysp, function(x) match(x[,2], alltp))

  cumhaz <- lapply(BHsp, function(x) x$hazard)

  # Idea is to have this at ALL event time points, not only those in the strata
  cumhazsp <- lapply(BHsp, function(bh) {
    approx(x = bh$time, y = bh$hazard, xout = time, method = "constant", yleft = 0)$y
  })

  # First the path of the - cumulative hazard for each row
  Lij_path <- mapply(function(x,y,cumhaz) {
    mapply(function(pos_left, pos_right) {
      if(pos_left == 0) {
        cumhaz[pos_right:length(cumhaz)] <- cumhaz[pos_right] # cumulative hazard doesn't increase after tstop
        return(-cumhaz)
      }

      cumhaz[(pos_left+1):length(cumhaz)] <- cumhaz[(pos_left+1):length(cumhaz)] - cumhaz[pos_left] # substract the hazard before tstart
      cumhaz[1:(pos_left+1)] <- 0 # make it 0 before tstop
      cumhaz[pos_right:length(cumhaz)] <- cumhaz[pos_right] # cumulative hazard doesn't increase after tstop

      return(-cumhaz) # I THINK
    }, x, y, SIMPLIFY = FALSE)
  }, pos_leftsp, pos_rightsp, cumhazsp, SIMPLIFY = FALSE)

  Lij_elp_path <- mapply(function(a,b) {
    mapply(function(x,y) x * y, a, b$elp, SIMPLIFY = FALSE)
  }, Lij_path, Ysp, SIMPLIFY = FALSE)

  Mij_path <- mapply(function(x,y,z) {
    mapply(function(a, pos, delta) {
      a[pos:length(a)] <- delta + a[pos]
      a
    }, x,y,z, SIMPLIFY = FALSE)

  }, Lij_elp_path, time_to_stop, deathsp, SIMPLIFY = FALSE)

  # Now to get Mi so sum these according to clusters
  # Each row is th Mi path of a cluster; the cluster is the name of the row

  Mit <- do.call(c, Mij_path) %>%
    do.call(rbind, .) %>%
    rowsum(., group = do.call(c, order_id), reorder = FALSE)

  # Elements of Hih(t)

  Mit_all <-
    Mit[order(as.numeric(rownames(Mit))),]
  # clusters by time points

  # pi_ht is clusters by time points
  # but I walso clusters that are missing within each strata
  # and add them there with zeros.
  pi_ht_all <- pi_ht %>%
    lapply(function(x) {
      rnames <- rownames(x)
      allnames <- as.character(unique(do.call(c, order_id)))
      missnames <- which(!(allnames %in% rnames))
      missmat <- matrix(0, nrow = length(missnames),
                        ncol = ncol(x),
                        dimnames = allnames[missnames])
      wholemat <- rbind(x, missmat)
      wholemat[order(as.numeric(rownames(wholemat))),]
    })


  # First element

  # terms 1 and 3
  term13 <- lapply(pi_ht_all, function(x) Mit_all - x)

  term2 <- lapply(pi_ht_all, function(x) x * Mit_all) %>%
    lapply(function(x) apply(x, 2, sum))

  term4 <- lapply(pi_ht_all, function(x) x^2) %>%
    lapply(function(x) apply(x, 2, sum))

  term24 <- mapply(function(a,b) b - a, term2, term4, SIMPLIFY = FALSE)

  H_hi_t <- mapply(function(a,b) {
    2 * (t(t(a) + b))
  }, term13, term24, SIMPLIFY = FALSE)

  # Ihat

  Ihat <- H_hi_t %>%
    lapply(function(x) x^2) %>%
    mapply(function(a,b) a * b,
           ., pi_ht_all, SIMPLIFY = FALSE) %>%
    mapply(function(a,b) {
      t(t(a) * b)
    }, ., nevent, SIMPLIFY = FALSE) %>%
    lapply(function(x) apply(x, 1, sum)) %>%
    lapply(sum) %>%
    do.call("+", .)

  # Jhat

  X <- object$model

  xsp <- split(as.data.frame(object$x), Y$strata)



  # lapply(function(x) {
  #   rnames <- rownames(x)
  #   allnames <- as.character(unique(id))
  #   missnames <- which(!(allnames %in% rnames))
  #   missmat <- matrix(0, nrow = length(missnames), ncol = ncol(x), dimnames = list(allnames[missnames]))
  #   wholemat <- rbind(x, missmat)
  #   wholemat[order(as.numeric(rownames(wholemat))),]
  # })


  zipit <- mapply(function(a,b) {
    lapply(a, function(vec) {
      vec * b
    })
  }, xsp, pij_hrowh_t, SIMPLIFY = FALSE) %>%
    mapply(function(a,b) {
      lapply(a, function(x) rowsum(x, b, reorder = FALSE))
    }, ., order_id, SIMPLIFY = FALSE)

  # zipit is a list(strata)(covariate)(clusters within strata X total times)
  # Add the missing cluster from all the matrices in zipit

  zipit_all <- lapply(zipit, function(y) {
    lapply(y, function(x) {
      rnames <- rownames(x)
      allnames <- unique(as.character(do.call(c, order_id)))
      missnames <- which(!(allnames %in% rnames))
      missmat <- matrix(0, nrow = length(missnames), ncol = ncol(x), dimnames = list(allnames[missnames]))
      wholemat <- rbind(x, missmat)
      wholemat[order(as.numeric(rownames(wholemat))),]
    })
  })

  J <- mapply(function(a,b) {
    lapply(b, function(x) {
      sum(apply(a * x, 1, sum))
    })
  }, H_hi_t, zipit_all, SIMPLIFY = FALSE) %>%
    lapply(function(x) do.call(c, x)) %>%
    do.call(rbind, .) %>%
    apply(., 2, sum)

  denominator <-  Ihat - t(J) %*% object$naive.var %*% J

  c(tstat = T_stat, var = denominator, pval = pchisq(T_stat^2 / denominator, 1, lower.tail = FALSE))

}
