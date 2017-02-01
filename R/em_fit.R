# This does a gamma frailty EM for fixed frailtypar


em_fit <- function(logfrailtypar, dist, pvfm,
                   Y, Xmat, # id,  # this is some data stuff
                   atrisk, # a list with a shit load of things that will not change with the EM
                   basehaz_line,  # need for log-likelihood
                   mcox = list(),
                   Cvec, lt = FALSE, Cvec_lt, # we need to start somewhere with the Cvec (E step comes first)
                   .control,
                   return_loglik = TRUE
) {

  # no events/time point, needed for the likelihood calculation
  #nev_tp <- tapply(X = Y[,3], INDEX = Y[,2], sum)

  .pars <- dist_to_pars(dist, logfrailtypar, pvfm)

  if(logfrailtypar < -100) stop("frailtypar virtually 0; try another starting value")

  if (isTRUE(.control$verbose)) {
    print(paste0(#"dist=", .pars$dist,
      "logfrailtypar= ", logfrailtypar,
      " / alpha=", .pars$alpha,
      " / bbeta=", .pars$bbeta))
  }
  #print("hello im in em_fit")

  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  # if the logfrailtypar is large, i.e. frailtypar is large, i.e. fr. variance close to 0, then
  if(!(dist %in% c("stable", "stable2")) &logfrailtypar > log(1/.control$zerotol)) {
    message("Frailty parameter very large, frailty variance close to 0")
    loglik <- mcox$loglik[length(mcox$loglik)]
    # loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
    #    sum(Y[,3]) - sum(nev_tp * log(nev_tp))

    if(isTRUE(return_loglik)) {
      if(isTRUE(.control$verbose)) print(paste("loglik = ",loglik))
      return(-loglik)
    }

  }

  # some things for the hazard calculation


  loglik_old = -Inf
  ncycles <- 0



  convergence <- FALSE
  while(!isTRUE(convergence)) {

    if(isTRUE(.control$fast_fit)) {
      e_step_val <- fast_Estep(Cvec, Cvec_lt, atrisk$nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
    } else {
      e_step_val <- Estep(Cvec, Cvec_lt, atrisk$nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
    }

    # a1 <- fast_Estep(Cvec + Cvec_lt, rep(0, length(Cvec)), nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
    # a2 <- fast_Estep(Cvec_lt, rep(0, length(Cvec)), rep(0, length(Cvec)), alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
    #
    # if(!isTRUE(all.equal(a1[,3] - a2[,3], e_step_val[,3]))) stop("sum ting wong")
    #
    # if(!isTRUE(all.equal(e_step_val[,1] / e_step_val[,2], a1[,1] / a1[,2]))) stop("e step not the same")

    # BAD idea:
    # rle(id)$lengths
    #
    # match(1:10, rep(1:10, each = 5))
    #
    # length(unique(id))
    # rep(1:278, rle(id)$lengths)
    #
     logz <- log((e_step_val[,1] / e_step_val[,2])[atrisk$order_id])
    # something only for the gamma:
    # logz <- log(rep((.pars$alpha + nev_id )/ (.pars$alpha + Cvec),   rle(id)$lengths))


    loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
     sum(e_step_val[,3]) + sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])# +  sum(nev_id * lp_individual)

    #
    # this is actually identical value:
    # loglik <- sum((log(basehaz_line) + t(mcox$coefficients %*% t(Xmat)))[Y[,3] == 1]) +
    #   sum(.pars$alpha * log(.pars$alpha) + lgamma(.pars$alpha + nev_id) - lgamma(.pars$alpha) -
    #         (.pars$alpha + nev_id) * log(.pars$alpha + Cvec)) +
    # sum(Y[,3]) - sum(nev_tp * log(nev_tp))# +  sum(nev_id * lp_individual)
    #
    if(loglik - loglik_old < 0) warning(paste0("likelihood decrease of ", loglik - loglik_old ))
    if((loglik - loglik_old) < .control$eps) break

    loglik_old <- loglik

    # print(paste0("loglik is ", loglik, " coef are ", paste0(mcox$coefficients, collapse = " ")))


    mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = NULL, offset = logz, init = NULL,
                                control = survival::coxph.control(), weights = NULL,
                                method = "breslow", rownames = NULL)


    #cc1 <- coxph(Surv(tstart, tstop, status) ~ x + offset(logz), dat1, method = "breslow")

    # NOTE: this ids what linear.predictors actually is:
    # exp(mcox$coefficients * (Xmat - mean(Xmat)) + logz)

    #cur <- survfit(cc1, newdata= data.frame(x = 0, logz = 0))

    #cur$cumhaz

    # How I calculate the cumulative hazard corresponding to each line in the data set...

    if(length(Xmat)==0) {
      lp <- mcox$linear.predictors
      g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))

    } else {
      lp <- mcox$linear.predictors + t(mcox$coefficients) %*% mcox$means
      g_x <- t(mcox$coefficients %*% t(Xmat))
    }

    explp <- exp(lp)



    # this is not really identical. Probably because shit is not scaled !
    # Cvec <- (nev_id - as.vector(rowsum(mcox$residuals, id))) / (e_step_val[,1] / e_step_val[,2])
# well funny enough this is fucking wrong

    # this is really wrong for no real reason?
    # cumhaz_line <- (Y[,3] - mcox$residuals) / exp(logz)# with covariates!



    # for the baseline hazard how the fuck is that gonna happen?
    # Idea: nrisk has the sum of elp who leave later at every tstop
    # esum has the sum of elp who enter at every tstart
    # indx groups which esum is right after each nrisk;
    # the difference between the two is the sum of elp really at risk at that time point.


    nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))


    nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]
    haz <- atrisk$nevent/nrisk # * newrisk
    cumhaz <- cumsum(haz)

    # baseline hazard for each tstop
    basehaz_line <- haz[atrisk$time_to_stop]
    cumhaz_0_line <- cumhaz[atrisk$time_to_stop]

    #cumhaz_tstop <- cumsum(haz)

    # for every tstop, this is the cumulative hazard at the following entry time.
    # indx2 <- findInterval(Y[,1], atrisk$time)


    # c(0, atrisk$time)[indx2[1:4]+1]

    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
    cumhaz_line <- cumhaz_0_line - cumhaz_tstart

# cumhaz_line[1:4]
# cumhaz_line_b[1:4]


    # finally, the cumulative hazard on each line is the difference
    # the trick used in emfrail() at the first place (with the residuals) does not work here
    # because agreg does something strange about scaling with offset.
#
#     hh <- getchz(Y, 1,  explp = exp(lp))
# #
#     hh
#     haz
#     cumhaz_line_b <- sapply(X = apply(as.matrix(Y[,c(1,2)]), 1, as.list),
#                           FUN = function(x)  sum(hh$haz_tev[x$start < hh$tev & hh$tev <= x$stop]))
#     #




    if(isTRUE(lt)) {
      Cvec_lt <- rowsum(x = cumhaz_tstart * exp(g_x), atrisk$order_id )
      # Cvec_lt <- tapply(X = cumhaz_tstart * exp(g_x),
      #                   INDEX = id,
      #                   FUN = sum)
    } else {
      Cvec_lt <- 0 * Cvec
    }

    Cvec <- rowsum( cumhaz_line * exp(g_x), atrisk$order_id)
#
#     Cvec_b <- tapply(X = cumhaz_line_b * exp(g_x), # * exp(g_x),
#                    INDEX = atrisk$order_id,
#                    FUN = sum)
#



    # .distribution does not carry around.



    ncycles <- ncycles + 1
    if(ncycles > .control$maxit) {
      warning(paste("did not converge in ", .control$maxit," iterations." ))
      break
    }

  }
  if(isTRUE(return_loglik)) {
    if(isTRUE(.control$verbose)) print(paste("loglik = ",loglik))
    return(-loglik)
  }  # for when maximizing

  # Standard error calculation
  # First part, second derivatives


  # this is a residual thing from the old way of calculating cumulative hazards

  tev <- atrisk$time[haz > 0]
  haz_tev = haz[haz > 0]
  #

  nev_tp <- atrisk$nevent[atrisk$nevent!=0]

  z_elp = exp(lp)
  elp = exp(lp)  / exp(logz)

  # message("calculating Information Matrix...")
  # by line !
  if(length(Xmat)>0) {
    x <- lapply(apply(Xmat, 1, list), function(x) x[[1]])
    x_z_elp <- mapply(function(a,b) a*b, x, z_elp, SIMPLIFY = FALSE)
    x_z_elp_H0 <- mapply(function(a,b,c) a*b*c, x, z_elp, cumhaz_line, SIMPLIFY = FALSE)
    x_elp_H0 <- mapply(function(a,b,c) a*b*c, x, z_elp / exp(logz), cumhaz_line, SIMPLIFY = FALSE)

    xx <- lapply(x, function(x) x %*% t(x) )
    xx_z_elp_H0 <- mapply(function(a,b, c) a * b * c, xx, z_elp, cumhaz_line, SIMPLIFY = FALSE)
    m_d2l_dgdg <- Reduce("+", xx_z_elp_H0)

    if(any(m_d2l_dgdg<0)) warning("negative eigen in dgdg")

    m_d2l_dhdg <- tev %>%
      lapply(function(tk) which(Y[,1] < tk & tk <= Y[,2])) %>%
      lapply(function(x) x_z_elp[x]) %>%
      lapply(function(...) Reduce("+", ...)) %>%
      do.call(rbind, .)
    # this is the most R piece of code I have ever written
  } else {
    m_d2l_dgdg <- NULL
    m_d2l_dhdg <- NULL
  }




  m_d2l_dhdh <- diag(nev_tp/haz_tev^2)
  if(any(m_d2l_dhdh<0)) warning("negative eigen in dhdh")



  Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))

  Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- m_d2l_dgdg

  Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- m_d2l_dhdh

  Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(m_d2l_dhdg)
  Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- m_d2l_dhdg


  if(any(eigen(Imat)$values<0)) warning("Imat naive negative eigenvalues")

 #  sqrt(diag(solve(I_full))) # this are the SE's, before adjusting for the frailty

  if(isTRUE(.control$fast_fit)) {
      estep_again <- fast_Estep(Cvec, Cvec_lt, atrisk$nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
      z <- estep_again[,1] / estep_again[,2]
      zz <- estep_again[,4]
    } else {
      estep_plusone <- Estep(Cvec, Cvec_lt, atrisk$nev_id+1, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
      estep_again <- Estep(Cvec, Cvec_lt, atrisk$nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
      zz <- estep_plusone[,1] /estep_again[,2]
      z <- estep_again[,1] / estep_again[,2]
    }




  dl1_dh <- nev_tp / haz_tev

  dl2_dh <- tev %>%
    lapply(function(tk) which(Y[,1] < tk & tk <= Y[,2])) %>%
    lapply(function(lin) sum(z_elp[lin])) %>%
    do.call(c, .)


  if(length(Xmat) > 0) {


    #sum(delta_ij * x_ij)
    dl1_dg <- Reduce("+", mapply(function(a,b) a*b, Y[,3], x, SIMPLIFY = FALSE))
    # sum z x H
    dl2_dg <- Reduce("+", x_z_elp_H0)


    # this one to add; removes the part with (EZ)^2 and adds part with E(Z^2)
    cor_dg <- x_elp_H0 %>%
      tapply(atrisk$order_id, function(...) Reduce("+", ...)) %>%
      lapply(function(x) x %*% t(x)) %>%
      mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
      Reduce("+", .)

    I_gg_loss <- (dl1_dg  - dl2_dg) %*% t(dl1_dg - dl2_dg) + cor_dg

    cor_dg_dh <- split(data.frame(elp, y1 = Y[,1], y2 = Y[,2]), atrisk$order_id) %>%
      lapply(function(dat) lapply(tev, function(tk) sum(dat$elp[dat$y1 < tk & tk <= dat$y2]))) %>%
      lapply(function(...) do.call(c, ...)) %>%
      mapply(function(a, b) a %*% t(b), tapply(x_elp_H0, atrisk$order_id, function(...) Reduce("+", ...))  ,. , SIMPLIFY = FALSE) %>%
      mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
      Reduce("+",.)

    I_gh_loss <- (dl1_dg  - dl2_dg) %*% t(dl1_dh - dl2_dh) + cor_dg_dh

  } else {
    I_gg_loss <- NULL
    I_gh_loss <- NULL
  }



  # also:
  # dl2_dh <- split(data.frame(elp, y1 = Y[,1], y2 = Y[,2]), id) %>%
  #   lapply(function(dat) lapply(hh$tev, function(tk) sum(dat$elp[dat$y1 < tk & tk <= dat$y2]))) %>%
  #   lapply(function(...) do.call(c, ...)) %>%  # these are the c_ik without the z man.
  #   mapply(function(a,b) a*b, ., z, SIMPLIFY = FALSE) %>%
  #   do.call(rbind, .) %>%
  #   apply(2,sum)

  # correction
  cor_dh <- split(data.frame(elp, y1 = Y[,1], y2 = Y[,2]), atrisk$order_id) %>%
    lapply(function(dat) lapply(tev, function(tk) sum(dat$elp[dat$y1 < tk & tk <= dat$y2]))) %>%
    lapply(function(...) do.call(c, ...)) %>%  # these are the c_ik without the z man.
    lapply(function(x) x %*% t(x)) %>%
    mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
    Reduce("+",.)

  I_hh_loss <- (dl1_dh - dl2_dh) %*% t(dl1_dh - dl2_dh) + cor_dh

  I_hh <- m_d2l_dhdh - I_hh_loss

  if(length(Xmat)>0) {
    I_gg <- m_d2l_dgdg - I_gg_loss
    I_hg <- m_d2l_dhdg - t(I_gh_loss)

    Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))

    Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- I_gg

    Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- I_hh

    Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(I_hg)
    Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- I_hg

  } else Imat <- I_hh


  #Imat %>% solve %>% diag %>% sqrt

  Vcov = solve(Imat)


  # with this one we will also need SE estimates and all the stuff
  if(!isTRUE(return_loglik)) {
    res = list(loglik = loglik,
               dist = dist,
               frailtypar = exp(logfrailtypar),
               haz = list(tev = tev, haz_tev = haz_tev),
               logz = logz,
               Cvec = Cvec,
               estep = e_step_val,
               coef = mcox$coefficients,
               Vcov = Vcov,
               pvfm = pvfm)

    res
  }


}


