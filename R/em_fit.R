# This does a gamma frailty EM for fixed frailtypar


em_fit <- function(logfrailtypar, dist, pvfm,
                   Y, Xmat, id,  # this is some data stuff
                   nev_id, newrisk,
                   basehaz_line,
                   mcox = list(),
                   explp, Cvec, lt = FALSE, Cvec_lt,
                   .control,
                   return_loglik = TRUE
) {

  # no events/time point, needed for the likelihood calculation
  nev_tp <- tapply(X = Y[,3], INDEX = Y[,2], sum)
  nev_tp <- nev_tp[nev_tp!=0]

  #print("hello im in em_fit")

  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  # if the logfrailtypar is large, i.e. frailtypar is large, i.e. fr. variance close to 0, then
  if(logfrailtypar > log(1/.control$zerotol)) {
    warning("Frailty parameter very large, frailty variance close to 0")
    loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
       sum(Y[,3]) - sum(nev_tp * log(nev_tp))

    if(isTRUE(return_loglik)) {
      if(isTRUE(.control$verbose)) print(paste("loglik = ",loglik))
      return(-loglik)
    }

  }

  .pars <- dist_to_pars(dist, logfrailtypar, pvfm)


  if (isTRUE(.control$verbose)) {
    print(paste0(#"dist=", .pars$dist,
                 "logfrailtypar= ", logfrailtypar,
                 " / alpha=", .pars$alpha,
                 " / bbeta=", .pars$bbeta))
  }



  loglik_old = -Inf
  ncycles <- 0


  convergence <- FALSE
  while(!isTRUE(convergence)) {

    if(dist=="gamma" & isTRUE(.control$fast_fit)) {
      e_step_val <- fast_Estep(Cvec, Cvec_lt, nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
    } else {
      e_step_val <- Estep(Cvec, Cvec_lt, nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
    }

    logz <- log(rep(e_step_val[,1] / e_step_val[,2],   rle(id)$lengths))
    # something only for the gamma:
    # logz <- log(rep((.pars$alpha + nev_id )/ (.pars$alpha + Cvec),   rle(id)$lengths))


    loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
     sum(e_step_val[,3]) + sum(Y[,3]) - sum(nev_tp * log(nev_tp))# +  sum(nev_id * lp_individual)

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
    } else {
      lp <- mcox$linear.predictors + t(mcox$coefficients) %*% mcox$means
    }

    hh <- getchz(Y = Y, newrisk = 1, explp = exp(lp) )

    # hh$tev
    # hh$haz_tev
    # plot(hh$time, hh$haz)
    # points(Y[,2], basehaz_line1, col = 2)

    # this is the baseline cumulative hazard for each line.
    # the idea is that this is only changes within an individual at the end of the line, and not in between. That's why it's correct.
    cumhaz_line <- sapply(X = apply(as.matrix(Y[,c(1,2)]), 1, as.list),
                          FUN = function(x)  sum(hh$haz_tev[x$start < hh$tev & hh$tev <= x$stop]))

    #
    basehaz_line <- hh$haz_tev[match(Y[,2], hh$tev)]

    if(length(Xmat)==0) {
      g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))
    } else {
      g_x <- t(mcox$coefficients %*% t(Xmat))
    }


    Cvec <- tapply(X = cumhaz_line * exp(g_x),
                   INDEX = id,
                   FUN = sum)

    # .distribution does not carry around.
    if(isTRUE(lt)) {
      cumhaz_lt_line <- sapply(X = apply(as.matrix(Y[,c(1,2)]), 1, as.list),
                               FUN = function(x) sum(hh$haz_tev[hh$tev <= x$start]))
      Cvec_lt <- tapply(X = cumhaz_lt_line * exp(g_x),
                        INDEX = id,
                        FUN = sum)
    } #else Cvec_lt <- 0 * Cvec


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

    m_d2l_dhdg <- hh$tev %>%
      lapply(function(tk) which(Y[,1] < tk & tk <= Y[,2])) %>%
      lapply(function(x) x_z_elp[x]) %>%
      lapply(function(...) Reduce("+", ...)) %>%
      do.call(rbind, .)
    # this is the most R piece of code I have ever written
  } else {
    m_d2l_dgdg <- NULL
    m_d2l_dhdg <- NULL
  }




  m_d2l_dhdh <- diag(nev_tp/hh$haz_tev^2)


#


 #  sqrt(diag(solve(I_full))) # this are the SE's, before adjusting for the frailty


  # now for the hell of the second one.

  estep_plusone <- Estep(Cvec, Cvec_lt, nev_id+1, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)
  estep_again <- Estep(Cvec, Cvec_lt, nev_id, alpha = .pars$alpha, bbeta = .pars$bbeta, pvfm = pvfm, dist = .pars$dist)

  zz <- estep_plusone[,1] /estep_again[,2]
  z <- e_step_val[,1] / e_step_val[,2]



  dl1_dh <- nev_tp / hh$haz_tev

  dl2_dh <- hh$tev %>%
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
      tapply(id, function(...) Reduce("+", ...)) %>%
      lapply(function(x) x %*% t(x)) %>%
      mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
      Reduce("+", .)

    I_gg_loss <- (dl1_dg  - dl2_dg) %*% t(dl1_dg - dl2_dg) + cor_dg

    cor_dg_dh <- split(data.frame(elp, y1 = Y[,1], y2 = Y[,2]), id) %>%
      lapply(function(dat) lapply(hh$tev, function(tk) sum(dat$elp[dat$y1 < tk & tk <= dat$y2]))) %>%
      lapply(function(...) do.call(c, ...)) %>%
      mapply(function(a, b) a %*% t(b), tapply(x_elp_H0, id, function(...) Reduce("+", ...))  ,. , SIMPLIFY = FALSE) %>%
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
  cor_dh <- split(data.frame(elp, y1 = Y[,1], y2 = Y[,2]), id) %>%
    lapply(function(dat) lapply(hh$tev, function(tk) sum(dat$elp[dat$y1 < tk & tk <= dat$y2]))) %>%
    lapply(function(...) do.call(c, ...)) %>%  # these are the c_ik without the z man.
    lapply(function(x) x %*% t(x)) %>%
    mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
    Reduce("+",.)

  I_hh_loss <- (dl1_dh - dl2_dh) %*% t(dl1_dh - dl2_dh) + cor_dh

  I_hh <- m_d2l_dhdh - I_hh_loss

  if(length(Xmat)>0) {
    I_gg <- m_d2l_dgdg - I_gg_loss
    I_hg <- m_d2l_dhdg - t(I_gh_loss)

    Imat <- matrix(0, ncol(Xmat) + length(hh$tev), ncol(Xmat) + length(hh$tev))

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
               haz = hh,
               z = exp(logz),
               Cvec = Cvec,
               estep = e_step_val,
               coef = mcox$coefficients,
               Vcov = Vcov,
               pvfm = pvfm)

    res
  }


}


