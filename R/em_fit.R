em_fit <- function(logfrailtypar, dist, pvfm,
                   Y, Xmat, # id,  # this is some data stuff
                   atrisk, # a list with stuff that won't change in the EM
                   basehaz_line,  # need for log-likelihood
                   mcox = list(),
                   Cvec, lt = FALSE, Cvec_lt, # we need to start somewhere with the Cvec (E step comes first)
                   inner_control, se,
                   return_loglik = TRUE
) {


  pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  if (isTRUE(inner_control$verbose)) {
    print(paste0(#"dist=", pars$dist,
      "logfrailtypar= ", logfrailtypar,
      " / alpha=", pars$alpha,
      " / bbeta=", pars$bbeta))
  }

  if(logfrailtypar < -100) warning("theta virtually 0; try another starting value")

  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  # if the logfrailtypar is large (i.e. frailty variance 0) then just return the Cox likelihood
  if(logfrailtypar > inner_control$lower_tol) {
    #message("Frailty parameter very large, frailty variance close to 0")
    loglik <- mcox$loglik[length(mcox$loglik)]

    if(isTRUE(return_loglik)) {
      if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
      return(-loglik)
    }

  }


  loglik_old = -Inf
  ncycles <- 0



  convergence <- FALSE
  while(!isTRUE(convergence)) {

    if(isTRUE(inner_control$fast_fit)) {
      e_step_val <- fast_Estep(Cvec, Cvec_lt, atrisk$nev_id, alpha = pars$alpha, bbeta = pars$bbeta, pvfm = pvfm, dist = pars$dist)
    } else {
      e_step_val <- Estep(Cvec, Cvec_lt, atrisk$nev_id, alpha = pars$alpha, bbeta = pars$bbeta, pvfm = pvfm, dist = pars$dist)
    }

     logz <- log((e_step_val[,1] / e_step_val[,2])[atrisk$order_id])

     # the last part of this is so that the log-lik is comparable with that from coxph
    loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) + sum(e_step_val[,3]) +
      sum(Y[,3]) - sum((atrisk$nevent * log(atrisk$nevent))[atrisk$nevent > 0])

    # if this happens, then something is going very wrong
    if(loglik < loglik_old - inner_control$lik_tol)
      warning(paste0("likelihood decrease of ", loglik - loglik_old ))

    if((loglik - loglik_old) < inner_control$eps) break

    loglik_old <- loglik

    # print(paste0("loglik is ", loglik, " coef are ", paste0(mcox$coefficients, collapse = " ")))

    mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = NULL, offset = logz, init = NULL,
                                control = survival::coxph.control(), weights = NULL,
                                method = "breslow", rownames = NULL)

    # NOTE: this is what linear.predictors actually is:
    # exp(mcox$coefficients * (Xmat - mean(Xmat)) + logz)


    if(length(Xmat)==0) {
      lp <- mcox$linear.predictors
      g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))
    } else {
      lp <- mcox$linear.predictors + as.numeric(t(mcox$coefficients) %*% mcox$means)
      g_x <- t(mcox$coefficients %*% t(Xmat))
    }

    explp <- exp(lp)

    # Calculation of the baseline hazard - inspired from what survfit() does

    nrisk <- rev(cumsum(rev(rowsum(explp, Y[, ncol(Y) - 1]))))
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))


    nrisk <- nrisk - c(esum, 0,0)[atrisk$indx]
    haz <- atrisk$nevent/nrisk
    cumhaz <- cumsum(haz)

    basehaz_line <- haz[atrisk$time_to_stop]

    # cumulative hazard for [0, tstop]
    cumhaz_0_line <- cumhaz[atrisk$time_to_stop]
    # cumulative hazard for [0, tstart]
    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]

    cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)

    if(isTRUE(lt)) {
      Cvec_lt <- rowsum(x = cumhaz_tstart * exp(g_x), atrisk$order_id , reorder = FALSE)
    } else {
      Cvec_lt <- 0 * Cvec
    }

    Cvec <- rowsum(cumhaz_line * exp(g_x), atrisk$order_id, reorder = FALSE)

    ncycles <- ncycles + 1

    if(ncycles > inner_control$maxit) {
      warning(paste("did not converge in ", inner_control$maxit," iterations." ))
      break
    }
  }

  if(isTRUE(return_loglik)) {
    if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
    return(-loglik)
  }


  # From this point on, the standard errors & return object

  tev <- atrisk$time[haz > 0]
  haz_tev = haz[haz > 0]


  if(!isTRUE(se)) {
    if(length(Xmat) == 0) {
      Vcov <- matrix(NA, length(tev), length(tev))
    } else {
      Vcov <- matrix(NA, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))
    }


    res = list(loglik = loglik, # this we need
               tev = tev, # event time points
               haz = haz_tev, # the Breslow estimator for ech tev
               nev_id = atrisk$nev_id,
               Cvec = Cvec, #the Lambdatildei
               estep = e_step_val, # the E step object, just keep it like that.
               coef = mcox$coefficients, # the maximized coefficients. I need this.
               Vcov = Vcov)
    return(res)
  }

  # Standard error calculation


  nev_tp <- atrisk$nevent[atrisk$nevent!=0]

  z_elp = exp(lp)
  elp = exp(lp)  / exp(logz)

  # Building E[d2l/dx^2 | Z]

  if(length(Xmat)>0) {
    x <- lapply(apply(Xmat, 1, list), function(x) x[[1]])
    x_z_elp <- Map(function(a,b) a*b, x, z_elp)
    x_z_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp, cumhaz_line)
    x_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp / exp(logz), cumhaz_line)

    xx <- lapply(x, function(x) x %*% t(x) )
    xx_z_elp_H0 <- Map(function(a,b, c) a * b * c, xx, z_elp, cumhaz_line)
    m_d2l_dgdg <- Reduce("+", xx_z_elp_H0)

    m_d2l_dhdg <-
      do.call(rbind,
              lapply(lapply(
               lapply(tev, function(tk) which(Y[,1] < tk & tk <= Y[,2])),
               function(x) x_z_elp[x]),
             function(...) Reduce("+", ...))
      )

    # m_d2l_dhdg_old <- tev %>%
    #   lapply(function(tk) which(Y[,1] < tk & tk <= Y[,2])) %>%
    #   lapply(function(x) x_z_elp[x]) %>%
    #   lapply(function(...) Reduce("+", ...)) %>% # instead of sum because this could be a matrix
    #   do.call(rbind, .)
    #
    # all.equal(m_d2l_dhdg, m_d2l_dhdg_old)
  } else {
    m_d2l_dgdg <- NULL
    m_d2l_dhdg <- NULL
  }

  m_d2l_dhdh <- diag(nev_tp/haz_tev^2)

  # Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))
  #
  # if(!is.null(mcox$coefficients)) {
  #   Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- m_d2l_dgdg
  #   Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(m_d2l_dhdg)
  #   Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- m_d2l_dhdg
  # }
  #
  #
  # Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- m_d2l_dhdh

  # Building E[dl/dx (dl/dx)' | Z]

  if(isTRUE(inner_control$fast_fit)) {
      estep_again <- fast_Estep(Cvec,
                                Cvec_lt,
                                atrisk$nev_id,
                                alpha = pars$alpha,
                                bbeta = pars$bbeta,
                                pvfm = pvfm,
                                dist = pars$dist)
      z <- estep_again[,1] / estep_again[,2]
      zz <- estep_again[,4]
    } else {
      estep_plusone <- Estep(Cvec,
                             Cvec_lt,
                             atrisk$nev_id+1,
                             alpha = pars$alpha,
                             bbeta = pars$bbeta,
                             pvfm = pvfm,
                             dist = pars$dist)
      estep_again <- Estep(Cvec,
                           Cvec_lt,
                           atrisk$nev_id,
                           alpha = pars$alpha,
                           bbeta = pars$bbeta,
                           pvfm = pvfm,
                           dist = pars$dist)
      zz <- estep_plusone[,1] /estep_again[,2]
      z <- estep_again[,1] / estep_again[,2]
    }



  dl1_dh <- nev_tp / haz_tev

  tl_ord <- findInterval(Y[,1], tev)
  tr_ord <- findInterval(Y[,2], tev, left.open = FALSE, rightmost.closed = FALSE)

  dl2_dh <- try(inf_mat_match(
    tl_ord,
    tr_ord,
    z_elp,
    length(tev)
  ))

  # this is a list of data frames - for each individual - in each one a vector of length tev
  # each thing is the sum of elp at risk at that tev from each cluster
  elp_to_tev <-  lapply(
    split.data.frame(data.frame(elp,
                                y1 = findInterval(Y[,1], tev),
                                y2 = findInterval(Y[,2], tev, left.open = FALSE, rightmost.closed = FALSE)),
                     atrisk$order_id),
    function(dat) inf_mat_match(dat$y1, dat$y2, dat$elp, length(tev))
  )


  if(length(Xmat) > 0) {


    tmp1 <- rowsum(do.call(rbind, x_elp_H0), atrisk$order_id, reorder = FALSE) * sqrt(zz - z^2)
    cor_dg <- Reduce("+",lapply(split(tmp1, 1:nrow(tmp1)), function(x) x %*% t(x)))

    # I_gg_loss <- cor_dg

    cor_dg_dh <- t(Reduce("+",
           Map(function(a,b) a %*% t(b),
        elp_to_tev,
        split(
          tmp1 *  sqrt(zz - z^2),
          1:nrow(tmp1))
        )
        )
    )

    # cor_dg_dh <- elp_to_tev %>%
    #   Map(function(a, b) a %*% t(b),
    #          tapply(x_elp_H0, atrisk$order_id,
    #                 function(...) Reduce("+", ...))  ,.) %>%
    #   mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
    #   Reduce("+",.)

    I_gh_loss <-  cor_dg_dh

  } else {
    cor_dg_dh <- NULL
    cor_dg <- NULL
  }


  a <- Map(function(a,b) a * b,
      elp_to_tev,
      sqrt(zz - z^2))
  m <- matrix(0, length(a[[1]]),  length(a[[1]]))
  m[upper.tri(m, diag = TRUE)] <- sumxxt(a, length(a[[1]]))
  cor_dh <- m + t(m) - diag(diag(m))

  # cor_dh <- elp_to_tev %>%  # these are the c_ik without the z.
  #   lapply(function(x) x %*% t(x)) %>%
  #   mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
  #   Reduce("+",.)

  I_hh <- m_d2l_dhdh - cor_dh

  if(length(Xmat)>0) {
    I_gg <- m_d2l_dgdg - cor_dg
    I_hg <- m_d2l_dhdg - t(cor_dg_dh)

    Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))

    Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- I_gg

    Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- I_hh

    Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(I_hg)
    Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- I_hg

  } else Imat <- I_hh

  Vcov = try(solve(Imat), silent = TRUE)

  if(!isTRUE(return_loglik)) {
    res = list(loglik = loglik,
               tev = tev,
               haz = haz_tev,
               nev_id = atrisk$nev_id,
               Cvec = Cvec,
               frail = e_step_val[,1] / e_step_val[,2],
               coef = mcox$coefficients,
               Vcov = Vcov,
               fitted = g_x + logz,
               cumhaz_line = cumhaz_line)

    res
  }


}


