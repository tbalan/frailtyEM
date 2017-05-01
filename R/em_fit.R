em_fit <- function(logfrailtypar, dist, pvfm,
                   Y, Xmat, # id,  # this is some data stuff
                   atrisk, # a list with a shit load of things that will not change with the EM
                   basehaz_line,  # need for log-likelihood
                   mcox = list(),
                   Cvec, lt = FALSE, Cvec_lt, # we need to start somewhere with the Cvec (E step comes first)
                   inner_control, se,
                   return_loglik = TRUE
) {

  # no events/time point, needed for the likelihood calculation
  #nev_tp <- tapply(X = Y[,3], INDEX = Y[,2], sum)


  .pars <- dist_to_pars(dist, logfrailtypar, pvfm)
  if (isTRUE(inner_control$verbose)) {
    print(paste0(#"dist=", .pars$dist,
      "logfrailtypar= ", logfrailtypar,
      " / alpha=", .pars$alpha,
      " / bbeta=", .pars$bbeta))
  }

  if(logfrailtypar < -100) warning("theta virtually 0; try another starting value")


  #print("hello im in em_fit")

  if(length(Xmat)==0) {
    g_x <- matrix(rep(0, nrow(Y)),ncol = 1)
  } else {
    g_x <- t(mcox$coefficients %*% t(Xmat))
  }

  # if the logfrailtypar is large (i.e. frailty variance 0) then just return the Cox likelihood
  if(logfrailtypar > inner_control$lower_tol) {
    #message("Frailty parameter very large, frailty variance close to 0")
    loglik <- mcox$loglik[length(mcox$loglik)]
    # loglik <- sum((log(basehaz_line) + g_x)[Y[,3] == 1]) +
    #    sum(Y[,3]) - sum(nev_tp * log(nev_tp))

    if(isTRUE(return_loglik)) {
      if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
      return(-loglik)
    }

  }

  # some things for the hazard calculation


  loglik_old = -Inf
  ncycles <- 0



  convergence <- FALSE
  while(!isTRUE(convergence)) {

    if(isTRUE(inner_control$fast_fit)) {
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

    if(loglik < loglik_old - inner_control$lik_tol)
      warning(paste0("likelihood decrease of ", loglik - loglik_old ))

    if((loglik - loglik_old) < inner_control$eps) break

    loglik_old <- loglik

    # print(paste0("loglik is ", loglik, " coef are ", paste0(mcox$coefficients, collapse = " ")))


    mcox <- survival::agreg.fit(x = Xmat, y = Y, strata = NULL, offset = logz, init = NULL,
                                control = survival::coxph.control(), weights = NULL,
                                method = "breslow", rownames = NULL)


    #cc1 <- coxph(Surv(tstart, tstop, status) ~ x + offset(logz), dat1, method = "breslow")

    # NOTE: this ids what linear.predictors actually is:
    # exp(mcox$coefficients * (Xmat - mean(Xmat)) + logz)

    # How I calculate the cumulative hazard corresponding to each line in the data set...

    if(length(Xmat)==0) {
      lp <- mcox$linear.predictors
      g_x <- t(matrix(rep(0, length(mcox$linear.predictors)), nrow = 1))
    } else {
      lp <- mcox$linear.predictors + as.numeric(t(mcox$coefficients) %*% mcox$means)
      g_x <- t(mcox$coefficients %*% t(Xmat))
    }

    explp <- exp(lp)

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

    cumhaz_tstart <- c(0, cumhaz)[atrisk$indx2 + 1]
    cumhaz_line <- (cumhaz_0_line - cumhaz_tstart)  #* explp #/ newrisk

    if(isTRUE(lt)) {
      Cvec_lt <- rowsum(x = cumhaz_tstart * exp(g_x), atrisk$order_id )
    } else {
      Cvec_lt <- 0 * Cvec
    }

    Cvec <- rowsum( cumhaz_line * exp(g_x), atrisk$order_id)

    ncycles <- ncycles + 1
    if(ncycles > inner_control$maxit) {
      warning(paste("did not converge in ", inner_control$maxit," iterations." ))
      break
    }
  }

  if(isTRUE(return_loglik)) {
    if(isTRUE(inner_control$verbose)) print(paste("loglik = ",loglik))
    return(-loglik)
  }  # for when maximizing


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
               Cvec = Cvec, #the Lambdatildei, I don't think I need that. But maybe I do?
               estep = e_step_val, # the E step object, just keep it like that.
               coef = mcox$coefficients, # the maximized coefficients. I need this.
               Vcov = Vcov) # the Vcov matrix

    return(res)
  }

  # Standard error calculation
  # First part, second derivatives


  # this is a residual thing from the old way of calculating cumulative hazards



  nev_tp <- atrisk$nevent[atrisk$nevent!=0]

  z_elp = exp(lp)
  elp = exp(lp)  / exp(logz)

  # message("calculating Information Matrix...")
  # by line !
  if(length(Xmat)>0) {
    x <- lapply(apply(Xmat, 1, list), function(x) x[[1]])
    x_z_elp <- Map(function(a,b) a*b, x, z_elp)
    x_z_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp, cumhaz_line)
    x_elp_H0 <- Map(function(a,b,c) a*b*c, x, z_elp / exp(logz), cumhaz_line)

    xx <- lapply(x, function(x) x %*% t(x) )
    xx_z_elp_H0 <- Map(function(a,b, c) a * b * c, xx, z_elp, cumhaz_line)
    m_d2l_dgdg <- Reduce("+", xx_z_elp_H0)

    # if(any(m_d2l_dgdg<0)) warning("negative eigen in dgdg")

    m_d2l_dhdg <-
      do.call(rbind,
              lapply(lapply(
               lapply(tev, function(tk) which(Y[,1] < tk & tk <= Y[,2])),
               function(x) x_z_elp[x]),
             function(...) Reduce("+", ...))
      )
    # function()

    # m_d2l_dhdg_old <- tev %>%
    #   lapply(function(tk) which(Y[,1] < tk & tk <= Y[,2])) %>%
    #   lapply(function(x) x_z_elp[x]) %>%
    #   lapply(function(...) Reduce("+", ...)) %>% # instead of sum because this could be a matrix
    #   do.call(rbind, .)
    #
    #all.equal(m_d2l_dhdg, m_d2l_dhdg_old)
    # this is the most R piece of code I have ever written
  } else {
    m_d2l_dgdg <- NULL
    m_d2l_dhdg <- NULL
  }

  m_d2l_dhdh <- diag(nev_tp/haz_tev^2)
  # if(any(m_d2l_dhdh<0)) warning("negative eigen in dhdh")



  Imat <- matrix(0, ncol(Xmat) + length(tev), ncol(Xmat) + length(tev))

  Imat[1:length(mcox$coefficients), 1:length(mcox$coefficients)] <- m_d2l_dgdg

  Imat[(length(mcox$coefficients)+1):nrow(Imat), (length(mcox$coefficients)+1):nrow(Imat)] <- m_d2l_dhdh

  Imat[1:length(mcox$coefficients), (length(mcox$coefficients)+1):nrow(Imat) ] <- t(m_d2l_dhdg)
  Imat[(length(mcox$coefficients)+1):nrow(Imat), 1:length(mcox$coefficients) ] <- m_d2l_dhdg

  # if(any(eigen(Imat)$values<0)) warning("Imat naive negative eigenvalues")

 #  sqrt(diag(solve(I_full))) # this are the SE's, before adjusting for the frailty

  if(isTRUE(inner_control$fast_fit)) {
      estep_again <- fast_Estep(Cvec,
                                Cvec_lt,
                                atrisk$nev_id,
                                alpha = .pars$alpha,
                                bbeta = .pars$bbeta,
                                pvfm = pvfm,
                                dist = .pars$dist)
      z <- estep_again[,1] / estep_again[,2]
      zz <- estep_again[,4]
    } else {
      estep_plusone <- Estep(Cvec,
                             Cvec_lt,
                             atrisk$nev_id+1,
                             alpha = .pars$alpha,
                             bbeta = .pars$bbeta,
                             pvfm = pvfm,
                             dist = .pars$dist)
      estep_again <- Estep(Cvec,
                           Cvec_lt,
                           atrisk$nev_id,
                           alpha = .pars$alpha,
                           bbeta = .pars$bbeta,
                           pvfm = pvfm,
                           dist = .pars$dist)
      zz <- estep_plusone[,1] /estep_again[,2]
      z <- estep_again[,1] / estep_again[,2]
    }




  dl1_dh <- nev_tp / haz_tev



  # dl2_dh_old <- tev %>%
  #   lapply(function(tk) which(Y[,1] < tk & tk <= Y[,2])) %>%
  #   lapply(function(lin) sum(z_elp[lin])) %>% # no Reduce or something because this ain't a matrix
  #   do.call(c, .)


  tl_ord <- findInterval(Y[,1], tev)
  tr_ord <- findInterval(Y[,2], tev, left.open = FALSE, rightmost.closed = FALSE)

  dl2_dh <- tryCatch(inf_mat_match(
    tl_ord,
    tr_ord,
    z_elp, #this is one dimensional!
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


    #sum(delta_ij * x_ij)
    dl1_dg <- Reduce("+", mapply(function(a,b) a*b, Y[,3], x, SIMPLIFY = FALSE))
    # sum z x H
    dl2_dg <- Reduce("+", x_z_elp_H0)


    # this one to add; removes the part with (EZ)^2 and adds part with E(Z^2)
    # cor_dg <- x_elp_H0 %>%
    #   tapply(atrisk$order_id, function(...) Reduce("+", ...)) %>%
    #   lapply(function(x) x %*% t(x)) %>%
    #   mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
    #   Reduce("+", .)

    # x_elp_H0 could be just a data frame for the purpose here.

    tmp1 <- rowsum(do.call(rbind, x_elp_H0), atrisk$order_id) * sqrt(zz - z^2)
    cor_dg <- Reduce("+",lapply(split(tmp1, 1:nrow(tmp1)), function(x) x %*% t(x)))

    I_gg_loss <- (dl1_dg  - dl2_dg) %*% t(dl1_dg - dl2_dg) + cor_dg


    cor_dg_dh <- t(Reduce("+",
           Map(function(a,b) a %*% t(b),
        elp_to_tev,
        split(
          tmp1 *  sqrt(zz - z^2),
          1:nrow(tmp1))
        )
        )
    )


    #
    # cor_dg_dh <- elp_to_tev %>%
    #   Map(function(a, b) a %*% t(b),
    #          tapply(x_elp_H0, atrisk$order_id,
    #                 function(...) Reduce("+", ...))  ,.) %>%
    #   mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
    #   Reduce("+",.)

    I_gh_loss <- (dl1_dg  - dl2_dg) %*% t(dl1_dh - dl2_dh) + cor_dg_dh

  } else {
    I_gg_loss <- NULL
    I_gh_loss <- NULL
  }




  cor_dh <- Reduce("+", lapply(
    Map(function(a,b) a * b,
        elp_to_tev,
        sqrt(zz - z^2)),
    function(x) x %*% t(x)
  )
  )





  # cor_dh <- elp_to_tev %>%  # these are the c_ik without the z man.
  #   lapply(function(x) x %*% t(x)) %>%
  #   mapply(function(a,b) a * b, ., zz - z^2, SIMPLIFY = FALSE) %>%
  #   Reduce("+",.)

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
    res = list(loglik = loglik, # this we need
               tev = tev, # event time points
               haz = haz_tev, # the Breslow estimator for ech tev
               nev_id = atrisk$nev_id,
               Cvec = Cvec, #the Lambdatildei, I don't think I need that. But maybe I do?
               estep = e_step_val, # the E step object, just keep it like that.
               coef = mcox$coefficients, # the maximized coefficients. I need this.
               Vcov = Vcov) # the Vcov matrix

    res
  }


}


