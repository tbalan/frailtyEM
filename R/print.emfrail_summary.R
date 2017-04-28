#' @export
#' @keywords internal
print.emfrail_summary <- function(x, print_opts = list(coef = TRUE,
                                                       dist = TRUE,
                                                       fit = TRUE,
                                                       frailty = TRUE),
                                  ...) {

  obj <- x
  # cat("Summary of emfrail fit\n")
  cat("Call: \n")
  dput(attr(obj, "call"))
  cat("\n")

  if(isTRUE(attr(x, "print_opts")$coef)) {
    if(!is.null(obj$coefmat)) {
      cat("Regression coefficients:\n")
      printCoefmat(obj$coefmat, ...)
    }
  }

  if(isTRUE(attr(x, "print_opts")$dist)) {
    cat("Estimated distribution:", obj$est_dist$dist, "/ left truncation:", obj$est_dist$left_truncation,"\n")
    if(obj$est_dist$dist == "pvf") {
      cat("PVF m =", obj$est_dist$pvfm," ")
      if(obj$est_dist$pvfm == -0.5) cat("(Inverse Gaussian)")
      cat("\n")
    }
    cat("\n")
  }


  if(isTRUE(attr(x, "print_opts")$fit)) {
    cat("Fit summary:\n")
    if(!is.null(x$ca_test))
      cat("Commenges-Andersen test for heterogeneity: p-val ", format(x$ca_test[3], digits = 3), "\n")
    cat("(marginal) no-frailty Log-likelihood:", round(obj$loglik[1], digits = 3), "\n")
    cat("(marginal) Log-likelihood:", round(obj$loglik[2], digits = 3), "\n")
    cat("LRT: 1/2 * pchisq(", format(obj$loglik[3], digits = 3),"), p-val ",
        format(obj$loglik[4], digits = 3), "\n", sep = "")

    cat("\n")
  }


  if(isTRUE(attr(x, "print_opts")$frailty)) {
    cat("Frailty summary:\n")
    cat("theta = ",
        round(obj$theta[1], digits = 3),
        " (",
        round(obj$theta[2], digits = 2),
        ") / 95% CI: [",
        round(obj$theta[3], digits = 3),
        ", ",
        round(obj$theta[4], digits = 3),
        "]\n", sep = "")

    # gamma and pvf have this
    if(!is.null(obj$fr_var))
      cat("variance = ",
          round(obj$fr_var[1], digits = 3),
          # " (",
          # round(obj$fr_var[2], digits = 2),
          # ") / 95% CI: [",
          " / 95% CI: [",
          round(obj$fr_var[3], digits = 3),
          ", ",
          round(obj$fr_var[4], digits = 3),
          "]\n", sep = "")

    # gamma-specific
    if(!is.null(obj$gamma_pars))
      with(obj$gamma_pars, {
        cat("Kendall's tau: ",
            round(tau[[1]], digits = 3),
            # " (",
            # round(tau[[2]], digits = 2),
            # ") / 95% CI: [",
            " / 95% CI: [",
            round(tau[[3]], digits = 3),
            ", ",
            round(tau[[4]], digits = 3),
            "]\n", sep = "")

        cat("Median concordance: ",
            round(kappa[[1]], digits = 3),
            # " (",
            # round(kappa[[2]], digits = 2),
            # ") / 95% CI: [",
            " / 95% CI: [",
            round(kappa[[3]], digits = 3),
            ", ",
            round(kappa[[4]], digits = 3),
            "]\n", sep = "")

        cat("E[log Z]: ",
            round(e_log_z[[1]], digits = 3),
            # " (",
            # round(e_log_z[[2]], digits = 2),
            # ") / 95% CI: [",
            " / 95% CI: [",
            round(e_log_z[[3]], digits = 3),
            ", ",
            round(e_log_z[[4]], digits = 3),
            "]\n", sep = "")

        cat("Var[log Z]: ",
            round(var_log_z[[1]], digits = 3),
            # " (",
            # round(var_log_z[[2]], digits = 2),
            # ") / 95% CI: [",
            " / 95% CI: [",
            round(var_log_z[[3]], digits = 3),
            ", ",
            round(var_log_z[[4]], digits = 3),
            "]\n", sep = "")

      }

           )

    if(!is.null(obj$stable_pars))
    with(obj$stable_pars, {
      cat("Kendall's tau: ",
          round(tau[[1]], digits = 3),
          # " (",
          # round(tau[[2]], digits = 2),
          # ") / 95% CI: [",
          " / 95% CI: [",
          round(tau[[3]], digits = 3),
          ", ",
          round(tau[[4]], digits = 3),
          "]\n", sep = "")

      cat("Median concordance: ",
          round(kappa[[1]], digits = 3),
          # " (",
          # round(kappa[[2]], digits = 2),
          # ") / 95% CI: [",
          " / 95% CI: [",
          round(kappa[[3]], digits = 3),
          ", ",
          round(kappa[[4]], digits = 3),
          "]\n", sep = "")

      cat("E[log Z]: ",
          round(e_log_z[[1]], digits = 3),
          # " (",
          # round(e_log_z[[2]], digits = 2),
          # ") / 95% CI: [",
          " / 95% CI: [",
          round(e_log_z[[3]], digits = 3),
          ", ",
          round(e_log_z[[4]], digits = 3),
          "]\n", sep = "")

      cat("Var[log Z]: ",
          round(var_log_z[[1]], digits = 3),
          # " (",
          # round(var_log_z[[2]], digits = 2),
          # ") / 95% CI: [",
          " / 95% CI: [",
          round(var_log_z[[3]], digits = 3),
          ", ",
          round(var_log_z[[4]], digits = 3),
          "]\n", sep = "")

      cat("Attenuation factor: ",
          round(attenuation[[1]], digits = 3),
          # " (",
          # round(attenuation[[2]], digits = 2),
          # ") / 95% CI: [",
          " / 95% CI: [",
          round(attenuation[[3]], digits = 3),
          ", ",
          round(attenuation[[4]], digits = 3),
          "]\n", sep = "")

    }

    )

    if(!is.null(obj$pvf_pars))
      with(obj$pvf_pars, {
        if(obj$est_dist$pvfm < 0) {
          cat("Kendall's tau: ",
              round(tau[[1]], digits = 3),
              #" (",
              #round(tau[[2]], digits = 2),
              #") / 95% CI: [",
              "/ 95% CI: [",
              round(tau[[3]], digits = 3),
              ", ",
              round(tau[[4]], digits = 3),
              "]\n", sep = "")

          cat("Median concordance: ",
              round(kappa[[1]], digits = 3),
              # " (",
              # round(kappa[[2]], digits = 2),
              # ") / 95% CI: [",
              " / 95% CI: [",
              round(kappa[[3]], digits = 3),
              ", ",
              round(kappa[[4]], digits = 3),
              "]\n", sep = "")


        }

        if(obj$est_dist$pvfm == -1/2) {
          # E log
          # Var log z...
        }


        if(obj$est_dist$pvfm > 0)
            cat("Estimated mass at 0:",
                obj$pvf_pars[[1]],
                "\n")
        #
        # cat("E[log Z]: ",
        #     round(e_log_z[[1]], digits = 3),
        #     " (",
        #     round(e_log_z[[2]], digits = 2),
        #     ") / 95% CI: [",
        #     round(e_log_z[[3]], digits = 3),
        #     ", ",
        #     round(e_log_z[[4]], digits = 3),
        #     "]\n", sep = "")
        #
        # cat("Var[log Z]: ",
        #     round(var_log_z[[1]], digits = 3),
        #     " (",
        #     round(var_log_z[[2]], digits = 2),
        #     ") / 95% CI: [",
        #     round(var_log_z[[3]], digits = 3),
        #     ", ",
        #     round(var_log_z[[4]], digits = 3),
        #     "]\n", sep = "")

      }

      )

    if(!is.null(obj$cens_test)) {
      cat("\n")
      cat("Score test for dependent censoring: p-val", format(obj$cens_test[2], digits = 3))
    }

    # if(!is.null(obj$pvf_pars))
    #   cat("Estimated mass at 0:",
    #       obj$pvf_pars[[1]],
    #       "\n")

    if(isTRUE(x$lik_ci)) cat("Confidence intervals based on the likelihood function") else
      cat("Confidence intervals based on the delta method")

    cat("\n")
  }



}
