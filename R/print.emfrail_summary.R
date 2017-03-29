#' @export
#' @keywords internal
print.emfrail_summary <- function(x,...) {

  obj <- x
  cat("Summary of emfrail fit\n")

  if(!is.null(obj$coefmat)) {
    cat("Regression coefficients:\n")
    printCoefmat(obj$coefmat, ...)
  }

  cat("Estimated distribution:", obj$est_dist$dist, "/ left truncation:", obj$est_dist$left_truncation,"\n")
  if(obj$est_dist$dist == "pvf") {
    cat("PVF m =", obj$est_dist$pvfm," ")
    if(obj$est_dist$pvfm == -0.5) cat("(Inverse Gaussian)")
    cat("\n")
  }
  cat("\n")

  cat("Fit summary:\n")
  if(!is.null(x$ca_test))
    cat("Commenges-Andersen test for heterogeneity: p-val ", format(x$ca_test[3], digits = 3), "\n")
  cat("(marginal) no-frailty Log-likelihood:", round(obj$loglik[1], digits = 3), "\n")
  cat("(marginal) Log-likelihood:", round(obj$loglik[2], digits = 3), "\n")
  cat("LRT: 1/2 * pchisq(", format(obj$loglik[3], digits = 3),"), p-val ",
      format(obj$loglik[4], digits = 3), "\n", sep = "")

  cat("\n")

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
        " (",
        round(obj$fr_var[2], digits = 2),
        ") / 95% CI: [",
        round(obj$fr_var[3], digits = 3),
        ", ",
        round(obj$fr_var[4], digits = 3),
        "]\n", sep = "")

  # gamma-specific
  if(!is.null(obj$gamma_pars))
  cat("Kendall's tau: ",
      round(obj$gamma_pars[1], digits = 3),
      " (",
      round(obj$gamma_pars[2], digits = 2),
      ") / 95% CI: [",
      round(obj$gamma_pars[3], digits = 3),
      ", ",
      round(obj$gamma_pars[4], digits = 3),
      "]\n", sep = "")
    # cat("Kendall's tau:",
    #     obj$gamma_pars[[1]],
    #     "\n")

  # pvf-specific
  if(!is.null(obj$pvf_pars))
    cat("Estimated mass at 0:",
        obj$pvf_pars[[1]],
        "\n")

  # stable-specific
  if(!is.null(obj$stable_pars)) {
    cat("Kendall's tau: ",
        round(obj$stable_pars[1], digits = 3),
        " (",
        round(obj$stable_pars[2], digits = 2),
        ") / 95% CI: [",
        round(obj$stable_pars[3], digits = 3),
        ", ",
        round(obj$stable_pars[4], digits = 3),
        "]\n", sep = "")
    cat("Attenuation factor: ",
        round(obj$stable_pars[5], digits = 2),
        " / Var[log(Z)] = ",
        round(obj$stable_pars[6], digits = 3),
        "\n",
        sep = "")
  }


  cat("\n")


}
