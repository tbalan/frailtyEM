#' @export
#' @keywords internal
print.emfrail <- function(x, ...) {

  obj <- x
  cat("Call: \n",
      "emfrail(",
      format(obj$.formula),
      ")\n",
      sep = "")

  names(obj$outer_m)[1] <- "theta"
  obj$outer_m[1] <- exp(obj$outer_m[1])

  names(obj$outer_m)[2] <- "loglik"
  obj$outer_m[2] <- (-1) * obj$outer_m[2]

  print(obj$outer_m)

  cat("\n")
  if(length(obj$inner_m$coef) > 0) {
    coefmat <- list(
      coef = obj$inner_m$coef,
      "exp(coef)" = exp(obj$inner_m$coef),
      "se(coef)" = sqrt(diag(obj$inner_m$Vcov)[seq_along(obj$inner_m$coef)]),
      "adjusted se" = sqrt(diag(obj$vcov_adj)[seq_along(obj$inner_m$coef)] ))

    coefmat$z <- coefmat$coef / coefmat$`se(coef)`
    coefmat$p <-  1 - pchisq(coefmat$z^2, df = 1)

    coefmat <- do.call(cbind, coefmat)

    printCoefmat(coefmat)
  }

  if(!is.null(obj$ca_test)) {
    cat("\n")
    cat("Commenges-Andersen test for heterogeneity: p-val", format(obj$ca_test[3], digits = 3))
  }

  invisible(obj)

}




