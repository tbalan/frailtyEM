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
  invisible(obj)

}




