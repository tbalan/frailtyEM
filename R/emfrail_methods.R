#' @export
#' @keywords internal
coef.emfrail <- function(object, ...) {
  object$coefficients
}

#' @export
#' @keywords internal
vcov.emfrail <- function(object, type = c("regular", "adjusted"), ...) {

  if("adjusted" %in% type & "regular" %in% type)
    return(list(regular = object$var,
         adjusted = object$var_adj))
  if(type == "adjusted")
    return(object$var_adj)
  if(type == "regular")
    return(object$var)
}

#' Residuals for frailty models
#'
#' @param object An \code{emfrail} object
#' @param type One of \code{cluster} or \code{individual}
#' @param ... Other arguments
#' @return A vector corresponding to the Martingale residuals, either for each cluster or for each individual (row of the data).
#'
#' @details For cluster \eqn{i}, individual \eqn{j} and observation row \eqn{k}, we write the cumulative hazard contribution as
#' \deqn{\Lambda_{ijk} = \exp(\beta^\top \mathbf{x}_{ijk}) \Lambda_{0, ijk}}
#' where \eqn{\Lambda_{0, ijk}} is the baseline cumulative hazard correspinding to the row \eqn{(i,j,k)}.
#'
#' When \code{type == "individual"}, the returned residuals are equal to \eqn{z_i \Lambda_{ijk}} where \eqn{z_i} is the estimated frailty in cluster \eqn{i}.
#' When \code{type == "cluster"}, the returned residuals are equal to \eqn{\sum_{j,k} \Lambda_{ijk}},
#' @export
#'
residuals.emfrail <- function(object, type = "group", ...) {
  object$residuals
}

#' @export
#' @method model.matrix emfrail
model.matrix.emfrail <- function(object, ...) {
  if(is.null(object$mm)) stop("emfrail must be called with model.matrix = TRUE in order to return the model matrix")
  object$mm
}


#' @export
#' @method model.frame emfrail
model.frame.emfrail <- function(formula, ...) {
  if(is.null(formula$mf)) stop("emfrail must be called with model = TRUE in order to return the model frame")
  formula$mf
}
