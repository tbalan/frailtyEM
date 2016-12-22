#' Print method for emfrail objects
#'
#' This function is included in case someone wants to access this output quickly,
#' and can look into the code.
#' @param obj An emfrail object
#'
#' @return Nothing
#' @export
print.emfrail <- function(obj) {
    outer_info <- obj[[1]]
    inner_info <- obj[[2]]
    cox_info <- obj[[3]]

    # cat('Object of class', class(obj), '\n')
    cat("Shared frailty model with frailty distribution:", inner_info$dist, "\n", "\n")

    if(inner_info$dist ==  "pvf") {
      cat("pvf m =", inner_info$pvfm, " ")
      if(inner_info$pvfm == -0.5) cat("(Inverse Gaussian)")
      if(inner_info$pvfm > 0) cat(" //e stimated mass at 0:", exp(-(inner_info$pvfm+1) / inner_info$pvfm *inner_info$theta))
      cat("\n")
    }
      cat()

    cat("Outer loops:", outer_info$fevals, "// optimizer:", rownames(outer_info), "\n")
    cat("(marginal) log-likelihood: ", -outer_info$value, "\n")
    cat("(no-frailty) log-likelihood:", cox_info$loglik[2], "\n\n")
    cat("LRT for frailty,",
        2 * (-outer_info$value - cox_info$loglik[2]),
        "on (chisq(1) + chisq(0))/2 => p =",
        pchisq(2 * (-outer_info$value - cox_info$loglik[2]), df = 1, lower.tail = FALSE)/2,
        "\n\n")


    cat("Frailty parameter:", inner_info$theta, "se: ", msm::deltamethod(~exp(x1), mean = outer_info$p1,
        cov = 1/attr(outer_info, "details")[[3]]), "\n")

    if (inner_info$dist %in% c("gamma", "pvf")) {
        low <- outer_info$p1 - 1.96 * sqrt(1/attr(outer_info, "details")[[3]])
        high <- outer_info$p1 + 1.96 * sqrt(1/attr(outer_info, "details")[[3]])

        lower_bound <- 1/exp(high)
        upper_bound <- 1/exp(low)

        cat("Frailty variance:", (1/inner_info$theta) %>% round(digits = 3) %>% format(nsmall = 2), "se:", msm::deltamethod(~1/exp(x1),
            mean = outer_info$p1, cov = 1/attr(outer_info, "details")[[3]]) %>% round(digits = 3) %>% format(nsmall = 2))
        cat(" // 95% CI: [", lower_bound %>% round(digits = 3) %>% format(nsmall = 2), ",", upper_bound %>% round(digits = 3) %>%
            format(nsmall = 2), "]")
        cat("\n")
    }




    cat("Regression coefficients:\n")



    tmp <- inner_info %>% with(cbind(coef, exp(coef), se_coef, se_coef_adj, coef/se_coef, 1 - pchisq((coef/se_coef)^2,
        1)))

    dimnames(tmp) <- list(names(inner_info$coef), c("coef", "exp(coef)", "se(coef)", "adjusted se", "z", "p"))

    printCoefmat(tmp, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE)
}




