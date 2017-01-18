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
      if(inner_info$pvfm > 0) cat(" //estimated mass at 0:", exp(-(inner_info$pvfm+1) / inner_info$pvfm *inner_info$theta))
      cat("\n")
    }
      cat()

    cat("Outer loops:", outer_info$fevals, "// optimizer:", rownames(outer_info), "\n")
    cat("(marginal) log-likelihood: ", -outer_info$value, "\n")
    cat("(no-frailty) log-likelihood:", cox_info$loglik[length(cox_info$loglik)], "\n\n")
    cat("LRT for frailty,",
        2 * (-outer_info$value - cox_info$loglik[length(cox_info$loglik)]),
        "on (chisq(1) + chisq(0))/2 => p =",
        pchisq(2 * (-outer_info$value - cox_info$loglik[length(cox_info$loglik)]), df = 1, lower.tail = FALSE)/2,
        "\n\n")


    # we always estimate the log of the frailty parameter so we need the delta method for the SE



    if (inner_info$dist %in% c("gamma", "pvf")) {

      cat("Frailty parameter:", inner_info$theta, "se: ", msm::deltamethod(~exp(x1), mean = outer_info$p1,
                                                                           cov = 1/attr(outer_info, "details")[[3]]), "\n")
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

    # for the stable; there are several parametrizations; here theta is a sort of thetatilde = theta + 1
    # if (inner_info$dist == "stable") {
    #   cat("pvf bbeta =",  1 - 1/inner_info$theta, "pvf alpha = ", 1/(1 - 1/inner_info$theta), "\n")
    #   cat("Frailty parameter where bbeta = 1 - (1 / (theta)), theta>1:", inner_info$theta, "se: ", msm::deltamethod(~exp(x1) , mean = outer_info$p1,
    #                                                                        cov = 1/attr(outer_info, "details")[[3]]), "\n")

      # bbeta <- 1 - 1/inner_info$theta
      # alphaprime <- bbeta - log(bbeta)
      #
      # cat("Frailty parameter (alpha = 1):", alphaprime,
      #     "se: ", msm::deltamethod(~exp(x1), mean = alphaprime ,
      #                              cov = 1/attr(outer_info, "details")[[3]]), "\n")
      #
      # cat("Kendall's tau")
    # }

    if (inner_info$dist == "stable") {
      theta <- exp(outer_info$p1) / (1 + exp(outer_info$p1))

      cat("pvf bbeta (attenuation factor) =",  1 - theta , "/ pvf alpha = 1\n")


      low <- outer_info$p1 - 1.96 * sqrt(1/attr(outer_info, "details")[[3]])
      high <- outer_info$p1 + 1.96 * sqrt(1/attr(outer_info, "details")[[3]])

      upper_bound <- exp(high) / (1 + exp(high))
      lower_bound <- exp(low) / (1 + exp(low))

      cat("E[logY] =", -1* (1 / (1-theta) - 1) * digamma(1),"\n")
      cat("theta = 1 - bbeta, 0<theta<1:",
          theta %>% round(digits = 3) %>% format(nsmall = 2),
          "se: ",
          msm::deltamethod(~exp(x1) / (exp(x1) + 1),
                           mean = outer_info$p1, # logfrailtypar
                           cov = 1/attr(outer_info, "details")[[3]]) %>% round(digits = 3) %>% format(nsmall = 2),
          "\n")

      cat(" // 95% CI: [", lower_bound %>% round(digits = 3) %>% format(nsmall = 2), ",", upper_bound %>% round(digits = 3) %>%
        format(nsmall = 2), "]")

      cat("\n")
      cat("Kendall's tau", (inner_info$theta / (1 + inner_info$theta)) %>% round(digits = 3) %>% format(nsmall = 2),"\n")

      }

    if(!is.null(inner_info$coef)) {
      cat("Regression coefficients:\n")



      tmp <- inner_info %>% with(cbind(coef, exp(coef), se_coef, se_coef_adj, coef/se_coef, 1 - pchisq((coef/se_coef)^2,
                                                                                                       1)))

      dimnames(tmp) <- list(names(inner_info$coef), c("coef", "exp(coef)", "se(coef)", "adjusted se", "z", "p"))

      printCoefmat(tmp, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE)
    }

    cat("\n\n")
    invisible(obj)

}




