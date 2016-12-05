#' Title
#'
#' @param obj
#'
#' @return muie
#' @export
print.emfrail <- function(obj) {
    outer_info <- obj[[1]]
    inner_info <- obj[[2]]
    
    # cat('Object of class', class(obj), '\n')
    cat("Shared frailty model with frailty distribution:", inner_info$dist, "\n", "\n")
    
    cat("Outer loops:", outer_info$fevals, "// optimizer:", rownames(outer_info), "\n")
    cat("(marginal) log-likelihood: ", -outer_info$value, "\n\n")
    
    cat("Frailty parameter (see parametrization):", inner_info$frailtypar, "se: ", msm::deltamethod(~exp(x1), mean = outer_info$p1, 
        cov = 1/attr(outer_info, "details")[[3]]), "\n")
    
    if (inner_info$dist %in% c("gamma", "pvf")) {
        low <- outer_info$p1 - 1.96 * sqrt(1/attr(outer_info, "details")[[3]])
        high <- outer_info$p1 + 1.96 * sqrt(1/attr(outer_info, "details")[[3]])
        
        lower_bound <- 1/exp(high)
        upper_bound <- 1/exp(low)
        
        cat("Frailty variance:", (1/inner_info$frailtypar) %>% round(digits = 3) %>% format(nsmall = 2), "se:", msm::deltamethod(~1/exp(x1), 
            mean = outer_info$p1, cov = 1/attr(outer_info, "details")[[3]]) %>% round(digits = 3) %>% format(nsmall = 2))
        cat(" // 95% CI: [", lower_bound %>% round(digits = 3) %>% format(nsmall = 2), ",", upper_bound %>% round(digits = 3) %>% 
            format(nsmall = 2), "]")
        cat("\n")
    }
    
    
    cat("Regression coefficients:\n")
    
    
    
    tmp <- inner_info %>% with(cbind(coef, exp(coef), se[1:length(coef)], coef/se[1:length(coef)], 1 - pchisq((coef/se[1:length(coef)])^2, 
        1)))
    
    dimnames(tmp) <- list(names(inner_info$coef), c("coef", "exp(coef)", "se(coef)", "z", "p"))
    
    printCoefmat(tmp, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE)
}

