summary.emfrail <- function(fit) {

  if (fit$est_dist$dist == "gamma") {
     var_est <- (fit$est$dist$frailtypar + fit$res$z$nev) / (fit$est$dist$frailtypar + fit$res$z$Lambda)^2

     low_z <- fit$res$z$z -
  }
}


(fit$est_dist$frailtypar + fit$res$z$nev) / (fit$est_dist$frailtypar + fit$res$z$Lambda)

qgamma(0.025, shape = fit$est_dist$frailtypar + fit$res$z$nev, scale = fit$est_dist$frailtypar + fit$res$z$Lambda)
lower_q <- qgamma(0.025, shape = fit$est_dist$frailtypar + fit$res$z$nev, rate = fit$est_dist$frailtypar + fit$res$z$Lambda)
upper_q <- qgamma(0.975, shape = fit$est_dist$frailtypar + fit$res$z$nev, rate = fit$est_dist$frailtypar + fit$res$z$Lambda)
lower_q

library(tidyverse)
library(plotly)

ppl1 <- data.frame(id = fit$res$z$id,
           z = fit$res$z$z,
           shape = fit$est_dist$frailtypar + fit$res$z$nev,
           rate = fit$est_dist$frailtypar + fit$res$z$Lambda) %>%
  arrange(z) %>%
  ggplot(aes(x = seq_along(id), y = z)) + geom_point() +
  geom_errorbar(aes(ymin = qgamma(0.025, shape = shape, rate = rate), ymax = qgamma(0.975, shape = shape, rate = rate),
                    id = id, gamma_shape = shape, gamma_rate = rate))


ggplotly(ppl1)
