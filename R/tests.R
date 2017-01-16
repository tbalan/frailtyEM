# # Rats --------------------------------------------------------------------
#
#
# dat <- survival::rats
# m1 <- emfrail(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#               .distribution = emfrail_distribution(dist = "gamma"))
# m1
#
#
# m2 <- emfrail(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#               .distribution = emfrail_distribution(dist = "pvf"))
# m2
#
# m3 <- emfrail(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#               .distribution = emfrail_distribution(dist = "stable"),
#               .control = emfrail_control(verbose = FALSE))
# m3
#
# m4 <- emfrail(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#               .distribution = emfrail_distribution(dist = "stable2"),
#               .control = emfrail_control(verbose = FALSE))
# m4
#
# # parfm doesn't work with this data set
# parfm(Surv(rep(0, nrow(dat)), time, status) ~ rx + sex, cluster = "litter",dist = "weibull", frailty = "gamma", data = dat)
#
# # profile likelihoods
# par(mfrow = c(2,2))
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#                         .distribution = emfrail_distribution(dist = "gamma"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "gamma")
#
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#                         .distribution = emfrail_distribution(dist = "pvf"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "IG")
#
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1/2),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1/2")
#
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1")
#
#
# par(mfrow = c(1,1))
#
#
#
# stabpar <- seq(from =1, to = 30, by = 1)
# plot(stabpar, emfrail_pll(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#                           .distribution = emfrail_distribution(dist = "stable"),
#                           .values = stabpar),
#      ylab = "log-likelihood",
#      xlab = "1-beta",
#      type = "l",
#      main = "stable")
#
# lines(stabpar, emfrail_pll(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
#                             .distribution = emfrail_distribution(dist = "stable2"),
#                             .values = stabpar),
#      ylab = "log-likelihood",
#      xlab = "kendall's tau (1-beta)",
#      type = "l",
#      main = "stable2", col = 2)
#
#
#
# # Kidney ------------------------------------------------------------------
#
# data(kidney)
# # type 'help(kidney)' for a description of the data set
# kidney$sex <- kidney$sex - 1
#
#
#
# m1 <- emfrail(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#               .distribution = emfrail_distribution(dist = "gamma"))
# m1
#
# parfm(Surv(time,status) ~ sex + age, cluster="id",
#       data=kidney, dist="weibull", frailty="gamma")
#
# coxph(Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + frailty(id), data = kidney, ties = "breslow")
#
#
# m2 <- emfrail(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#               .distribution = emfrail_distribution(dist = "pvf"))
# m2
#
# parfm(Surv(time,status) ~ sex + age, cluster="id",
#       data=kidney, dist="weibull", frailty="ingau")
#
# # this fails
# m3 <- emfrail(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#               .distribution = emfrail_distribution(dist = "stable"))
# m3
#
# # this fails
# m4 <- emfrail(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#               .distribution = emfrail_distribution(dist = "stable2"))
# m4
#
#
# # profile likelihoods
# par(mfrow = c(2,2))
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#                         .distribution = emfrail_distribution(dist = "gamma"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "gamma")
#
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "IG")
#
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1/2),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1/2")
#
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1")
#
#
#
# par(mfrow = c(1,1))
#
#
#
# stabpar <- seq(from =1, to = 50, by = 1)
# plot(stabpar, emfrail_pll(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#                           .distribution = emfrail_distribution(dist = "stable"),
#                           .values = stabpar),
#      ylab = "log-likelihood",
#      xlab = "1-beta",
#      type = "l",
#      main = "stable")
#
# lines(stabpar, emfrail_pll(.data =  kidney, .formula = Surv(rep(0, nrow(kidney)), time, status) ~  sex + age + cluster(id),
#                            .distribution = emfrail_distribution(dist = "stable2"),
#                            .values = stabpar),
#       ylab = "log-likelihood",
#       xlab = "kendall's tau (1-beta)",
#       type = "l",
#       main = "stable2", col = 2)
#
#
#
#
#
# # simulated data ----------------------------------------------------------
#
# set.seed(1)
# x <- sample(c(0,1), 300, TRUE)
# z <- rep(rgamma(100, 1, 1), each = 3)
#
# time <- rexp(300, rate = z * exp(0.5*x) )
#
# censtime <- 5
# status <- rep(1, 300)
#
# status[time >= censtime] <- 0
# time[status == 0] <- censtime
# time0 <- rep(0, 300)
#
# dd <- data.frame(id = rep(1:100, each = 3), x = x, time0 = time0, time = time, status = status)
#
#
# m1 <- emfrail(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "gamma"))
# m1
#
# #
# parfm(Surv(time, status)~x, cluster="id",
#       data=dd, dist="weibull", frailty="gamma")
#
# coxph(Surv(rep(0, nrow(dd)), time, status) ~  x + frailty(id), data = dd, ties = "breslow")
#
#
# m2 <- emfrail(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "pvf"))
# m2
#
# parfm(Surv(time,status) ~ x, cluster="id",
#       data=dd, dist="weibull", frailty="ingau")
#
# # this works.
# m3 <- emfrail(.data = dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "stable"))
# m3
#
# # this also works.
# m4 <- emfrail(.data =  dd, .formula =  Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "stable2"))
# m4
#
#
# # profile likelihoods
# par(mfrow = c(2,2))
#
# frvar <- seq(from = 0.1, to = 2.9, by = 0.18)
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "gamma"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "gamma")
#
#
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "IG")
#
#
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1/2),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1/2")
#
#
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1")
#
#
#
# # This has a very strange shape for the log-likelihood in theta.
# par(mfrow = c(1,1))
#
# # l(theta) here means with L(c) = exp(-alpha c^(1 - 1/theta))
# # x axis - 1/stabpar means L(c) = exp(-alpha c^(1 - x))
#
#
# # emfrail_pll(.data =  dd,
# #             .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
# #             .distribution = emfrail_distribution(dist = "stable"),
# #             .values = 3)
#
#
# stabpar <- seq(from =1, to = 5, by = .1)
# plot(1/stabpar, emfrail_pll(.data =  dd,
#                           .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                           .distribution = emfrail_distribution(dist = "stable"),
#                           .values = stabpar),
#      ylab = "log-likelihood",
#      xlab = "1-beta",
#      type = "l",
#      main = "stable")
#
# lines(1/stabpar, emfrail_pll(.data =  dd,
#                            .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                            .distribution = emfrail_distribution(dist = "stable2"),
#                            .values = stabpar),
#       ylab = "log-likelihood",
#       xlab = "kendall's tau (1-beta)",
#       type = "l",
#       main = "stable2", col = 2)
#
#
# # simulated data with IG frailty ------------------------------------------
# library(statmod)
#
# set.seed(1)
# x <- sample(c(0,1), 300, TRUE)
# z <- rep(rinvgauss(100), each = 3)
#
# time <- rexp(300, rate = z * exp(0.5*x) )
#
# censtime <- 5
# status <- rep(1, 300)
#
# status[time >= censtime] <- 0
# time[status == 0] <- censtime
# time0 <- rep(0, 300)
#
# dd <- data.frame(id = rep(1:100, each = 3), x = x, time0 = time0, time = time, status = status)
#
# # this is very different with the coxph
# m1 <- emfrail(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "gamma"))
# m1
#
# #
# parfm(Surv(time, status)~x, cluster="id",
#       data=dd, dist="weibull", frailty="gamma")
#
# coxph(Surv(rep(0, nrow(dd)), time, status) ~  x + frailty(id), data = dd, ties = "breslow")
#
#
# m2 <- emfrail(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "pvf"))
# m2
#
# parfm(Surv(time,status) ~ x, cluster="id",
#       data=dd, dist="weibull", frailty="ingau")
#
# # this works.
# m3 <- emfrail(.data = dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "stable"),
#               .control =emfrail_control( verbose = TRUE))
# m3
#
# # this also works.
# m4 <- emfrail(.data =  dd, .formula =  Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#               .distribution = emfrail_distribution(dist = "stable2"),
#               .control = emfrail_control(verbose = TRUE))
# m4
# # beta is 0.8055, so the marginal effect would be 0.481 x 0.8055 = 0.387
#
# coxph(Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id), data = dd)
# # 0.354 funny
#
#
# parfm(Surv(rep(0, nrow(dd)), time, status) ~  x,
#       cluster = "id", data = dd,
#       dist = "exponential", frailty = "possta")
#
# # profile likelihoods
# par(mfrow = c(2,2))
#
# frvar <- seq(from = 0.1, to = 1.5, by = 0.1)
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "gamma"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "gamma")
#
#
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf"),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "IG")
#
#
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1/2),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1/2")
#
#
# plot(frvar, emfrail_pll(.data =  dd, .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                         .distribution = emfrail_distribution(dist = "pvf",
#                                                              pvfm = 1),
#                         .values = 1/frvar),
#      ylab = "log-likelihood",
#      xlab = "frailty variance",
#      type = "l",
#      main = "pvf m=1")
#
#
#
# par(mfrow = c(1,1))
# stabpar <- seq(from =1, to = 10, by = .1)
# plot(1/stabpar, emfrail_pll(.data =  dd,
#                             .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                             .distribution = emfrail_distribution(dist = "stable"),
#                             .values = stabpar),
#      ylab = "log-likelihood",
#      xlab = "1-beta",
#      type = "l",
#      main = "stable")
#
# lines(1/stabpar, emfrail_pll(.data =  dd,
#                              .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#                              .distribution = emfrail_distribution(dist = "stable2"),
#                              .values = stabpar),
#       ylab = "log-likelihood",
#       xlab = "kendall's tau (1-beta)",
#       type = "l",
#       main = "stable2", col = 2)
#
#
# emfrail(.data =  dd,
#         .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#         .distribution = emfrail_distribution(dist = "stable2"),
#         .control = emfrail_control(opt_fit = FALSE, verbose = TRUE))
#
# emfrail(.data =  dd,
#         .formula = Surv(rep(0, nrow(dd)), time, status) ~  x + cluster(id),
#         .distribution = emfrail_distribution(dist = "stable"),
#         .control = emfrail_control(opt_fit = FALSE, verbose = TRUE))
#
#
#
#
# # Left truncation
# # simulate 300 clusters of size 5
# set.seed(17)
# x <- sample(c(0,1), 5 * 300, TRUE)
# u <- rep(rgamma(300,1,1), each = 5)
# stime <- rexp(5*300, rate = u * exp(x))
# ltime <- runif(5 * 300)
#
# library(tidyverse)
# d <- data.frame(id = rep(1:300, each = 5),
#                 x = x,
#                 stime = stime,
#                 u = u,
#                 ltime = ltime,
#                 status = 1)
#
#
# d1 <- d %>% filter(stime > ltime)
# # this is the same as the cph (naive):
# mymod <- d1 %>%
#   emfrail(Surv(ltime, stime, status)~ x + cluster(id), .control = emfrail_control(verbose = TRUE))
#
# # this is the correct way here:
# mymod <- d1 %>%
#   emfrail(Surv(ltime, stime, status)~ x + cluster(id), .control = emfrail_control(verbose = TRUE),
#           .distribution = emfrail_distribution(left_truncation = TRUE))
