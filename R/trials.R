
parfm::asthma
data("asthma")



parfm(Surv(Begin, End, Status)~Drug, cluster="Patid",
      data=asthma,
      dist="weibull", frailty="gamma")


coxph(Surv(Begin, End, Status)~Drug + frailty(Patid),
      data=asthma,
      ties = "breslow")

ml<- emfrail(asthma, Surv(Begin, End, Status)~Drug + cluster(Patid))


dat <- survival::rats
m1 <- emfrail(.data =  dat, .formula = Surv(rep(0, nrow(dat)), time, status) ~  rx + sex + cluster(litter),
              .distribution = emfrail_distribution(dist = "gamma"),
              .control = emfrail_control(opt_fit = TRUE))
m1

ddd <- boot::channing
head(ddd)


m2 <- emfrail(.data =  ddd, .formula = Surv(rep(0, nrow(ddd)), exit, cens) ~  sex + cluster(1:nrow(ddd)),
              .distribution = emfrail_distribution(dist = "gamma"),
              .control = emfrail_control(opt_fit = TRUE))

m2 <- emfrail(.data =  ddd, .formula = Surv(rep(0, nrow(ddd)), exit, cens) ~  sex + cluster(1:nrow(ddd)),
              .distribution = emfrail_distribution(dist = "gamma", frailtypar = 100000),
              .control = emfrail_control(opt_fit = FALSE))
m2


set.seed(2018)
frail <- rep(rgamma(100, 1, 1), each = 5)
x <- sample(c(0,1), 500, TRUE)
time <- rexp(exp(1/2 * x) * frail)

d1 <- data.frame(id = rep(1:100, each = 5),
           time,
           x,
           status = rep(1, 500))
d1_lt <- d1 %>%
  mutate(tstart = rexp(500)) %>%
  filter(tstart < time)


coxph(Surv(tstart, time, status) ~ x + frailty(id), ties = "breslow", data = d1_lt)



parfm(Surv(tstart, time, status) ~ x, cluster = "id", frailty = "gamma", data = d1_lt)

m3 <- emfrail(.data =  d1_lt, .formula = Surv(tstart, time, status) ~  x + cluster(id),
              .control = emfrail_control(opt_fit = TRUE))
)

m4 <- emfrail(.data =  d1_lt, .formula = Surv(tstart, time, status) ~  x + cluster(id),
              .control = emfrail_control(opt_fit = TRUE))
)
m4

m5 <- emfrail(.data =  d1_lt, .formula = Surv(tstart, time, status) ~  x + cluster(id),
              .control = emfrail_control(opt_fit = FALSE, verbose = TRUE),
              .distribution = emfrail_distribution(left_truncation = TRUE, frailtypar = 1/0.04))
)

pl1 <- lapply(c(0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.01),
       function(th) {
         emfrail(.data =  d1_lt, .formula = Surv(tstart, time, status) ~  x + cluster(id),
                 .control = emfrail_control(opt_fit = FALSE, verbose = TRUE),
                 .distribution = emfrail_distribution(left_truncation = TRUE, frailtypar = 1/th))

       })

emfrail(.data =  d1_lt, .formula = Surv(tstart, time, status) ~  x + cluster(id),
        .control = emfrail_control(opt_fit = FALSE, verbose = TRUE, eps = 0.00001, fast_fit = FALSE),
        .distribution = emfrail_distribution(left_truncation = TRUE, frailtypar = 1/0.005))



m5 <- emfrail(.data =  d1_lt, .formula = Surv(tstart, time, status) ~  x + cluster(id),
                    .control = emfrail_control(opt_fit = TRUE, fast_fit = FALSE),
                    .distribution = emfrail_distribution(left_truncation = TRUE, frailtypar = 1/0.06))
)

m5
