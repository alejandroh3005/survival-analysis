library(survival)
library(flexsurv)
library(msm)
source("fitparametric.R")

set.seed(23)

#number of individuals, true lambda
n <- 1000
lambda <- 5 

#generate survival and censoring times
survival_times <- rexp(n, rate = lambda)
censoring_times <- rexp(n, rate = 1/2*lambda)

#find observed times and event indicator, create dataframe
observed_times = pmin(survival_times, censoring_times)
event <- as.numeric(survival_times <= censoring_times)

dat <- data.frame(
  observed_times = observed_times,
  event = event
)

#create surv object - necessary for fitting survival models in R
dat.surv <- Surv(time = dat$observed_times, event = dat$event, type = "right")

#fit exponential model using flexsurv
model <- flexsurvreg(dat.surv ~ 1, dist = "exponential")

#see how the model looks
model

#look just at the coefficient estimate/CI
model$res

#WARNING: using model$cov gives you the SE estimate for log(lambda), not lambda!
#In general, be careful with using estimates from R output
model$cov

#estimate of survival function at t = 0.5
surv_est <- exp(-.5*model$res["rate", "est"])

#finding the CI manually
#find by calculus that h'(lambda) = -0.5 exp(-0.5*lambda)
#using delta method
surv_se <- sqrt((-0.5*exp(-0.5*model$res["rate", "est"]))^2*model$res["rate", "se"]^2)

#construct 95% confidence interval
surv_ci_lower <- surv_est - 1.96*surv_se
surv_ci_upper <- surv_est + 1.96*surv_se


#finding the CI using deltamethod from msm
surv_se_msm <- deltamethod(g=~exp(-0.5*x1), mean = model$res["rate", "est"],
            cov = model$res["rate", "se"]^2, ses=TRUE)

#note that surv_se (calculated manually) and surv_se_msm (using the package) agree!
surv_ci_lower_msm <- surv_est - 1.96*surv_se_msm
surv_ci_upper_msm <- surv_est + 1.96*surv_se_msm

#inference on mean using fitparametric
fitparametric(dat.surv, dist = "exp", feature = "mean")

#inference on median (0.5 quantile) using fitparametric
#to get another quantile, e.g. if you want 0.3 quantile, use pi = 0.3
fitparametric(dat.surv, dist = "exp", feature = "quantile", pi = 0.5)

#inference on survival probability at t = 0.5, i.e. P(T > 0.5)
#can change t to do the same for different time points
#note that this matches what we got earlier!
fitparametric(dat.surv, dist = "exp", feature = "survival", t = 0.5)

#inference on prob. of surviving longer than 1, given survival longer than 0.5 
#i.e. P(T > 1 | T > 0.5)
#can change t and t0 to calculate different probabilities
fitparametric(dat.surv, dist = "exp", feature = "condsurvival", t = 1, t0 = 0.5)

#What about a Weibull model?
fitparametric(dat.surv, dist = "weibull", feature = "condsurvival", t = 1, t0 = 0.5)
