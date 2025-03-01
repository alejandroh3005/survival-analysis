set.seed(123)

n <- 1000  # Number of individuals
true_mean <- 50  # True mean survival time

# Generate true survival times from an exponential distribution
survival_times <- rexp(n, rate = 1 / true_mean)

# Generate censoring times
flip <- rbinom(n, size = 1, prob = 0.3)
censoring_times <- ifelse(flip == 1, Inf, 5)

# Observed times are the minimum of survival and censoring times
observed_times <- pmin(survival_times, censoring_times)

# Create a censoring indicator: 1 = event, 0 = censored
event <- as.numeric(survival_times <= censoring_times)

# Combine into a data frame
dat <- data.frame(
  survival_times = survival_times,
  observed_times = observed_times,
  event = event
)

#finds proportion of those in the risk set at t = 5,
#i.e. uncensored and alive at time 5, who survived until t = 10
prop <- mean(dat[dat$observed_times > 5 & dat$event == 1,]$observed_times > 10)

#estimates the number of censored individuals who survived until t = 10
est_censored_survived <- sum(dat$event == 0)*prop

#counts number of uncensored individuals who survived until time t = 10
num_uncensored_survived <- sum(dat[dat$event == 1,]$observed_times > 10)

#final estimate!
(est_censored_survived + num_uncensored_survived)/n

#compare to the truth:
1 - pexp(q = 10, rate = 1/true_mean)
mean(dat$survival_times > 10)
