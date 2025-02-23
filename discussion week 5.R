library(survival)
library(survminer)
library(survMisc)
library(ggplot2)

set.seed(23)

# Simulate survival times (Exponential distribution with rate = 0.1)
n <- 100
survival_times <- rexp(n, rate = 0.1)

# Simulate censoring times (Uniform distribution)
censoring_times <- runif(n, min = 5, max = 20)

# Determine observed times and censoring status
observed_times <- pmin(survival_times, censoring_times)
event <- (survival_times <= censoring_times) 

dat <- data.frame(observed_times = observed_times,
                  event = event)

# Create survival object
dat_surv <- Surv(time = observed_times, event = event)

# Fit Kaplan-Meier estimate
km_fit <- survfit(dat_surv ~ 1, data = dat, conf.type = "log-log")

# Plot using ggsurvplot, no confidence interval (conf.int = F)
ggsurvplot(km_fit, data = dat_surv, ggtheme = theme_minimal(), conf.int = F,
           title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability")

#Add confidence intervals (conf.int = T)
ggsurvplot(km_fit, data = dat_surv, conf.int = T,
           ggtheme = theme_minimal(),
           title = "Kaplan-Meier Survival Curve", xlab = "Time", ylab = "Survival Probability")

#Add a risk table (risk.table = T)
ggsurvplot(km_fit, data = dat_surv, conf.int = T, risk.table = T,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability")

#Add a legend (legend.title, legend.labs, legend = "right")
#legend argument allows you to place the legend
ggsurvplot(km_fit, data = dat_surv, conf.int = T, risk.table = T,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability", 
           legend.title = "Survival Probability",
           legend.labs = c("Survival Estimate"),
           legend = "right") 

#Position legend so that it's inside the graph
#first coordinate in legend = c(x,y) is horizontal, second is vertical 
#e.g. c(0,0) would be bottom left, c(1,0) would be bottom right, 
#c(0,1) top left, etc
#something like c(0.8, 0.8) usually looks good (top right, but not too far)
ggsurvplot(km_fit, data = dat_surv, conf.int = T, risk.table = T,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability", 
           legend.title = "Survival Probability",
           legend.labs = c("Survival Estimate"),
           legend = c(0.8, 0.8)) 

#play around with colors
#conf.int.fill changes colors of the confidence interval
#palette changes colors of the estimator and the text for the risk set
ggsurvplot(km_fit, data = dat_surv, conf.int = T, conf.int.fill = "lightblue",
           risk.table = T,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability", 
           legend.title = "Survival Probability",
           legend.labs = c("Survival Estimate"),
           legend = c(0.8, 0.8),
           palette = "purple") 

#remove censoring marks (censor = F)
ggsurvplot(km_fit, data = dat_surv, conf.int = T,
           risk.table = T,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability", 
           legend.title = "Survival Probability",
           legend.labs = c("Survival Estimate"),
           legend = c(0.8, 0.8),
           palette = "purple",
           censor = F) 

#change censoring marks (censor.shape = 18; diff. numbers for diff. shapes)
ggsurvplot(km_fit, data = dat_surv, conf.int = T,
           risk.table = T,
           ggtheme = theme_minimal(), title = "Kaplan-Meier Survival Curve", 
           xlab = "Time", ylab = "Survival Probability", 
           legend.title = "Survival Probability",
           legend.labs = c("Survival Estimate"),
           legend = c(0.8, 0.8),
           palette = "purple",
           censor.shape = 18) 

##############
##############
##############

#Get point estimate and confidence interval for survival at t = 3
summary(km_fit, times = 3)

#Get point estimate and confidence interval for survival at t = 3 and t = 4
summary(km_fit, times = c(3, 4))

#Get the same for a sequence of times using the seq function
summary(km_fit, times = seq(from = 0, to = 5, by = 0.1))

#What if you want a different significance confidence interval?
#create a new survfit object and use the conf.int parameter to set confidence
km_fit_90 <- survfit(dat_surv ~ 1, data = dat, conf.type = "log-log", 
                     conf.int = 0.99)
summary(km_fit_90, times = 3)

#how can we get the median along with confidence interval?
quantile(km_fit, prob = 0.5)

#arbitrary quantiles?
quantile(km_fit, prob = 0.3)

#again, you can get multiple quantiles at once
quantile(km_fit, prob = seq(from = 0.1, to = 0.8, by = 0.1))

#and, again, by using another survfit object, you can change CI significance
quantile(km_fit_90, prob = 0.5)

##############
##############
##############

#simulate data with treatment and control groups

n <- 200
#randomly assign to treatment (1) and control (0)
group <- sample(c(1, 0), n, replace = TRUE)

# Generate survival times with different rates for each group
survival_times <- ifelse(group == 1, rexp(n, rate = 0.1), 
                         rexp(n, rate = 0.2))

#censoring times
censoring_times <- runif(n, min = 5, max = 20)

#observed times and event status
observed_times <- pmin(survival_times, censoring_times)
event <- as.numeric(survival_times <= censoring_times)

dat <- data.frame(observed_times = observed_times,
                  event = event,
                  group = group)

# Create a survival object
surv_obj <- Surv(time = observed_times, event = event)

# Fit Kaplan-Meier model stratified by group
km_fit <- survfit(surv_obj ~ group, data = dat, conf.type = "log-log")

# Plot survival curves for both groups together
ggsurvplot(km_fit, data = dat,
           palette = c("blue", "red"), 
           legend.title = "Group",
           legend.labs = c("Control", "Treatment"),
           title = "Kaplan-Meier Survival Curves",
           xlab = "Time", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal())

#same as above, but with confidence intervals, risk set, and p-value
ggsurvplot(km_fit, data = dat,
           conf.int = TRUE, 
           risk.table = TRUE, 
           pval = TRUE,
           palette = c("blue", "red"), 
           legend.title = "Group",
           legend.labs = c("Control", "Treatment"),
           title = "Kaplan-Meier Survival Curves",
           xlab = "Time", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal())

#we can "facet" this to display the plots side-by-side
p <- ggsurvplot(km_fit, data = dat,
           conf.int = T,
           risk.table = T,
           palette = c("blue", "red"), 
           legend.title = "Group",
           legend.labs = c("Control", "Treatment"),
           title = "Kaplan-Meier Survival Curves",
           xlab = "Time", 
           ylab = "Survival Probability",
           ggtheme = theme_minimal())

p$plot + theme_bw() + facet_wrap(~group)

#if we want to change the names of the facets (from 0 and 1 to ctrl/trt)
group_names <- c(`0` = "Control", `1` = "Treatment")
p$plot + theme_bw() + facet_wrap(~group, labeller = as_labeller(group_names))

##############
##############
##############

#compare survival probabilities at t = 5 using a Wald statistic
summary_km <- summary(km_fit, times = 5)
summary_km

#find survival probabilities and standard errors for both groups
surv_control <- summary_km$surv[1]
surv_treatment <- summary_km$surv[2]
se_control <- summary_km$std.err[1]
se_treatment <- summary_km$std.err[2]

#calculate Wald statistic
wald <- (surv_control - surv_treatment)/sqrt(se_control^2 + se_treatment^2)

#calculate two-sided p-value
p_value <- 2*(1 - pnorm(abs(wald)))

#other ways to calculate the p-value
2*pnorm(abs(wald), lower.tail = F)
2*pnorm(-abs(wald))

##############
##############
##############

#logrank test to compare overall curves
survdiff(surv_obj ~ group, data = dat)
surv_pvalue(km_fit, data = dat, method = "survdiff")

#WGB test to do a weighted log-rank test
#do this by setting rho = 1 for survdiff, or using "GB" for surv_pvalue
survdiff(surv_obj ~ group, data = dat, rho = 1)
surv_pvalue(km_fit, data = dat, method = "GB")

#set rho = 1.5 to do a Tarone-Ware test for survdiff, or "TB" for surv_pvalue
survdiff(surv_obj ~ group, data = dat, rho = 1.5)
surv_pvalue(km_fit, data = dat, method = "TW")

#Peto-Peto-Prentice
surv_pvalue(km_fit, data = dat, method =  "PP")

#Fleming-Harrington test
surv_pvalue(km_fit, data = dat, method = "FH")

#comp from survMisc seems quite buggy for me (and others online)
comp(ten(km_fit))

#rather than giving p-values, it just gives ranks.
#However, you can get the p-values back by using the Z statistic.
#if the Z statistic for a particular test is given by Z, the p-value is:
#2*(1 - pnorm(abs(Z)))

##############
##############
##############

set.seed(23)
#generate data
n <- 300  

#age will be a confounder
age_group <- sample(c("Young", "Old"), n, replace = TRUE, prob = c(0.5, 0.5))

#age affects group assignment
#0 = control, 1 = treatment
group <- ifelse(runif(n) < ifelse(age_group == "Young", 0.7, 0.3), 0, 1)  

#age also affects survival, thus it is a confounder
base_rate_young <- ifelse(group == 1, 0.04, 0.04)
base_rate_old <- ifelse(group == 1, 0.04, 0.14)

#assign hazards
rate <- ifelse(age_group == "Young", base_rate_young, base_rate_old)

#generate survival times using an exponential distribution
survival_times <- rexp(n, rate = rate)

#generate random censoring times
censoring_times <- runif(n, min = 5, max = 25)

#observed times and event status
observed_times <- pmin(survival_times, censoring_times)
event_status <- as.numeric(survival_times <= censoring_times)

#create data frame
surv_data <- data.frame(time = observed_times, 
                        status = event_status, group = factor(group), 
                        age_group = factor(age_group))

#KM
km_fit <- survfit(Surv(time, status) ~ group, data = surv_data)

#log-rank, ignoring confounder
survdiff(Surv(time, status) ~ group, data = surv_data)

#stratified log-rank
survdiff(Surv(time, status) ~ group + strata(age_group), data = surv_data)
      