library(survival)
library(ggplot2)
library(survminer)
library(msm)

set.seed(23)

#sample size
n <- 1000  

#income, minimum value is 0
income <- pmax(rnorm(n, mean = 50000, sd = 15000), 0)

#treatment assignment biased by income 
#wealthier patients more likely to get surgery
treatment_prob <- plogis((income - 50000) / 10000)
treatment <- rbinom(n, 1, treatment_prob) 

#doctor performing the surgery
doctors <- factor(sample(1:5, n, replace = TRUE))

#generate survival times using weibull
#different baseline hazard for each doctor.
baseline_hazard <- c(0.01, 0.02, 0.015, 0.025, 0.02)
doctor_effect <- baseline_hazard[as.numeric(doctors)]  

#parameters for the Weibull model
treatment_effect <- 2  
income_effect <- 0.000005 
shape <- 0.35
scale <- exp(doctor_effect + treatment * treatment_effect 
             + income * income_effect)

time_to_death <- rweibull(n, shape = shape, scale = scale)

#censoring time (study ends after 15 years)
censor_time <- runif(n, min = 0, max = 4)

#generate observed times and event indicators
observed_time <- pmin(time_to_death, censor_time)
event <- as.numeric(time_to_death <= censor_time)

#create data frame
dat <- data.frame(observed_time, event, treatment, 
                  income, doctor = as.factor(doctors))

#######################################################
#######################################################
#######################################################

#first attempt to answer the question
cox1 <- coxph(Surv(observed_time, event) ~ treatment, data = dat)

#summary
summary(cox1)

#get confidence intervals of other levels using the level argument
#e.g. getting 90%, 99%, and 99.9% CIs, respectively:
confint(cox1, level = 0.9)
confint(cox1, level = 0.99)
confint(cox1, level = 0.999)

#second attempt to answer the question
#now adjusting for income and doctor
cox2 <- coxph(Surv(observed_time, event) ~ treatment + income + doctor, 
              data = dat)

summary(cox2)

#adjusting for income and doctor, but with re-scaled income for interpretability
dat$income_scaled <- dat$income/10000
cox3 <- coxph(Surv(observed_time, event) ~ treatment + income_scaled + doctor, 
              data = dat)

summary(cox3)

#adjusting for income and doctor, but dichotomizing income for interpretability
dat$income_binary <- ifelse(dat$income < 50000, 0, 1)
cox4 <- coxph(Surv(observed_time, event) ~ treatment + income_binary + doctor, 
              data = dat)

summary(cox4)

#adjusting for income, stratifying for doctor (scaled income)
cox5 <- coxph(Surv(observed_time, event) ~ treatment + income_scaled 
              + strata(doctor), data = dat)

summary(cox5)

#first attempt to include interaction for high/low income
cox6 <- coxph(Surv(observed_time, event) ~ treatment + treatment*(income >50000) 
              + doctor, data = dat)

summary(cox6)

#second attempt
cox7 <- coxph(Surv(observed_time, event) ~ treatment*(income >50000) + treatment*(income >50000) 
              + doctor + income, data = dat)

summary(cox7)

#third attempt (the best one!)
dat$rxhigh <- dat$treatment*(dat$income > 50000)
dat$rxlow <- dat$treatment*(dat$income <= 50000)

cox8 <- coxph(Surv(observed_time, event) ~ rxhigh + rxlow + income 
              + doctor, data = dat)

summary(cox8)

#interaction with stratification
cox9 <- coxph(Surv(observed_time, event) ~ rxhigh + rxlow + income_scaled
              + strata(doctor), data = dat)

summary(cox9)

#use ANOVA to see whether this interaction is significant
#compare models 9 and 5
anova(cox5, cox9)

#compare models 3 and 8 with log likelihood
l_3 <- logLik(cox3)
l_8 <- logLik(cox8)

-2*(l_3 - l_8)




#graphing survival curves for individuals given a cox model
new_data <- data.frame(
  treatment = c(1, 0, 1), #first, second, third patients did, didn't, did receive experimental
  income = c(30000, 70000, 80000),
  doctor = factor(c(1, 3, 4), levels = c(1, 2, 3, 4, 5))
)

surv_pred <- survfit(cox2, newdata = new_data)

ggsurvplot(surv_pred, 
           data = new_data,
           conf.int = TRUE,  #CIs
           risk.table = TRUE,  #risk table
           ggtheme = theme_minimal(),
           legend.labs = c("Treated, Income $30k, Doctor 1", 
                           "Placebo, Income $70k, Doctor 3",
                           "Treated, Income $80k, Doctor 4"),
           title = "Estimated Survival Curves for Selected Individuals")




#constructing a confidence interval around a hazard ratio
pe_hr <- exp(coef(cox5)["treatment"] + 2*coef(cox5)["income_scaled"])

se_hr <- deltamethod(g =~ exp((x1+2*x2)),
            mean = coef(cox5)[c("treatment", "income_scaled")],
            cov = vcov(cox5)[c("treatment", "income_scaled"),
                               c("treatment", "income_scaled")])

lower_ci <- pe_hr - 1.96*se_hr
upper_ci <- pe_hr + 1.96*se_hr