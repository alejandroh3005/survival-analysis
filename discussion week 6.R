# Load necessary libraries
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)

# Set seed for reproducibility
set.seed(123)

#sample size
n <- 1000

#0 = low pollution, 1 = high
#0 = nonsmoker, 1 = smoker
pollution <- rbinom(n, 1, 0.5)
age <- rnorm(n, mean = 50, sd = 10)
smoking <- rbinom(n, 1, 0.3)

#generate survival times
linear_predictor <- 0.5 * pollution + 0.02 * age + 0.8 * smoking + 0.03 * age*smoking
hazard_rate <- 0.01 * exp(linear_predictor)
survival_time <- rexp(n, rate = hazard_rate)

#censoring times
censoring_time <- rexp(n, rate = 0.01)

#observed times and event indicator
observed_time <- pmin(survival_time, censoring_time)
event <- ifelse(survival_time <= censoring_time, 1, 0)

#dataset
dat <- data.frame(
  id = 1:n,
  pollution = pollution,
  age = age,
  smoking = smoking,
  survival_time = observed_time,
  event = event
)

surv_obj <- Surv(dat$survival_time, dat$event)

#fit KM estimators
km_fit <- survfit(surv_obj ~ pollution + smoking, data = dat)

#KM curve
ggsurvplot(
  km_fit,
  data = dat,
  pval = TRUE,
  conf.int = TRUE,
  palette = "jco",
  legend.labs = c("Low Pollution, Non-Smoker", 
                  "High Pollution, Non-Smoker", 
                  "Low Pollution, Smoker", 
                  "High Pollution, Smoker"),
  legend = c(0.8, 0.8),
  xlab = "Time (days)", 
  ylab = "Survival Probability",
  break.time.by = 20,  #break x-axis by 20 units
  risk.table = TRUE,
  risk.table.height = 0.25, #adjust risk table height
  ggtheme = theme_minimal(),
  title = "KM curves by pollution exposure and smoking status",
  pval.coord = c(200, 0.5)
)

#Cox model, no interaction
cox_model <- coxph(surv_obj ~ pollution + age + smoking, data = dat)
summary(cox_model)

#Cox model, interaction
cox_model_interaction <- coxph(surv_obj ~ pollution + age + smoking + age:smoking, data = dat)
summary(cox_model_interaction)

#Use ANOVA to check whether the interaction is significant
anova(cox_model, cox_model_interaction)
