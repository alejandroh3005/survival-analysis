# Load required libraries
library(ggplot2)  # For visualization

# 1. Simulate Survival Data
set.seed(123)

n <- 1000  # Number of individuals
true_mean <- 50  # True mean survival time

# Generate true survival times from an exponential distribution
survival_times <- rexp(n, rate = 1 / true_mean)

# Generate censoring times
censoring_times <- rexp(n, rate = 1 / 60)

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

# Distribution of event times vs. observed times
dist_plot <- ggplot(dat, aes(x = survival_times, fill = "Event Times")) +
  geom_density(alpha = 0.5, color = NA) +
  geom_density(aes(x = observed_times, fill = "Observed Times"), 
               alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("Event Times" = "red", 
                               "Observed Times" = "blue")) +
  labs(title = "Distribution of Event Times vs. Observed Times",
       x = "Time",
       y = "Density",
       fill = "Legend") +
  theme_minimal()

# Display the plots
print(dist_plot)

#pretending observed times are survival times
mean(dat$observed_times)

#throwing out all censored times, just taking the average of uncensored obs
mean(dat[dat$event == 1,]$survival_times)

#mean of the actual survival times (not possible to do in real life!)
mean(dat$survival_times)
