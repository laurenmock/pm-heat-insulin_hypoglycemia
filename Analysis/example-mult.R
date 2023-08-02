library(dlnm)
library(mgcv)
library(data.table)
library(ggplot2)
library(viridis)

source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_basis.R")
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_pen.R")
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_pred.R")
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_plot.R")

# Create lagged data
data("chicagoNMMAPS")
setDT(chicagoNMMAPS)
lags <- 15
chicagoNMMAPS[, paste0("temp.", 0:lags) := shift(temp, 0:lags, type = "lag")]
chicagoNMMAPS[, paste0("rhum.", 0:lags) := shift(rhum, 0:lags, type = "lag")]
chicagoNMMAPS[, paste0("o3.", 0:lags) := shift(o3, 0:lags, type = "lag")]
chicagoNMMAPS[, paste0("pm10.", 0:lags) := shift(pm10, 0:lags, type = "lag")]

# separate summer/winter data
data <- chicagoNMMAPS[complete.cases(chicagoNMMAPS)]
summer_data <- data[month > 4 & month < 10]
winter_data <- data[month < 5 | month > 9]

# Set year/month/dow as factors for control in model
summer_data$year <- factor(summer_data$year)
summer_data$month <- factor(summer_data$month)
summer_data$dow <- factor(summer_data$dow)

# Summer model
# summer_deaths <- tdlnm(death ~ as.factor(year) + as.factor(month) + as.factor(dow),
#                        data = summer_data,
#                        exposure.data = as.matrix(summer_data[, paste0("temp.", 0:15)]),
#                        monotone = TRUE)
# (s_summer_deaths <- summary(summer_deaths))
# plot(s_summer_deaths, ylab = "degrees C", flab = "Additional\ndeaths")
# plot(s_summer_deaths, "slice", time = 1, xlab = "degrees C", ylab = "Additional deaths")
# plot(s_summer_deaths, "slice", val = 30, main = "30 degrees C", ylab = "Additional deaths")
#
# # Winter model
# winter_deaths <- tdlnm(death ~ as.factor(year) + as.factor(month) + as.factor(dow),
#                        data = winter_data,
#                        exposure.data = -as.matrix(winter_data[, paste0("temp.", 0:15)]),
#                        monotone = TRUE)
# (s_winter_deaths <- summary(winter_deaths, cenval = -20))
# plot(s_winter_deaths, ylab = "degrees C below zero", flab = "Additional deaths")
# plot(s_winter_deaths, "slice", time = 4, xlab = "degrees C below zero", ylab = "Additional deaths")
# plot(s_winter_deaths, "slice", val = 15, main = "15 degrees C below zero", ylab = "Additional deaths")

# list of exposures
x <- list(
  temp = as.matrix(summer_data[, paste0("temp.", 0:15)]),
  # rhum = as.matrix(summer_data[, paste0("rhum.", 0:15)]),
  pm10 = as.matrix(summer_data[, paste0("pm10.", 0:15)])
  # o3 = as.matrix(summer_data[, paste0("o3.", 0:15)])
)

# set up basis and cross basis
cb <- dlmm_basis(x, basis_all = list(fun = "ps", df = 4), interactions = "noself")
pen <- dlmm_pen(cb, shrink = T)

# fit model
sum_death_gam <- bam(death ~ as.factor(year) * as.factor(month) + as.factor(dow) + cb,
                     data = summer_data, family = poisson())#, paraPen = list(cb = pen))

# predict multiplicative effect holding one variable constant
q <- c(0.1, 0.25, 0.5, 0.75, 0.9)

pred <- lapply(1:5, function(p) {
  dlmm_pred(cb, sum_death_gam, xscale = sapply(x, sd), bylag = rep(0.5, 4),
            xmarg = sapply(x, quantile, probs = q[p]))
  })

exp <- "temp"
plot_dat <- rbindlist(lapply(1:5, function(p) {
  data.table(fit = pred[[p]]$marginal[[exp]]$fit,
             lag = pred[[p]]$lag_vals[[exp]],
             co_perc = q[p])
  }))

# plot independent and joint effects (heat map)
ggplot(plot_dat, aes(x = lag, y = fit.V1, color = factor(co_perc), linetype = factor(co_perc))) +
  geom_line(size = 2) +
  theme_bw() +
  labs(x = "Lag", y = "Additional deaths for sd unit increase in exposure", title = exp,
       color = "Percentiles of\nco-exposures", linetype = "Percentiles of\nco-exposures")

dlmm_plot(pred[[3]], "temp") + labs(title = "temp")
dlmm_plot(pred[[3]], "pm10") + labs(title = "pm10")
dlmm_plot(pred[[3]], "pm10", "temp") + labs(title = "temp/pm10")
