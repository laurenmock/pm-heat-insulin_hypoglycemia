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

# list of exposures
x <- list(
  temp = as.matrix(summer_data[, paste0("temp.", 0:15)]),
  pm10 = as.matrix(summer_data[, paste0("pm10.", 0:15)])
)

# set up basis and cross basis
cb <- dlmm_basis(x, basis_all = list(fun = "ps", df = 4), interactions = "noself")

# fit model
formula <- as.formula(paste("death ~ year * month + dow +",
                            paste(colnames(cb), collapse = "+")))
model_for_predictions <- bam(formula, data = cbind.data.frame(summer_data, cb), family = poisson())

# Function for counterfactual exposure predictions
cf_predict <- function(model, data, cross_basis, cf_exposure) {
  
  cf_mat <- lapply(cf_exposure, function(x) {
    matrix(x, nrow(data), length(x), byrow = TRUE)
  })
  
  cf_cb <- dlmm_basis(cf_mat, 
                      basis_each = attributes(cross_basis)$basis,
                      interactions = attributes(cross_basis)$interactions)
  
  return(predict(model, cbind(data, cf_cb), type = "response"))
  
}

w1x1 <- cf_predict(model = model_for_predictions, data = summer_data, cross_basis = cb,
                   cf_exposure = list(temp = rep(30,16), pm10 = rep(50, 16)))

w1x0 <- cf_predict(model = model_for_predictions, data = summer_data, cross_basis = cb,
                   cf_exposure = list(temp = rep(30,16), pm10 = rep(40, 16)))

w0x1 <- cf_predict(model = model_for_predictions, data = summer_data, cross_basis = cb,
                   cf_exposure = list(temp = rep(20,16), pm10 = rep(50, 16)))

w0x0 <- cf_predict(model = model_for_predictions, data = summer_data, cross_basis = cb,
                   cf_exposure = list(temp = rep(20,16), pm10 = rep(40, 16)))

mean(w1x1) - mean(w1x0) - mean(w0x1) + mean(w0x0)
