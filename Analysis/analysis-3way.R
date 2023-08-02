require(magrittr)
require(haven)
library(readr)
require(lubridate)
require(dplyr)
library(fst)
library(data.table)
library(survival)
library(dlnm)
library(mgcv)
library(ggplot2)

# load modules
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_basis.R")
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_pen.R")
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_pred.R")
source("M:/External Users/KevinJos/code/dlm/Functions/dlmm_plot.R")

## Fit DLMM Interaction Model

load("M:/External Users/KevinJos/data/dlm/cco.RData")

cco <- cco[sample(1:nrow(cco), 20000),] # sample for testing

x <- list(
  drug = as.matrix(cco[, paste0("drug_lag", 0:6)]),
  temp = as.matrix(cco[, paste0("lag", 0:6)]),
  pm25 = as.matrix(cco[, paste0("pm25_lag", 0:6)])
)

cb <- dlmm_basis(x, basis_all = list(fun = "cr", df = 4), interactions = "noself") # all, noself, none
formula <- as.formula(paste("case ~ year + month + dow + cluster(bene_id) + ",
                            paste(colnames(cb), collapse = "+")))
# formula <- as.formula(paste("case ~ year + month + dow + cluster(bene_id) + ",
#                             "pm25_lag0*lag0*drug_lag0")) # simple alternative (no dlm)
sum_clogit <- clogit(formula, data = cbind.data.frame(cco, cb), method = "efron")

# Function for counterfactual exposure predictions
cf_predict <- function(model, data, cross_basis, cf_exposure) {
  
  
  cf_mat <- lapply(cf_exposure, function(x) {
    matrix(x, nrow(data), length(x), byrow = TRUE)
  })
  
  cf_cb <- dlmm_basis(cf_mat, 
                      basis_each = attributes(cross_basis)$basis,
                      interactions = attributes(cross_basis)$interactions)
  
  return(predict(model, cbind.data.frame(data, cf_cb), type = "risk"))
  
}

# simple alternative (no dlm)
# cf_predict <- function(model, data, cross_basis, cf_exposure) {
#   
#   cf <- cbind(pm25_lag0 = rep(cf_exposure$pm25_lag0, nrow(data)),
#               lag0 = rep(cf_exposure$lag0, nrow(data)),
#               drug_lag0 = rep(cf_exposure$lag0, cf_exposure$drug_lag0))
#   
#   return(predict(model, cbind.data.frame(subset(data, select = -c(pm25_lag0, lag0, drug_lag0)), cf), type = "risk"))
#   
# }

### Three-way RERI (confusing)

u1w1x1 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(drug = rep(1, 7),
                                      temp = rep(0.9, 7),
                                      pm25 = rep(12, 7)))

u1w1x0 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(drug = rep(1, 7),
                                      temp = rep(0.9, 7), 
                                      pm25 = rep(8, 7)))

u1w0x1 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(drug = rep(1, 7),
                                      temp = rep(0.7, 7), 
                                      pm25 = rep(12, 7)))

u1w0x0 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                     cf_exposure = list(drug = rep(1, 7),
                                        temp = rep(0.7, 7), 
                                        pm25 = rep(8, 7)))

u0w1x1 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                     cf_exposure = list(drug = rep(0, 7),
                                        temp = rep(0.9, 7), 
                                        pm25 = rep(12, 7)))

u0w1x0 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                     cf_exposure = list(drug = rep(0, 7),
                                        temp = rep(0.9,7), 
                                        pm25 = rep(8, 7)))

u0w0x1 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                     cf_exposure = list(drug = rep(0, 7),
                                        temp = rep(0.7,7), 
                                        pm25 = rep(12, 7)))

u0w0x0 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                     cf_exposure = list(drug = rep(0, 7), 
                                        temp = rep(0.7,7), 
                                        pm25 = rep(8, 7)))

mean(u1w1x1) - mean(u1w1x0) - mean(u1w0x1) - mean(u0w1x1) +
  mean(u1w0x0) + mean(u0w1x0) + mean(u0w0x1) - mean(u0w0x0) # calculates 3-way RERI
