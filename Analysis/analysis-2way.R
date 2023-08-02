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
library(viridis)

# load modules
source("M:/External Users/LaurenMoc/dlm/Functions/dlmm_basis.R")
source("M:/External Users/LaurenMoc/dlm/Functions/dlmm_pen.R")
source("M:/External Users/LaurenMoc/dlm/Functions/dlmm_pred.R")
source("M:/External Users/LaurenMoc/dlm/Functions/dlmm_plot.R")

## Fit DLMM Interaction Model

load("M:/External Users/LaurenMoc/data/cco.RData")

# FOR TESTING - reduce sample size
uidx <- unique(cco$bene_id)
usamp <- sample(uidx, 10000, replace = FALSE)
cco <- subset(cco, bene_id %in% usamp)

x <- list(
  temp = as.matrix(cco[, paste0("lag", 0:14)]),
  pm25 = as.matrix(cco[, paste0("pm25_lag", 0:14)])
)

cb <- dlmm_basis(x, basis_all = list(fun = "cr", df = 4), interactions = "noself") # all, noself, none
formula <- as.formula(paste("case ~ year + month + dow + cluster(bene_id) + ",
                            paste(colnames(cb), collapse = "+")))
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

### RERI

w1x1 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(temp = rep(0.9,15), pm25 = rep(12, 15)))

w1x0 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(temp = rep(0.9,15), pm25 = rep(8, 15)))

w0x1 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(temp = rep(0.7,15), pm25 = rep(12, 15)))

w0x0 <- cf_predict(model = sum_clogit, data = cco, cross_basis = cb,
                   cf_exposure = list(temp = rep(0.7,15), pm25 = rep(8, 15)))

mean(w1x1) - mean(w1x0) - mean(w0x1) + mean(w0x0) # calculate RERI

### Multiplicative Interaction Plots

# predict multiplicative effect holding one variable constant
q <- c(0.1, 0.25, 0.5, 0.75, 0.9)

pred <- lapply(1:length(q), function(p) {
  dlmm_pred(cb, sum_clogit, xscale = sapply(x, sd), bylag = rep(0.5, 4),
            xmarg = sapply(x, quantile, probs = q[p]))
  })

exp <- "pm25"
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
dlmm_plot(pred[[3]], "pm25") + labs(title = "pm25")
dlmm_plot(pred[[3]], "pm25", "temp") + labs(title = "temp/pm25")
