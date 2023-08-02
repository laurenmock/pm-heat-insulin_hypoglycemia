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

## Collect PM2.5 into a single DF

# setwd("M:/External Users/Shared/To KevinJos/aqdh-pm2-5-o3-no2-concentrations-zipcode-contiguous-us-2000-2016-pm2-5-csv/PM2.5.csv/daily")
# 
# pm_list <- list.files(pattern = ".csv")
# pm_df <- data.frame()
# 
# for (i in 1:length(pm_list)) {
#   
#   pm_tmp <- read_csv(pm_list[[i]])
#   pm_tmp$date <- ymd(substr(pm_list[[i]], 1, 8))
#   pm_df <- rbind(pm_df, pm_tmp)
#   rm(pm_tmp)
#   
# }
# 
# colnames(pm_df) <- tolower(colnames(pm_df))
# write_csv(pm_df, "M:/External Users/Shared/To KevinJos/daily-pm25.csv")

### Begin Analysis ###

# load data
setwd("M:/External Users/KevinJos/code/dlm")
temp_df <- read_sas("M:/External Users/Shared/To KevinJos/DLNM Insulin/BiDirectional/dlnm_insulin_hi_pcnt.sas7bdat")
temp_df <- temp_df[,c(1:2,17:ncol(temp_df))]
pm_df <- read_csv("M:/External Users/Shared/To KevinJos/daily-pm25.csv")
df <- read_sas("M:/External Users/Shared/To KevinJos/DLNM Insulin/BiDirectional/cc_hypo_bi.sas7bdat")

colnames(df) <- tolower(colnames(df))
colnames(temp_df) <- tolower(colnames(temp_df))
colnames(pm_df) <- tolower(colnames(pm_df))

cases <- df[df$cases == 1,] # subset cases

#create control days 
for (i in 1:4) { 
  
  tmp_forw <- df[(df$dates == df$casedate + 7*i) & month(df$dates) == month(df$casedate),] #create copy of cases dataframe to have "forward" controls
  tmp_back <- df[(df$dates == df$casedate - 7*i) & month(df$dates) == month(df$casedate),] #create copy of cases dataframe to have "backwards" controls
  
  # Combine eligible controls 
  tmp = bind_rows(tmp_back, tmp_forw) 
  if (i == 1){
    controls = tmp
  }else{
    controls = bind_rows(controls, tmp)
  }
  
  rm(tmp_back); rm(tmp_forw); rm(tmp); gc()
  
}

#create case/control variables 
controls$case <- 0
cases$case <- 1

# Combine cases and controls
cco <- bind_rows(cases, controls)

## Merging

setDT(cco)
setDT(temp_df)
setDT(pm_df)

# Merge in Temperatures
cco <- temp_df[cco, on = .(bene_id, dates)]

# Merge in PM2.5
pm_df <- select(pm_df, zip, pm25, date)
setnames(pm_df, c("zipcode", "pm25_lag0", "dates"))
setkey(pm_df, zipcode, dates)
pm_df[, paste0("pm25_lag", 1:14) := shift(.SD, 1:14, type = "lag"),
      by = "zipcode", .SDcols = "pm25_lag0"]
cco <- pm_df[cco, on = .(zipcode, dates)]

cco <- cco[!is.na(cco$pm25_lag14),]
cco <- subset(cco, month(dates) %in% c(5,6,7,8,9)) # use summer months

# factor dow, month, and year
cco$year <- factor(year(cco$dates))
cco$month <- factor(month(cco$dates))
cco$dow <- factor(weekdays(cco$dates))

save(cco, file = "M:/External Users/KevinJos/data/dlm/cco.RData")

rm(cases, controls, pm_df, temp_df, df); gc()
