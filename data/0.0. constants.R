#-----------------
# Library
#-----------------

setwd("~/nsaph_projects/causal-recurrent-point-exposure/data/")
.libPaths("/n/User/R/ifxrstudio/RELEASE_3_18")

library(data.table)
library(arrow) # C++ dplyr
library(fst) # fast serialization of dataframe
library(lubridate) # process time

#-----------------
# Directory
#-----------------

scratch <- "~/nsaph_scratch/causal-recurrent-point-exposure/" 
dir_confounders <- "~/nsaph_projects/analytic/confounders/"
dir_temps <- "~/nsaph_projects/analytic/temperature_seasonal_zipcode/"
dir_events <- "~/nsaph_projects/analytic/recurrent_events/"
dir_dmork <- "~/nsaph_projects/pm25-adrd-cvd-recurrent_events/CVD_recur/"
dir <- "~/nsaph_projects/causal-recurrent-point-exposure/data/"

#-----------------
# Some other constants
#-----------------

#some other constants
tau <- 48.0 # Administrative censoring time # Why?
start <- lags <- 24 # Study start (# mo after age 65 for baseline) 

# state_code <- 
#   c("AL", "AZ", "AR", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "ID", "IL", "IN", "IA",
#     "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ",
#     "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT",
#     "VT", "VA", "WA", "WV", "WI", "WY")
# loc_name <- "USA"

#-----------------
# Confounders names
#-----------------

#ZIP code level??
ps_confounders <- 
  c("poverty", 
    "popdensity", 
    "medianhousevalue", 
    "medhouseholdincome", 
    "pct_owner_occ",      # ???
    "pct_blk", 
    "hispanic",           # pct_hispanic?
    "education",     
    "smoke_rate",    
    "mean_bmi", 
    "ozone",      
    "no2",        
    "max_temp",           # Why not min?
    "min_humid")          # Why not max?
ps_confounders_sl <- c(ps_confounders, 
                       "state", 
                       "year", 
                       "month", 
                       "lat",          # Of ZIP code?
                       "lng")          # Of ZIP code?
rho_confounders_sl <- c(ps_confounders_sl, 
                        "sex", 
                        "race", 
                        "dual",        # ???
                        "prev_events") # ???

#-----------------
# Analysis head
#-----------------

# cat("Running analyses for:", loc_name)
# cat("\nLags =", lags, "\nStart =", start, "\nEnd = ", tau, "\n")
# 
# if (!dir.exists(paste0(scratch, loc_name))) {
#   dir.create(paste0(scratch, loc_name))
#   dir.create(paste0(getwd(), "/results/", loc_name))
# }