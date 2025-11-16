################################
# Extracting confounder and exposure data
################################

setwd("~/nsaph_projects/causal-recurrent-point-exposure/data/")
source("0.0. constants.R")

#--------------------
# Load confounders based on ZIP code
#--------------------

# (Permanent) ZIP code data
load(paste0(dir_dmork, "zip_code_db.rda")) # ZIP code data based on 2010 census
setDT(zip_code_db) # convert lists and data frames to data.table
zip_code_db[, zip := as.integer(zipcode)] # create a zip column as integer

# Yearly confounders at ZIP code
zip_yr_dat <- 
  rbindlist(lapply(2000:2016, function(y) {
    confounders <- readRDS(paste0(dir_confounders, "aggregate_data_", y, ".rds"))
    setDT(confounders)
    confounders
  })) # extract confounders at ZIP code level from 2000-2016
zip_yr_dat <- zip_yr_dat[, .(zip, year, poverty, popdensity,
                             medianhousevalue, medhouseholdincome, pct_owner_occ,
                             pct_blk, hispanic, education, 
                             smoke_rate, mean_bmi)] # select the columns
zip_yr_dat[, zip := as.integer(zip)] # create a zip column as integer
conf_dat <- merge(zip_yr_dat, zip_code_db[, .(zip, state, lat, lng)], by = "zip") # merge zip data with confounders

# create monthly zip confounders from yearly ones 
year_mo <- data.table(year = rep(2000:2016, each = 12), month = 1:12)
conf_dat <- merge(conf_dat, year_mo, by = "year", allow.cartesian = T)

#--------------------
# Load exposure data
#--------------------

# Load exposure data
for (e in c("pm25", "no2", "ozone", "max_temp", "min_humid")) {
  exp_dat <- read_feather(paste0(scratch, e, "_mo.feather"))
  exp_dat[, year := as.integer(as.character(year))] # convert from factor
  conf_dat <- merge(conf_dat, exp_dat, by = c("year", "month", "zip"), all.x = T)
}
setkey(conf_dat, zip, year, month)
conf_dat <- conf_dat[complete.cases(conf_dat)]
conf_dat[, covars := TRUE]

#--------------------
# save into a data file
#--------------------

write_feather(conf_dat, paste0(scratch, "conf_exposure.feather"))