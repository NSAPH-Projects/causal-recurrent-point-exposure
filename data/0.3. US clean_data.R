################################
# Clean data by states
###############################

setwd("~/nsaph_projects/causal-recurrent-point-exposure/data/")
source("0.0. constants.R")

if (!dir.exists(paste0(scratch, "US"))) {
  dir.create(paste0(scratch, "US"))
}

#--------------------
# Load confounders data
#--------------------

conf_dat <- read_feather(paste0(scratch, "conf_exposure.feather"))

# Loop years --------------------------------
files <- list.files(dir_events)
for (year in 2000:2016) {
  cat("\n", year)
  
  # Denominator --------------------------------
  cat("...loading")
  denom <- read_feather(paste0(dir_events, "bene_track_all_", year, ".feather"),
                        col_select = c("zip", "year", "month", "bene_id", "t",
                                       "age_dob", "dob", "dod", "sex",
                                       "race", "dual", "HMO")) # loading beneficiary data
  setDT(denom)
  setkey(denom, bene_id, year, month, zip)
  cat("\n")
  print(denom[, .N, by = year], "\n")
  denom <- denom[t >= 0 & t <= tau] # between ages 65 and admin censoring
  
  
  # + Censoring time  --------------------------------
  cat("...censoring: ")
  # drop data pre-age 65 and non-unique zip
  denom[, n_zip := uniqueN(zip), by = bene_id] # people who have more than 1 zip code, i.e. who may have moved
  cat(denom[n_zip > 1, .N], ">1 zip ") # .N is the number of observations satisfying the filter condition
  denom <- denom[n_zip == 1] # only retain people with only 1 zip code
  
  # drop non FFS
  denom[, FFS := (HMO == 0)] # HMO is Medicare Advantage program that may not have hospitalization info
  cat(denom[!(FFS),.N], "non-FFS ")
  denom <- denom[FFS == TRUE]
  denom <- merge(denom, conf_dat[, .(zip, year, month, covars)],
                 by = c("year", "month", "zip"), all.x = TRUE) # merge ZIP code confounders data
  denom[is.na(covars), covars := FALSE]
  cat(denom[!(covars),.N], "miss. covars ")
  denom <- denom[covars == TRUE]
  denom[, min_t := min(t), by = bene_id]
  cat(denom[min_t != 0,.N], "no t=0 ")
  denom <- denom[min_t == 0] # only include people who start at age 65, i.e., having min_t = 0
  
  # drop times after any gaps in FFS or covariates
  setkey(denom, bene_id, year, month)
  denom[, ordered_t := 0:(.N - 1), by = bene_id]
  cat(denom[ordered_t != t,.N], "missing times \n")
  denom <- denom[ordered_t == t]
  denom[, min_t := min(t), by = bene_id][, max_t := max(t) + 1 - 1e-6, by = bene_id]
  denom[, n_obs := .N, by = bene_id]
  denom[, censor_t := max_t] # last time observed in the data
  
  
  # + Death --------------------------------
  cat("...death...")
  setkey(denom, bene_id, year, month)
  denom[!is.na(dod), dod_yr := year(first(dod)), by = bene_id]
  denom[, dod_yr := mean(dod_yr, na.rm = T), by = bene_id] # year of death
  denom[!is.na(dod), dod_mo := month(first(dod)), by = bene_id]
  denom[, dod_mo := mean(dod_mo, na.rm = T), by = bene_id] # month of death
  denom[!is.na(dod), dod_day := day(first(dod)), by = bene_id]
  denom[, dod_day := mean(dod_day, na.rm = T), by = bene_id] # date of death
  denom[, rm := FALSE]
  denom[dod_yr < year, rm := TRUE] # remove people who are already dead at the time
  denom[dod_yr == year & dod_mo < month, rm := TRUE]
  cat(denom[rm == TRUE,.N], "after death ")
  denom <- denom[rm == FALSE] # remove rows of people who are already dead at the time 
  denom[dod_yr == year & dod_mo == month,
        death_t := t + (dod_day - 0.5) / 
          days_in_month(ymd(paste0(dod_yr, "-", dod_mo, "-", dod_day)))] # calculate time at death
  denom[, death_t := mean(death_t, na.rm = T), by = bene_id] # time of death (t) for each individual 
  denom[, dod_yr := NULL][, dod_mo := NULL][, dod_day := NULL][, rm := NULL] # remove these columns
  
  
  # Drop denom data after censoring/death
  denom[!is.nan(death_t), 
        max_t := ceiling(min(censor_t, death_t)) - 1e-6, by = bene_id]
  denom <- denom[(t < max_t) & (t >= 0)]
  yr_bene_id <- unique(denom$bene_id)
  
  cat(denom[, .N], "obs", length(yr_bene_id), "bene \n")
  
  # write the results
  write_feather(denom, paste0(scratch, "US/denom_yr", year, ".feather"))
  
  # Hospitalization events --------------------------------
  cat("...hospitalizations...")
  hosp <- read_feather(paste0(dir_events, "/hosp_track_all_", year, ".feather")) # hospitalization events
  setDT(hosp)
  setkey(hosp, bene_id, year, month)
  hosp <- hosp[bene_id %in% yr_bene_id]
  # merge with denominator data
  hosp <- merge(hosp, 
                denom[, .(bene_id, censor_t, death_t, year, month)], 
                by = c("bene_id", "year", "month"))
  
  cat(hosp[, .N], "hosp", hosp[, uniqueN(bene_id)], "bene \n")
  
  write_feather(hosp, paste0(scratch, "US/hosp_yr", year, ".feather"))
  
  
} # end loop years
