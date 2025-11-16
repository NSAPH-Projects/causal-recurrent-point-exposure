################################
# Combine data
################################

setwd("~/nsaph_projects/causal-recurrent-point-exposure/data/")
source("0.0. constants.R")

# state_code <- "AL" #one of AL, AZ, CA, GA, ID, MT, NV, NY, UT, WA
state_code <- as.character(Sys.getenv('STATE_CODE'))
exposure_code <- "pm25" #one of pm25, no2, ozone, min_humid, max_temp
end_year <- floor(2016 - (tau/12))

if (!dir.exists(paste0(scratch, "states/"))) {
  dir.create(paste0(scratch, "states/"))
}

loc_name <- paste0("states/", state_code, "/")
if (!dir.exists(paste0(scratch, loc_name))) {
  dir.create(paste0(scratch, loc_name))
  dir.create(paste0(scratch, loc_name, "time_dat/"))
}

#-----------------------
# Loop over cohorts
#----------------------- 

conf_dat <- read_feather(paste0(scratch, "conf_exposure.feather"))
setDT(conf_dat)
setkey(conf_dat, year, month, zip)
conf_dat <- conf_dat[state %in% state_code]
zip_list <- conf_dat[, zip]

for (y in 2000:end_year) {
  cat("\n", y, "")
  
  # Load cohort data
  dat <- read_feather(paste0(scratch, "US/denom_yr", y, ".feather"))
  setDT(dat)
  setkey(dat, zip, year, month, bene_id)
  
  # # Subsample data (10%)
  # all_bene <- unique(dat$bene_id)
  # set.seed(y)
  # bene_samp <- sample(all_bene, ceiling(length(all_bene) * .1))
  # dat <- dat[bene_id %in% bene_samp]
  
  # Check living in state/location of analysis
  dat <- dat[, cohort := y]
  dat <- dat[zip %in% zip_list]
  dat <- dat[max_t > start, .(zip, year, month, bene_id, dob, sex, race, dual, cohort, t,
                death_t, censor_t, max_t)]
  dat <- merge(dat, conf_dat, by = c("year", "month", "zip"), all.x = TRUE)
  setkey(dat, bene_id, t)
  
  # Loop over times, save each time point as separate file
  for (time in 0:tau) {
    cat(time, "")
    time_dat <- dat[t == time]
    write_feather(time_dat, 
                  paste0(scratch, loc_name, "time_dat/time_", time, "_year_", y, ".feather"))
  }
}

#-----------------------
# Combine all years for each time t
#-----------------------

files <- list.files(paste0(scratch, loc_name, "time_dat/"))
for (time in 0:tau) {
  cat("\nTime", time, "")
  dat <- rbindlist(lapply(2000:end_year, function(y) {
    read_feather(paste0(scratch, loc_name, "time_dat/time_", time, "_year_", y, ".feather"))
  }))
  setkey(dat, bene_id)
  cat(dat[, .N], "obs")
  write_feather(dat, paste0(scratch, loc_name, "time_dat/time_", time, "_combined.feather"))
  unlink(paste0(scratch, loc_name, "time_dat/", 
                files[grep(paste0("time_", time, "_year_"), files, fixed = TRUE)]))
  if (time == 0) {
    all_bene <- dat$bene_id
    save(all_bene, file = paste0(scratch, loc_name, "all_bene.rda"))
    rm(all_bene)
  }
}
rm(dat,time_dat)

#-----------------------
# Hospitalization data
#-----------------------

# Load list of bene_id
load(paste0(scratch, loc_name, "all_bene.rda"))

# Load events
event_dat <- rbindlist(lapply(2000:end_year, function(y) {
  read_feather(paste0(scratch, "US/hosp_yr", y, ".feather"))
}))
event_dat <- event_dat[(bene_id %in% all_bene) & (t >= 0) & (t < tau + 1)]
setkey(event_dat, bene_id, t)

# CVD any diagnosis code 
cvd_any <- 
  event_dat[((atrial_fibdx_25) | (cardiac_arrestdx_25) | (acute_myocard_infarcdx_25)), 
            .(bene_id, t)]

# drop events within 2 days of prev hosp
for (j in 1:3) { 
  cvd_any[, t_lag := shift(t, 1, type = "lag"), by = bene_id] # NA returned because we are subsetting by bene_id
  cvd_any[, t_diff := t - t_lag]
  cvd_any <- cvd_any[is.na(t_diff) | t_diff > 2/30]
}
cvd_any[, t_lag := NULL][, t_diff := NULL]

# 
cvd_any <- rbind(cvd_any[, .(bene_id, t)], 
                 data.table(bene_id = all_bene, t = -1.0))
cvd_any <- cvd_any[!duplicated(cvd_any)]
setnames(cvd_any, "t", "time")
setkey(cvd_any, bene_id, time)
cvd_any[, t := floor(time)]
cvd_any[, event_num := 0]
cvd_any[t >= 0, event_num := 1:.N, by = bene_id]
write_feather(cvd_any[t != -1], paste0(scratch, loc_name, "cvd_any_event_times.feather"))

N_events <- cvd_any[, .(N = sum(time >= 0)), by = bene_id]
write_feather(N_events, paste0(scratch, loc_name, "cvd_any_N_events.feather"))
