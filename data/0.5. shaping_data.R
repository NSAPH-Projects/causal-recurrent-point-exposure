################################
# Shaping data to the desired format
################################

setwd("~/nsaph_projects/causal-recurrent-point-exposure/data/")
source("0.0. constants.R")

# state_code <- "AL" #one of AL, AZ, CA, GA, ID, MT, NV, NY, UT, WA
state_code <- as.character(Sys.getenv('STATE_CODE')) 
exposure_code <- "pm25" #one of pm25, no2, ozone, min_humid, max_temp
loc_name <- paste0("states/", state_code, "/")

#-----------------------
# Combine all data
#-----------------------

# combine all data until start date
for (time in 0:(start-1)) {
  
  tmp_dat <- read_feather(paste0(scratch, loc_name, "time_dat/time_", time, "_combined.feather"))
  tmp_dat <- tmp_dat[, -c("year", "month", "dual", "max_t", "t", "covars", "lat", "lng")]
  
  if (time == 0) {
    dat <- tmp_dat
  }
  
  dat <- rbind(dat, tmp_dat) 
  dat <- dat[, .(zip, dob, sex, race, cohort, death_t, censor_t, state, 
                 poverty = sum(poverty), 
                 popdensity = sum(popdensity), 
                 medianhousevalue = sum(medianhousevalue), 
                 medhouseholdincome = sum(medhouseholdincome), 
                 pct_owner_occ = sum(pct_owner_occ), 
                 pct_blk = sum(pct_blk), 
                 hispanic = sum(hispanic), 
                 education = sum(education), 
                 smoke_rate = sum(smoke_rate), 
                 mean_bmi = sum(mean_bmi), 
                 pm25 = sum(pm25), 
                 no2 = sum(no2), 
                 ozone = sum(ozone), 
                 max_temp = sum(max_temp), 
                 min_humid = sum(min_humid)), by = bene_id]
  dat <- unique(dat, by = "bene_id")
}
rm(tmp_dat)

# take the average
dat <- dat[, .(zip, dob, sex, race, cohort, death_t, censor_t, state, 
               poverty = poverty/start, 
               popdensity = popdensity/start, 
               medianhousevalue = medianhousevalue/start, 
               medhouseholdincome = medhouseholdincome/start, 
               pct_owner_occ = pct_owner_occ/start, 
               pct_blk = pct_blk/start, 
               hispanic = hispanic/start, 
               education = education/start, 
               smoke_rate = smoke_rate/start, 
               mean_bmi = mean_bmi/start, 
               pm25 = pm25/start, 
               no2 = no2/start, 
               ozone = ozone/start, 
               max_temp = max_temp/start, 
               min_humid = min_humid/start), by = bene_id]

#change names according to the mock data format
dat[, x := min(death_t, censor_t, na.rm = TRUE) - start, by = bene_id] 
dat[, delta := ifelse(!is.na(death_t), as.numeric(death_t <= censor_t), 0), by = bene_id]

write_feather(dat, paste0(scratch, loc_name, "indv_data.feather"))

#-----------------------
# Creating events data
#-----------------------

dat <- read_feather(paste0(scratch, loc_name, "indv_data.feather"))
eval_times <- seq(0.5, tau-start+1, by = 0.5)
dat[, paste0("NE_", eval_times) := 0]

events <- read_feather(paste0(scratch, loc_name, "cvd_any_event_times.feather"))
events <- events[time >= start]
events[, time_to_start := time - start]

progress <- 0
for (i in 1:nrow(events)) {
  tmp <- (eval_times >= as.numeric(events[i,"time_to_start"]))
  nm <- paste0("NE_", eval_times[eval_times >= as.numeric(events[i,"time_to_start"])])
  dat[bene_id == events[i,"bene_id"], (nm) := eval(as.symbol(nm)) + 1]
  
  # cat(i, "/", nrow(events), "\r")
  if ((i/nrow(events)*10) >= (progress+1)) {
    progress <- progress + 1
    cat(progress*10, "% completed \n")
  }
}

summary(dat[,(NE_25)])

write_feather(dat, paste0(scratch, loc_name, "data_for_analysis.feather"))
