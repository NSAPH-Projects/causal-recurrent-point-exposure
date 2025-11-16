################################
# Monthly exposure data by zip code
################################

setwd("~/nsaph_projects/causal-recurrent-point-exposure/data/")
source("0.0. constants.R")

# Exposures directory
dir_exposures <- "~/nsaph_data/exposures/exposure/"
dir_temp <- "~/nsaph_data/gridmet_flat/"

# Daily air pollution data, edit inputs to capture specific range of days
year_range <- 2000:2016
month_range <- 1:12 #formatC(1:5, width = 2, format = "d", flag = "0")

# PM25
file <- "pm25/PM25_v2" # pm25/PM25_v2, no2/NO2_v2, ozone/O3_v2
pm25_dat <- rbindlist(lapply(year_range, function(y) {
  cat(y, "")
  rbindlist(lapply(month_range, function(m) { # rbindlist: rbind data.table or data.frame from a list 
    rbindlist(lapply(1:days_in_month(ymd(paste0(y, "/", m, "/1"))), function(d) { # days_in_month: number of days in the month component of the date-time object, ymd: transform dates stored in character and numeric vectors to Date objects 
      dat <- readRDS(paste0(dir_exposures, file, "/daily/",
                            y, sprintf("%02d", m), sprintf("%02d", d), ".rds")) # load the daily data file
      setDT(dat)
      dat$year <- y # add year and month to the data
      dat$month <- m
      dat
    }))[, .(pm25 = mean(pm25)), by = .(year, month, ZIP)] # group means
  }))
}))
setkey(pm25_dat, year, month, ZIP) # kinda sorting by columns
pm25_dat[, zip := as.integer(ZIP)][, ZIP := NULL] # change column name from ZIP -> zip
write_feather(pm25_dat, paste0(scratch, "pm25_mo.feather"))

# no2
file <- "no2/NO2_v2" # pm25/PM25_v2, no2/NO2_v2, ozone/O3_v2
no2_dat <- rbindlist(lapply(year_range, function(y) { # similar processing
  cat(y, "")
  rbindlist(lapply(month_range, function(m) {
    rbindlist(lapply(1:days_in_month(ymd(paste0(y, "/", m, "/1"))), function(d) {
      dat <- readRDS(paste0(dir_exposures, file, "/daily/",
                            y, sprintf("%02d", m), sprintf("%02d", d), ".rds"))
      setDT(dat)
      dat$year <- y
      dat$month <- m
      
      dat
    }))[, .(no2 = mean(no2)), by = .(year, month, ZIP)]
  }))
}))
setkey(no2_dat, year, month, ZIP)
no2_dat[, zip := as.integer(ZIP)][, ZIP := NULL]
write_feather(no2_dat, paste0(scratch, "no2_mo.feather"))

# ozone
file <- "ozone/O3_v2" # pm25/PM25_v2, no2/NO2_v2, ozone/O3_v2
ozone_dat <- rbindlist(lapply(year_range, function(y) { # similar processing
  cat(y, "")
  rbindlist(lapply(month_range, function(m) {
    rbindlist(lapply(1:days_in_month(ymd(paste0(y, "/", m, "/1"))), function(d) {
      dat <- readRDS(paste0(dir_exposures, file, "/daily/",
                            y, sprintf("%02d", m), sprintf("%02d", d), ".rds"))
      setDT(dat)
      dat$year <- y
      dat$month <- m
      dat
    }))[, .(ozone = mean(ozone)), by = .(year, month, ZIP)]
  }))
}))
setkey(ozone_dat, year, month, ZIP)
ozone_dat[, zip := as.integer(ZIP)][, ZIP := NULL]
write_feather(ozone_dat, paste0(scratch, "ozone_mo.feather"))

# Daily max temperature, min humidity, heat index data
year_range <- 2000:2016

# load gridmet max temp
max_temp <- rbindlist(lapply(year_range, function(y) {
  cat(y, "")
  load(paste0(dir_temp, 
              "maximum_air_temperature/", y,
              "_maximum_air_temperature_by_zip.RData"))
  df$zip <- as.integer(rownames(df))
  setDT(df)
  max_temp <- melt(df, id.vars = "zip", variable.name = "date", value.name = "max_temp") # wide to long table
  max_temp[, date := ymd(date)]
  max_temp[, year := year(date)]
  max_temp[, month := month(date)]
  max_temp[, .(max_temp = mean(max_temp)), by = .(zip, year, month)]
}))
write_feather(max_temp, paste0(scratch, "max_temp_mo.feather"))

# load gridmet min humidity
min_humid <- rbindlist(lapply(year_range, function(y) { # similar to max_temp 
  cat(y, "")
  load(paste0(dir_temp, 
              "minimum_relative_humidity/", y,
              "_minimum_relative_humidity_by_zip.RData"))
  df$zip <- as.integer(rownames(df))
  setDT(df)
  min_humid <- melt(df, id.vars = "zip", variable.name = "date", value.name = "min_humid")
  min_humid[, date := ymd(date)]
  min_humid[, year := year(date)]
  min_humid[, month := month(date)]
  min_humid[, .(min_humid = mean(min_humid)), by = .(zip, year, month)]
}))
write_feather(min_humid, paste0(scratch, "min_humid_mo.feather"))
