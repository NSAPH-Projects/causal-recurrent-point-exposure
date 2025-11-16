source("0.0. constants.R")
library(xtable)

# Load zip covars
conf_dat <- read_feather(paste0(scratch, loc_name, "/conf_dat.feather"))
setDT(conf_dat)

# Load t = 0
dat <- read_feather(paste0(scratch, loc_name, "/time_dat/time_0_combined.feather"))
dat[, year := as.integer(as.character(year))]
dat[, month := as.integer(as.character(month))]
dat <- merge(dat, conf_dat, by = c("zip", "year", "month"), all.x = TRUE)

# Individual characteristics
dat[, .N]
dat[, sum(max_t)]
xtable(t(dat[, .N, by = year][order(year)]))
out <- rbind(
dcast(dat[, .(.N, perc = .N / dat[, .N]), by = .(year, sex)][order(year, sex)],
      sex ~ year, value.var = "N"),
dcast(dat[race %in% c(1, 2, 4, 5), .(.N, perc = .N / dat[, .N]), by = .(year, race)][order(year, race)],
      race ~ year, value.var = "N"),
dcast(dat[!(race %in% c(1, 2, 4, 5)), .(.N, perc = .N / dat[, .N], race = 1), by = .(year)][order(year)],
       race ~ year, value.var = "N"),
dcast(dat[, .(.N, perc = .N / dat[, .N]), by = .(year, dual)][order(year, dual)],
      dual ~ year, value.var = "N"), fill = TRUE
)
cols <- as.character(2000:2010)
table <- xtable(out[, ..cols])
print(table, include.rownames = F)

summary(dat[, pm25]); hist(dat$pm25)
summary(dat[max_t < tau, max_t])
dat[death_t < Inf, .(.N, mean(death_t))]
dat[censor_t < tau & censor_t < death_t, .(.N, mean(censor_t), quantile(censor_t, c(.25, .75)))]
dat[, sum(max_t)] # person-months of data
paste0(dcast(dat[death_t < tau & death_t <= censor_t, .N, by = .(year)][order(year)],
      . ~ year, value.var = "N"), collapse = " & ")
paste0(dcast(dat[censor_t < tau, .N, by = .(year)][order(year)],
      . ~ year, value.var = "N"), collapse = " & ")

# Survival plots
library(survival)
sfit <- survfit(Surv(dat$max_t - start, dat$death_t == dat$max_t) ~ 1)
plot(sfit, ylim = c(.9, 1), conf.int = T)
cfit <- survfit(Surv(dat$max_t - start, dat$censor_t == dat$max_t & dat$max_t < tau) ~ 1)
plot(cfit, ylim = c(.8, 1), conf.int = T)

# Neighborhood characteristics
dat[, .(mean = mean(poverty), sd = sd(poverty))]
dat[, .(mean = mean(popdensity), sd = sd(popdensity))]
dat[, .(mean = mean(medianhousevalue), sd = sd(medianhousevalue))]
dat[, .(mean = mean(medhouseholdincome), sd = sd(medhouseholdincome))]
dat[, .(mean = mean(pct_owner_occ), sd = sd(pct_owner_occ))]
dat[, .(mean = mean(pct_blk), sd = sd(pct_blk))]
dat[, .(mean = mean(pct_asian), sd = sd(pct_asian))]
dat[, .(mean = mean(pct_native), sd = sd(pct_native))]
dat[, .(mean = mean(pct_white), sd = sd(pct_white))]
dat[, .(mean = mean(hispanic), sd = sd(hispanic))]
dat[, .(mean = mean(education), sd = sd(education))]
dat[, .(mean = mean(smoke_rate), sd = sd(smoke_rate))]
dat[, .(mean = mean(mean_bmi), sd = sd(mean_bmi))]
dat[, .(.N, perc = .N/dat[, .N]), by = state][order(perc, decreasing = T)]

# Visualize entry/exit
library(ggplot2)
dat[, floor_max_t := floor(max_t)]
dat[, entry := .N, by = .(year, month)]
dat[, exit := .N, by = .(year, month, floor_max_t)]
dat[order(floor_max_t), exit_cum := cumsum(exit), by = .(year, month, floor_max_t)]
dat[, .(.N), by = .(year, month, floor_max_t)][order(year, month, floor_max_t)]
ggplot(dat[, .(.N), by = .(year, month)
           ][order(year, month)],
       aes(x = year, y = N)) +
  geom_point()


# Neighborhood confounding avg corr
zip_yr_dat <- 
  rbindlist(lapply(2000:2016, function(y) {
    confounders <- readRDS(paste0(dir_confounders, "aggregate_data_", y, ".rds"))
    setDT(confounders)
    confounders
  }))
zip_yr_dat <- zip_yr_dat[, .(zip, year, poverty, popdensity,
                             medianhousevalue, medhouseholdincome, pct_owner_occ,
                             pct_blk, hispanic, education, 
                             smoke_rate, mean_bmi)]
zip_yr_dat[, zip := as.integer(zip)]
mean(sapply(names(zip_yr_dat)[-c(1:2)], function(n) {
  var_by_year <- dcast(zip_yr_dat, zip ~ year, value.var = n)[, -1]
  diag(cor(var_by_year[complete.cases(var_by_year)])[-1, -17])
}))
