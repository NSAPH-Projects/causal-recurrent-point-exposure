#################################
# Visualization
#################################

#------------------------
# Processing the output, giving variance
#------------------------

# setting directories
setwd("~/nsaph_projects/causal-recurrent-point-exposure/")
source("code/0. general functions.R")

# Directories
scratch <- "~/nsaph_scratch/causal-recurrent-point-exposure/" 
dir_results <- "~/nsaph_projects/causal-recurrent-point-exposure/results/"

# load data
fn <- "AZ_buffer30"
dat <- read_feather(paste0(scratch, "results/", fn, ".feather"))

# some parameters
t_fits <- 1:24
nsamples <- nrow(dat)
cutoff <- 0

# calculate full variance matrix
vars <- c()
EU <- c()
mu <- c()
for (t in 1:length(t_fits)) {
  # extract info
  clnm <- paste0(c("mu_1_", "mu_0_", "eta_1_", "eta_0_"), t_fits[t])
  est <- colSums(dat[, ..clnm]) / nsamples
  mu <- c(mu, est)
  
  # get variance
  vars <- cbind(vars, 
                dat[, get(paste0("mu_1_", t_fits[t]))] - est[1],
                dat[, get(paste0("mu_0_", t_fits[t]))] - est[2],
                dat[, get(paste0("eta_1_", t_fits[t]))] - est[3],
                dat[, get(paste0("eta_0_", t_fits[t]))] - est[4])
  EU <- c(EU, -nsamples, -nsamples, -nsamples, -nsamples)
}

EU <- diag(EU)
vars <- solve(EU) %*% t(vars) %*% vars %*% t(solve(EU))

# output object
out <- data.frame(size = double(),
                  method = character(),
                  cutoff = character(),
                  time = double(),
                  estimand = character(),
                  est = double(),
                  sd = double())

# calculate estimates and their variances
cnt <- 0
for (t in 1:length(t_fits)) {
  
  # eta and mu
  out[cnt+1:4,] <- list(size = nsamples,
                        method = "One-step AIPW",
                        cutoff = cutoff,
                        time = t_fits[t],
                        estimand = c("mu_1", "mu_0", "eta_1", "eta_0"),
                        est = mu[4*(t-1)+1:4],
                        sd = sqrt(diag(vars[4*(t-1)+1:4, 4*(t-1)+1:4])))
  
  # log of while alive strategy
  dG <- rep(0, nrow(vars))
  dG[1+4*(t-1)] <- 1/(mu[1+4*(t-1)]) #dmu1
  dG[2+4*(t-1)] <- -1/(mu[1+4*(t-1)]) #dmu0
  dG[3+4*(0:(t-1))] <- -1/(1+sum(mu[3+4*(0:(t-1))])) #deta1
  dG[4*(1:t)] <- 1/(1+sum(mu[4*(1:t)])) #deta0
  
  out[cnt+5,] <- list(size = nsamples,
                      method = "One-step AIPW",
                      cutoff = cutoff,
                      time = t_fits[t],
                      estimand = "log causal",
                      est = log((mu[1+4*(t-1)]*(1+sum(mu[4*(1:t)])))/(mu[2+4*(t-1)]*(1+sum(mu[3+4*(0:(t-1))])))),
                      sd = sqrt(t(as.matrix(dG)) %*% vars %*% as.matrix(dG)))
  cnt <- cnt+5
}

save(list = c("out"),
     file = paste0(dir_results, "results_", fn, ".RData"))

#------------------------
# visualization
#------------------------

# setting directories
setwd("~/nsaph_projects/causal-recurrent-point-exposure/")
dir_results <- "results/"
dir_plots <- "figures/"

# Load libraries
library(dplyr)
library(ggplot2)

# load analysis results
fn <- "AZ_buffer30"
load(paste0(dir_results, "results_", fn, ".RData"))

# calculate confidence intervals
out$varG <- 0
out$lbG <- 0
out$ubG <- 0
out$lb <- 0
out$ub <- 0

# log transformation for mu
out$varG[out$estimand %in% c("mu_1", "mu_0")] <- ((1/out$est[out$estimand %in% c("mu_1", "mu_0")])^2)*((out$sd[out$estimand %in% c("mu_1", "mu_0")])^2)
out$lbG[out$estimand %in% c("mu_1", "mu_0")] <- log(out$est[out$estimand %in% c("mu_1", "mu_0")]) - 1.96 * sqrt(out$varG[out$estimand %in% c("mu_1", "mu_0")]) 
out$ubG[out$estimand %in% c("mu_1", "mu_0")] <- log(out$est[out$estimand %in% c("mu_1", "mu_0")]) + 1.96 * sqrt(out$varG[out$estimand %in% c("mu_1", "mu_0")]) 
out$lb[out$estimand %in% c("mu_1", "mu_0")] <- exp(out$lbG[out$estimand %in% c("mu_1", "mu_0")])
out$ub[out$estimand %in% c("mu_1", "mu_0")] <- exp(out$ubG[out$estimand %in% c("mu_1", "mu_0")])

# log-log transformation for eta
out$varG[out$estimand %in% c("eta_1", "eta_0")] <- ((1/(out$est[out$estimand %in% c("eta_1", "eta_0")]*log(out$est[out$estimand %in% c("eta_1", "eta_0")])))^2)*((out$sd[out$estimand %in% c("eta_1", "eta_0")])^2)
out$lbG[out$estimand %in% c("eta_1", "eta_0")] <- log(-log(out$est[out$estimand %in% c("eta_1", "eta_0")])) - 1.96 * sqrt(out$varG[out$estimand %in% c("eta_1", "eta_0")]) 
out$ubG[out$estimand %in% c("eta_1", "eta_0")] <- log(-log(out$est[out$estimand %in% c("eta_1", "eta_0")])) + 1.96 * sqrt(out$varG[out$estimand %in% c("eta_1", "eta_0")]) 
out$lb[out$estimand %in% c("eta_1", "eta_0")] <- exp(-exp(out$lbG[out$estimand %in% c("eta_1", "eta_0")]))
out$ub[out$estimand %in% c("eta_1", "eta_0")] <- exp(-exp(out$ubG[out$estimand %in% c("eta_1", "eta_0")]))

# while alive strategy
out$varG[out$estimand %in% c("log causal")] <- (out$sd[out$estimand %in% c("log causal")])^2
out$lbG[out$estimand %in% c("log causal")] <- out$est[out$estimand %in% c("log causal")] - 1.96 * sqrt(out$varG[out$estimand %in% c("log causal")]) 
out$ubG[out$estimand %in% c("log causal")] <- out$est[out$estimand %in% c("log causal")] + 1.96 * sqrt(out$varG[out$estimand %in% c("log causal")]) 
out$lb[out$estimand %in% c("log causal")] <- exp(out$lbG[out$estimand %in% c("log causal")])
out$ub[out$estimand %in% c("log causal")] <- exp(out$ubG[out$estimand %in% c("log causal")])
out$est[out$estimand %in% c("log causal")] <- exp(out$est[out$estimand %in% c("log causal")])

# adding some categories
# out$which <- ifelse(out$estimand %in% c("eta_1", "eta_0"), "eta(t)", "mu(t)")
out$which <- ifelse(out$estimand %in% c("eta_1", "eta_0"), "Survival Probability", "Expected # of CVD-related Hospitalizations")
out$which[out$estimand %in% c("log causal")] <- "while alive causal"
# out$counterfactual <- ifelse(out$estimand %in% c("eta_1", "mu_1"), "1", "0")
out$counterfactual <- ifelse(out$estimand %in% c("eta_1", "mu_1"), "High", "Low")
out$counterfactual[out$estimand %in% c("log causal")] <- NA

# eta and mu plot
png(filename = paste0(dir_plots, fn, "_mu_eta.png"), 
    width = 900, height = 500, res = 100)
out %>%
  filter(estimand %in% c("mu_1", "mu_0", "eta_1", "eta_0")) %>%
  ggplot(aes(x = time, y = est, color = counterfactual)) +
  geom_line() +
  geom_line(aes(x = time, y = ub), linetype = "dashed") +
  geom_line(aes(x = time, y = lb), linetype = "dashed") + 
  facet_wrap(vars(which), scales = "free_y") + 
  labs(x = "Month (t)", y = "Estimates", color = "Exposure to PM2.5") + 
  theme(legend.position="top")
dev.off()

# while alive strategy plot
png(filename = paste0(dir_plots, fn, "_causal.png"), 
    width = 900, height = 500, res = 100)
out %>%
  filter(estimand %in% c("log causal")) %>%
  ggplot(aes(x = time, y = est)) +
  geom_hline(aes(yintercept = 1), linetype = "solid", color = "darkgrey") +
  geom_line() +
  geom_line(aes(x = time, y = ub), linetype = "dashed") +
  geom_line(aes(x = time, y = lb), linetype = "dashed") +
  labs(x = "Month (t)", y = "Estimates", title = "While-alive Causal Composite Estimates") + 
  theme(legend.position="top")
dev.off()

#################################
# Exposure on a map
#################################

library(ggplot2)
library(maps)
library(mapdata)
library(scales)
#library(tigris)

setwd("~/nsaph_projects/causal-recurrent-point-exposure/")
.libPaths("/n/User/R/ifxrstudio/RELEASE_3_18")

library(data.table)
library(arrow) # C++ dplyr
library(fst) # fast serialization of dataframe
library(lubridate) # process time

state <- map_data("state")
arizona <- subset(state, region=="arizona")

counties <- map_data("county")
arizona_county <- subset(counties, region=="arizona")

scratch <- "~/nsaph_scratch/causal-recurrent-point-exposure/" 
dir_results <- "results/"
dir_plots <- "figures/"

state_code <- "AZ" #one of CA, AZ, NY
loc_name <- paste0("states/", state_code, "/")
dat <- read_feather(paste0(scratch, loc_name, "data_for_analysis.feather"))
setDT(dat)
setkey(dat, cohort)
quantile(dat$pm25, probs = c(0.25, 0.35, 0.5, 0.65, 0.75))

# summaries
as.vector(dat[, .N, by = cohort][,2]) #N
as.vector(dat[sex == 1, .N, by = cohort][,2]) #male
as.vector(dat[sex == 2, .N, by = cohort][,2]) #female
as.vector(dat[sex == 0, .N, by = cohort][,2]) #other sex

as.vector(dat[race == 1, .N, by = cohort][,2]) #white
as.vector(dat[race == 2, .N, by = cohort][,2]) #black
as.vector(dat[race == 4, .N, by = cohort][,2]) #asian
as.vector(dat[race == 5, .N, by = cohort][,2]) #asian
as.vector(dat[race == 6, .N, by = cohort][,2]) #north native american
as.vector(dat[race %in% c(0,3), .N, by = cohort][,2]) #unknown/other

as.vector(dat[death_t < 48.99, .N, by = cohort][,2])
as.vector(dat[censor_t < 48.99, .N, by = cohort][,2])

conf_dat <- read_feather(paste0(scratch, "conf_exposure.feather"))
setDT(dat)
conf_dat <- conf_dat[state == "AZ"]
conf_dat <- conf_dat[, .(zip, lng, lat)]
conf_dat <- conf_dat[!duplicated(conf_dat)]

dat <- merge(dat, conf_dat,
             by = c("zip"), all.x = TRUE)
dat[, lng_jitter := lng + rnorm(n = nrow(dat), mean = 0, sd = 0.1)]
dat[, lat_jitter := lat + rnorm(n = nrow(dat), mean = 0, sd = 0.1)]
dat <- as.data.frame(dat)

arizona_cities <- data.frame(city = c("Phoenix", "Tucson", "Yuma"), 
                             lng = c(-112.074036, -110.911789, -114.650398), 
                             lat = c(33.448376, 32.253460, 32.698437))

AZ_map <- ggplot(data=arizona, mapping=aes(x=long, y=lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black", fill="gray") + 
  geom_point(aes(x=lng_jitter, y=lat_jitter, col = pm25), 
             data = dat, size = 1, alpha=.1, 
             inherit.aes = FALSE) +
  scale_color_gradientn(colours = c("red", "red", "yellow", "green", "green"), 
                        values = rescale(c(quantile(dat$pm25, probs = c(1, 0.75, 0.5, 0.25, 0))), from = range(dat$pm25))) +
  geom_polygon(data=arizona_county, fill=NA, color="white") + 
  geom_polygon(color="black", fill=NA) +
  geom_point(data = arizona_cities, aes(x=lng, y=lat), inherit.aes = FALSE) +
  geom_text(data = arizona_cities, aes(x=lng, y=lat, label=city), size = 3, hjust=0, vjust=-1, inherit.aes = FALSE) +
  ggtitle('Arizona map and PM2.5 exposure') + 
  labs(col = bquote("PM2.5 ("*mu*"g/m"^.(3)*")")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
AZ_map

png(filename = paste0(dir_plots, "AZ_PM25.png"), 
    width = 600, height = 500, res = 100)
AZ_map
dev.off()
