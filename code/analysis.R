#################################
# Analysis for data
#################################

#-------------------
# Directories
#-------------------

# setting directories
setwd("~/nsaph_projects/causal-recurrent-point-exposure/code")

# Load source code files
source("0. general functions.R")

# Directories
scratch <- "~/nsaph_scratch/causal-recurrent-point-exposure/" 
dir <- "~/nsaph_projects/causal-recurrent-point-exposure/code"
dir_results <- "~/nsaph_projects/causal-recurrent-point-exposure/results/"

# s <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # simulation replicate (1-1000)
state_code <- "AZ" #one of CA, AZ, NY
# state_code <- as.character(Sys.getenv('STATE_CODE'))
exposure_code <- "pm25" #one of pm25, no2, ozone, min_humid, max_temp
loc_name <- paste0("states/", state_code, "/")

# if (!dir.exists(paste0(dir_results, state_code))) {
#   dir.create(paste0(dir_results, state_code))
# }

#-------------------
# Prepare data
#-------------------

# load the dataset
dat <- read_feather(paste0(scratch, loc_name, "data_for_analysis.feather"))
setDT(dat)

summary(dat$pm25)
# Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.377   4.753   6.737   7.011   9.106  22.232 
hist(dat$pm25, xlab = "pm25", main = "Histogram of pm25")
abline(v = median(dat$pm25), col = "red", lty = 2)

summary(dat$no2)
summary(dat$ozone)
summary(dat$max_temp)
summary(dat$min_humid)

dat[, a := as.numeric(pm25 >= median(pm25))]

# dat[, a := ifelse(pm25 <= quantile(pm25, 0.35), 0, ifelse(pm25 >= quantile(pm25, 0.65), 1, NA))]
# dat <- dat[!is.na(a),]

# dat[, a := ifelse(pm25 <= quantile(pm25, 0.25), 0, ifelse(pm25 >= quantile(pm25, 0.75), 1, NA))]
# dat <- dat[!is.na(a),]

sum(dat$a == 1)
sum(dat$a == 0)

# some descriptives
summary(dat$pm25)
mean(dat$a) #0.5
#                                AZ        
mean(dat$delta)                  # 0.033   
mean(dat$delta[dat$a == 1])      # 0.032   
mean(dat$delta[dat$a == 0])      # 0.034   

mean(dat$NE_25 > 0)              # 0.038   
mean(dat$NE_25[dat$a == 1] > 0)  # 0.040   
mean(dat$NE_25[dat$a == 0] > 0)  # 0.036   

mean(dat$NE_25)                  # 0.058   
mean(dat$NE_25[dat$a == 1])      # 0.062   
mean(dat$NE_25[dat$a == 0])      # 0.054   
sum(dat$NE_25)

# data preparation
summary(dat$censor_t)
dat[, censor_t_p := ifelse(censor_t < 48.99, censor_t, NaN)]
dat[, death_t_p := ifelse(!is.nan(death_t), death_t, 48.99999)]
dat[, x := ifelse(!is.nan(censor_t_p), pmin(death_t_p, censor_t_p), death_t_p) - 24]
dat[, delta := ifelse(!is.nan(censor_t_p), as.numeric(death_t_p <= censor_t_p), 1)]

# factor covariates
dat[, cohort := as.factor(cohort)]
dat[, sex := as.factor(sex)]
dat[, race := as.factor(race)]
# str(dat[, cohort])
# str(dat[, sex])
# str(dat[, race])

#-------------------
# Parameters
#-------------------

# t_fits <- c(3, 6, 9, 12, 15, 18, 21, 24)
t_fits <- 1:24
kfolds <- 5
tau <- 24
cutoff <- 0
incr <- 0.01
pi.library <- c("SL.glm", "SL.lgb")
# event.library <- c("survSL.rfsrc", "survSL.coxph", "survSL.weibreg", "survSL.km") 
# cens.library <- c("survSL.rfsrc", "survSL.coxph", "survSL.weibreg", "survSL.km")
c.library <- c("SL.glm", "SL.lgb") 
d.library <- c("SL.glm", "SL.lgb")
covnames <- c("cohort", "sex", "race", "poverty", "popdensity", 
             "medianhousevalue", "medhouseholdincome", 
             "pct_owner_occ", "pct_blk", "hispanic", 
             "education", "smoke_rate", "mean_bmi")

#-------------------
# Analysis
#-------------------

# extract information from the data
start <- Sys.time()
nsamples <- nrow(dat)
acovnames <- c("a", covnames)

# estimation grid
eval.times <- seq(0, tau, by = incr)
fit.times <- quantile(eval.times, probs = seq(0, 0.95, by = 0.05))
nu <- length(eval.times)
nf <- length(fit.times)

#result holders
dat[, paste0("mu_1_", t_fits) := 0]
dat[, paste0("mu_0_", t_fits) := 0]
dat[, paste0("eta_1_", t_fits) := 0]
dat[, paste0("eta_0_", t_fits) := 0]

# creating folds so that delta 1 and 0 appear in all folds
ae.00 <- which((dat[, delta] == 0) & (dat[, a] == 0))
ae.01 <- which((dat[, delta] == 0) & (dat[, a] == 1))
ae.10 <- which((dat[, delta] == 1) & (dat[, a] == 0))
ae.11 <- which((dat[, delta] == 1) & (dat[, a] == 1))

folds.00 <- sample(rep(1:kfolds, length = length(ae.00)))
folds.01 <- sample(rep(1:kfolds, length = length(ae.01)))
folds.10 <- sample(rep(1:kfolds, length = length(ae.10)))
folds.11 <- sample(rep(1:kfolds, length = length(ae.11)))

idx <- rep(NA, nsamples)
idx[ae.00] <- folds.00
idx[ae.01] <- folds.01
idx[ae.10] <- folds.10
idx[ae.11] <- folds.11
rm(ae.00, ae.01, ae.10, ae.11, folds.00, folds.01, folds.10, folds.11)

cat("Start fitting ... \n \n")

# cross fitting
for (fold in 1:kfolds) {
  
  start_crsfit <- Sys.time() 
  
  # training and test indexes
  idx_tst <- idx == fold
  ntst <- sum(idx_tst)
  
  #---------- PROPENSITY SCORE
  
  pi.hat <- SuperLearner(
    Y = dat[!idx_tst, a],
    X = dat[!idx_tst, ..covnames],
    newX = dat[idx_tst, ..covnames],
    family='binomial',
    SL.library = pi.library
  )$SL.predict
  pi.hat <- pmax(pi.hat, cutoff)
  pi.hat <- pmin(pi.hat, 1-cutoff)
  pi.hat <- ifelse(dat[idx_tst, a] == 1, pi.hat, 1-pi.hat)
  
  # print progress
  end <- Sys.time()
  cat("Fitted propensity model in", hms_span(start_crsfit, end), "\n")
  
  #---------- SURVIVAL FUNCTION
  
  # Fit survival estimator for terminal event and censoring k(t; a, l), h(t; a, l)
  newX <- copy(dat[idx_tst, ..acovnames])
  setDT(newX)
  newX <- newX[, a := 1]
  newX.0 <- copy(newX)
  setDT(newX.0)
  newX.0 <- newX.0[, a:= 0]
  newX <- rbind(newX.0, newX)
  rm(newX.0)
  
  # Fitting H model using survival random forest
  fit.H <- survSL.rfsrc(time = dat[!idx_tst, x],
                        event = dat[!idx_tst, delta],
                        X = dat[!idx_tst, ..acovnames],
                        newX = newX,
                        new.times = eval.times, 
                        obsWeights = rep(1,sum(!idx_tst)),
                        id = NULL)
  
  # return predictions for A = 0
  h.hat.0 <- fit.H$pred[1:ntst,]
  h.hat.0 <- pmax(h.hat.0, cutoff) # truncation
  if(any(h.hat.0 == 0)) h.hat.0[h.hat.0 == 0] <- min(h.hat.0[h.hat.0 > 0]) # prevent 0
  
  # return predictions for A = 1
  h.hat.1 <- fit.H$pred[-(1:ntst),]
  h.hat.1 <- pmax(h.hat.1, cutoff) # truncation
  if(any(h.hat.1 == 0)) h.hat.1[h.hat.1 == 0] <- min(h.hat.1[h.hat.1 > 0]) # prevent 0
  
  # return predictions at X
  h.x <- sapply(1:ntst, function(i) {
    if (dat[idx_tst, a][i] == 1) {
      out <- stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(dat[idx_tst, x][i])
    } else {
      out <- stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(dat[idx_tst, x][i])
    }
    out
  })
  
  # print progress
  rm(fit.H)
  end2 <- Sys.time()
  cat("Fitted death model in", hms_span(end, end2), "\n")
  
  # Fitting K model using survival random forest
  fit.K <- survSL.rfsrc(time = dat[!idx_tst, x],
                        event = 1 - dat[!idx_tst, delta],
                        X = dat[!idx_tst, ..acovnames],
                        newX = rbind(newX, dat[!idx_tst, ..acovnames]),
                        new.times = eval.times, 
                        obsWeights = rep(1,sum(!idx_tst)),
                        id = NULL)
  
  # return predictions for A = 0
  k.hat.0 <- fit.K$pred[1:ntst,]
  k.hat.0 <- pmax(k.hat.0, cutoff) # truncation
  if(any(k.hat.0 == 0)) k.hat.0[k.hat.0 == 0] <- min(k.hat.0[k.hat.0 > 0]) # prevent 0
  
  # return predictions for A = 1
  k.hat.1 <- fit.K$pred[(ntst+1):(2*ntst),]
  k.hat.1 <- pmax(k.hat.1, cutoff) # truncation
  if(any(k.hat.1 == 0)) k.hat.1[k.hat.1 == 0] <- min(k.hat.1[k.hat.1 > 0]) # prevent 0
  
  # return predictions for training data
  k.fitted <- fit.K$pred[-(1:(2*ntst)),]
  k.fitted <- pmax(k.fitted, cutoff) # truncation
  if(any(k.fitted == 0)) k.fitted[k.fitted == 0] <- min(k.fitted[k.fitted > 0]) # prevent 0
  k.x.trn <- sapply(1:sum(!idx_tst), function(i) stepfun(eval.times, c(1,k.fitted[i,]), right = FALSE)(dat[!idx_tst, x][i]))
  
  # return predictions at X for test data
  k.x <- sapply(1:ntst, function(i) {
    if (dat[idx_tst, a][i] == 1) {
      out <- stepfun(eval.times, c(1,k.hat.1[i,]), right = FALSE)(dat[idx_tst, x][i])
    } else {
      out <- stepfun(eval.times, c(1,k.hat.0[i,]), right = FALSE)(dat[idx_tst, x][i])
    }
    out
  })
  
  # print progress
  rm(fit.K, k.fitted)
  end <- Sys.time()
  cat("Fitted censoring model in", hms_span(end2, end), "\n")
  
  #---------- FITTING b FUNCTION

  for (t in 1:length(t_fits)) {
    
    start_tfits <- Sys.time()
    
    b.tmp.1 <- matrix(0, nrow = ntst, ncol = nf)
    b.tmp.0 <- matrix(0, nrow = ntst, ncol = nf)
    c.pred.1 <- rep(0, ntst)
    c.pred.0 <- rep(0, ntst)
    d.pred.1 <- rep(0, ntst)
    d.pred.0 <- rep(0, ntst)
    
    idx.d <- dat[!idx_tst, delta] == 1
    y <- dat[!idx_tst, delta] * dat[!idx_tst, get(paste0("NE_", t_fits[t]))] / k.x.trn
    
    for (u in 1:nf) {
      
      idx.c <- dat[!idx_tst, x] > fit.times[u]
      
      crit.c <- sum(idx.c & idx.d) < 10
      crit.d <- (length(unique(dat[!idx_tst, delta][idx.c])) == 1) | (sum(idx.c) < 10)
      
      if (!crit.c) {
        
        # fit the c model
        fit.c <- SuperLearner(
          Y = y[idx.c & idx.d],          # I(X>u) = 1 and Delta = 1
          X = dat[!idx_tst, ..acovnames][idx.c & idx.d,],
          family = gaussian(),
          newX = newX,
          method = "method.NNLS",
          SL.library = c.library
        )$SL.predict
        c.pred.0 <- fit.c[1:ntst] 
        c.pred.1 <- fit.c[-c(1:ntst)]
        rm(fit.c)
      }
      
      if (!crit.d) {
        # fit the d model
        fit.d <- SuperLearner(
          Y = dat[!idx_tst, delta][idx.c],              # I(X>u) = 1
          X = dat[!idx_tst, ..acovnames][idx.c,],
          family = binomial(),
          newX = newX,
          method = "method.CC_LS",
          SL.library = d.library
        )$SL.predict
        d.pred.0 <- fit.d[1:ntst]
        d.pred.1 <- fit.d[-c(1:ntst)]
        rm(fit.d)
      }
      
      b.tmp.1[,u] <- c.pred.1 * d.pred.1
      b.tmp.0[,u] <- c.pred.0 * d.pred.0
    }
    
    b.hat.1 <- sapply(1:ntst, function(i) {
      approx(x = fit.times, y = b.tmp.1[i,], xout = eval.times, rule = 2)$y
    })
    b.hat.0 <- sapply(1:ntst, function(i) {
      approx(x = fit.times, y = b.tmp.0[i,], xout = eval.times, rule = 2)$y
    })
    
    rm(b.tmp.0, b.tmp.1, c.pred.0, c.pred.1, d.pred.0, d.pred.1)
    
    #---------- CALCULATE INPUTVECS
    
    h.t.1 <- sapply(1:ntst, function(i) stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(t_fits[t]))
    h.t.0 <- sapply(1:ntst, function(i) stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(t_fits[t]))
    h.x.v.t <- sapply(1:ntst, function(i) {
      if (dat[idx_tst, a][i] == 1) {
        out <- stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(max(t_fits[t], dat[idx_tst, x][i]))
      } else {
        out <- stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(max(t_fits[t], dat[idx_tst, x][i]))
      }
      out
    })
    
    f.0.1 <- b.hat.1[1,] #because HK = 1 when u = 0
    f.0.0 <- b.hat.0[1,]
    b.x <- sapply(1:ntst, function(i) {
      if (dat[idx_tst, a][i] == 1) {
        out <- stepfun(eval.times, c(b.hat.1[1,i],b.hat.1[,i]), right = FALSE)(dat[idx_tst, x][i])
      } else {
        out <- stepfun(eval.times, c(b.hat.0[1,i],b.hat.0[,i]), right = FALSE)(dat[idx_tst, x][i])
      }
      out
    })
    
    int.f <- sapply(1:ntst, function(i) {
      if (dat[idx_tst, a][i] == 1) {
        vals <- diff(-log(k.hat.1[i,])) * ((b.hat.1[-nu,i] + b.hat.1[-1,i])/2) #trapezoidal rule
        if(any(eval.times[-1] > dat[idx_tst, x][i])) vals[eval.times[-1] > dat[idx_tst, x][i]] <- 0 #only integrate until x
        out <- sum(vals)
      } else {
        vals <- diff(-log(k.hat.0[i,])) * ((b.hat.0[-nu,i] + b.hat.0[-1,i])/2) #trapezoidal rule
        if(any(eval.times[-1] > dat[idx_tst, x][i])) vals[eval.times[-1] > dat[idx_tst, x][i]] <- 0 #only integrate until x
        out <- sum(vals)
      }
      out
    })
    
    uvt <- sapply(1:nu, function(k) {max(eval.times[k], t_fits[t])})
    int.h <- sapply(1:ntst, function(i) {
      if (dat[idx_tst, a][i] == 1) {
        huvt <- stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(uvt)
        vals <- diff(1/k.hat.1[i,]) * ((huvt[-nu] / h.hat.1[i,-nu] + huvt[-1] / h.hat.1[i,-1])/2)
        if(any(eval.times[-1] > dat[idx_tst, x][i])) vals[eval.times[-1] > dat[idx_tst, x][i]] <- 0
        out <- sum(vals)
      } else {
        huvt <- stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(uvt)
        vals <- diff(1/k.hat.0[i,]) * ((huvt[-nu] / h.hat.0[i,-nu] + huvt[-1] / h.hat.0[i,-1])/2)
        if(any(eval.times[-1] > dat[idx_tst, x][i])) vals[eval.times[-1] > dat[idx_tst, x][i]] <- 0
        out <- sum(vals)
      }
      out
    })
    
    #---------- CALCULATE ESTIMATORS
    
    rm(b.hat.0, b.hat.1)
    ipwts_1 <- (as.numeric(dat[idx_tst, a] == 1) * dat[idx_tst, delta]) / (pi.hat * k.x) 
    ipwts_0 <- (as.numeric(dat[idx_tst, a] == 0) * dat[idx_tst, delta]) / (pi.hat * k.x)
    
    # IPW part
    ipwest <- cbind(ipwts_1 * dat[idx_tst, get(paste0("NE_", t_fits[t]))], 
                    ipwts_0 * dat[idx_tst, get(paste0("NE_", t_fits[t]))], 
                    ipwts_1 * as.numeric(dat[idx_tst, x] > t_fits[t]), 
                    ipwts_0 * as.numeric(dat[idx_tst, x] > t_fits[t]))
    
    aug_1 <- cbind(((as.numeric(dat[idx_tst, a] == 1) - pi.hat) * f.0.1) / pi.hat, 
                   ((as.numeric(dat[idx_tst, a] == 0) - pi.hat) * f.0.0) / pi.hat,
                   ((as.numeric(dat[idx_tst, a] == 1) - pi.hat) * h.t.1) / pi.hat,
                   ((as.numeric(dat[idx_tst, a] == 0) - pi.hat) * h.t.0) / pi.hat)
    
    aug_2 <- cbind((as.numeric(dat[idx_tst, a] == 1) * (1-dat[idx_tst, delta]) * b.x) / pi.hat,
                   (as.numeric(dat[idx_tst, a] == 0) * (1-dat[idx_tst, delta]) * b.x) / pi.hat,
                   (as.numeric(dat[idx_tst, a] == 1) * (1-dat[idx_tst, delta]) * h.x.v.t) / (pi.hat * h.x * k.x),
                   (as.numeric(dat[idx_tst, a] == 0) * (1-dat[idx_tst, delta]) * h.x.v.t) / (pi.hat * h.x * k.x))
    
    aug_3 <- cbind((as.numeric(dat[idx_tst, a] == 1) * int.f) / pi.hat, 
                   (as.numeric(dat[idx_tst, a] == 0) * int.f) / pi.hat, 
                   (as.numeric(dat[idx_tst, a] == 1) * int.h) / pi.hat, 
                   (as.numeric(dat[idx_tst, a] == 0) * int.h) / pi.hat)
    
    infl <- ipwest - aug_1 + aug_2 - aug_3
    dat[idx_tst, paste0("mu_1_", t_fits[t]) := infl[,1]]
    dat[idx_tst, paste0("mu_0_", t_fits[t]) := infl[,2]]
    dat[idx_tst, paste0("eta_1_", t_fits[t]) := infl[,3]]
    dat[idx_tst, paste0("eta_0_", t_fits[t]) := infl[,4]]
    
    end <- Sys.time()
    cat("Time t =", t_fits[t], "finished in", hms_span(start_tfits, end), "\n")
  }
  
  rm(k.hat.0, k.hat.1, h.hat.0, h.hat.1, newX)
  
  # print progress
  end <- Sys.time()
  cat("Fold ", fold, " finished in", hms_span(start_crsfit, end), "\n\n")
  
} # end of cross fitting

#---------- SUMMARIZE RESULTS

# # output object
# out <- data.frame(size = double(),
#                   method = character(), 
#                   cutoff = character(), 
#                   time = double(),
#                   estimand = character(), 
#                   est = double(), 
#                   sd = double())
# 
# cnt <- 0
# for (t in 1:length(t_fits)) {
#   clnm <- paste0(c("mu_1_", "mu_0_", "eta_1_", "eta_0_"), t_fits[t])
#   est <- colSums(dat[, ..clnm]) / nsamples
#   vars <- cbind(dat[, get(paste0("mu_1_", t_fits[t]))] - est[1],
#                 dat[, get(paste0("mu_0_", t_fits[t]))] - est[2],
#                 dat[, get(paste0("eta_1_", t_fits[t]))] - est[3],
#                 dat[, get(paste0("eta_0_", t_fits[t]))] - est[4])
#   EU <- diag(c(-nsamples, -nsamples,
#                -nsamples, -nsamples))
#   vars <- solve(EU) %*% t(vars) %*% vars %*% t(solve(EU))
#   assign(paste0("vars.", t_fits[t]), vars)
#   
#   out[cnt+1:4,] <- list(size = nsamples,
#                         method = "One-step AIPW",
#                         cutoff = cutoff,
#                         time = t_fits[t],
#                         estimand = c("mu_1", "mu_0", "eta_1", "eta_0"),
#                         est = est,
#                         sd = sqrt(diag(vars)))
#   cnt <- cnt+4
# }

# print progress
end <- Sys.time()
cat("Function fitting finished in", hms_span(start, end), "\n\n")

#-------------------
# Save results
#-------------------

# save(list = c("out", paste0("vars.", t_fits)), 
#      file = paste0(dir_results, "results_", state_code, "_buffer30.RData"))

write_feather(dat, paste0(scratch, "results/", state_code, "_buffer30.feather"))
