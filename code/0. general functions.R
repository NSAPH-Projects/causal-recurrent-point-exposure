#################################
#Loading libraries
#################################

setwd("~/nsaph_projects/causal-recurrent-point-exposure/code")
.libPaths("/n/User/R/ifxrstudio/RELEASE_3_18")

library(reda) # generate recurrent event data
library(survSuperLearner) # Westling - ensemble learning to estimate conditional survival functions from right-censored time-to-event data
library(SuperLearner) # super learner
library(randomForestSRC) # fast unified random forests for survival, regression, and classification
library(survival) # survival analysis

library(CFsurvival) # Westling2023 paper
library(lightgbm)
library(quadprog)

library(data.table)
library(arrow) # C++ dplyr
library(fst) # fast serialization of dataframe
library(lubridate) # process time

# library(devtools)
# install.packages("reda")
# install.packages("SuperLearner")
# install.packages("randomForestSRC")
# install.packages("survival")
# install.packages("xgboost")
# install.packages("randomForest")
# install.packages("lightgbm")
# install.packages("quadprog")
# devtools::install_github("tedwestling/CFsurvival")
# devtools::install_github("tedwestling/survSuperLearner")

#################################
#Helper functions
#################################

#' Time printing
hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

#################################
#Functions to generate data
#################################

#' Weibull rate function: rho(t) = f(t)/S(t-) for baseline hazard
rho_weibull <- function(t, lambda, k) {
  # lambda: scale
  # k: shape
  out <- k * t ^ (k - 1) * lambda ^ (-k)
  return(out)
}

#-------------------
#SuperLearner wrappers
#-------------------

SL.lgb <- function(Y, X, newX, family, obsWeights, id,
                   depth = 3, rate = 0.01, rounds = 200) {
  stopifnot(require("lightgbm"))
  
  params <- list(task = "train",
                 max_depth = depth,
                 learning_rate = rate)
  
  if (family$family == 'gaussian') {
    params$objective <- "regression"
  } else if (family$family == 'binomial') {
    params$objective <- "binary"
  }
  
  m <- lgb.train(params,
                 lgb.Dataset(data = as.matrix(X), label = Y, weight = obsWeights),
                 nrounds = rounds, verbose = -1)
  
  pred <- predict(m, as.matrix(newX))
  
  fit = list(object = m)
  class(fit) <- 'SL.lgb'
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.lgb <- function(object, newdata, ...) {
  stopifnot(require("lightgbm"))
  pred <- predict(object$object, as.matrix(newdata))
  return(pred)
}

SL.lgb.1 <- function(..., depth=3, rate=0.01, rounds=50) {
  SL.lgb(..., depth=depth, rate=rate, rounds=rounds) }
SL.lgb.2 <- function(..., depth=3, rate=0.1, rounds=50) {
  SL.lgb(..., depth=depth, rate=rate, rounds=rounds) }
SL.lgb.3 <- function(..., depth=3, rate=0.01, rounds=200) {
  SL.lgb(..., depth=depth, rate=rate, rounds=rounds) }
SL.lgb.4 <- function(..., depth=3, rate=0.1, rounds=200) {
  SL.lgb(..., depth=depth, rate=rate, rounds=rounds) }


# Define HAL function for SuperLearner
SL.hal <- function(Y, X, newX, family, obsWeights, id,
                   num_knots = 10, max_degree = 1) {
  stopifnot(require("hal9001"))
  
  # if (family$family == 'gaussian') {
  #   params$objective <- "regression"
  # } else if (family$family == 'binomial') {
  #   params$objective <- "binary"
  # }
  
  m <- fit_hal(as.matrix(X), Y,
               num_knots = num_knots,
               max_degree = max_degree)
  
  pred <- predict(m, as.matrix(newX))
  
  fit = list(object = m)
  class(fit) <- 'SL.hal'
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.hal <- function(object, newdata, ...) {
  stopifnot(require("hal9001"))
  pred <- predict(object$object, as.matrix(newdata))
  return(pred)
}



