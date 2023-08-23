require(rstan)       # 2.17.3
require(StanHeaders) # 2.17.2
require(RcppEigen)   # 0.3.3.3.0
require(BH)          # 1.65.0-1
require(purrr)       # 0.3.4
# require(future)    # ?
# require(furrr)     # ?
# require(loo)       # ?
# require(dplyr)     # ?
# require(tidyr)     # ?
# The parallel version 
library(foreach)     # 1.4.4
library(doParallel)  # 1.0.11
library(parallel)    # 3.5.1
library(Rmpi)        # for TACC; 0.6.7


rm(list = ls())
source("sim_functions_small.R")

# Experimental factors
# param_list <- list(
#   mobility = c(.2,.1),
#   J = 100,#c(100, 50),
#   rho = FALSE,#c(FALSE, TRUE),
#   tau = "same"# c("same", "different")
# )
# 
# params <- param_list %>% cross_df()
# params$reps <- 1 # Number of replications?
# params$seed <- 1:nrow(params) + sample(2^30, size = 1) # set seed
# print(paste0(nrow(params), " fully crossed experimental conditions have been generated."))
# params

# run simulations ---------------------------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(error=recover)

# for desktop
# n_cores <- detectCores() - 2
# cl <- makeCluster(n_cores)
# for TACC
cl <- getMPIcluster()
registerDoParallel(cl)

n_tasks <- 3
data <- UDF_sample(mobility = .1, J = 100, rho= FALSE, tau="same")

system.time(
results <-
  foreach (task =  1:n_tasks, .combine = rbind, 
           .packages=c("rstan"), .export=c("optimizing")) %dopar% {
    set.seed(task)
    tmp <- run_mmrem(data, mobility = .1, J = 100, rho= FALSE, tau="same")
    
    c(task, tmp)
  } 
)

# save results ------------------------------------------------------------
session_info <- sessionInfo()
run_date <- date()
save(results, session_info, run_date, file = "sim_1rep_small.Rdata")
