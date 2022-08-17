library(rstan)
library(loo)
library(purrr)
library(future)
library(furrr)
library(dplyr)

rm(list = ls())
source("sim_functions.R")

# Experimental factors
param_list <- list(
  mobility = c(.2,.1),
  J = c(100, 50),
  rho = c(FALSE, TRUE),
  tau = c("same", "different")
)

params <- param_list %>% cross_df()
params$reps <- 2 # Number of replications?
params$seed <- 1:nrow(params) + sample(2^30, size = 1) # set seed
print(paste0(nrow(params), " fully crossed experimental conditions have been generated."))
params

# run simulations ---------------------------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 2)
options(error=recover)
plan("multisession") 

system.time(
  results <-
    params %>%
    mutate(res = future_pmap(., .f = runSim, .options = furrr_options(seed=NULL))) %>% 
    dplyr::select(-c(mobility, J, rho, tau)) %>% 
    tidyr::unnest(cols = res)
) 


# save results ------------------------------------------------------------
session_info <- sessionInfo()
run_date <- date()
save(params, results, session_info, run_date, file = "sim_2rep.Rdata")
