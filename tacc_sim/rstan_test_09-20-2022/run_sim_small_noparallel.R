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
# plan("multisession") 

system.time(
  results <- runSim(reps = 1, mobility = .1, J = 100, rho = FALSE, tau = "same")
    # params %>%
    # mutate(res = future_pmap(., .f = runSim, .options = furrr_options(seed=NULL))) %>% 
    # dplyr::select(-c(mobility, J, rho, tau)) %>% 
    # tidyr::unnest(cols = res)
) 


# save results ------------------------------------------------------------
session_info <- sessionInfo()
run_date <- date()
save(results, session_info, run_date, file = "sim_1rep_small.Rdata")
