library(dplyr)    # 3.5.2
library(tidyr)    # 3.5.3
library(lme4)     # 3.5.3
library(sandwich)
library(lmtest)   # 3.5.2
library(lmerTest) # 3.5.2
library(lfe)      # 3.5.2
library(truncnorm)# 3.5.3
library(purrr)    # 3.5.3
library(psych)    # 3.5.3
library(future)   # 3.5.3
library(furrr)    # 3.5.2

rm(list = ls())
source("QP_simfunctions.R")
source("Pusto_functions.R")

#--------------------------------------------------------
# simulation parameters 
#--------------------------------------------------------
design_factors <- list(
  gamma000 = 0, gamma100 = c(0.5), 
  gamma010 = c(0.5), gamma002 = c(0.5), 
  G = c(70, 245), H = c(20, 70), # H = school; G = neighborhood
  ICC_g = c(.05), ICC_h = c(0.05), 
  sparse = .1, J = c(100),
  L1cov_m = 0, L1cov_sd = 10,  L2cov_m = 0, L2cov_sd = 1,
  assumption = c("met")
)

params <- 
  cross_df(design_factors) %>% 
  mutate(iterations = 50, seed = 20220525 + 1:n())

#--------------------------------------------------------
# run simulations in parallel - future + furrr workflow
#--------------------------------------------------------
source_obj <- ls()
cluster <- start_parallel(source_obj = source_obj, 
                          setup = "register")
system.time(
  results <-
    plyr::mdply(params, .f = run_sim, .parallel = TRUE)
)

stop_parallel(cluster)

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------
session_info <- sessionInfo()
run_date <- date()
save(results, params, session_info, run_date, file = "sim_results.Rdata")
