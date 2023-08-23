rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

library(dplyr, quietly = TRUE, warn.conflicts = FALSE) # 1.0.10
library(tidyr)       # 0.8.1
require(purrr)       # 0.3.4
require(rstan)       # 2.17.3
require(StanHeaders) # 2.17.2
require(RcppEigen)   # 0.3.3.3.0
require(BH)          # 1.65.0-1
require(loo)         # 2.0.0
require(doParallel)
require(Rmpi)

# Parse command line arguments --------------------------------------------
extract_int_arg <- function(arg_string, arg_name, default) {
  res <- 
    arg_string %>%
    stringr::str_extract(paste(arg_name, "[0-9]+")) %>%
    stringr::str_sub(stringr::str_length(arg_name) + 2, -1) %>%
    as.integer()
  
  if (is.na(res)) res <- default
  res
}

arg_string <- paste(args, collapse = " ")

cores_n <- extract_int_arg(arg_string, "-cores", 48L)
# batch <- extract_int_arg(arg_string, "-batch", 1L)
reps <- extract_int_arg(arg_string, "-reps", 5L)


# Source the function -----------------------------------------------------
source("sim_functions.R")
source_obj <- ls()


# Experimental Design -----------------------------------------------------
set.seed(20220926)

# Experimental factors
design_factors <- list(
  mobility = c(.2,.1),
  J = c(100, 50),
  rho = c(FALSE, TRUE),
  tau = c("same", "different")
)

# combine into a design set
params <- 
  purrr::cross_df(design_factors) %>%
  dplyr::mutate(
    iterations = reps,
    cores_val = cores_n,
    seed = 710724962 + 1:n()#round(runif(1) * 2^30) + 1:n() 
  )

range(params$seed)
params


# Run simulations in parallel ---------------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(error=recover)

tm <- system.time(
  results <- plyr::mdply(params, .fun = runSim, .parallel = TRUE)
)

tm

write.csv(results, "results.csv", row.names = F)
# Save results and details ------------------------------------------------
session_info <- sessionInfo()
run_date <- date()

results_name <- paste0("Simulation-results.Rdata")
save(results, tm, session_info, run_date, file = results_name)
