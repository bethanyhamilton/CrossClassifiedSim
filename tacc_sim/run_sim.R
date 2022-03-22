source("sim_functions.R")


library(dplyr)
library(purrr)
library(Pusto)


set.seed(20222203)
# Experimental factors
mobility <- c(.2,.1)
J <- c(100, 50)
rho <- c(FALSE, TRUE)
tau <- c("same", "different" )


# Number of replications How many??
reps <- 2


param_list <- list(
  mobility=mobility, 
  J=J, 
  rho = rho, 
  tau = tau,
  reps = reps
)


#not_right..
params <- 
  param_list %>%
  cross_df() %>%
  mutate(
    seed = 1:n() + sample(2^30, size = 1),
  )

rownames(params) <- NULL

print(paste0(nrow(params), " fully crossed experimental conditions have been generated."))



source_obj <- ls()




cluster <- Pusto::start_parallel(source_obj = source_obj, setup = "register")

tm <- system.time(
  results <- plyr::mdply(params, .f = runSim, 
                         .parallel = TRUE)
)

tm 
parallel::stopCluster(cluster)

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "simulation_results.Rdata")

