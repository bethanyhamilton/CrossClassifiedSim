# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#installation-of-rstan
remove.packages(c("StanHeaders", "rstan"))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
if (file.exists(".RData")) file.remove(".RData")

# To verify your installation, you can run the RStan example/test model:
example(stan_model, package = "rstan", run.dontrun = TRUE)

# Loading the package
library(rstan)
library(loo)
library(dplyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# my model --------------------------------------------------------
rm(list = ls())

# design
mobility <- c(.2)
J <- c(100)
rho <- c(FALSE)
tau <- c("same")

# Data Generation
data <- UDF_sample(mobility = mobility, J = J, rho = rho, tau = tau) 

# Make data list
ccrem_mmrem_data <- getdata_ccrem_mmrem(data)
rccrem_data <- getdata_rccrem(data)

# model fit ------------------------------------------------------
## MM ------------------------------------------------------------
mm_stan_fit <-
  stan(
    model_code = mmrem_model,
    data = ccrem_mmrem_data,
    iter = 2000,
    warmup = 1000,
    thin = 10,
    chains = 4,
    cores = 6,
    control = list(adapt_delta = 0.999, max_treedepth = 25)
  )

library(bayesplot)
# Extract posterior draws for later use
posterior_cp <- as.array(mm_stan_fit)
posterior_ncp <- as.array(mm_stan_fit)
lp_cp <- log_posterior(mm_stan_fit)
head(lp_cp)
np_cp <- nuts_params(mm_stan_fit)
head(np_cp)

# trace plot
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, 
           pars = c("sigma_u2", "sigma_y", "popint", "popslope"), 
           np = np_cp) 

# mcmc_nuts_divergence
color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp)

# R-hat (should be < 1.1)
print(mm_stan_fit)
rhats <- rhat(mm_stan_fit)[1:4]
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

# n_eff (a useful heuristic is to worry about any neff/N less than 0.1)
ratios_cp <- neff_ratio(mm_stan_fit)[1:4]
print(ratios_cp)
mcmc_neff(ratios_cp, size = 2)

## CC ------------------------------------------------------------
cc_stan_fit <-
  stan(
    model_code = ccrem_model,
    data = ccrem_mmrem_data,
    iter = 2000,
    warmup = 1000,
    thin = 10,
    chains = 4,
    cores = 6,
    control = list(adapt_delta = 0.999, max_treedepth = 25)
  )

# Extract posterior draws for later use
fit <- cc_stan_fit
posterior_cp <- as.array(fit)
posterior_ncp <- as.array(fit)
lp_cp <- log_posterior(fit)
head(lp_cp)
np_cp <- nuts_params(fit)
head(np_cp)

# trace plot
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, 
           pars = c("sigma_u2[1]", "sigma_u2[2]", "sigma_u2[3]", "sigma_u2[4]",
                    "sigma_y", "popint", "popslope"), 
           np = np_cp) 

# mcmc_nuts_divergence
color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp)

# R-hat (should be < 1.1)
print(fit)
rhats <- rhat(fit)[1:7]
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

# n_eff (a useful heuristic is to worry about any neff/N less than 0.1)
ratios_cp <- neff_ratio(fit)[1:7]
print(ratios_cp)
mcmc_neff(ratios_cp, size = 2)

## rCC -----------------------------------------------------------
rcc_stan_fit <-
  stan(
    model_code = rccrem_model_new,
    data = rccrem_data,
    iter = 3000,
    warmup = 1500,
    thin = 10,
    chains = 4,
    cores = 6,
    control = list(adapt_delta = 0.999, max_treedepth = 25)
  )

# Extract posterior draws for later use
fit <- rcc_stan_fit
posterior_cp <- as.array(fit)
posterior_ncp <- as.array(fit)
lp_cp <- log_posterior(fit)
head(lp_cp)
np_cp <- nuts_params(fit)
head(np_cp)

# trace plot
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, 
           pars = c("sigma_u2[1]", "sigma_u2[2]", "sigma_u2[3]", "sigma_u2[4]",
                    "rho_12", "rho_13", "rho_14", "rho_23", "rho_24", "rho_34",
                    "sigma_y", "popint", "popslope"), 
           np = np_cp) 

# mcmc_nuts_divergence
color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp)

# R-hat (should be < 1.1)
print(fit)
fit_df <- summary(fit)$summary
fit_df <- fit_df[row.names(fit_df) %in% c("sigma_u2[1]", "sigma_u2[2]", "sigma_u2[3]", "sigma_u2[4]",
                                          "rho_12", "rho_13", "rho_14", "rho_23", "rho_24", "rho_34",
                                          "sigma_y", "popint", "popslope"),]
rhats <- fit_df[,10]
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

# n_eff (a useful heuristic is to worry about any neff/N less than 0.1)
ratios_cp <- fit_df[,9]
print(ratios_cp)
mcmc_neff(ratios_cp, size = 2)

# Run Models Funcs -----------------------------------------------
res_mmrem <- run_mmrem(data,
                       mobility = mobility, J = J, rho = rho, tau = tau)


res_ccrem <- run_ccrem(data,
                       mobility = mobility, J = J, rho = rho, tau = tau)

res_rccrem <- run_rccrem(data,
                         mobility = mobility, J = J, rho = rho, tau = tau) 

# Run Sims ----------------------------------------------------------------

results <- runmodels(mobility = mobility, J = J, rho=rho, tau=tau,
                     mmrem=TRUE, ccrem=TRUE, rccrem=TRUE)


runSim <- function(reps, mobility, J, rho, tau, mmrem=TRUE, ccrem=TRUE, rccrem=TRUE, seed = NULL, ...){
  
  suppressPackageStartupMessages(require(purrr, quietly = TRUE, warn.conflicts = FALSE))
  suppressPackageStartupMessages(require(dplyr, quietly = TRUE, warn.conflicts = FALSE))
  
  if (!is.null(seed)) set.seed(seed)
  
  replicates <- purrr::rerun(reps, {
    runmodels(mobility = mobility, J = J, rho=rho, tau=tau,
              mmrem = mmrem, ccrem = ccrem, rccrem= rccrem)
  })
  
  result <- do.call(rbind, replicates)
  result
  
}


