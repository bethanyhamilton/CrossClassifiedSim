# Data Generating Model----------------------------------------------------
generate_dat <- function(gamma000, gamma100, gamma010, gamma002,
                         G, H, ICC_g, ICC_h, sparse, J,
                         L1cov_m, L1cov_sd, L2cov_m, L2cov_sd,
                         assumption) {
  
  # set sigma, tau_G00 and tau_H00 based on ICC
  tau_G00 = ICC_g # Neighborhood
  tau_H00 = ICC_h # School
  sigma = sqrt(1-tau_G00-tau_H00)
  
  # data assignment
  dat <- 
    rerun(.n = H, sample(1:(sparse * G), J, replace = TRUE)) %>% 
    map_df(~ tibble(neighid = .x), .id = "schid") %>% 
    as.data.frame() %>% 
    mutate(
      schid = as.numeric(schid),
      neighid = round(schid * (G / H)) + neighid,
      neighid = ifelse(neighid > G, neighid - G, neighid)
    ) 
  
  # neighborhood data
  # create between-neighborhood variance of X 
  X_bw_neigh <- 
    dat %>% 
    group_by(neighid) %>% 
    summarise() %>% 
    mutate(X_bw_neigh = rnorm(nrow(.), mean = L1cov_m, sd = sqrt(.2 * L1cov_sd^2)))
  
  tau_G10 <- if (assumption == "random slopes") 0.05 else 0.00
  
  if (assumption == "exogeneity") {
    # when exogeneity assumption is violated
    r_g <- 0.4 # correlation between X_bw_neigh and b_0g0
    
    neighbordata <-
      X_bw_neigh %>% 
      mutate(
        W = rnorm(nrow(.), mean = L2cov_m, sd = L2cov_sd),
        v_0g0 = rnorm(nrow(.), mean = 0, sd = sqrt((1 - r_g ^ 2) * tau_G00)),
        b_0g0 = r_g * sqrt(tau_G00 / (.2 * L1cov_sd ^ 2)) * X_bw_neigh + v_0g0, # neighborhood random effect
        
        v_1g0 = rnorm(nrow(.), mean = 0, sd = sqrt((1 - r_g ^ 2) * tau_G10)),
        b_1g0 = r_g * sqrt(tau_G10 / (.2 * L1cov_sd ^ 2)) * X_bw_neigh + v_1g0 # random slope for X
      ) %>% 
      dplyr::select(-c(v_0g0, v_1g0))
    
  } else {
    
    # when all assumptions are met:
    neighbordata <-
      X_bw_neigh %>% 
      mutate(
        W = rnorm(nrow(.), mean = L2cov_m, sd = L2cov_sd),
        b_0g0 = rnorm(nrow(.), mean = 0, sd = sqrt(tau_G00)), # neighborhood random effect
        b_1g0 = rnorm(nrow(.), mean = 0, sd = sqrt(tau_G10)), # random slope for X
      ) 
  }
  
  dat <- dat %>% left_join(neighbordata, by = "neighid") 
  
  # school data
  # create between-school variance of X
  schooldata <- 
    dat %>% 
    group_by(schid) %>%
    summarise() %>%
    mutate(
      Z = rnorm(H, mean = L2cov_m, sd = L2cov_sd), # neighborhood-level Z
      c_00h = rnorm(H, mean = 0, sd = sqrt(tau_H00)),
      X_bw_school = rnorm(nrow(.), mean = L1cov_m, sd = sqrt(.2*L1cov_sd^2))
    )
  
  dat <- dat %>% left_join(schooldata, by = "schid")
  
  # student data
  dat <-
    dat %>%
    mutate(
      # student ID
      stuid = 1:nrow(.),
      
      # student-level X
      X_within = rnorm(nrow(.), mean = L1cov_m, sd = sqrt(.6*L1cov_sd^2)),
      X = X_bw_neigh + X_bw_school + X_within,
      
      # student-level residuals u
      u_sd = if (assumption == "heterosced") sigma * sqrt(exp((2 * 15 * (X - L1cov_m) - L1cov_sd^2) / (2 * 15^2))) else sigma,
      u = rnorm(nrow(.), mean = 0, sd = u_sd),
      
      # outcome variable
      y = gamma000 + gamma100 * X + gamma010 * W + gamma002 * Z + b_1g0 * X_within + b_0g0 + c_00h + u
    )

  return(dat)
}

# Model-fitting/Estimation-------------------------------------------------

## CCREM
estimate_ccrem <- function(dat) {
  
  # estimation
  model <- lmer(y ~ 1 + X + W + Z + (1 | schid) + (1 | neighid), 
                data = dat)
  summary <- summary(model)
  
  fixed_est <- 
    summary$coefficients %>% 
    as_tibble(rownames = "cov") %>% 
    dplyr::select(cov, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
    mutate(
      method = "CCREM",
      # convergence
      converged = (is.na(is.na(model@optinfo$conv$lme4)[1]))
    ) %>% 
    dplyr::select(cov, method, everything()) %>% 
    filter(cov %in% c("X", "W", "Z"))
  
  return(fixed_est)
}

## OLS-CRVE
estimate_ols <- function(dat) {
  
  # estimation
  model_ols <- felm(y ~ X + W + Z | 0 | 0 | schid + neighid, 
                    data = dat, psdef = FALSE) # see the vignette
  summary <- coeftest(model_ols)
  
  # fixed effects
  fixed_est <- 
    summary[2:4, c(1, 2, 4)] %>% 
    as_tibble(rownames = "cov") %>% 
    mutate(method = "OLS") %>% 
    dplyr::select(cov, method, est = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) 
  
  return(fixed_est)
}

## FE-CRVE
estimate_fe <- function(dat) {
  
  # estimation
  ## felm(equation | fixed effect | 0 | clustering)
  model_fem <- felm(y ~ X | schid + neighid | 0 | schid + neighid, 
                    data = dat, psdef = FALSE)
  summary <- summary(model_fem)
  
  # fixed effects
  fixed_est <- 
    summary$coefficients[, c(1, 2, 4),drop = FALSE] %>% 
    as_tibble(rownames = "cov") %>% 
    mutate(method = "FE") %>% 
    dplyr::select(cov, method, est = Estimate, se = `Cluster s.e.`, pval = `Pr(>|t|)`)
  
  return(fixed_est)
}


# bind_results

estimate <- function(dat, CI_level = .95) {
  
  crit <- qnorm((1 + CI_level) / 2)
  
  results <- bind_rows(
    estimate_ccrem(dat),
    estimate_ols(dat),
    estimate_fe(dat)
  ) %>% 
    mutate(
      var = se^2,
      lower_bound = est - crit * se,
      upper_bound = est + crit * se,
    ) %>%
    as_tibble()

  return(results)
}


# Performance calculations ------------------------------------------------
calc_performance <- function(results, CI_level = .95) {
  
  abs_crit <- results %>%
    group_by(method, cov) %>%
    group_modify(~ calc_absolute(.x, estimates = est, true_param = param))
  
  rel_crit <- results %>%
    group_by(method, cov) %>%
    group_modify(~ calc_relative(.x, estimates = est, true_param = param)) 
  
  # Relative Criteria for Variance Estimators
  rel_crit_val <- results %>% 
    group_by(method, cov) %>%
    group_modify(~ calc_relative_var(.x, estimates = est, var_estimates = var))
  
  # Hypothesis Testing
  crit <- qnorm((1 + CI_level) / 2)
  rejection_rate <- results %>% 
    dplyr::select(-var) %>% 
    mutate(rej_rate = ifelse(abs(est - param)/se >= crit, 1, 0)) %>% 
    group_by(method, cov) %>% 
    summarise(rej_rate = mean(rej_rate), .groups = "drop")
  
  power <- results %>%
    group_by(method, cov) %>%
    group_modify(~ calc_rejection(.x, p_values = pval)) %>% 
    rename(power = rej_rate) 
  
  #  Confidence Intervals
  conf_int <- results %>%
    group_by(method, cov) %>%
    group_modify(~ calc_coverage(.x, lower_bound = lower_bound, 
                                 upper_bound = upper_bound, 
                                 true_param = param))
  
  # Convergence Rate
  convergence <- results %>% 
    group_by(method, cov) %>% 
    summarise(convergence_rate = sum(converged)/n(), .groups = "drop")
  
  performance_measures <- rejection_rate %>% 
    left_join(abs_crit, by = c("method", "cov")) %>% 
    left_join(rel_crit, by = c("method", "cov", "K")) %>% 
    left_join(rel_crit_val, by = c("method", "cov", "K")) %>% 
    left_join(power, by = c("method", "cov", "K")) %>% 
    left_join(conf_int, by = c("method", "cov", "K")) %>% 
    left_join(convergence, by = c("method", "cov"))
  
  return(performance_measures)
}

# Simulation driver -------------------------------------------------------
run_sim <- function(iterations, gamma000, gamma100, gamma010, gamma002, 
                    G, H, ICC_g, ICC_h, sparse, J, 
                    L1cov_m, L1cov_sd, L2cov_m, L2cov_sd, assumption,
                    seed = NULL) {
  require(dplyr)    # 3.5.2
  require(tidyr)    # 3.5.3
  require(lme4)     # 3.5.3
  require(sandwich)
  require(lmtest)   # 3.5.2
  require(lmerTest) # 3.5.2
  require(lfe)      # 3.5.2
  require(truncnorm)# 3.5.3
  require(purrr)    # 3.5.3
  # require(simhelpers) # citation("simhelpers")
  require(psych)    # 3.5.3
  require(future)   # 3.5.3
  require(furrr)    # 3.5.2
  
  if (!is.null(seed)) set.seed(seed)
  
  results <-
    rerun(iterations, {
      generate_dat(
        gamma000 = gamma000, gamma100 = gamma100,
        gamma010 = gamma010, gamma002 = gamma002,
        G = G, H = H, ICC_g = ICC_g, ICC_h = ICC_h, 
        sparse = sparse, J = J,
        L1cov_m = L1cov_m, L1cov_sd = L1cov_sd,
        L2cov_m = L2cov_m, L2cov_sd = L2cov_sd, assumption = assumption
      ) %>%
      estimate()
    }) %>%
    bind_rows() %>%
    mutate(
      param = recode(cov, X = gamma100, W = gamma010, Z = gamma002)
    )
  
  results
  # calc_performance(results)

}
