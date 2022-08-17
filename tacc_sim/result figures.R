# calc_performance
library(dplyr)
library(simhelpers)
library(knitr)

rm(list =ls())
load("sim_10rep.Rdata")
names(results)

# Table 4. rel_bias: intercept --------------------------------------------
results %>%
  mutate(params = 100) %>% 
  group_by(model, mobility, J, rho, tau) %>%
  group_modify(~ calc_relative(.x, estimates = int_mean, true_param = params,
                               perfm_criteria = c("relative bias"))) 

# Table 5. rel_bias: slope ------------------------------------------------
results %>%
  mutate(params = 10) %>% 
  group_by(model, mobility, J, rho, tau) %>%
  group_modify(~ calc_relative(.x, estimates = slope_mean, true_param = params,
                               perfm_criteria = c("relative bias"))) 

# Table 6. rel_bias: L1 residuals' SD -------------------------------------
results %>%
  mutate(params = 100) %>% 
  group_by(model, mobility, J, rho, tau) %>%
  group_modify(~ calc_relative(.x, estimates = sigma_y_mean, true_param = params,
                               perfm_criteria = c("relative bias")))

# Table 7. avg. rel_bias: L2 residuals' SD --------------------------
results %>%
  # dplyr::select(model, starts_with("sigma_u"), starts_with("rho")) %>%
  # dplyr::select(model, ends_with("median")) %>% 
  rowwise() %>% 
  mutate(median = median(c(sigma_u2_1_median, sigma_u2_2_median, 
                           sigma_u2_3_median, sigma_u2_4_median))) %>% 
  mutate(params = 100,
         sigma_u2_median = ifelse(model %in% c("ccrem", "rccrem"),
                                  median, sigma_u2_median)) %>% 
  group_by(model, mobility, J, rho, tau) %>%
  group_modify(~ calc_relative(.x, estimates = sigma_u2_median, true_param = params,
                               perfm_criteria = c("relative bias")))

# Table 8. avg. rel_bias: L2 residuals' total SD --------------------------
results %>%
  # dplyr::select(model, starts_with("sigma_u"), starts_with("rho")) %>% 
  # dplyr::select(model, ends_with("mean")) %>% 
  mutate(params = 100,
         sigma_u2_median = ifelse(model == "ccrem",
                                sqrt(sigma_u2_1_median^2 + sigma_u2_2_median^2 + sigma_u2_3_median^2 + sigma_u2_4_median^2),
                                sigma_u2_median)) %>% 
  mutate(sigma_u2_median = ifelse(model == "rccrem",
                                sqrt(sigma_u2_1_median^2 + sigma_u2_2_median^2 + sigma_u2_3_median^2 + sigma_u2_4_median^2 + 2*(rho_12_median + rho_13_median + rho_14_median + rho_23_median + rho_24_median + rho_34_median)),
                                sigma_u2_median)) %>% 
  group_by(model, mobility, J, rho, tau) %>%
  group_modify(~ calc_relative(.x, estimates = sigma_u2_median, true_param = params,
                               perfm_criteria = c("relative bias")))

# Table 9. median coverage: intercept -------------------------------------
# using group_modify()
welch_res %>%
  mutate(params = mean_diff) %>%
  group_by(n1, n2, mean_diff, method) %>%
  group_modify(~ calc_coverage(.x, lower_bound = lower_bound, upper_bound = upper_bound, true_param = params)) %>%
  kable(digits = 5)

results %>%
  # dplyr::select(model, starts_with("sigma_u"), starts_with("rho")) %>% 
  # dplyr::select(model, ends_with("mean")) %>% 
  mutate(params = 100,
         CI_level = .95,
         crit = qnorm((1 + CI_level) / 2),
         lower_bound = est - crit * se,
         upper_bound = est + crit * se,
         sigma_u2_median = ifelse(model == "ccrem",
                                  sqrt(sigma_u2_1_median^2 + sigma_u2_2_median^2 + sigma_u2_3_median^2 + sigma_u2_4_median^2),
                                  sigma_u2_median)) %>% 
  mutate(sigma_u2_median = ifelse(model == "rccrem",
                                  sqrt(sigma_u2_1_median^2 + sigma_u2_2_median^2 + sigma_u2_3_median^2 + sigma_u2_4_median^2 + 2*(rho_12_median + rho_13_median + rho_14_median + rho_23_median + rho_24_median + rho_34_median)),
                                  sigma_u2_median)) %>% 
  group_by(model, mobility, J, rho, tau) %>%
  group_modify(~ calc_relative(.x, estimates = sigma_u2_median, true_param = params,
                               perfm_criteria = c("relative bias")))
names(results)
head(results)
