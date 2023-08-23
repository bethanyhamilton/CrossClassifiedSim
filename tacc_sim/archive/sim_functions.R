# Models ---------------------------------------------------------
## mmrem_model ---------------------------------------------------
mmrem_model<- '
data {
  int<lower=1> n;  // number of observations
  real y[n];       // outcome
  int<lower=1> J;  // number of group level coefficients; parameters
  real x[n];
  matrix[n, 4] w;
  int group1[n];
  int group2[n];
  int group3[n];
  int group4[n];
}


parameters {
  real<lower=0> sigma_u2;
  real<lower=0> sigma_y;
  real popint;
  real popslope;
  vector[J] u2;
}


model {

  //precision is 1/variance. variance - 10,000
  popint ~ normal(100,100);
  popslope ~ normal(10,100);

  for (j in 1:J){
    u2[j] ~ normal(0,sigma_u2);
  }


  sigma_u2 ~ uniform(0,100);
  sigma_y ~ uniform(0, 100);

  for (i in 1:n){
  y[i] ~ normal(popint + popslope*x[i] + w[i,1] * u2[group1[i]] + w[i,2] * u2[group2[i]] + w[i,3] * u2[group3[i]] + w[i,4] * u2[group4[i]], sigma_y);
  } 
}


generated quantities {
vector[n] log_lik;
for (i in 1:n) {
log_lik[i] = normal_lpdf(y[i]| popint + popslope*x[i] + w[i,1] * u2[group1[i]] + w[i,2] * u2[group2[i]] + w[i,3] * u2[group3[i]] + w[i,4] * u2[group4[i]], sigma_y); 
// loop to run model prior to running Loo. 
}
}


'

## ccrem_model ---------------------------------------------------
ccrem_model<- '
data {
int<lower=1> n;  // number of observations
real y[n]; // observations
int<lower=1> J; // number of group level coefficients; parameters
real x[n];
int group1[n];
int group2[n];
int group3[n];
int group4[n];
}


parameters {
vector<lower=0>[4] sigma_u2;
real<lower=0> sigma_y;
real popint;
real popslope;
matrix[J, 4] u2;
//vector[4] u2[J]; // use this if you are doing it the other way.
}

model {

popint ~ normal(100,100);
popslope ~ normal(10,100);
  
  
for (j in 1:J){
for (k in 1:4){
u2[j,k] ~ normal(0,sigma_u2[k]);
}
}

sigma_u2 ~ uniform(0,100);
sigma_y ~ uniform(0, 100);


for (i in 1:n){
y[i] ~ normal(popint + popslope*x[i] + u2[group1[i],1] + u2[group2[i],2] + u2[group3[i],3] + u2[group4[i],4], sigma_y);
}


}

generated quantities {
vector[n] log_lik;
for (i in 1:n) {
log_lik[i] = normal_lpdf(y[i]| popint + popslope*x[i] + u2[group1[i],1] + u2[group2[i],2] + u2[group3[i],3] + u2[group4[i],4], sigma_y); 
// loop to run model prior to running Loo. 
}
}
'
## rccrem_model --------------------------------------------------
rccrem_model_new <- 
  'data {
  int<lower=1> n; // number of observations
  int<lower=1> J; //the number of levels of each type of school
  real x[n];
  int group1[n];
  int group2[n];
  int group3[n];
  int group4[n];
  matrix[4, 4] W;
  real y[n];      // observations
  vector<lower=0>[4] zeros;
}

parameters {
  real popint;
  real popslope;
  real<lower=0> sigma_y;
  vector<lower=0>[4] xi_u2;
  matrix[J,4] B_raw;
  cov_matrix[4] Tau_B_raw; // Sigma -> Tau_B_raw 
}

transformed parameters {
  matrix[J, 4] u2;
  cov_matrix[4] Sigma_B_raw;
  
  for (j in 1:J){
        for (k in 1:4){
             u2[j,k] = xi_u2[k]*B_raw[j,k];
        }
    }
  Sigma_B_raw = inverse(Tau_B_raw);
}



model {
  popint ~ normal(100,100);
  popslope ~ normal(10,100);
  sigma_y ~  uniform(0, 100);
  
  for (k in 1:4){
    xi_u2[k] ~ uniform(0, 100);
  }
  
  for (j in 1:J){
    B_raw[j,1:4] ~ multi_normal(zeros, Tau_B_raw); //u2 -> B_raw
  }
  
  Tau_B_raw ~ wishart(5, W);

  for (i in 1:n){
    y[i] ~ normal(popint + popslope*x[i] + u2[group1[i],1] + u2[group2[i],2] + u2[group3[i],3] + u2[group4[i],4], sigma_y);
  }
  
}

generated quantities {
  
  vector<lower=0>[4] sigma_u2;
  real<lower=-1,upper=1> rho_12;
  real<lower=-1,upper=1> rho_13;
  real<lower=-1,upper=1> rho_14;
  real<lower=-1,upper=1> rho_23;
  real<lower=-1,upper=1> rho_24;
  real<lower=-1,upper=1> rho_34;
  vector[n] log_lik;
  
  for (k in 1:4){
    sigma_u2[k] = xi_u2[k]*sqrt(Sigma_B_raw[k,k]);
  } 
  rho_12 = Sigma_B_raw[1,2]/sqrt(Sigma_B_raw[1,1]*Sigma_B_raw[2,2]);
  rho_13 = Sigma_B_raw[1,3]/sqrt(Sigma_B_raw[1,1]*Sigma_B_raw[3,3]);
  rho_14 = Sigma_B_raw[1,4]/sqrt(Sigma_B_raw[1,1]*Sigma_B_raw[4,4]);
  rho_23 = Sigma_B_raw[2,3]/sqrt(Sigma_B_raw[2,2]*Sigma_B_raw[3,3]);
  rho_24 = Sigma_B_raw[2,4]/sqrt(Sigma_B_raw[2,2]*Sigma_B_raw[4,4]);
  rho_34 = Sigma_B_raw[3,4]/sqrt(Sigma_B_raw[3,3]*Sigma_B_raw[4,4]);
    for (i in 1:n) {
    log_lik[i] = normal_lpdf(y[i]|popint + u2[group1[i],1] + u2[group2[i],2] + u2[group3[i],3] + u2[group4[i],4] + popslope*x[i], sigma_y); 
  }
}
'
# Data Generation ------------------------------------------------
UDF_sample <- function(mobility, J, rho, tau) {
  require(MASS)
  
  m <- 10
  mu.X_init <- 0
  sigma2.X_init <- 1
  mu.u2_init <- rep(0,4)
  sigma2_init <- 10
  beta_init <- c(100,10)
  n<-m*J
  
  #use rho and tau
  
  if(tau == "same" & rho== FALSE){
    sigma2.u2_init <- matrix(0, 4, 4)  # Overwrite the old matrix
    diag(sigma2.u2_init) <-c(rep(10, 4))
    
  }else if(tau == "different" & rho== FALSE){
    
    sigma2.u2_init <- matrix(0, 4, 4)  # Overwrite the old matrix
    diag(sigma2.u2_init) <- c(10, 12, 14, 16)
    
    
  }else if(tau == "same" & rho== TRUE){
    
    rho <- rbind(c(10,6,5,4),c(6,10,6,5),c(5,6,10,6),c(4,5,6,10))*.1
    tau1 <- c(rep(10, 4))
    sigma2.u2_init <- sqrt(tau1)%*%t(sqrt(tau1))*rho
    
    
  }else{
    
    rho <- rbind(c(10,6,5,4),c(6,10,6,5),c(5,6,10,6),c(4,5,6,10))*.1
    tau2 <- c(10, 12, 14, 16)
    
    sigma2.u2_init <- sqrt(tau2)%*%t(sqrt(tau2))*rho
    
  }
  
  
  
  
  #Step 1: Generate data through sampling
  u2_samp <- mvrnorm(J, mu.u2_init, Sigma=sigma2.u2_init)
  l2_samp <- matrix(nrow=J,ncol=length(mu.u2_init)+1); colnames(l2_samp) <- c("ID",paste("u",1:length(mu.u2_init),sep='')); l2_samp <- data.frame(l2_samp)
  l2_samp[,"ID"]<-1:J; l2_samp[,-1]<-u2_samp
  
  X_samp <- rnorm(n, mu.X_init, sqrt(sigma2.X_init))
  e_samp <- rnorm(n,0, sqrt(sigma2_init))
  
  l1_samp <- matrix(nrow=n, ncol=3); colnames(l1_samp) <- c("ID", "X", "e");l1_samp <- data.frame(l1_samp)
  l1_samp[,"ID"]<-1:n; l1_samp[,"X"]<-X_samp; l1_samp[,"e"]<-e_samp
  
  l1_assign <- matrix(trunc(runif(n*length(mu.u2_init), min=1, max=J+1)),nrow=n)
  l1_l2_xref <-matrix(nrow=n, ncol=length(mu.u2_init)+1); colnames(l1_l2_xref) <- c("l1_ID",paste("l2_ID",1:length(mu.u2_init),sep='')); l1_l2_xref <- data.frame(l1_l2_xref)
  l1_l2_xref[,"l1_ID"]<-1:n; l1_l2_xref[,2]<-l1_assign[,1];
  l1_l2_xref[1:trunc(mobility*n),-c(1,2)]<-l1_assign[1:trunc(mobility*n),-1];l1_l2_xref[(trunc(mobility*n)+1):n,-c(1,2)]<-l1_assign[(trunc(mobility*n)+1):n,1]
  
  data_samp <- matrix(nrow=n,ncol=length(mu.u2_init)+2); colnames(data_samp) <- c(paste("group",1:length(mu.u2_init),sep=''),"y","X"); data_samp <- data.frame(data_samp)
  data_samp[,1:length(mu.u2_init)]<-l1_l2_xref[,-1]
  data_samp[,"X"]<-l1_samp$X
  
  data_samp[,"y"] <- beta_init[1] + l2_samp$u1[data_samp[,"group1"]] + l2_samp$u2[data_samp[,"group2"]] + l2_samp$u3[data_samp[,"group3"]] + l2_samp$u4[data_samp[,"group4"]] +
    beta_init[2]*data_samp[,"X"] + l1_samp$e
  
  return(data_samp)
}


#tm <- system.time(mydata <- UDF_sample(mobility = c(.2), J = c(100),rho = c(FALSE), tau = c("same")))


# Make data list ----------------------------------------------------------
getdata_ccrem_mmrem <- function(data){
  # pre_data
  k = 4
  id <- c(paste("group",1:k,sep=""))
  gm = 1
  predictor="X"
  
  # data 
  y <- data$y
  J <- length(levels(as.factor(data[,id[1]])))
  n <- dim(data)[1]
  
  if(gm == 1)
  {
    x <- data[,predictor] - mean(data[,predictor])
  } else
  {
    x <- data[,predictor]
  }
  
  W <- diag(k)
  w <- matrix(rep(1/k,n*k),ncol=k, nrow=n)
  
  for(i in 1:k) {assign(paste("group",i,sep=''),data[,id[i]])}
  
  for (i in 1:n) {
    
    if (group1[i]==group2[i]){
      w[i,1]=w[i,1]+0.25
      w[i,2]=0
    } 
    
    if (group1[i]==group3[i]){
      w[i,1]=w[i,1]+0.25
      w[i,3]=0
    } 
    else 
      if (group2[i]==group3[i]){
        w[i,2]=w[i,2]+0.25
        w[i,3]=0
      } 
    
    
    if (group1[i]==group4[i]){
      w[i,1]=w[i,1]+0.25
      w[i,4]=0
    } 
    else {
      if (group2[i]==group4[i]){
        w[i,2]=w[i,2]+0.25
        w[i,4]=0
      } 
      
      else 
        if (group3[i]==group4[i]){
          w[i,3]=w[i,3]+0.25
          w[i,4]=0
        }
    }
    
  }
  
  group1 <- data$group1
  group2 <- data$group2
  group3 <- data$group3
  group4 <- data$group4
  
  ccrem_mmrem_data <- list(y = y, J = J, n = n, w = w, x = x, 
                           group1 = group1, group2 = group2, 
                           group3 = group3, group4 = group4)
  
  return(ccrem_mmrem_data)
  
}

#system.time(listdata <- getdata_ccrem_mmrem(mydata))
getdata_rccrem <- function(data){
  #Check on this
  k = 4
  id <- c(paste("group",1:k,sep=""))
  gm = 1
  predictor="X"
  
  ####
  
  y <- data$y
  
  J <- length(levels(as.factor(data[,id[1]])))
  n <- dim(data)[1]
  
  if(gm == 1)
  {
    x <- data[,predictor] - mean(data[,predictor])
  } else
  {
    x <- data[,predictor]
  }
  
  W <- diag (k)
  
  
  for(i in 1:k) {assign(paste("group",i,sep=''),data[,id[i]])}
  
  group1 <- data$group1
  group2 <- data$group2
  group3 <- data$group3
  group4 <- data$group4
  
  zeros = rep(0, 4)
  
  rccrem_data <- list(y=y, J=J, n=n, group1=group1, group2=group2, group3=group3, group4=group4, x= x, W=W, zeros = zeros) 
  
  return(rccrem_data)
  
}
#system.time(listdata <- getdata_rccrem(mydata))

# Run Models Funcs --------------------------------------------------------
run_mmrem <-
  function(data,
           mobility, J, rho, tau,
           iterations = 1000,
           warmup_val = 500,
           thin_val = 10,
           chains_val = 4,
           cores_val = 6,
           adapt_delta = 0.999,
           max_treedepth = 25) {
    
    
    ccrem_mmrem_data <- getdata_ccrem_mmrem(data)
    
    mm_stan_fit <-
      stan(
        model_code = mmrem_model,
        data = ccrem_mmrem_data,
        iter = iterations,
        warmup = warmup_val,
        thin = thin_val,
        chains = chains_val,
        cores = cores_val,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
      )
    
    
    # LOO 
    mm_loo <- loo(mm_stan_fit, k_threshold = 0.7, save_psis = TRUE)
    
    
    LL_mm <- extract_log_lik(mm_stan_fit)
    mm_waic <- waic(LL_mm)
    
    summary <- summary(mm_stan_fit)$summary
    summary <- summary[1:4,]
    summary <- summary %>% as.data.frame() %>% tibble::rownames_to_column("paramtype") 
      
    return(
    summary %>% 
          tidyr::pivot_longer(!paramtype, names_to = "name", values_to = "value") %>% 
          mutate(paramtype = paste(paramtype, name, sep = "_")) %>% 
          tidyr::pivot_wider(!name, names_from = "paramtype", values_from = "value") %>% 
          mutate(model = "mmrem",
                 mobility = mobility, J = J, rho = rho, tau = tau) %>% 
          dplyr::select(model, mobility, J, rho, tau, starts_with("popint"),
                        starts_with("popslope"), starts_with("sigma_y"),
                        starts_with("sigma_u")) %>% 
          mutate(loo_est = mm_loo$estimates["looic", "Estimate"],
                 waic_est = mm_waic$estimates["waic", "Estimate"]))
  }


run_ccrem <-
  function(data,
           mobility, J, rho, tau,
           iterations = 1000,
           warmup_val = 500,
           thin_val = 10,
           chains_val = 4,
           cores_val = 6,
           adapt_delta = 0.999,
           max_treedepth = 25) {
    
    
    ccrem_mmrem_data <- getdata_ccrem_mmrem(data)
    
    
    cc_stan_fit <-
      stan(
        model_code = ccrem_model,
        data = ccrem_mmrem_data,
        iter = iterations,
        warmup = warmup_val,
        thin = thin_val,
        chains = chains_val,
        cores = cores_val,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
      )
    
    # LOO 
    cc_loo <- loo(cc_stan_fit, k_threshold = 0.7, save_psis = TRUE)
    
    
    LL_cc <- extract_log_lik(cc_stan_fit)
    cc_waic <- waic(LL_cc)
    
    summary <- summary(cc_stan_fit)$summary
    summary <- summary[1:7,]
    summary <- summary %>% as.data.frame() %>% tibble::rownames_to_column("paramtype") 
    
    return(
      summary %>% 
        tidyr::pivot_longer(!paramtype, names_to = "name", values_to = "value") %>% 
        mutate(paramtype = paste(paramtype, name, sep = "_")) %>% 
        tidyr::pivot_wider(!name, names_from = "paramtype", values_from = "value") %>% 
        mutate(model = "ccrem",
               mobility = mobility, J = J, rho = rho, tau = tau) %>% 
        dplyr::select(model, mobility, J, rho, tau, starts_with("popint"),
                      starts_with("popslope"), starts_with("sigma_y"),
                      starts_with("sigma_u")) %>% 
        mutate(loo_est = cc_loo$estimates["looic", "Estimate"],
               waic_est = cc_waic$estimates["waic", "Estimate"]))
  }

run_rccrem <-
  function(data,
           mobility, J, rho, tau,
           iterations = 1000,
           warmup_val = 500,
           thin_val = 10,
           chains_val = 4,
           cores_val = 6,
           adapt_delta = 0.999,
           max_treedepth = 25) {
    
    rccrem_data <- getdata_rccrem(data)
    
    
    rcc_stan_fit <-
      stan(
        model_code = rccrem_model_new,
        data = rccrem_data,
        iter = iterations,
        warmup = warmup_val,
        thin = thin_val,
        chains = chains_val,
        cores = cores_val,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
      )
    
    
    # LOO 
    rcc_loo <- loo(rcc_stan_fit, k_threshold = 0.7, save_psis = TRUE)
    
    LL_rcc <- extract_log_lik(rcc_stan_fit)
    rcc_waic <- waic(LL_rcc)
    
    summary <- summary(rcc_stan_fit)$summary
    summary <- 
      bind_rows(summary["popint",], summary["popslope",], summary["sigma_y",],
              summary["sigma_u2[1]",], summary["sigma_u2[2]",],
              summary["sigma_u2[3]",], summary["sigma_u2[4]",],
              summary["rho_12",], summary["rho_13",], summary["rho_13",],
              summary["rho_23",], summary["rho_24",], summary["rho_34",]) %>% 
      mutate(paramtype = c("popint", "popslope", "sigma_y", 
                           "sigma_u2[1]", "sigma_u2[2]", "sigma_u2[3]",
                           "sigma_u2[4]", "rho_12", "rho_13", "rho_14",
                           "rho_23", "rho_24", "rho_34")) %>% 
      dplyr::select(paramtype, everything())
    
    
    return(
      summary %>% 
        tidyr::pivot_longer(!paramtype, names_to = "name", values_to = "value") %>% 
        mutate(paramtype = paste(paramtype, name, sep = "_")) %>% 
        tidyr::pivot_wider(!name, names_from = "paramtype", values_from = "value") %>% 
        mutate(model = "rccrem",
               mobility = mobility, J = J, rho = rho, tau = tau) %>% 
        dplyr::select(model, mobility, J, rho, tau, starts_with("popint"),
                      starts_with("popslope"), starts_with("sigma_y"),
                      starts_with("sigma_u"), starts_with("rho")) %>% 
        mutate(loo_est = rcc_loo$estimates["looic", "Estimate"],
               waic_est = rcc_waic$estimates["waic", "Estimate"]))
    
  }



# Run Sims ----------------------------------------------------------------
runmodels <- function(mobility = mobility, J = J, rho=rho, tau=tau, 
                      mmrem=TRUE, ccrem=TRUE, rccrem=TRUE){
  suppressPackageStartupMessages(require(dplyr, quietly = TRUE, warn.conflicts = FALSE))
  
  suppressPackageStartupMessages(require(rstan, quietly = TRUE, warn.conflicts = FALSE))
  
  suppressPackageStartupMessages(require(loo, quietly = TRUE, warn.conflicts = FALSE))
  
  data <- UDF_sample(mobility = mobility, J = J, rho=rho, tau=tau)
  
  res <- data.frame()
  
  
  if (mmrem) {
    res_mmrem <- run_mmrem(data, mobility = mobility, J = J, rho = rho, tau = tau)
    
    res <- bind_rows(list(res, res_mmrem))
    
  }
  
  if (ccrem) {
    res_ccrem <- run_ccrem(data, mobility = mobility, J = J, rho = rho, tau = tau)
    
    res <- bind_rows(list(res, res_ccrem))
    
  }
  
  if (rccrem) {
    res_rccrem <- run_rccrem(data, mobility = mobility, J = J, rho = rho, tau = tau)
    
    res <- bind_rows(list(res, res_rccrem))
    
    
  }
  
  res
}



#res <- runmodels(data= mydata)

runSim <- function(reps, mobility, J, rho, tau, mmrem=TRUE, ccrem=TRUE, rccrem=TRUE, seed = NULL, ...){
  
  suppressPackageStartupMessages(require(purrr, quietly = TRUE, warn.conflicts = FALSE))
  suppressPackageStartupMessages(require(dplyr, quietly = TRUE, warn.conflicts = FALSE))
  
  if (!is.null(seed)) set.seed(seed)
  
  replicates <- rerun(reps, {
    runmodels(mobility = mobility, J = J, rho=rho, tau=tau,
              mmrem = mmrem, ccrem = ccrem, rccrem= rccrem)
  }) %>% 
    dplyr::bind_rows() 
  
  replicates
}

# system.time(
# results <- runSim(reps=1, mobility = c(.2), J = c(100),rho = c(FALSE), tau = c("same"))
# )

# convergence_diagnostics
# results %>% dplyr::select(ends_with("Rhat"))
# 
# results %>% dplyr::select(ends_with("sd"), ends_with("eff")) %>% 
#   mutate(popint_mcse = popint_sd/sqrt(popint_n_eff),
#          popslope_mcse = popslope_sd/sqrt(popslope_n_eff),
#          sigma_y_mcse = sigma_y_sd/sqrt(sigma_y_n_eff))
