rm(list=ls())
set.seed(123)
library(tidyverse)
########## simulation for no effect modification from X or H.

n <- 10000 #total number of observations
n_sim <- 1000
scale_factor <- -2.1 # a factor to adjust the model for d such that P(d=1) = 0.1 (with value -2.1) or 0.45 (with value -0.6)

library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)
tictoc::tic()

sim_result <- foreach(k = 1:n_sim, .combine = cbind, .packages = c("tidyverse")) %dopar% {
  x1 <- runif(n,0.5,1)
  x2 <- rbinom(n,1,0.5)
  hsb <- rbinom(n, 1, x1)
  v <- map_dbl(1:n, ~rbinom(1, 1, 0.6*x1[.]+0.1*x2[.]))
  case_tracing <- rbinom(n, 1, 0.05 + 0.05 * x2)
  d <- map_dbl(1:n, ~rbinom(1, 1, exp(scale_factor-0.5*v[.]+0.5*x1[.]-x2[.]+0.1*case_tracing[.])))
  # d <- map_dbl(1:n, ~ifelse(rbinom(1,1,0.1), 1-d[.], d[.])) # misclasification rate of 10%

  reasons_T <- map_chr(1:n, function(i){
    if(case_tracing[i] == 1){
      TT <- ifelse(rbinom(1,1, 0.9+0.1*x2[i]), "case-tracing", NA)
    } else {
      u <- runif(1,0,1)
      temp <- (1- exp(scale_factor+0.5*x1[i]-x2[i]+0.1*case_tracing[i]))/(1- exp(scale_factor-0.5+0.5*x1[i]-x2[i]+0.1*case_tracing[i]))
      symptom_factor <- (1-d[i]) * (1-v[i] + v[i]*temp) + d[i] * (1- exp(scale_factor+0.5*x1[i]-x2[i]+0.1*case_tracing[i]))
      if (u < 0.4 * hsb[i] * symptom_factor * x1[i]) {TT <- "symptoms"}
      else if(u < 0.4 * hsb[i] * symptom_factor * x1[i] + 0.2 * v[i] + 0.2 * x1[i]) {TT <- "disease-unrelated"}
      else if(u < 0.4 * hsb[i] * symptom_factor * x1[i] + 0.2 * v[i] + 0.2 * x1[i] + 0.4 * d[i] * (1-v[i])) {TT <- "other"}
      else {TT <- NA}
    }
    TT
  }) # the prob of T=1 & reasons = "a reason"
  # reasons_T <- map_chr(1:n, ~ifelse(rbinom(1,1,0.2), NA, reasons_T[.])) # 20% missing reasons for testing
  
  ############### data analysis for testing due to symptoms
  data1 <- data.frame(x1, x2, v, d, reasons_T) %>% filter(reasons_T == "symptoms")
  glm.fit1 <- glm(d~v+x1+x2, data = data1, family = "binomial")
  beta_v1 <- summary(glm.fit1)$coefficient[2,1]
  sd_beta_v1 <- summary(glm.fit1)$coefficient[2,2]
  ve_symptom <- 1 - exp(beta_v1)
  sd_ve_symptom <- exp(beta_v1) * sd_beta_v1
  
  ############### data analysis for testing due to disease-unrelated reasons
  data2 <- data.frame(x1, x2, v, d, reasons_T) %>% filter(reasons_T == "disease-unrelated")
  glm.fit2 <- glm(d~v+x1+x2, data = data2, family = binomial(link = "log"), start = rep(-1,4))
  beta_v2 <- summary(glm.fit2)$coefficient[2,1]
  sd_beta_v2 <- summary(glm.fit2)$coefficient[2,2]
  ve_unrelated <- 1 - exp(beta_v2)
  sd_ve_unrelated <- exp(beta_v2) * sd_beta_v2
  
  ############### data analysis for testing due to case-tracing
  data3 <- data.frame(x1, x2, v, d, reasons_T) %>% filter(reasons_T == "case-tracing")
  glm.fit3 <- glm(d~v+x1+x2, data = data3, family = binomial(link = "log"), start = rep(-1,4))
  beta_v3 <- summary(glm.fit3)$coefficient[2,1]
  sd_beta_v3 <- summary(glm.fit3)$coefficient[2,2]
  ve_casetracing <- 1 - exp(beta_v3)
  sd_ve_casetracing <- exp(beta_v3) * sd_beta_v3
  
  ############## data analysis pooling all reasons
  data_all <- data.frame(x1, x2, v, d, reasons_T) %>% filter(!is.na(reasons_T))
  glm.fitall <- glm(d~v+x1+x2, data = data_all, family = "binomial")
  beta_vall <- summary(glm.fitall)$coefficient[2,1]
  sd_beta_vall <- summary(glm.fitall)$coefficient[2,2]
  ve_all <- 1 - exp(beta_vall)
  sd_ve_all <- exp(beta_vall) * sd_beta_vall
  
  ############## stratified analysis combining disease-unrelated reasons and case-tracing
  ve_stratified1 <- (ve_unrelated/sd_ve_unrelated^2 + ve_casetracing/sd_ve_casetracing^2)/(1/sd_ve_unrelated^2 + 1/sd_ve_casetracing^2)
  sd_ve_stratified1 <- sqrt(1/(1/sd_ve_unrelated^2 + 1/sd_ve_casetracing^2))
  
  ############# stratified analysis combining all three reasons
  ve_stratifiedall <- (ve_symptom/sd_ve_symptom^2 + ve_unrelated/sd_ve_unrelated^2 + ve_casetracing/sd_ve_casetracing^2)/(1/sd_ve_symptom^2 + 1/sd_ve_unrelated^2 + 1/sd_ve_casetracing^2)
  sd_ve_stratifiedall <- sqrt(1/(1/sd_ve_symptom^2 + 1/sd_ve_unrelated^2 + 1/sd_ve_casetracing^2))
  
  
  ############ output 4 estimators
  c(ve_symptom, sd_ve_symptom, ve_all, sd_ve_all, ve_stratified1, sd_ve_stratified1, ve_stratifiedall, sd_ve_stratifiedall)
}
tictoc::toc()
stopCluster(cl)

est_index <- c(1,3,5,7)
sd_index <- c(2,4,6,8)
truth <- 1-exp(-0.5)
summary_table <- data.frame(bias = apply(sim_result[est_index,], 1, mean) - truth,
                            ese = apply(sim_result[est_index,], 1, sd),
                            ase = apply(sim_result[sd_index,], 1, mean),
                            cp = apply(abs((sim_result[est_index,]-truth)/sim_result[sd_index,]) <= qnorm(0.975), 1, mean)
)
rownames(summary_table) <- c("symptoms", "pooling", "stratifeid1", "stratified_all")

round(summary_table,2)
xtable::xtable(summary_table)
