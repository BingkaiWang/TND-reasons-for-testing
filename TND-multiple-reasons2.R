rm(list=ls())
set.seed(123)
library(tidyverse)
setwd("~/Downloads/")
########## simulation for effect modification from X1 and H.

n <- 1000 #total number of observations
n_sim <- 1000
disease_popularity <- 0.1
if(disease_popularity == 0.1) {scale_factor <- -2.1} else {scale_factor <- -0.6}

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
  d <- map_dbl(1:n, ~rbinom(1, 1, exp(scale_factor-0.5*v[.]*hsb[.] - 0.1*v[.]*x1[.]+0.5*x1[.]-x2[.]+0.1*case_tracing[.])))
  mean(d)
  reasons_T <- map_chr(1:n, function(i){
    if(case_tracing[i] == 1){
      TT <- ifelse(rbinom(1,1, 0.9+0.1*x2[i]), "case-tracing", NA)
    } else {
      u <- runif(1,0,1)
      temp <- (1- exp(scale_factor+0.5*x1[i]-x2[i]+0.1*case_tracing[i]))/(1- exp(scale_factor-0.5*hsb[i]+0.4*x1[i]-x2[i]+0.1*case_tracing[i]))
      symptom_factor <- (1-d[i]) * (1-v[i] + v[i]*temp) + d[i] * (1- exp(scale_factor+0.5*x1[i]-x2[i]+0.1*case_tracing[i]))
      if (u < 0.4 * hsb[i] * symptom_factor * x1[i]) {TT <- "symptoms"}
      else if(u < 0.4 * hsb[i] * symptom_factor * x1[i] + 0.2 * v[i] + 0.2 * x1[i]) {TT <- "disease-unrelated"}
      else if(u < 0.4 * hsb[i] * symptom_factor * x1[i] + 0.2 * v[i] + 0.2 * x1[i] + 0.4 * d[i] * (1-v[i])) {TT <- "other"}
      else {TT <- NA}
    }
    TT
  }) # the prob of T=1 & reasons = "a reason"
  x_eval <- seq(0.5, 0.99, by = 0.01)
  
  ############### data analysis for testing due to symptoms
  data1 <- data.frame(x1, x2, v, d, reasons_T) %>% filter(reasons_T == "symptoms")
  glm.fit1 <- glm(d~v*x1+x2, data = data1, family = "binomial")
  beta_v1 <- summary(glm.fit1)$coefficient[2,1]
  beta_vx1 <-  summary(glm.fit1)$coefficient[5,1]
  sd_beta_v1 <- summary(glm.fit1)$coefficient[2,2]
  re1 <- 1-exp(beta_v1 + beta_vx1 * x_eval)
  vcov1 <- vcov(glm.fit1)[c(2,5), c(2,5)]
  sd_re1 <- exp(beta_v1 + beta_vx1 * x_eval) * map_dbl(1: length(x_eval), ~sqrt(t(c(1, x_eval[.])) %*% vcov1 %*% c(1, x_eval[.])))
  
  ############### data analysis for testing due to disease-unrelated reasons
  data2 <- data.frame(x1, x2, v, d, reasons_T) %>% filter(reasons_T == "disease-unrelated")
  glm.fit2 <- glm(d~v*x1+x2, data = data2, family = binomial(link = "log"), start = rep(-1,5))
  beta_v2 <- summary(glm.fit2)$coefficient[2,1]
  beta_vx2 <-  summary(glm.fit2)$coefficient[5,1]
  sd_beta_v2 <- summary(glm.fit2)$coefficient[2,2]
  re2 <- 1-exp(beta_v2 + beta_vx2 * x_eval)
  vcov2 <- vcov(glm.fit2)[c(2,5), c(2,5)]
  sd_re2 <- exp(beta_v2 + beta_vx2 * x_eval) * map_dbl(1: length(x_eval), ~sqrt(t(c(1, x_eval[.])) %*% vcov2 %*% c(1, x_eval[.])))
  
  ############### data analysis for testing due to case-tracing
  data3 <- data.frame(x1, x2, v, d, reasons_T) %>% filter(reasons_T == "case-tracing")
  glm.fit3 <- glm(d~v*x1+x2, data = data3, family = binomial(link = "log"), start = rep(-1,5))
  beta_v3 <- summary(glm.fit3)$coefficient[2,1]
  beta_vx3 <-  summary(glm.fit3)$coefficient[5,1]
  sd_beta_v3 <- summary(glm.fit3)$coefficient[2,2]
  re3 <- 1-exp(beta_v3 + beta_vx3 * x_eval)
  vcov3 <- vcov(glm.fit3)[c(2,5), c(2,5)]
  sd_re3 <- exp(beta_v2 + beta_vx2 * x_eval) * map_dbl(1: length(x_eval), ~sqrt(t(c(1, x_eval[.])) %*% vcov3 %*% c(1, x_eval[.])))
  
  ############## data analysis pooling all reasons
  data_all <- data.frame(x1, x2, v, d, reasons_T) %>% filter(!is.na(reasons_T))
  glm.fitall <- glm(d~v*x1+x2, data = data_all, family = "binomial")
  beta_vall <- summary(glm.fitall)$coefficient[2,1]
  beta_vxall <-  summary(glm.fitall)$coefficient[5,1]
  sd_beta_vall <- summary(glm.fitall)$coefficient[2,2]
  reall <- 1-exp(beta_vall + beta_vxall * x_eval)
  vcovall <- vcov(glm.fitall)[c(2,5), c(2,5)]
  sd_reall <- exp(beta_vall + beta_vxall * x_eval) * map_dbl(1: length(x_eval), ~sqrt(t(c(1, x_eval[.])) %*% vcovall %*% c(1, x_eval[.])))
  
  ############## stratified analysis combining disease-unrelated reasons and case-tracing
  m2 <- nrow(data2)
  m3 <- nrow(data3)
  w2 <- m2/(m2+m3)
  w3 <- 1-w2
  ve_stratified1 <- w2 * re2 + w3 * re3
  sd_ve_stratified1 <- sqrt(w2^2 * sd_re2^2 + w3^2 * sd_re3^2)
  
  ############# stratified analysis combining all three reasons
  ve_stratifiedall <- (re1/sd_re1^2 + re2/sd_re2^2 + re3/sd_re3^2)/(1/sd_re1^2 + 1/sd_re2^2 + 1/sd_re2^2)
  sd_ve_stratifiedall <- sqrt(1/(1/sd_re1^2 + 1/sd_re2^2 + 1/sd_re3^2))
  
  
  ############ output 4 estimators
  c(re1, sd_re1, reall, sd_reall, ve_stratified1, sd_ve_stratified1, ve_stratifiedall, sd_ve_stratifiedall)
}
tictoc::toc()
stopCluster(cl)


plot_data <- data.frame(x = seq(0.5, 0.99, by = 0.01)) %>% 
  mutate(ve1_truth = 1-exp(-0.5-0.1 * x), ve_truth = 1-exp(-0.5-0.1*x) * x - exp(-0.1*x) * (1-x)) %>%
  mutate(ve_symp = apply(sim_result[1:50,],1, mean)) %>%
  mutate(ve_str1 = apply(sim_result[201:250,],1,mean)) %>%
  mutate(ve_symp_upper = map_dbl(1:50, ~mean(sim_result[.,] + qnorm(0.975, 0, sim_result[.+50,])))) %>%
  mutate(ve_symp_lower = map_dbl(1:50, ~mean(sim_result[.,] - qnorm(0.975, 0, sim_result[.+50,])))) %>%
  mutate(ve_str1_upper = map_dbl(201:250, ~mean(sim_result[.,] + qnorm(0.975, 0, sim_result[.+50,])))) %>%
  mutate(ve_str1_lower = map_dbl(201:250, ~mean(sim_result[.,] - qnorm(0.975, 0, sim_result[.+50,])))) 

ggplot(plot_data) +
  geom_path(aes(x = x, y = ve1_truth, linetype = "Truth VE(X,1)", color ="Truth VE(X,1)"), size = 1) +
  geom_path(aes(x = x, y = ve_symp, linetype = "Mean of estimates", color = "Mean of estimates"),size  =1) + 
  geom_ribbon(aes(x = x, ymin = ve_symp_lower, ymax=ve_symp_upper, fill = "95% C.I."), alpha = 0.3) +
  ylim(c(-0.5,1)) + xlim(c(0.5,1)) +
  theme_bw() + labs(colour = "Line", linetype = "Line", fill = "Fill") +
  theme(text = element_text(size=16), legend.position = c(0.65,0.15), legend.box = "horizontal") + labs(x = "X", y = "VE(X,1)") +
  scale_fill_manual(values = c("darkblue")) + scale_colour_manual(values = c("blue", "black")) 
ggsave(paste0("ve1-d", disease_popularity,".png"), dpi = 600, width = 7, height = 5)

ggplot(plot_data) +
  geom_path(aes(x = x, y = ve_truth, linetype = "Truth VE(X)", color ="Truth VE(X)"), size = 1) +
  geom_path(aes(x = x, y = ve_str1, linetype = "Mean of estimates", color = "Mean of estimates"),size  =1) + 
  geom_ribbon(aes(x = x, ymin = ve_str1_lower, ymax=ve_str1_upper, fill = "95% C.I."), alpha = 0.3) +
  ylim(c(-0.5,1)) + xlim(c(0.5,1)) +
  theme_bw() +
  theme(text = element_text(size=16), legend.position = c(0.65,0.15), legend.box = "horizontal") + labs(x = "X", y = "VE(X)") +
  labs(colour = "Line", linetype = "Line", fill = "Fill") +
  scale_fill_manual(values = c("red")) + scale_colour_manual(values = c("red", "black")) 
ggsave(paste0("ve-d", disease_popularity,".png"), dpi = 600, width = 7, height = 5)





