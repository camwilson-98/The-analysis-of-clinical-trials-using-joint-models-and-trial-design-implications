#### Assumption 7: follow-up schedule ####
# Simulate the power for different follow-up schedules (vary length and number of visits)
# Simulate a sufficiently complex longitudinal trajectory. i.e. a quartic polynomial

# Efficient way of installing previously not installed packages
list.of.packages <- c("tidyverse","survival","survminer","lme4",
                      "ggplot2","JM","joineR","simsurv","flexsurv",
                      "MASS","nlme","msm","tictoc","foreach","doParallel",
                      "parallel", "rlecuyer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))suppressPackageStartupMessages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

# Set Working Directory
setwd("~ /OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R")

# Load libraries
library(simsurv); library(flexsurv); library(survival); library(survminer);
library(dplyr); library(MASS); library(nlme); library(JM); 
library(msm); library(tictoc); library(MASS)

haz <- function(t, x, betas, ...) {
  betas[["shape"]] * exp(
    betas[["betaEvent_intercept"]] +
      betas[["betaEvent_binary"]] * x[["Z1"]] +
      betas[["betaEvent_assoc"]] * (
        betas[["betaLong_intercept"]] +
          betas[["betaLong_slope"]] * t +
          betas[["betaLong_quadratic"]] * (t^2) +
          betas[["betaLong_cubic"]] * (t^3) +
          betas[["betaLong_quartic"]] * (t^4) +
          betas[["betaLong_binary"]] * x[["Z1"]]
      )
  )
}

# Function to carry out the simulations
sim4_run <- function(N, p1, beta, gamma, alpha, lambda, a, f,
                     mu_theta, covar_theta, sigmasq_e) {
  # Then we construct data frames with the true parameter
  # values and the covariate data for each individual
  
  # Population (fixed effect) parameters
  betas <- data.frame(
    shape                = rep(lambda, N),
    betaEvent_intercept  = rep(0, N),
    betaEvent_binary     = rep(alpha, N),
    betaEvent_assoc      = rep(beta, N),
    betaLong_binary      = rep(gamma, N),
    betaLong_intercept   = rep(0, N),
    betaLong_slope       = rep(0, N),
    betaLong_quadratic   = rep(0, N),
    betaLong_cubic       = rep(0, N),
    betaLong_quartic     = rep(0, N)
  )
  
  # Individual-specific (random effect) parameters
  p <- 4
  b_corrmat <- matrix(c(rep(NA, (p+1)^2)), nrow = p+1, ncol = p+1)
  for (i in c(1:length(mu_theta))) {
    for (j in c(1:length(mu_theta))) {
      b_corrmat[i,j] <- covar_theta[i,j]/sqrt(covar_theta[i,i]*covar_theta[j,j])
    }
  }
  b_sds     <- sqrt(diag(covar_theta))
  b_means   <- mu_theta
  b_z       <- MASS::mvrnorm(n = N, mu = b_means, Sigma = b_corrmat)
  b         <- sapply(1:length(b_sds), function(x) b_sds[x] * b_z[,x])
  betas$betaLong_intercept <- betas$betaLong_intercept + b[,1]
  betas$betaLong_slope     <- betas$betaLong_slope     + b[,2]
  betas$betaLong_quadratic <- betas$betaLong_quadratic + b[,3]
  betas$betaLong_cubic     <- betas$betaLong_cubic     + b[,4]
  betas$betaLong_quartic   <- betas$betaLong_quartic   + b[,5]
  
  # Covariate data
  covdat <- data.frame(
    Z1 = as.integer(c(rep(0,N*p1), rep(1,N*(1-p1)))) # a binary covariate
  )
  
  # Then we simulate the survival times based on the
  # hazard function, covariates, and true parameter values
  times <- simsurv(hazard = haz, x = covdat, betas = betas, maxt = 2,
                   nodes = 7, interval = c(0,2.1), rootsolver = "uniroot")
  times$eventtime <- ifelse(times$eventtime == 0, 0.0001, times$eventtime)
  
  # Simulate observed biomarker data at every month
  sim_dat <- cbind(times, covdat, betas) %>% 
    dplyr::select(-shape, -betaEvent_intercept, -betaEvent_binary,-betaEvent_assoc, -betaLong_binary) %>%
    dplyr::slice(rep(1:n(), 25)) %>%
    dplyr::arrange(id) %>%
    dplyr::rename(trt = Z1) %>%
    dplyr::mutate(obstime = as.numeric(rep(seq(0,2,length.out = 25), length(unique(id)))),
                  p = rep(p,N*25))
  
  sim_dat <- sim_dat %>% dplyr::mutate(theta0 = betaLong_intercept,
                                       theta1 = betaLong_slope,
                                       theta2 = betaLong_quadratic,
                                       theta3 = betaLong_cubic,
                                       theta4 = betaLong_quartic,
                                       biomarker = theta0 + theta1*obstime + theta2*(obstime^2) + theta3*(obstime^3) + theta4*(obstime^4) + 
                                         gamma*trt + rnorm(n = N*25, mean = 0, sd = sqrt(sigmasq_e))) %>%
    dplyr::select(id, trt, obstime, biomarker, eventtime, status) %>%
    dplyr::ungroup() 
  
  # Apply uniform censoring
  sim_dat <- sim_dat %>% group_by(id) %>%
    dplyr::mutate(cnsrtime = runif(n = 1, min = f, max = a+f),
                  status = case_when(
                    (status == 0) | (cnsrtime < eventtime) ~ 0,
                    (status == 1) & (cnsrtime >= eventtime) ~ 1),
                  eventtime = min(eventtime, cnsrtime)) %>%
    dplyr::filter(obstime <= eventtime) %>% # remove values after event time, as biomarker is endogenous
    dplyr::ungroup()
    
  # Split into data sets with varying spacing of longitudinal follow-up visit 
  monthly.intervals <- c(1,3,6,12)
  outdata <- c()
  
  for (mi in monthly.intervals){
    
    num.intervals <- (12*(a+f)/mi) + 1
    sim_dat_split <- sim_dat %>% dplyr::filter(obstime %in% seq(0,2,length.out = num.intervals))
    
    # Fit a mixed-effects longitudinal model
    library(nlme)
    if (mi == 12) {
      
      lmeFit_split <- nlme::lme(biomarker ~ obstime + trt, data = sim_dat_split, random = ~ obstime | id,
                                na.action = na.omit,
                                control = lmeControl(maxIter = 50, msMaxIter = 50, niterEM = 25, 
                                                     tolerance = 1e-08, opt = "optim"))
      
    } else {
      
      lmeFit_split <- nlme::lme(biomarker ~ obstime + I(obstime^2) + I(obstime^3) + trt, data = sim_dat_split, random = ~ obstime | id,
                                na.action = na.omit,
                                control = lmeControl(maxIter = 50, msMaxIter = 50, niterEM = 25, 
                                                     tolerance = 1e-08, opt = "optim"))
      
    }
    
    # Cox model
    library(survival)
    dat.id_split <- sim_dat_split %>% dplyr::distinct(id, .keep_all = TRUE)
    survFit_split <- survival::coxph(Surv(eventtime, status) ~ trt, data = dat.id_split, 
                                     x = TRUE, na.action = na.omit)
    
    # Joint model - Maximum likelihood framework
    library(JM)
    jointFit <- jointModel(lmeFit_split, survFit_split, timeVar = "obstime", 
                           method = "Cox-PH-GH",
                           control = list(iter.qN = 5, iter.EM = 5, 
                                          tol1 = 1e-08, tol2 = 1e-08, 
                                          GHk = 5, GKk = 7))
    
    # Obtain estimates, standard errors and Wald test p-value for overall treatment effect
    alpha_est <- summary(jointFit)[["CoefTable-Event"]]["trt","Value"] 
    beta_est <- summary(jointFit)[["CoefTable-Event"]]["Assoct","Value"]
    gamma_est <- summary(jointFit)[["CoefTable-Long"]]["trt","Value"]
    
    est <- alpha_est + beta_est*gamma_est
    
    mean <- c(alpha_est, beta_est, gamma_est)
    cov <- vcov(jointFit)[c("T.trt","T.alpha","Y.trt"), c("T.trt","T.alpha","Y.trt")]
    names(mean) <- c("alpha_est", "beta_est", "gamma_est")
    colnames(cov) <- c("alpha_est", "beta_est", "gamma_est")
    rownames(cov) <- c("alpha_est", "beta_est", "gamma_est")
    ses <- deltamethod(~ x1 + x2*x3, mean = mean, cov = cov, ses = TRUE) # calculate SE for non-linear combination of parameters using delta method
    
    cil <- est + qnorm(.025) * ses
    ciu <- est + qnorm(.975) * ses
    wald <- est/ses
    #pval <- 2*pnorm(-abs(wald))
    pval <- pnorm(wald)
    
    bias <- est - (alpha + beta*gamma)
    power <- (as.numeric(pval) < 0.05)
    
    outdata <- c(outdata, mi, est, ses, bias, power)
    
  }
  
  outdata
  
}

## Produce results 
# Define model parameters
N = 200 # sample size
p1 <- 0.5 # proportion of patients allocated treatment 1 
T1ER <- 0.05 # alpha/sign level
T2ER <- 0.2 # beta/1-power
beta <- 0.3 # relation of longitudinal model to survival model
gamma <- -0.1 # fixed treatment effect in longitudinal model
alpha <- -0.3 # treatment effect in survival model
lambda <- 0.85 # baseline hazard rate in survival model
a <- 1.25 # accrual period in years 
f <- 0.75 # follow-up period in years
sigmasq_e = 0.16 # variance of measurement error
mu_theta <- c(0,2,-1,-0.5,0.5)
covar_theta <- matrix(c(1.2, 0.2, 0.4, 0.1, 0.4, 
                         0.2, 0.7, 0.1, 0.2, 0.2,
                         0.4, 0.1, 0.9, 0.3, 0.3,
                         0.1, 0.2, 0.3, 0.2, 0.1,
                         0.4, 0.2, 0.3, 0.1, 0.4), nrow = 5, ncol = 5)

set.seed(1234)
sim4_dat <- replicate(400, sim4_run(N = 200, p1, beta, gamma, alpha, lambda, a, f,
                                     mu_theta, covar_theta, sigmasq_e))

#### Output results ####
write.csv(sim4_dat, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4_results.csv", row.names=FALSE, quote=FALSE) 
########################

library(readr)
sim4_dat <- read_csv("Simulations/Output Tables/sim4_results.csv")

## Summary table
sim4_dat_summary <- rowMeans(sim4_dat, na.rm = T)
sim4_dat_summary2 <- rbind(sim4_dat_summary[1:5],
                           sim4_dat_summary[6:10],
                           sim4_dat_summary[11:15],
                           sim4_dat_summary[16:20])
colnames(sim4_dat_summary2) <- c("Follow-up visit spacing (months)",
                                 "Est", "SE", "Bias", "power")

sim4_dat_summary2 <- as.data.frame(sim4_dat_summary2) %>%
  dplyr::mutate(power.mce = sqrt(power*(1-power)/400),
                power.cil = ifelse(power - 1.96*power.mce < 0, 0, power - 1.96*power.mce),
                power.ciu = ifelse(power + 1.96*power.mce > 1, 1, power + 1.96*power.mce))

colnames(sim4_dat_summary2) <- c("Follow-up visit spacing (months)",
                                 "Est. Trt Effect", "SE", "Bias", "Power",
                                 "Power MCE", "Power CI Lower", "Power CI Upper")

sim4_dat_summary2$`Est. Trt Effect` <- as.numeric(format(round(sim4_dat_summary2$`Est. Trt Effect`, 3), nsmall = 3))
sim4_dat_summary2$SE <- as.numeric(format(round(sim4_dat_summary2$SE, 3), nsmall = 3))
sim4_dat_summary2$Bias <- as.numeric(format(round(sim4_dat_summary2$Bias, 3), nsmall = 3))
sim4_dat_summary2$Power <- as.numeric(format(round(sim4_dat_summary2$Power, 3), nsmall = 3))
sim4_dat_summary2$`Power MCE` <- as.numeric(format(round(sim4_dat_summary2$`Power MCE`, 3), nsmall = 3))
sim4_dat_summary2$`Power CI Lower` <- as.numeric(format(round(sim4_dat_summary2$`Power CI Lower`, 3), nsmall = 3))
sim4_dat_summary2$`Power CI Upper` <- as.numeric(format(round(sim4_dat_summary2$`Power CI Upper`, 3), nsmall = 3))

write.csv(sim4_dat_summary2, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4.csv", row.names=FALSE, quote=FALSE)

## Plot the estimated treatment effect under each follow-up schedule
plot_dat <- as.data.frame(t(sim4_dat)) # transpose, rows are iterations, columns are statistics
colnames(plot_dat)[c(2,7,12,17)] <- c("Est. Trt Effect (monthly)",
                                      "Est. Trt Effect (3-monthly)",
                                      "Est. Trt Effect (6-monthly)",
                                      "Est. Trt Effect (12-monthly)")

# monthly vs. 3-monthly
library(ggplot2)
p1 <- ggplot(data = plot_dat, aes(x = `Est. Trt Effect (monthly)`,
                                  y = `Est. Trt Effect (3-monthly)`)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Comparison of monthly with 3-monthly follow-up visits")
p1
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4_p1.pdf")

# monthly vs. 6-monthly
library(ggplot2)
p2 <- ggplot(data = plot_dat, aes(x = `Est. Trt Effect (monthly)`,
                                  y = `Est. Trt Effect (6-monthly)`)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Comparison of monthly with 6-monthly follow-up visits")
p2
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4_p2.pdf")

# monthly vs. 12-monthly
library(ggplot2)
p3 <- ggplot(data = plot_dat, aes(x = `Est. Trt Effect (monthly)`,
                                  y = `Est. Trt Effect (12-monthly)`)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  labs(title = "Comparison of monthly with 12-monthly follow-up visits")
p3
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4_p3.pdf")

## Plot power by sample size
outdata <- rep(NA,20)

# Simulate using varying N
set.seed(1234)
for (i in seq(from = 100, to = 400, length.out = 7)) {

  x <- replicate(400, sim4_run(N = as.numeric(i), p1, beta, gamma, alpha, lambda, a, f,
                                      mu_theta, covar_theta, sigmasq_e))

  x_summary <- rowMeans(x)

  outdata <- cbind(outdata, x_summary)

}

outdata <- outdata[,-1]
outdata <- as.data.frame(outdata)
colnames(outdata) <- c("N=100", "N=150", "N=200", "N=250", "N=300", "N=350", "N=400")

#### Output results ####
write.csv(outdata, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4_results2.csv", row.names=FALSE, quote=FALSE) 
########################

library(readr)
outdata <- read_csv("Simulations/Output Tables/FINAL/sim4_results2.csv")

plotdata <- t(as.data.frame(outdata))
plotdata2 <- rbind(plotdata[,1:5], plotdata[,6:10], plotdata[,11:15], plotdata[,16:20])
plotdata2 <- as.data.frame(cbind(plotdata2,
                   rep(seq(from = 100, to = 400, length.out = 7), 4))) %>%
  dplyr::select(V1, V5, V6)
colnames(plotdata2) <- c("Follow-up visit spacing (months)",
                         "Power", "N")

# Plot
library(ggplot2)
p4 <- ggplot(data = plotdata2, aes(x = N, y = Power, linetype = factor(`Follow-up visit spacing (months)`))) +
  geom_point(aes(x = N, y = Power)) +
  geom_line(aes(x = N, y = Power, linetype = factor(`Follow-up visit spacing (months)`)), se = FALSE, col = "black") +
  geom_hline(yintercept = 0.8, linetype = "solid", col = "black") +
  theme_classic() +
  labs(title = "Power by sample size for varying follow-up frequency",
       x = "Sample size (N)",
       y = "Power (1-$\\beta$)",
       linetype = "Follow-up frequency") +
  xlim(100, 400)
p4
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim4_p4.pdf")

# Greater power in follow-up schedules with more visits and evenly spaced visits.
# For trial design this means by including an extra X visits we can reduce the
# sample size by X. Study team should discuss what would be a greater burden
# on patients and study conduct: extra X patients or extra X follow-up visits.
