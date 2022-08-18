################################################
# Title: Power by simulation                   #
# Author: Cameron Wilson (219006159)           #
################################################

#### Efficient way of installing previously not installed packages ####
list.of.packages <- c("tidyverse","survival","survminer","lme4",
                      "ggplot2","JM","joineR","simsurv","flexsurv",
                      "MASS","nlme","msm","tictoc","foreach","doParallel",
                      "parallel", "rlecuyer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))suppressPackageStartupMessages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

#### Set Working Directory ####
setwd("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations")

#### Set the seed for reproducibility ####
set.seed(1234)

#### Replicate Chen et al Section 5 ####
## Using Schoenfelds extended formula
# Define model parameters 
T1ER <- 0.05 # alpha/sign level
T2ER <- 0.2 # beta/1-power
p1 <- 0.5 # proportion of patients allocated treatment 1 
beta <- 0.3 # relation of longitudinal model to survival model
gamma <- -0.1 # fixed treatment effect in longitudinal model
alpha <- -0.3 # treatment effect in survival model
lambda <- 0.85 # baseline hazard rate in survival model
a <- 1.25 # accrual period in years 
f <- 0.75 # follow-up period in years
N <- 200 # sample size
mu_theta0 <- 0 # expected value of random intercept
mu_theta1 <- 3 # expected value of random slope
var_theta0 <- 1.2 # variance of random intercept
var_theta1 <- 0.7 # variance of random slope
cov_theta <- 0.2 # covariance of random intercept and slope
sigmasq_e <- 0.16 # variance of measurement error

# Define survival function
surv_fun <- function(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt, t) {
  exp(-lambda*exp(beta*mu_theta0 + (alpha+beta*gamma)*trt)*
        (exp(beta*mu_theta1*t)-1)/(beta*mu_theta1))
}

# Calculate the number of events required
req_events <- function(T1ER, T2ER, beta, gamma, alpha, p1) {
  ((qnorm(p = 1-(T1ER), mean = 0, sd = 1, lower.tail = TRUE) +
      qnorm(p = 1-T2ER, mean = 0, sd = 1, lower.tail = TRUE))^2)/
    (p1*(1-p1)*((beta*gamma + alpha)^2))
}

# Calculate the proportion of patients that will die - Schoenfeld method
prop_die <- function(a, f, beta, gamma, alpha, p1, survfun, lambda, mu_theta0, mu_theta1) {
  Delta <- exp(beta*gamma + alpha)
  d_B <- 1 - (1/6)*(survfun(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt = 1, t = f) + 
                      4*survfun(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt = 1, t = f + 0.5*a) + 
                      survfun(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt = 1, t = f + a))
  d_A <- 1 - (1/6)*(survfun(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt = 0, t = f) + 
                      4*survfun(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt = 0, t = f + 0.5*a) + 
                      survfun(lambda, beta, gamma, alpha, mu_theta0, mu_theta1, trt = 0, t = f + a))
  
  d <- p1*d_A + (1 - p1)*d_B
  d
}

## Sample size from power 
# Calculate the required sample size 
sample_size <- function(T1ER, T2ER, beta, gamma, alpha, p1, a, f, survfun, lambda, mu_theta0, mu_theta1) {
  req_events(T1ER, T2ER, beta, gamma, alpha, p1) / 
    prop_die(a, f, beta, gamma, alpha, p1, survfun, lambda, mu_theta0, mu_theta1)
}

N_est <- sample_size(T1ER, T2ER, beta, gamma, alpha, p1, a, f, survfun = surv_fun, lambda, mu_theta0, mu_theta1)

## Power from sample size 
power <- function(T1ER, beta, gamma, alpha, p1, a, f, survfun, lambda, N, mu_theta0, mu_theta1) {
  x <- sqrt(N*prop_die(a, f, beta, gamma, alpha, p1, survfun, lambda, mu_theta0, mu_theta1)*
              p1*(1 - p1)*((beta*gamma + alpha)^2)) - 
    qnorm(p = 1 - (T1ER), mean = 0, sd = 1, lower.tail = TRUE)
  beta_est <- pnorm(x, mean = 0, sd = 1, lower.tail = TRUE)
  beta_est
}

power(T1ER, beta, gamma, alpha, p1, a, f, survfun = surv_fun, lambda, N = N, mu_theta0, mu_theta1)

#### Power by simulation ####

## Replicate Chen et al. Section 5 now using simulation 
# Assumes an exponential distribution with time-dependent effect. Need to use root finding
# Load libraries
library(simsurv); library(flexsurv); library(survival); library(survminer);
library(dplyr); library(MASS); library(nlme); library(JM); 
library(msm); library(tictoc); library(MASS)

# Define exponential hazard function
haz <- function(t, x, betas, ...) {
  betas[["shape"]] * exp(
    betas[["betaEvent_intercept"]] +
      betas[["betaEvent_binary"]] * x[["Z1"]] +
      betas[["betaEvent_assoc"]] * (
        betas[["betaLong_intercept"]] +
          betas[["betaLong_slope"]] * t +
          betas[["betaLong_binary"]] * x[["Z1"]]
      )
  )
}

# Function to carry out the simulations
sim_run <- function(N, p1, beta, gamma, alpha, lambda, a, f,
                    mu_theta0, mu_theta1, var_theta0, var_theta1, cov_theta, sigmasq_e) {
  
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
    betaLong_slope       = rep(0, N)
  )
  
  # Individual-specific (random effect) parameters
  b_corrmat <- matrix(c(1, cov_theta/(sqrt(var_theta0*var_theta1)), cov_theta/(sqrt(var_theta0*var_theta1)), 1), 2, 2)
  b_sds     <- c(sqrt(var_theta0), sqrt(var_theta1))
  b_means   <- c(mu_theta0, mu_theta1)
  b_z       <- MASS::mvrnorm(n = N, mu = b_means, Sigma = b_corrmat)
  b         <- sapply(1:length(b_sds), function(x) b_sds[x] * b_z[,x])
  betas$betaLong_intercept <- betas$betaLong_intercept + b[,1]
  betas$betaLong_slope     <- betas$betaLong_slope     + b[,2]
  
  # Covariate data
  covdat <- data.frame(
    Z1 = as.integer(c(rep(0,N*p1), rep(1,N*(1-p1)))) # a binary covariate
  )
  
  # Then we simulate the survival times based on the
  # hazard function, covariates, and true parameter values
  times <- simsurv(hazard = haz, x = covdat, betas = betas, maxt = 2,
                   nodes = 7, interval = c(0,2.001), rootsolver = "uniroot")
  times$eventtime <- ifelse(times$eventtime == 0, 0.0001, times$eventtime)
  
  # Simulate observed biomarker data
  sim_dat <- cbind(times, covdat, betas) %>% 
    dplyr::select(id, Z1, eventtime, status, betaLong_intercept, betaLong_slope) %>%
    dplyr::slice(rep(1:n(), 5)) %>%
    dplyr::arrange(id) %>%
    dplyr::rename(trt = Z1) %>%
    dplyr::mutate(obstime = as.numeric(rep(seq(0,2,length.out = 5), length(unique(id)))),
                  theta0 = betaLong_intercept,
                  theta1 = betaLong_slope,
                  biomarker = theta0 + theta1*obstime + gamma*trt + rnorm(n = N*5, mean = 0, sd = sqrt(sigmasq_e))) %>%
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
  
  # Fit a mixed-effects longitudinal model
  library(nlme)
  lmeFit <- nlme::lme(biomarker ~ obstime + trt, data = sim_dat, random = ~ obstime | id,
                      na.action = na.omit,
                      control = lmeControl(maxIter = 50, msMaxIter = 50, niterEM = 25,
                                           opt = "optim"))
  
  # Cox proportional hazards survival model
  library(survival)
  dat.id <- sim_dat %>% dplyr::distinct(id, .keep_all = TRUE)
  survFit <- survival::coxph(Surv(eventtime, status) ~ trt, data = dat.id, 
                             x = TRUE, na.action = na.omit)
  
  
  # Joint model - Maximum likelihood framework
  library(JM)
  jointFit <- jointModel(lmeFit, survFit, timeVar = "obstime",
                         method = "Cox-PH-GH",
                         control = list(iter.EM = 5, 
                                        tol1 = 1e-08, tol2 = 1e-08, 
                                        GHk = 5, GKk = 7))
  
  # Obtain estimates, standard errors and Wald test p-value for overall treatment effect
  alpha_est <- summary(jointFit)[["CoefTable-Event"]]["trt","Value"] 
  beta_est <- summary(jointFit)[["CoefTable-Event"]]["Assoct","Value"]
  gamma_est <- summary(jointFit)[["CoefTable-Long"]]["trt","Value"]
  
  est <- alpha_est + beta_est*gamma_est # overall treatment effect estimate
  
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
  pval <- pnorm(wald) # one-sided p-value
  
  # Need this for power using formula 
  prop_die <- prop.table(table(dat.id$status))[2]
  
  # Return bias & power for the overall treatment effect
  c(bias = est - (alpha + beta*gamma), 
    coverage = ((cil < (alpha + beta*gamma)) & (ciu > (alpha + beta*gamma))),
    power = (as.numeric(pval) < 0.05),
    prop_die = prop_die,
    effect = est,
    SE = ses)
  
}

## Recreate Table 4 last column
beta <- c(rep(0.3, 4), 0.1, 0.4, 0.8, rep(0.3, 4), rep(0.4, 3))
gamma <- c(-0.1, -0.4, -0.8, -1.2, rep(-0.4,10))
var_theta0 <- rep(1.2, length(beta))
var_theta1 <- c(rep(0.7,7), 1, 1.5, 2, 4, rep(0.7, 3))
cov_theta <- c(rep(0.2,11), -0.8, -0.4, 0.4)
power1 <- power(T1ER = 0.05, beta = beta, gamma = gamma, alpha = -0.3, p1 = 0.5, 
                a = 1.25, f = 0.75, survfun = surv_fun, lambda = 0.85, N = 200, mu_theta0, mu_theta1) # calculate power using formula
bias <- rep(0, length(beta))
coverage <- rep(0, length(beta))
prop_die_sim <- rep(0, length(beta))
effect <- rep(0, length(beta))
SE <- rep(0, length(beta))
power2 <- rep(0, length(beta))
table4 <- data.frame(beta, gamma, var_theta0, var_theta1, cov_theta, power1, bias, coverage, prop_die_sim, effect, SE, power2)
table4 <- table4[c(1:14),]

# Require Monte-Carlo SE of power to be less than 0.025. Given that Monte Carlo SE for power = sqrt(power*(1-power)/n)
# and assuming power = 0.5 we need 
n <- (0.5*0.5)/(0.025^2)

# calculate power by simulation
library(parallel); library(foreach); library(doParallel);
library(rlecuyer); library(tictoc)
set.seed(1234, kind = "L'Ecuyer-CMRG")
numcores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  numcores,
  type = "FORK"
)
parallel::clusterSetRNGStream(my.cluster)
print(my.cluster)
req.packages <- list.of.packages
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

Nbeta <- length(beta)
Nreps <- 400

x <- foreach(i = 1:Nbeta,.packages=req.packages,.combine=rbind) %dopar% {
  rowMeans(replicate(Nreps, sim_run(N = 200, p1 = 0.5, beta = table4$beta[i],
                                    gamma = table4$gamma[i], alpha = -0.3, lambda = 0.85, a = 1.25, f = 0.75,
                                    mu_theta0 = 0, mu_theta1 = 3, var_theta0 = table4$var_theta0[i], var_theta1 = table4$var_theta1[i],
                                    cov_theta = table4$cov_theta[i], sigmasq_e = 0.16)), na.rm = T)
}
parallel::stopCluster(my.cluster)

#### Output results ####
write.csv(x, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim0_results.csv", row.names=FALSE, quote=FALSE) 
########################

library(readr)
x <- read_csv("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim0_results.csv")

for(i in 1:Nbeta) {
  table4$bias[i] <-         sprintf('%.3f', round(x[i,1], 3))
  table4$coverage[i] <-     sprintf('%.3f', round(x[i,2], 3))
  table4$power2[i] <-       sprintf('%.3f', round(x[i,3], 3))
  table4$prop_die_sim[i] <- sprintf('%.3f', x[i,4])
  table4$effect[i] <-       sprintf('%.3f', round(x[i,5], 3))
  table4$SE[i] <-           sprintf('%.3f', round(x[i,6], 3))
}

table4$power2 <- as.numeric(table4$power2)
table4$prop_die_sim <- as.numeric(table4$prop_die_sim)
table4 <- table4 %>% dplyr::mutate(simpower = paste0(power2, " (",
                                                     round(power2 - 1.96*sqrt((power2*(1-power2))/Nreps), 3),
                                                     ", ",
                                                     round(power2 + 1.96*sqrt((power2*(1-power2))/Nreps), 3), ")"),
                                   x = sqrt(N*prop_die_sim*p1*(1-p1)*((alpha + beta*gamma)^2)) -
                                     qnorm(p = 1 - T1ER, mean = 0, sd = 1, lower.tail = TRUE),
                                   power3 = pnorm(x, mean = 0, sd = 1, lower.tail = TRUE)) %>%
  dplyr::select(-x, -effect, -SE, -power2)

table4$power3 <- sprintf('%.3f', round(table4$power3, 3))

colnames(table4)[6:11] <- c("Formula: Power", "Sim: Bias", "Sim: Coverage",
                            "Proportion that die", "Sim: Power", "Formula: Power*")

# * indicates power was calculated by formula using the mean number of events from simulations
# rather than Schoenfeld's approximation using the survival function (Simpson's rule)

write.csv(table4, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim0.csv", row.names=FALSE, quote=FALSE)

# The power by simulation has replicated the power calculations using the
# formula from Chen et al. Now we can relax the assumptions made through the
# simulation and calculate power in more flexible scenarios. We can then assess
# whether the formula is robust to these assumptions or whether it is too rigid.