#### Assumption 5: no missing longitudinal data ####
# Simulate the power when longitudinal data is MAR and MNAR
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
setwd("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R")

cols <- c("#7D1D3F", "#3A506B", "#D5CAD6")

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
sim2_run <- function(N, p1, beta, gamma, alpha, lambda, a, f,
                     mu_theta, covar_theta, sigmasq_e,
                     outcome, mechanism, percentage) {
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
    dplyr::slice(rep(1:n(), 5)) %>%
    dplyr::arrange(id) %>%
    dplyr::rename(trt = Z1) %>%
    dplyr::mutate(obstime = as.numeric(rep(seq(0,2,length.out = 5), length(unique(id)))),
                  p = rep(p,N*5))
  
  sim_dat <- sim_dat %>% dplyr::mutate(theta0 = betaLong_intercept,
                                       theta1 = betaLong_slope,
                                       theta2 = betaLong_quadratic,
                                       theta3 = betaLong_cubic,
                                       theta4 = betaLong_quartic,
                                       biomarker = theta0 + theta1*obstime + theta2*(obstime^2) + theta3*(obstime^3) + theta4*(obstime^4) + 
                                         gamma*trt + rnorm(n = N*5, mean = 0, sd = sqrt(sigmasq_e))) %>%
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
    
  # Introducing missing data mechanism
  dat.id <- sim_dat %>% dplyr::distinct(id, .keep_all = TRUE)
  
  if (outcome == "longitudinal") {
    
    if (mechanism == "MCAR") {
      
      n <- length(sim_dat$biomarker)
      u <- runif(n, min = 0, max = 1)
      sim_dat$biomarker <- ifelse((u <= percentage/100) & (sim_dat$obstime != 0), NA, sim_dat$biomarker)
      
    } else if (mechanism == "MAR") {
      
      n <- length(sim_dat$biomarker[sim_dat$trt == 1])
      u <- runif(n, min = 0, max = 1)
      sim_dat$biomarker[sim_dat$trt == 1] <- ifelse((u <= percentage/100) & (sim_dat$obstime[sim_dat$trt == 1] != 0), NA, sim_dat$biomarker[sim_dat$trt == 1])
      
    }
    
  } else if (outcome == "survival") {
    
    if (mechanism == "MCAR") {
      
      n <- length(dat.id$eventtime)
      u <- runif(n, min = 0, max = 1)
      dat.id$eventtime <- ifelse(u <= percentage/100, NA, dat.id$eventtime)
      
    } else if (mechanism == "MAR") {
      
      n <- length(dat.id$eventtime[dat.id$trt == 1])
      u <- runif(n, min = 0, max = 1)
      dat.id$eventtime[dat.id$trt == 1] <- ifelse(u <= percentage/100, NA, dat.id$eventtime[dat.id$trt == 1])
      
    } else if (mechanism == "MNAR") {
      
      n <- length(dat.id$eventtime[dat.id$status == 1])
      u <- runif(n, min = 0, max = 1)
      dat.id$eventtime[dat.id$status == 1] <- ifelse(u <= percentage/100, NA, dat.id$eventtime[dat.id$status == 1])
      
    }
    
    dat.id$status <- ifelse(is.na(dat.id$eventtime), NA, dat.id$status)
    dat.id$biomarker <- ifelse(is.na(dat.id$eventtime), NA, dat.id$biomarker)
    dat.id$obstime <- ifelse(is.na(dat.id$eventtime), NA, dat.id$obstime)
    
    # If survival data missing, longitudinal data must logically also be missing
    ids <- dat.id$id[is.na(dat.id$eventtime)]
    sim_dat$biomarker[sim_dat$id %in% ids] <- NA
    sim_dat$obstime[sim_dat$id %in% ids] <- NA
    sim_dat$eventtime[sim_dat$id %in% ids] <- NA
    sim_dat$status[sim_dat$id %in% ids] <- NA
    
  } else if (outcome == "none") {
    
    sim_dat <- sim_dat
    
  }
  
  # Fit a mixed-effects longitudinal model
  library(nlme)
  lmeFit <- nlme::lme(biomarker ~ obstime + I(obstime^2) + I(obstime^3) + trt, data = sim_dat, random = ~ obstime | id,
                      na.action = na.omit,
                      control = lmeControl(maxIter = 50, msMaxIter = 50, niterEM = 25, 
                                           tolerance = 1e-08, opt = "optim"))
    
  # Cox model
  library(survival)
  survFit <- survival::coxph(Surv(eventtime, status) ~ trt, data = dat.id, 
                             x = TRUE, na.action = na.omit)
    
  # Joint model - Maximum likelihood framework
  library(JM)
  jointFit <- jointModel(lmeFit, survFit, timeVar = "obstime", 
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
  coverage <- ((cil < (alpha + beta*gamma)) & (ciu > (alpha + beta*gamma)))
    
  c(est, ses, cil, ciu, bias, power, coverage)
  
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

sim2_dat <- data.frame(ref = seq(from = 1, to = 21, by = 1),
                       outcome = c("none", rep("longitudinal", 8), rep("survival", 12)),
                       mechanism = c(NA, rep("MCAR", 4), rep("MAR", 4), rep("MCAR", 4), rep("MAR", 4), rep("MNAR", 4)),
                       percentage = c(0, rep(c(5, 10, 20, 50), 5)),
                       est = rep(NA, 21),
                       ses = rep(NA, 21),
                       cil = rep(NA, 21),
                       ciu = rep(NA, 21),
                       bias = rep(NA, 21), 
                       bias.se = rep(NA, 21), 
                       power = rep(NA, 21),
                       coverage = rep(NA, 21))

# calculate power by simulation
library(parallel); library(foreach); library(doParallel);
library(rlecuyer); library(tictoc)
set.seed(1234, kind = "L'Ecuyer-CMRG")
numcores = parallel::detectCores()
my.cluster <- parallel::makeCluster(
  numcores,
  type = "FORK"
)
print(my.cluster)
req.packages <- list.of.packages
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

Nbeta = length(sim2_dat$ref)
Nreps = 400

x <- foreach(i = 1:Nbeta,.packages=req.packages) %dopar% {
  set.seed(1234)
  replicate(Nreps, sim2_run(N = 200, p1, beta, gamma, alpha, lambda, a, f,
                            mu_theta, covar_theta, sigmasq_e,
                            outcome = sim2_dat$outcome[i],
                            mechanism = sim2_dat$mechanism[i],
                            percentage = sim2_dat$percentage[i]))
}
parallel::stopCluster(my.cluster)

#### Output results ####
write.csv(x, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_results.csv", row.names=FALSE, quote=FALSE) 
########################

library(readr)
x <- read_csv("Simulations/Output Tables/sim2_results.csv")

sim2_results <- as.data.frame(t(x)) %>%
  dplyr::mutate(ref = sort(rep(seq(from = 1, to = 21, by = 1), Nreps)),
                outcome = c(rep("None", Nreps),
                         rep("Longitudinal", Nreps*8),
                         rep("Survival", Nreps*12)),
                mechanism = c(rep(NA, Nreps), rep("MCAR", 4*Nreps), rep("MAR", 4*Nreps), rep("MCAR", 4*Nreps), rep("MAR", 4*Nreps), rep("MNAR", 4*Nreps)),
                percentage = c(rep(0, Nreps), rep(c(rep(5, Nreps), rep(10, Nreps), rep(20, Nreps), rep(50, Nreps)), 5))) %>%
  group_by(ref, outcome, mechanism, percentage) %>%
  summarise(Est = as.numeric(format(round(mean(V1, na.rm = TRUE), 3), nsmall = 3)),
            SE = as.numeric(format(round(mean(V2, na.rm = TRUE), 3), nsmall = 3)),
            CIL = as.numeric(format(round(mean(V3, na.rm = TRUE), 3), nsmall = 3)),
            CIU = as.numeric(format(round(mean(V4, na.rm = TRUE), 3), nsmall = 3)),
            Bias = as.numeric(format(round(mean(V5, na.rm = TRUE), 3), nsmall = 3)),
            Bias.mce = as.numeric(format(round(sqrt(sum((V1 -  (alpha + beta*gamma))^2)/(Nreps*(Nreps-1))), 5), nsmall = 5)),
            Bias.cil = as.numeric(format(round(Bias - 1.96*Bias.mce, 5), nsmall = 5)),
            Bias.ciu = as.numeric(format(round(Bias + 1.96*Bias.mce, 5), nsmall = 5)),
            Power = as.numeric(format(round(mean(V6, na.rm = TRUE), 3), nsmall = 3)),
            Coverage = as.numeric(format(round(mean(V7, na.rm = TRUE), 3), nsmall = 3))) %>%
  ungroup() %>%
  dplyr::mutate(Power.mce = as.numeric(format(round(sqrt(Power*(1-Power)/Nreps), 3), nsmall = 3)),
                Coverage.mce = as.numeric(format(round(sqrt(Coverage*(1-Coverage)/Nreps), 3), nsmall = 3)),
                Power.cil = as.numeric(format(round(Power - 1.96*Power.mce, 3), nsmall = 3)),
                Power.ciu = as.numeric(format(round(Power + 1.96*Power.mce, 3), nsmall = 3)),
                Coverage.cil = as.numeric(format(round(Coverage - 1.96*Coverage.mce, 3), nsmall = 3)),
                Coverage.ciu = as.numeric(format(round(Coverage + 1.96*Coverage.mce, 3), nsmall = 3))) %>%
  dplyr::select(ref, outcome, mechanism, percentage, Est, SE, CIL, CIU,
                Bias, Bias.mce, Bias.cil, Bias.ciu,
                Power, Power.mce, Power.cil, Power.ciu,
                Coverage, Coverage.mce, Coverage.cil, Coverage.ciu)

# Bound confidence intervals for power and coverage between 0-1
sim2_results$Power.cil <- ifelse(sim2_results$Power.cil < 0, 0, sim2_results$Power.cil)
sim2_results$Power.ciu <- ifelse(sim2_results$Power.ciu > 1, 1, sim2_results$Power.ciu)
sim2_results$Coverage.cil <- ifelse(sim2_results$Coverage.cil < 0, 0, sim2_results$Coverage.cil)
sim2_results$Coverage.ciu <- ifelse(sim2_results$Coverage.ciu > 1, 1, sim2_results$Coverage.ciu)

write.csv(sim2_results, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2.csv", row.names=FALSE, quote=FALSE)

# Plot the treatment effect estimates to show bias
library(ggplot2)
sim2_plot <- ggplot(data = sim2_results, aes(x = Est, y = 1, xmin = CIL, xmax = CIU, group = factor(ref))) +
  geom_point(size = 1) +
  geom_errorbar(aes(xmin = CIL, xmax = CIU), width = 0.2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  theme_classic() +
  labs(x = "Estimated Overall Treatment Effect") +
  facet_grid(rows = vars(factor(outcome), factor(mechanism), factor(percentage)), scales = "fixed") +
  theme(legend.position="right",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
sim2_plot
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_plot.pdf")
# Bias is the distance from the true overall treatment effect (alpha + beta*gamma)

for (i in c(1:length(sim2_results$ref))) {
  index <- ((Nreps*(i-1)) + 1):(Nreps*i)
  plot_data <- data.frame(Est = c(t(x[1, index])), y = c(1:Nreps))

  p <- ggplot(data = plot_data, aes(x = Est, y = y)) +
    geom_point(col = "black", alpha = 0.3) +
    geom_point(data = sim2_results[i,], aes(x = Est, y = 200, group = factor(ref)), size = 3) +
    geom_errorbar(data = sim2_results[i,], aes(x = Est, y = 200, xmin = CIL, xmax = CIU), width = 0.2, col = "black") +
    geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
    theme_classic() +
    labs(x = "Estimated Overall Treatment Effect",
         y = "Repetition number") +
    theme(legend.position="none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text.y = element_blank())
  p
  ggsave(paste0("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_plot_", i, ".png"))
}
# Bias is the distance from the true overall treatment effect (alpha + beta*gamma)

# Table
sim2_table <- data.frame(ref = sort(rep(1:21,3)),
                         outcome = c(rep("none", 3), rep("Longitudinal", 24), rep("Survival", 36)),
                         mechanism = c(rep(NA, 3), rep("MCAR", 12), rep("MAR", 12), rep("MCAR", 12), rep("MAR", 12), rep("MNAR", 12)),
                         percentage = c(rep(0, 3), rep(sort(rep(c(5,10,20,50), 3)), 5)),
                         `Missing data\n mechanism` = rep(c("Est", "Monte Carlo SE", "95% CI"), 21),
                         Est = rep(NA, 63),
                         Bias = rep(NA, 63),
                         Power = rep(NA, 63),
                         Coverage = rep(NA, 63))

for (i in seq(from = 1, to = 21, by = 1)) {
  sim2_table$Est[3*i-2] <- sim2_results$Est[i]
  sim2_table$Est[3*i-1] <- sim2_results$SE[i]
  sim2_table$Est[3*i] <- paste0("(", sim2_results$CIL[i], ", ", sim2_results$CIU[i], ")")
  sim2_table$Bias[3*i-2] <- sim2_results$Bias[i]
  sim2_table$Bias[3*i-1] <-  NA
  sim2_table$Bias[3*i] <- NA
  sim2_table$Power[3*i-2] <- sim2_results$Power[i]
  sim2_table$Power[3*i-1] <- sim2_results$Power.mce[i]
  sim2_table$Power[3*i] <- paste0("(", sim2_results$Power.cil[i], ", ", sim2_results$Power.ciu[i], ")")
  sim2_table$Coverage[3*i-2] <- sim2_results$Coverage[i]
  sim2_table$Coverage[3*i-1] <- sim2_results$Coverage.mce[i]
  sim2_table$Coverage[3*i] <- paste0("(", sim2_results$Coverage.cil[i], ", ", sim2_results$Coverage.ciu[i], ")")
}

library(knitr); library(kableExtra)
sim2_table <- sim2_table %>% dplyr::mutate(images = paste0("<","img src=","sim2_plot_", ref,".png",">")) %>%
  dplyr::select(-ref, -outcome, -mechanism, -percentage) %>%
  kbl(linesep = "", row.names=T,
      caption="Simulation 1: Summary table (Missing data)", align = "c", vline = "",
      col.names = c("Data generating\n mechanism", "Est. Trt\n Effect",
                    "Bias", "Power", "Coverage", "Plot"), escape = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options="scale_down") %>%
  pack_rows("No missing", 1, 3, color = cols[3], background = cols[1]) %>%
  pack_rows("Longitudinal", 4, 27, color = cols[3], background = cols[1]) %>%
  pack_rows("Survival", 28, 63, color = cols[3], background = cols[1]) %>%
  footnote(general = "1. Vertical axis of plots is the repetition number.
           2. Dashed line represents true treatment effect.
           3. Interval is 95% confidence interval.")

# Bivariate plot of estimated treatment effect against standard error
for (i in c(1:length(sim2_results$ref))) {

  index <- ((Nreps*(i-1)) + 1):(Nreps*i)
  plot_data <- data.frame(Est = c(t(x[1, index])),
                          SE = c(t(x[2, index])))

  p <- ggplot(data = plot_data, aes(x = Est, y = SE)) +
    geom_point(col = "black") +
    theme_classic() +
    labs(x = "Estimated Overall Treatment Effect",
         y = "Standard Error",
         title = paste("Bivariate plot for\n", sim2_results$outcome[i],
                       sim2_results$mechanism[i], paste0(sim2_results$percentage[i], "%"), "missing model")) +
    theme(legend.position="none")
  p
  ggsave(paste0("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_plot_bivariate_", i, ".png"))
}

# Bivariate plot of estimated treatment effect between pairs
for (i in seq(from = 5, to = 21, by = 4)) {

  for (j in seq(from = 5, to = 21, by = 4)) {

    index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
    index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
    plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                            SE.i = c(t(x[2, index.i])),
                            Est.j = c(t(x[1, index.j])),
                            SE.j = c(t(x[2, index.j])))

    if (i > j) {

      # Top diagonal is bivariate plot of estimated treatment effect
      x.min <- min(plot_data$Est.j, na.rm = T)
      x.max <- max(plot_data$Est.j, na.rm = T)
      y.min <- min(plot_data$Est.i, na.rm = T)
      y.max <- max(plot_data$Est.i, na.rm = T)
      min <- min(x.min, y.min)
      max <- max(x.max, y.max)
      p <- ggplot(data = plot_data, aes(x = Est.j, y = Est.i)) +
        geom_point() +
        geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
        geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
        geom_abline(intercept = 0, slope = 1) +
        theme_classic() +
        scale_x_continuous(limits = c(min, max)) +
        scale_y_continuous(limits = c(min, max)) +
        labs(x = paste("Est. Trt Effect in", sim2_results$outcome[j], sim2_results$mechanism[j], "model"),
             y = paste("Est. Trt Effect in", sim2_results$outcome[i], sim2_results$mechanism[i], "model"))
      p
      ggsave(paste0("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_plot_bivariate_pairs_", i, j, ".png"))

    } else if (i < j) {

      # Bottom diagonal is bivariate plot of standard error
      x.min <- min(plot_data$SE.j, na.rm = T)
      x.max <- max(plot_data$SE.j, na.rm = T)
      y.min <- min(plot_data$SE.i, na.rm = T)
      y.max <- max(plot_data$SE.i, na.rm = T)
      min <- min(x.min, y.min)
      max <- max(x.max, y.max)
      p <- ggplot(data = plot_data, aes(x = SE.j, y = SE.i)) +
        geom_point() +
        geom_abline(intercept = 0, slope = 1) +
        theme_classic() +
        scale_x_continuous(limits = c(min, max)) +
        scale_y_continuous(limits = c(min, max)) +
        labs(x = paste("Std. Err. in", sim2_results$outcome[j], sim2_results$mechanism[j], "model"),
             y = paste("Std. Err. in", sim2_results$outcome[i], sim2_results$mechanism[i], "model"))
      p
      ggsave(paste0("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_plot_bivariate_pairs_", i, j, ".png"))

    }

  }

}

# Zip plot: plot of confidence intervals fractionally ranked by mod(zi)
plot_data2 <- data.frame()
for (i in c(1:length(sim2_results$ref))) {

  index <- ((Nreps*(i-1)) + 1):(Nreps*i)
  plot_data <- data.frame(t(x[, index])) %>% dplyr::mutate(ref = rep(i, n()))
  colnames(plot_data) <- c("est", "ses", "cil", "ciu", "bias", "power", "coverage", "ref")
  plot_data <- plot_data %>% dplyr::mutate(modzi = abs(bias/ses)) %>%
    arrange(desc(modzi)) %>% dplyr::mutate(percentile = rev(seq(from = 0.25, to = 100, length.out = n())))
  plot_data$coverage.cil <- rep(sim2_results$Coverage.cil[sim2_results$ref == i], Nreps)
  plot_data$coverage.ciu <- rep(sim2_results$Coverage.ciu[sim2_results$ref == i], Nreps)
  plot_data2 <- rbind(plot_data2, plot_data)

}

p <- ggplot(data = plot_data2, aes(x = est, y = percentile, xmin = cil, xmax = ciu,
                                   colour = factor(coverage), group = factor(ref))) +
  geom_point(size = 0.5, col = "black") +
  geom_errorbar(aes(xmin = cil, xmax = ciu, color = factor(coverage)), width = 0.1, alpha = 0.5) +
  geom_vline(xintercept = (alpha + beta*gamma), lty = 1, col = cols[1]) +
  geom_hline(data = plot_data2, aes(yintercept = 100*coverage.cil), lty = 1, col = cols[1]) +
  geom_hline(data = plot_data2, aes(yintercept = 100*coverage.ciu), lty = 1, col = cols[1]) +
  theme_classic() +
  facet_wrap(facets = vars(factor(ref)),
             ncol = 3, nrow = 7, scales = "fixed", strip.position = "top") +
  labs(x = "95% Confidence Intervals for Estimated Treatment Effect",
       y = "Fractional centile of |z|") +
  scale_color_manual(values = cols[c(2,3)],
                     name = "CI Covers True Value",
                     labels = c("Non-coverer", "Coverer"),
                     guide = guide_legend(reverse=TRUE)) +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = c(5,50,95), limits = c(0,100)) +
  scale_x_continuous(breaks = c(-1.25, -1, -0.75, -0.5, -0.25, 0, 0.25))
p
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim2_zipplot.png")

# When the method has 95% coverage (also seen in the table) the colour of the intervals
# switches at the 95th percentile on the vertical axis.

# Horizontal lines are the Monte Carlo 95% confidence intervals for the percentage coverage.