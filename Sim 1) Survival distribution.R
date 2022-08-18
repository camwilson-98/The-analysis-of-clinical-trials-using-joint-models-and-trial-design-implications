#### Assumption 3: the survival density function follows an exponential distribution ####
# Simulate the power and bias when the survival distribution is something other than 
# the exponential 

# Efficient way of installing previously not installed packages
list.of.packages <- c("tidyverse","survival","survminer","lme4",
                      "ggplot2","JM","joineR","simsurv","flexsurv",
                      "MASS","nlme","msm","tictoc","foreach","doParallel",
                      "parallel", "rlecuyer","expint")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))suppressPackageStartupMessages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

# Set Working Directory
setwd("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R")

cols <- c("#7D1D3F", "#3A506B", "#D5CAD6")

# Generalised gamma hazard function 
haz <- function(t, x, betas, ...) {
  ((betas[["shape2"]] * betas[["scale"]] * (betas[["scale"]]*t)^(betas[["shape"]] - 1) * 
      exp(-(betas[["scale"]]*t)^betas[["shape2"]])) / 
     (expint::gammainc(betas[["shape"]]/betas[["shape2"]], (betas[["scale"]]*t)^betas[["shape2"]]))) *
    exp(betas[["betaEvent_intercept"]] +
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
sim1_run <- function(N, p1, beta, gamma, alpha, lambda, shape, shape2, 
                     a, f, mu_theta, covar_theta, sigmasq_e,
                     dist) {
  # Then we construct data frames with the true parameter
  # values and the covariate data for each individual
  
  # Set generalised gamma parameters according to desired survival distribution
  if (dist == "Gen Gamma") {
    
    shape_param  <- shape  # alpha
    shape2_param <- shape2  # p
    lambda_param <- lambda+0.25  # lambda 
    
  } else if (dist == "Gamma") {
    
    shape_param  <- shape        
    shape2_param <- 1           # Gen Gamma reduces to Gamma when shape2 = 1
    lambda_param <- lambda+1.10
    
  } else if (dist == "Weibull") {
    
    shape_param  <- shape
    shape2_param <- shape       # Gen Gamma reduces to Weibull when shape = shape2
    lambda_param <- lambda
    
  } else if (dist == "Exponential") {
    
    shape_param  <- 1
    shape2_param <- 1           # Gen Gamma reduces to Exponential when shape = shape2 = 1
    lambda_param <- lambda+0.05
    
  }
  
  # Population (fixed effect) parameters
  betas <- data.frame(
    shape                = rep(shape_param, N),
    scale                = rep(lambda_param, N),
    shape2               = rep(shape2_param, N), 
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
    
  # Fit a mixed-effects longitudinal model
  library(nlme)
  lmeFit <- nlme::lme(biomarker ~ obstime + I(obstime^2) + I(obstime^3) + trt, data = sim_dat, random = ~ obstime | id,
                      na.action = na.omit,
                      control = lmeControl(maxIter = 50, msMaxIter = 50, niterEM = 25, 
                                           tolerance = 1e-08, opt = "optim"))
  
  # Cox proportional hazards survival model
  library(survival)
  dat.id <- sim_dat %>% dplyr::distinct(id, .keep_all = TRUE)
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
  
  # Include the proportion of events to adjust survival distribution parameters if needed.
  # Need a similar number of events to ensure comparability between the distributions. 
  prop_die <- prop.table(table(dat.id$status))[2]
  
  bias <- est - (alpha + beta*gamma)
  power <- (as.numeric(pval) < 0.05)
  coverage <- ((cil < (alpha + beta*gamma)) & (ciu > (alpha + beta*gamma)))
  
  c(est, ses, cil, ciu, bias, power, coverage, prop_die)
  
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
shape <- 2
shape2 <- 1.5
a <- 1.25 # accrual period in years 
f <- 0.75 # follow-up period in years
sigmasq_e = 0.16 # variance of measurement error
mu_theta <- c(0,2,-1,-0.5,0.5)
covar_theta <- matrix(c(1.2, 0.2, 0.4, 0.1, 0.4, 
                        0.2, 0.7, 0.1, 0.2, 0.2,
                        0.4, 0.1, 0.9, 0.3, 0.3,
                        0.1, 0.2, 0.3, 0.2, 0.1,
                        0.4, 0.2, 0.3, 0.1, 0.4), nrow = 5, ncol = 5)

sim1_dat <- data.frame(ref = seq(from = 1, to = 4, by = 1),
                       dist = c("Gen Gamma", "Gamma", "Weibull", "Exponential"),
                       est = rep(NA, 4),
                       ses = rep(NA, 4),
                       cil = rep(NA, 4),
                       ciu = rep(NA, 4),
                       bias = rep(NA, 4), 
                       power = rep(NA, 4),
                       coverage = rep(NA, 4),
                       prop_die = rep(NA,4))

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

Nbeta = length(sim1_dat$ref)
Nreps = 400

x <- foreach(i = 1:Nbeta,.packages=req.packages) %dopar% {
  set.seed(1234)
  replicate(Nreps, sim1_run(N = 200, p1, beta, gamma, alpha, lambda, shape, shape2,
                                 a, f, mu_theta, covar_theta, sigmasq_e,
                                 dist = sim1_dat$dist[i]))
}
parallel::stopCluster(my.cluster)

#### Output results ####
write.csv(x, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim1_results.csv", row.names=FALSE, quote=FALSE) 
########################

library(readr)
x <- read_csv("Simulations/Output Tables/sim1_results.csv")

sim1_results <- as.data.frame(t(x)) %>%
  dplyr::mutate(dist = c(rep("Gen Gamma", Nreps),
                         rep("Gamma", Nreps),
                         rep("Weibull", Nreps),
                         rep("Exponential", Nreps))) %>%
  group_by(dist) %>%
  summarise(Est = as.numeric(format(round(mean(V1, na.rm = T), 3), nsmall = 3)),
            SE = as.numeric(format(round(mean(V2, na.rm = T), 3), nsmall = 3)),
            CIL = as.numeric(format(round(mean(V3, na.rm = T), 3), nsmall = 3)),
            CIU = as.numeric(format(round(mean(V4, na.rm = T), 3), nsmall = 3)),
            Bias = as.numeric(format(round(mean(V5, na.rm = T), 3), nsmall = 3)),
            Bias.mce = as.numeric(format(round(sqrt(sum((V1 -  (alpha + beta*gamma))^2)/(Nreps*(Nreps-1))), 5), nsmall = 5)),
            Bias.cil = as.numeric(format(round(Bias - 1.96*Bias.mce, 5), nsmall = 5)),
            Bias.ciu = as.numeric(format(round(Bias + 1.96*Bias.mce, 5), nsmall = 5)),
            Power = as.numeric(format(round(mean(V6, na.rm = T), 3), nsmall = 3)),
            Coverage = as.numeric(format(round(mean(V7, na.rm = T), 3), nsmall = 3)),
            Prop_die = as.numeric(format(round(mean(V8, na.rm = T), 3), nsmall = 3))) %>%
  ungroup() %>%
  dplyr::mutate(Power.mce = as.numeric(format(round(sqrt(Power*(1-Power)/Nreps), 3), nsmall = 3)),
                Coverage.mce = as.numeric(format(round(sqrt(Coverage*(1-Coverage)/Nreps), 3), nsmall = 3)),
                Power.cil = as.numeric(format(round(Power - 1.96*Power.mce, 3), nsmall = 3)),
                Power.ciu = as.numeric(format(round(Power + 1.96*Power.mce, 3), nsmall = 3)),
                Coverage.cil = as.numeric(format(round(Coverage - 1.96*Coverage.mce, 3), nsmall = 3)),
                Coverage.ciu = as.numeric(format(round(Coverage + 1.96*Coverage.mce, 3), nsmall = 3))) %>%
  dplyr::select(dist, Est, SE, CIL, CIU,
                Bias, Bias.mce, Bias.cil, Bias.ciu,
                Power, Power.mce, Power.cil, Power.ciu,
                Coverage, Coverage.mce, Coverage.cil, Coverage.ciu, Prop_die)

# Bound confidence intervals for power and coverage between 0-1
sim1_results$Power.cil <- ifelse(sim1_results$Power.cil < 0, 0, sim1_results$Power.cil)
sim1_results$Power.ciu <- ifelse(sim1_results$Power.ciu > 1, 1, sim1_results$Power.ciu)
sim1_results$Coverage.cil <- ifelse(sim1_results$Coverage.cil < 0, 0, sim1_results$Coverage.cil)
sim1_results$Coverage.ciu <- ifelse(sim1_results$Coverage.ciu > 1, 1, sim1_results$Coverage.ciu)

sim1_results <- sim1_results[c(3,2,4,1),]
write.csv(sim1_results, "~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim1.csv", row.names=FALSE, quote=FALSE)

# Plot the treatment effect estimates to show bias
sim1_results$dist <- factor(sim1_results$dist, 
                            levels = c("Gen Gamma", "Gamma", "Weibull", "Exponential"), 
                            labels = c("Gen Gamma", "Gamma", "Weibull", "Exponential"), ordered = T)
library(ggplot2)
sim1_plot <- ggplot(data = sim1_results, aes(x = Est, y = 1, xmin = CIL, xmax = CIU, group = factor(dist))) +
  geom_point(size = 1, col = "black") +
  geom_errorbar(aes(xmin = CIL, xmax = CIU), width = 0.2, col = "black") +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  theme_classic() +
  labs(x = "Estimated Overall Treatment Effect") +
  facet_grid(rows = vars(dist), scales = "fixed") +
  theme(legend.position="right",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())
sim1_plot
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim1_plot.pdf")

for (i in c(1:length(sim1_results$dist))) {
  index <- ((Nreps*(i-1)) + 1):(Nreps*i)
  plot_data <- data.frame(Est = c(t(x[1, index])), y = c(1:Nreps))

  p <- ggplot(data = plot_data, aes(x = Est, y = y)) +
    geom_point(col = "black", alpha = 0.3) +
    geom_point(data = sim1_results[i,], aes(x = Est, y = 200, group = factor(dist, labels = paste(sim1_results$dist[i]))), size = 3) +
    geom_errorbar(data = sim1_results[i,], aes(x = Est, y = 200, xmin = CIL, xmax = CIU), width = 0.2, col = "black") +
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
  ggsave(paste0("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim1_plot_", i, ".png"))
}
# Bias is the distance from the true overall treatment effect (alpha + beta*gamma)

# Table
sim1_table <- data.frame(ref = sort(rep(1:4,3)),
                         labels = c(rep("Gen Gamma", 3), rep("Gamma", 3), rep("Weibull", 3), rep("Exponential", 3)),
                         dist = rep(c("Est", "Monte Carlo SE", "95% CI"), 4),
                         `Est. Trt\n Effect`= c(c(sim1_results$Est[sim1_results$dist == "Gen Gamma"], sim1_results$SE[sim1_results$dist == "Gen Gamma"],
                                                  paste0("(", sim1_results$CIL[sim1_results$dist == "Gen Gamma"], ", ", sim1_results$CIU[sim1_results$dist == "Gen Gamma"], ")")),
                                                c(sim1_results$Est[sim1_results$dist == "Gamma"], sim1_results$SE[sim1_results$dist == "Gamma"],
                                                  paste0("(", sim1_results$CIL[sim1_results$dist == "Gamma"], ", ", sim1_results$CIU[sim1_results$dist == "Gamma"], ")")),
                                                c(sim1_results$Est[sim1_results$dist == "Weibull"], sim1_results$SE[sim1_results$dist == "Weibull"],
                                                  paste0("(", sim1_results$CIL[sim1_results$dist == "Weibull"], ", ", sim1_results$CIU[sim1_results$dist == "Weibull"], ")")),
                                                c(sim1_results$Est[sim1_results$dist == "Exponential"], sim1_results$SE[sim1_results$dist == "Exponential"],
                                                  paste0("(", sim1_results$CIL[sim1_results$dist == "Exponential"], ", ", sim1_results$CIU[sim1_results$dist == "Exponential"], ")"))),
                         Bias = c(c(sim1_results$Bias[sim1_results$dist == "Gen Gamma"], NA,NA),
                                  c(sim1_results$Bias[sim1_results$dist == "Gamma"], NA,NA),
                                  c(sim1_results$Bias[sim1_results$dist == "Weibull"], NA,NA),
                                  c(sim1_results$Bias[sim1_results$dist == "Exponential"], NA,NA)),
                         Power = c(c(sim1_results$Power[sim1_results$dist == "Gen Gamma"], sim1_results$Power.mce[sim1_results$dist == "Gen Gamma"],
                                     paste0("(", sim1_results$Power.cil[sim1_results$dist == "Gen Gamma"], ", ", sim1_results$Power.ciu[sim1_results$dist == "Gen Gamma"], ")")),
                                   c(sim1_results$Power[sim1_results$dist == "Gamma"], sim1_results$Power.mce[sim1_results$dist == "Gamma"],
                                     paste0("(", sim1_results$Power.cil[sim1_results$dist == "Gamma"], ", ", sim1_results$Power.ciu[sim1_results$dist == "Gamma"], ")")),
                                   c(sim1_results$Power[sim1_results$dist == "Weibull"], sim1_results$Power.mce[sim1_results$dist == "Weibull"],
                                     paste0("(", sim1_results$Power.cil[sim1_results$dist == "Weibull"], ", ", sim1_results$Power.ciu[sim1_results$dist == "Weibull"], ")")),
                                   c(sim1_results$Power[sim1_results$dist == "Exponential"], sim1_results$Power.mce[sim1_results$dist == "Exponential"],
                                     paste0("(", sim1_results$Power.cil[sim1_results$dist == "Exponential"], ", ", sim1_results$Power.ciu[sim1_results$dist == "Exponential"], ")"))),
                         Coverage = c(c(sim1_results$Coverage[sim1_results$dist == "Gen Gamma"], sim1_results$Coverage.mce[sim1_results$dist == "Gen Gamma"],
                                        paste0("(", sim1_results$Coverage.cil[sim1_results$dist == "Gen Gamma"], ", ", sim1_results$Coverage.ciu[sim1_results$dist == "Gen Gamma"], ")")),
                                      c(sim1_results$Coverage[sim1_results$dist == "Gamma"], sim1_results$Coverage.mce[sim1_results$dist == "Gamma"],
                                        paste0("(", sim1_results$Coverage.cil[sim1_results$dist == "Gamma"], ", ", sim1_results$Coverage.ciu[sim1_results$dist == "Gamma"], ")")),
                                      c(sim1_results$Coverage[sim1_results$dist == "Weibull"], sim1_results$Coverage.mce[sim1_results$dist == "Weibull"],
                                        paste0("(", sim1_results$Coverage.cil[sim1_results$dist == "Weibull"], ", ", sim1_results$Coverage.ciu[sim1_results$dist == "Weibull"], ")")),
                                      c(sim1_results$Coverage[sim1_results$dist == "Exponential"], sim1_results$Coverage.mce[sim1_results$dist == "Exponential"],
                                        paste0("(", sim1_results$Coverage.cil[sim1_results$dist == "Exponential"], ", ", sim1_results$Coverage.ciu[sim1_results$dist == "Exponential"], ")"))))

library(knitr); library(kableExtra)
sim1_table <- sim1_table %>% dplyr::mutate(images = paste0("<","img src=","sim1_plot_", ref,".png",">")) %>%
  dplyr::select(-ref, -labels) %>%
  kbl(linesep = "", row.names=T,
      caption="Simulation 1: Summary table (Survival distribution)", align = "c", vline = "",
      col.names = c("Data generating\n mechanism", "Est. Trt\n Effect",
                    "Bias", "Power", "Coverage", "Plot"), escape = FALSE, booktabs = TRUE) %>%
  kable_styling(latex_options="scale_down") %>%
  kable_paper("basic", full_width = F, position = "center") %>%
  pack_rows("Generalised Gamma", 1, 3, color = cols[3], background = cols[1]) %>%
  pack_rows("Gamma", 4, 6, color = cols[3], background = cols[1]) %>%
  pack_rows("Weibull", 7, 9, color = cols[3], background = cols[1]) %>%
  pack_rows("Exponential", 10, 12, color = cols[3], background = cols[1]) %>%
  collapse_rows(columns = c(6), valign = "top") %>%
  footnote(general = "1. Vertical axis of plots is the repetition number.
           2. Dashed line represents true treatment effect.")

# Bivariate plot of estimated treatment effect against standard error
i <- 3
index <- ((Nreps*(i-1)) + 1):(Nreps*i)
plot_data <- data.frame(Est = c(t(x[1, index])),
                        SE = c(t(x[2, index])))
p11 <- ggplot(data = plot_data, aes(x = Est, y = SE)) +
  geom_point(col = "black") +
  theme_classic() +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[i], "model")) +
  theme(legend.position="none")

i <- 2
index <- ((Nreps*(i-1)) + 1):(Nreps*i)
plot_data <- data.frame(Est = c(t(x[1, index])),
                        SE = c(t(x[2, index])))
p22 <- ggplot(data = plot_data, aes(x = Est, y = SE)) +
  geom_point(col = "black") +
  theme_classic() +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[i], "model")) +
  theme(legend.position="none")

i <- 4
index <- ((Nreps*(i-1)) + 1):(Nreps*i)
plot_data <- data.frame(Est = c(t(x[1, index])),
                        SE = c(t(x[2, index])))
p33 <- ggplot(data = plot_data, aes(x = Est, y = SE)) +
  geom_point(col = "black") +
  theme_classic() +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[i], "model")) +
  theme(legend.position="none")

i <- 1
index <- ((Nreps*(i-1)) + 1):(Nreps*i)
plot_data <- data.frame(Est = c(t(x[1, index])),
                        SE = c(t(x[2, index])))
p44 <- ggplot(data = plot_data, aes(x = Est, y = SE)) +
  geom_point(col = "black") +
  theme_classic() +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[i], "model")) +
  theme(legend.position="none")

# Bivariate plot of estimated treatment effect between pairs
i <- 3 
j <- 2
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$SE.j, na.rm = T)
x.max <- max(plot_data$SE.j, na.rm = T)
y.min <- min(plot_data$SE.i, na.rm = T)
y.max <- max(plot_data$SE.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p12 <- ggplot(data = plot_data, aes(x = SE.i, y = SE.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Std. Err. in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[j], "model"))

i <- 3 
j <- 4
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$SE.j, na.rm = T)
x.max <- max(plot_data$SE.j, na.rm = T)
y.min <- min(plot_data$SE.i, na.rm = T)
y.max <- max(plot_data$SE.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p13 <- ggplot(data = plot_data, aes(x = SE.i, y = SE.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Std. Err. in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[j], "model"))

i <- 3 
j <- 1
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$SE.j, na.rm = T)
x.max <- max(plot_data$SE.j, na.rm = T)
y.min <- min(plot_data$SE.i, na.rm = T)
y.max <- max(plot_data$SE.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p14 <- ggplot(data = plot_data, aes(x = SE.i, y = SE.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Std. Err. in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[j], "model"))

i <- 2 
j <- 4
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$SE.j, na.rm = T)
x.max <- max(plot_data$SE.j, na.rm = T)
y.min <- min(plot_data$SE.i, na.rm = T)
y.max <- max(plot_data$SE.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p23 <- ggplot(data = plot_data, aes(x = SE.i, y = SE.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Std. Err. in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[j], "model"))

i <- 2 
j <- 1
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$SE.j, na.rm = T)
x.max <- max(plot_data$SE.j, na.rm = T)
y.min <- min(plot_data$SE.i, na.rm = T)
y.max <- max(plot_data$SE.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p24 <- ggplot(data = plot_data, aes(x = SE.i, y = SE.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Std. Err. in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[j], "model"))

i <- 4 
j <- 1
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$SE.j, na.rm = T)
x.max <- max(plot_data$SE.j, na.rm = T)
y.min <- min(plot_data$SE.i, na.rm = T)
y.max <- max(plot_data$SE.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p34 <- ggplot(data = plot_data, aes(x = SE.i, y = SE.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Std. Err. in", sim1_results$dist[i], "model"),
       y = paste("Std. Err. in", sim1_results$dist[j], "model"))

i <- 2 
j <- 3
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$Est.j, na.rm = T)
x.max <- max(plot_data$Est.j, na.rm = T)
y.min <- min(plot_data$Est.i, na.rm = T)
y.max <- max(plot_data$Est.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p21 <- ggplot(data = plot_data, aes(x = Est.i, y = Est.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Est. Trt Effect in", sim1_results$dist[j], "model"))

i <- 4 
j <- 3
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$Est.j, na.rm = T)
x.max <- max(plot_data$Est.j, na.rm = T)
y.min <- min(plot_data$Est.i, na.rm = T)
y.max <- max(plot_data$Est.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p31 <- ggplot(data = plot_data, aes(x = Est.i, y = Est.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Est. Trt Effect in", sim1_results$dist[j], "model"))

i <- 1 
j <- 3
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$Est.j, na.rm = T)
x.max <- max(plot_data$Est.j, na.rm = T)
y.min <- min(plot_data$Est.i, na.rm = T)
y.max <- max(plot_data$Est.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p41 <- ggplot(data = plot_data, aes(x = Est.i, y = Est.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Est. Trt Effect in", sim1_results$dist[j], "model"))

i <- 4 
j <- 2
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$Est.j, na.rm = T)
x.max <- max(plot_data$Est.j, na.rm = T)
y.min <- min(plot_data$Est.i, na.rm = T)
y.max <- max(plot_data$Est.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p32 <- ggplot(data = plot_data, aes(x = Est.i, y = Est.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Est. Trt Effect in", sim1_results$dist[j], "model"))

i <- 1 
j <- 2
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$Est.j, na.rm = T)
x.max <- max(plot_data$Est.j, na.rm = T)
y.min <- min(plot_data$Est.i, na.rm = T)
y.max <- max(plot_data$Est.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p42 <- ggplot(data = plot_data, aes(x = Est.i, y = Est.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Est. Trt Effect in", sim1_results$dist[j], "model"))

i <- 1 
j <- 4
index.i <- ((Nreps*(i-1)) + 1):(Nreps*i)
index.j <- ((Nreps*(j-1)) + 1):(Nreps*j)
plot_data <- data.frame(Est.i = c(t(x[1, index.i])),
                        SE.i = c(t(x[2, index.i])),
                        Est.j = c(t(x[1, index.j])),
                        SE.j = c(t(x[2, index.j])))
x.min <- min(plot_data$Est.j, na.rm = T)
x.max <- max(plot_data$Est.j, na.rm = T)
y.min <- min(plot_data$Est.i, na.rm = T)
y.max <- max(plot_data$Est.i, na.rm = T)
min <- min(x.min, y.min)
max <- max(x.max, y.max)
p43 <- ggplot(data = plot_data, aes(x = Est.i, y = Est.j)) +
  geom_point() +
  geom_hline(yintercept = alpha + beta*gamma, lty = 2) +
  geom_vline(xintercept = alpha + beta*gamma, lty = 2) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(min,max)) +
  scale_y_continuous(limits = c(min,max)) +
  labs(x = paste("Est. Trt Effect in", sim1_results$dist[i], "model"),
       y = paste("Est. Trt Effect in", sim1_results$dist[j], "model"))

sim1plots <- ggpubr::ggarrange(p11, p21, p31, p41,
                              p12, p22, p32, p42,
                              p13, p23, p33, p43, 
                              p14, p24, p34, p44, nrow = 4, ncol = 4)

# Zip plot: plot of confidence intervals fractionally ranked by mod(zi)
plot_data2 <- data.frame()
for (i in c(1:length(sim1_results$dist))) {

  index <- ((Nreps*(i-1)) + 1):(Nreps*i)
  plot_data <- data.frame(t(x[, index])) %>% dplyr::mutate(ref = rep(i, n()))
  colnames(plot_data) <- c("est", "ses", "cil", "ciu", "bias", "power", "coverage", "prop_die", "ref")
  plot_data <- plot_data %>% dplyr::mutate(modzi = abs(bias/ses)) %>%
    arrange(desc(modzi)) %>% dplyr::mutate(percentile = rev(seq(from = 0.25, to = 100, length.out = n())))
  plot_data$coverage.cil <- rep(sim1_results$Coverage.cil[i], Nreps)
  plot_data$coverage.ciu <- rep(sim1_results$Coverage.ciu[i], Nreps)
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
  facet_wrap(facets = vars(factor(ref, labels = c("Gen Gamma", "Gamma", "Weibull", "Exponential"))),
             ncol = 2, nrow = 2, scales = "fixed", strip.position = "top") +
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
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim1_zipplot.png")

# When the method has 95% coverage (also seen in the table) the colour of the intervals
# switches at the 95th percentile on the vertical axis.

# Horizontal lines are the Monte Carlo 95% confidence intervals for the percentage coverage.