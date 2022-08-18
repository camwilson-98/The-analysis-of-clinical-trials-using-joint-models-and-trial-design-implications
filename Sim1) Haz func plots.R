cols <- c("#7D1D3F", "#3A506B", "#D5CAD6")

library(expint)

plotfun <- function(dist) { 
  
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
  
  if (dist == "Gen Gamma") {
    
    shape_param  <- shape  # alpha
    shape2_param <- shape2  # p
    lambda_param <- lambda+0.25  # lambda 
    
  } else if (dist == "Gamma") {
    
    shape_param  <- shape        
    shape2_param <- 1           # Gen Gamma reduces to Gamma when shape2 = 1
    lambda_param <- lambda+1.10
    
  } else if (dist == "Weibull") {
    
    # shape_param  <- shape
    # shape2_param <- shape       # Gen Gamma reduces to Weibull when shape = shape2
    # lambda_param <- lambda
    
    shape_param <- 0.5
    shape2_param <- 0.5
    lambda_param <- 2
    
  } else if (dist == "Exponential") {
    
    shape_param  <- 1
    shape2_param <- 1           # Gen Gamma reduces to Exponential when shape = shape2 = 1
    lambda_param <- lambda+0.15
    
  }
  
  betas <- data.frame(
    shape                = rep(shape_param, 100),
    scale                = rep(lambda_param, 100),
    shape2               = rep(shape2_param, 100), 
    betaEvent_intercept  = rep(0, 100),
    betaEvent_binary     = rep(alpha, 100),
    betaEvent_assoc      = rep(beta, 100),
    betaLong_binary      = rep(gamma, 100),
    betaLong_intercept   = rep(0, 100),
    betaLong_slope       = rep(0, 100),
    betaLong_quadratic   = rep(0, 100),
    betaLong_cubic       = rep(0, 100),
    betaLong_quartic     = rep(0, 100)
  )
  
  t <- seq(from = 0, to = 2, length.out = 100)
  covdata <- data.frame(Z1 = rep(1, 100))
  
  y <- haz(t, x = covdata, betas)
  y
  
}

lambda <- 0.85 # baseline hazard rate in survival model
shape <- 2
shape2 <- 1.5

plotdata1 <- data.frame(t = seq(from = 0, to = 2, length.out = 100),
                        dist = rep("Exponential", 100),
                        haz = plotfun(dist = "Exponential"))
plotdata2 <- data.frame(t = seq(from = 0, to = 2, length.out = 100),
                        dist = rep("Weibull", 100),
                        haz = plotfun(dist = "Weibull"))
plotdata3 <- data.frame(t = seq(from = 0, to = 2, length.out = 100),
                        dist = rep("Gamma", 100),
                        haz = plotfun(dist = "Gamma"))
plotdata4 <- data.frame(t = seq(from = 0, to = 2, length.out = 100),
                        dist = rep("Gen Gamma", 100),
                        haz = plotfun(dist = "Gen Gamma"))

plotdata <- rbind(plotdata1, plotdata2, plotdata3, plotdata4)

library(ggplot2)
hazplot <- ggplot(data = plotdata, aes(x = t, y = haz, colour = dist)) + 
  geom_smooth(se = FALSE) +
  theme_classic() +
  labs(title = "Simulation 1: Hazard functions of the \nparametric distributions",
       x = "Time (years)",
       y = "Hazard function",
       colour = "Distribution") +
  scale_colour_manual(values=c("black", cols))
hazplot
ggsave("~/OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R/Simulations/Output Tables/sim1_hazfunc_plot.png")

######### TEST TO BALANCE PROPORTION THAT DIE
# Function to carry out the simulations
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
    
    # shape_param  <- shape
    # shape2_param <- shape       # Gen Gamma reduces to Weibull when shape = shape2
    # lambda_param <- lambda
    
    shape_param <- 0.5
    shape2_param <- 0.5
    lambda_param <- 2
    
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
  
  # Cox proportional hazards survival model
  dat.id <- sim_dat %>% dplyr::distinct(id, .keep_all = TRUE)
  
  # Include the proportion of events to adjust survival distribution parameters if needed.
  # Need a similar number of events to ensure comparability between the distributions. 
  prop_die <- prop.table(table(dat.id$status))[2]
  
  prop_die
  
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
                       prop_die = rep(NA,4))

Nreps <- 5

for (i in c(1:length(sim1_dat$ref))) { 
  set.seed(1234)
  x <- replicate(Nreps, sim1_run(N = 200, p1, beta, gamma, alpha, lambda, shape, shape2, 
                            a, f, mu_theta, covar_theta, sigmasq_e,
                            dist = sim1_dat$dist[i]))
  sim1_dat$prop_die[i] <- mean(x)
}
sim1_dat
