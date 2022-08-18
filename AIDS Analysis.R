################################################
# Title: Joint Modelling Analysis on Aids Data #
# Author: Cameron Wilson (219006159)           #
################################################

# Program description: Applying the JMbayes2 package to the Aids data from joineR package

#### Efficient way of installing previously not installed packages ####
list.of.packages <- c("tidyverse","survival","survminer","lme4",
                      "ggplot2","JM","joineR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))suppressPackageStartupMessages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

#### Set Working Directory ####
setwd("~ /OneDrive/Documents/University/Masters/12. Dissertation Project/Background Reading/JM In R")

#### Set the seed for reproducibility ####
set.seed(123)

#### Read in the aids data from joineR package ####
library(joineR)
aids = as.data.frame(aids)

library(dplyr)
aids <- aids %>% dplyr::rename(id = patient,
                               time = Time)
# This dataset describes a randomized clinical trial (Goldman et al., 1996) 
# in which both survival and longitudinal data were collected to compare the 
# efficacy and safety of two antiretroviral drugs, namely ddI (didanosine) 
# and ddC (zalcitabine), in treating HIV-infected patients intolerant or 
# failing zidovudine (AZT) therapy.

# id - patient identifier. Integer
# time - time to death or censoring. Numeric
# death - event indicator. 0 = right-censoring, 1 = event. Integer
# obstime - measurement times for repeated CD4 count measurements. Numeric
# CD4 - CD4 cell counts measured at obstime. Numeric
# drug - drug indicator. Factor
# gender - gender indicator. Factor
# prevOI - opportunistic infection indicator. AIDS = AIDS diagnosis at entry,
# noAIDS = no previous infection. Factor
# AZT - AZT intolerance / failure indiciator. Factor

# Data is hierarchical as repeat measurements of CD4 cell count over time 
# for each patient. 

table(aids2$obstime) # CD4 cell count measured at baseline, 2 months, 6 months,
# 12 months and 18 months
table(aids2$drug[unique(aids2$id)]) # 239 patients on ddC
# 228 patients on ddI
aids2$id <- as.factor(aids2$id)
aids2$drug <- as.factor(aids2$drug)
aids2$gender <- as.factor(aids2$gender)
aids2$prevOI <- as.factor(aids2$prevOI)
aids2$AZT <- as.factor(aids2$AZT)

# No missing data present - this is a cleaned data set. 
# No real outliers in the continuous variables either.

# # Log transformation of CD4
# hist(aids2$CD4) # right-skewed
# hist(log(aids2$CD4)) # normal is a better approx to log CD4
# # Some CD4 values are 0. Set to 0.001 to allow log transformation. 
# # Perform sensitivity analyis setting to 0.01 & 0.00001 to see if it makes any difference
# aids2$CD4 <- ifelse(aids2$CD4 == 0, 0.001, aids2$CD4)
# aids2 <- aids2 %>% dplyr::mutate(logCD4 = log(CD4))

# Plot trajectory of CD4 cell count by death event or censored
death.labs <- c("Alive/censored","Died")
names(death.labs) <- c(0,1)

library(ggplot2)
p0 <- ggplot(aids2, aes(x = obstime, y = CD4, group = id)) +
  geom_line(lty = 2) +
  facet_grid( ~ death,
              labeller = labeller(death = death.labs)) +
  theme_classic() +
  labs(y = "CD4+ Lymphocyte Count",
       x = "Observation Time (months)")

p1 <- ggplot(aids2, aes(x = obstime, y = CD4, group = id)) +
  geom_line(lty = 2) +
  facet_grid( ~ prevOI) +
  theme_classic() +
  labs(y = "CD4+ Lymphocyte Count",
       x = "Observation Time (months)")

# Plot trajectory of CD4 cell count for a subset
set.seed(1234)
ids <- sample(unique(aids2$id), size = 10, replace = FALSE)
library(dplyr)
aids_sub <- aids2 %>% dplyr::filter(id %in% ids)
p2 <- ggplot(data = aids_sub, aes(x = obstime, y = CD4, group = id, colour = factor(death, labels = c("Alive/Censored", "Died")))) +
  geom_point() +
  geom_line(aes(colour = factor(death, labels = c("Alive/Censored", "Died")))) +
  theme_classic() +
  theme(legend.position = "top") +
  labs(y = "CD4+ Lymphocyte Count",
       x = "Observation Time (months)",
       colour = "Survival outcome") +
  scale_colour_manual(values = cols[1:2]) +
  xlim(0, 15) + 
  ylim(0, 20)

# Plot the trajectory of CD4 cell count over time with mean overlaid
library(dplyr)
aids_mean <- aids2 %>% dplyr::group_by(obstime) %>%
  summarise(mean = mean(CD4))
library(ggplot2)
p3 <- ggplot(data = aids_mean, aes(x = obstime, y = mean)) +
  geom_line(col = cols[1]) +
  geom_line(data = aids2, aes(x = obstime, y = CD4, group = id), lty = 2, alpha = 0.3) +
  theme_classic() +
  labs(y = "CD4+ Lymphocyte Count",
       x = "Observation Time (months)",
       title = "Informative Droput Bias")

# Group exploratory plots 
plots1 <- ggpubr::ggarrange(p0, p1, p2, p3, labels = c("A", "B", "C", "D"),
                            ncol = 2, nrow = 2)

# Summary statistics 
table(aids2$obstime) # summary of number of longitudinal visits
aids.id <- aids2 %>% dplyr::distinct(id, .keep_all = TRUE)
table(aids.id$death, aids.id$drug) # number of deaths vs. censored by treatment 
CD4.sum <- aids2 %>% dplyr::group_by(obstime, drug) %>%
  dplyr::summarise(n = n(),
                   mean = round(mean(CD4), 2),
                   sd = round(sd(CD4), 2)) %>% 
  ungroup()
# mean log-CD4 count over time by treatment, beware that this ignores the correlation within subjects

library(tidyr)
CD4.sum1 <- CD4.sum %>% dplyr::filter(drug == "ddC") %>%
  dplyr::select(-drug) %>% tidyr::gather("statistic", "value", 2:4) %>% 
  arrange(obstime) %>% dplyr::select(-obstime)
CD4.sum1$value[CD4.sum1$statistic == "n"] <- sprintf('%i', CD4.sum1$value[CD4.sum1$statistic == "n"])

CD4.sum2 <- CD4.sum %>% dplyr::filter(drug == "ddI") %>%
  dplyr::select(-drug) %>% gather("statistic", "value", 2:4) %>% 
  arrange(obstime) %>% dplyr::select(-obstime)
CD4.sum2$value[CD4.sum2$statistic == "n"] <- sprintf('%i', CD4.sum2$value[CD4.sum2$statistic == "n"])

CD4.sum.final <- cbind(CD4.sum1, CD4.sum2)[,-3]
colnames(CD4.sum.final) <- c("", "ddC (N=237)", "ddI (N=230)")

table(aids.id$drug)
table(aids.id$drug, aids.id$gender) # gender split by treatment 
table(aids.id$drug, aids.id$prevOI) # AIDS diagnosis split by treatment
table(aids.id$drug, aids.id$AZT) # AZT intolerance / failure split by treatment

# Plot of mean CD4+ count 
CD4.sum <- CD4.sum %>% dplyr::mutate(cil = mean - 1.96*sd,
                                     ciu = mean + 1.96*sd)
p4 <- ggplot(data = aids2, aes(x = factor(obstime), y = CD4)) +
  geom_boxplot(aes(colour = drug)) +
  theme_classic() +
  labs(x = "Follow-up time (months)",
       y = "CD4+ Count",
       title = "Mean CD4+ Count over Follow-up \nby Treatment Group",
       colour = "Treatment Group") +
  scale_colour_manual(values = cols[1:2])

# Survival summary
library(survival); library(survminer)
surv.fit <- survfit(Surv(time, death) ~ drug, data = aids.id) 
ggsurvplot(surv.fit, 
           xlab = "Time (months)",
           ylim = c(0,1)
)
quantile(surv.fit)
print(surv.fit)
survdiff(Surv(time, death) ~ drug, data = aids.id)

# Median follow-up
fit0 <- survfit(Surv(time,1-death) ~ drug, data = aids.id) # switch event and censoring variable to derive median follow-up
median1 <- round(summary(fit0)$table[1,7],1); median1
median2 <- round(summary(fit0)$table[2,7],1); median2

#### 1) Time-varying covariate approach ####
# Allow for delayed entry 
ids <- unique(aids2$id)
aids2$tstart <- aids2$obstime
aids2$tstop <- NA
aids2$death.int <- 0

for (i in c(1:length(aids2$id)-1)) {
  
  aids2$tstop[i] <- aids2$obstime[i+1]
  
}

for (i in c(1:length(ids))) {
  
  n <- length(aids2$obstime[aids2$id == ids[i]])
  aids2$tstop[aids2$id == ids[i]][n] <- aids2$time[aids2$id == ids[i]][n]
  aids2$death.int[aids2$id == ids[i]][n] <- aids2$death[aids2$id == ids[i]][n]
  
}

library(survival)
surv.fit1 <- coxph(Surv(time = tstart, time2 = tstop, event = death.int) ~ drug + CD4 + cluster(id), data = aids2)
summary(surv.fit1)

# Add gender, prevOI, and AZT; then, perform backwards automatic stepwise model selection based on AIC
surv.fit2 <- coxph(Surv(time = tstart, time2 = tstop, event = death.int) ~ drug + CD4 + cluster(id) +
                     gender + prevOI + AZT, data = aids2)
summary(surv.fit2)

library(MASS)
surv.fit3 <- stepAIC(surv.fit2, direction = "both", trace = TRUE,
                     scope = list(lower = ~ drug + CD4,
                                  upper = surv.fit2))

# Summarise results: HR, SE(logHR), 95% CI for HR
HR1 <- round(summary(surv.fit3)$coefficients[1,2], 3) # HR for treatment effect
SE1 <- round(summary(surv.fit3)$coefficients[1,4], 3) # SE for treatment effect coef 
CI1 <- paste0("(", 
              round(exp(confint(surv.fit3)[1,1]), 3), 
              ", ", round(exp(confint(surv.fit3)[1,2]), 3), 
              ")") # 95% CI for the HR 
summary_table <- data.frame(HR = paste0(HR1, " (", SE1, ")"),
                            CI = CI1)

#### 2) Two-stage approach ####
# Model selection based on LR test
# Fit a linear mixed-effects longitudinal model - random intercepts 
library(nlme)
lmeFit1 <- nlme::lme(CD4 ~ obstime + drug, data = aids2, random = ~ 1 | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "ML")
summary(lmeFit1)

# Fit a linear mixed-effects longitudinal model - random intercepts & slopes
lmeFit2 <- nlme::lme(CD4 ~ obstime + drug, data = aids2, random = ~ obstime | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "ML")
summary(lmeFit2)

library(lmtest)
lmtest::lrtest(lmeFit1, lmeFit2) # random slopes are sign to model fit 

# Fit a non-linear mixed-effects longitudinal model - random intercepts & slopes
lmeFit3 <- nlme::lme(CD4 ~ as.numeric(obstime) + drug, data = aids2, random = ~ obstime + I(as.numeric(obstime)^2)| id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "ML")
summary(lmeFit3)
lmtest::lrtest(lmeFit2, lmeFit3) # not significant to the model

# Fit a linear mixed-effects longitudinal model with gender
lmeFit4 <- nlme::lme(CD4 ~ obstime + drug + gender, data = aids2, random = ~ obstime | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "ML")
summary(lmeFit4)
lmtest::lrtest(lmeFit2, lmeFit4) # gender not sign to model fit 

# Fit a linear mixed-effects longitudinal model with prevOI
lmeFit5 <- nlme::lme(CD4 ~ obstime + drug + prevOI, data = aids2, random = ~ obstime | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "ML")
summary(lmeFit5)
lmtest::lrtest(lmeFit2, lmeFit5) # prevOI sign to model fit 

# Fit a linear mixed-effects longitudinal model with AZT
lmeFit6 <- nlme::lme(CD4 ~ obstime + drug + prevOI + AZT, data = aids2, random = ~ obstime | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim", method = "ML"))
summary(lmeFit6)
lmtest::lrtest(lmeFit5, lmeFit6) # AZT not sign to model fit 

# Fit a linear mixed-effects longitudinal model with interaction between drug and obstime - random intercepts & slopes
lmeFit7 <- nlme::lme(CD4 ~ as.numeric(obstime)*drug + prevOI, data = aids2, random = ~ obstime | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "ML")
summary(lmeFit7)
lmtest::lrtest(lmeFit5, lmeFit7) # interaction is not sign to the model fit 

# Final model is lmeFit5, now fit using restricted log-likelihood
lmeFit5 <- nlme::lme(CD4 ~ obstime + drug + prevOI, data = aids2, random = ~ obstime | id,
                     control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50,
                                          opt = "optim"), method = "REML")
summary(lmeFit5)

# Generate individualised predictions of CD4 profile
aids2 <- aids2 %>% dplyr::mutate(y_hat = predict(lmeFit5, newdata = .))

# Fit the Cox survival model with predicted longitudinal response 
# prevOI was sign in TVC approach so include and test using AIC
library(survival)
survFit1 <- coxph(Surv(time = tstart, time2 = tstop, event = death.int) ~ drug + y_hat + prevOI + AZT + gender + cluster(id), data = aids2, x = TRUE)
survFit2 <- stepAIC(survFit1, direction = "both", trace = TRUE,
                    scope = list(lower = ~ drug + y_hat,
                                 upper = survFit1))
summary(survFit2) # prevOI still sign. Has been included in longitudinal and survival submodels alike drug

# Summarise results
HR2 <- round(summary(survFit2)$coefficients[1,2], 3) # HR for treatment effect
SE2 <- round(summary(survFit2)$coefficients[1,4], 3) # SE for treatment effect coef 
CI2 <- paste0("(", 
              round(exp(confint(survFit2)[1,1]), 3), 
              ", ", 
              round(exp(confint(survFit2)[1,2]), 3), 
              ")") # 95% CI for the HR 
summary_table <- rbind(summary_table, c(paste0(HR2, " (", SE2, ")"), CI2))

#### Two-stage structured time points approach ####
outdata <- data.frame(rbind(NA,NA))
for (i in c(3,5,10,100)){
  # Estimate the longitudinal profile at structured timepoints 
  aids_structured <- aids2 %>% dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::select(-obstime, -y_hat, -tstart, -tstop, -death.int) %>%
    dplyr::slice(rep(1:n(), i)) %>%
    dplyr::arrange(id) %>%
    dplyr::mutate(obstime = as.numeric(rep(seq(0,18,length.out = i), length(unique(id))))) %>%
    dplyr::filter(obstime <= time) %>%
    dplyr::mutate(y_hat_structured = predict(lmeFit5, newdata = .))
  
  # Allow for delayed entry 
  ids <- unique(aids_structured$id)
  aids_structured$tstart <- aids_structured$obstime
  aids_structured$tstop <- NA
  aids_structured$death.int <- 0
  
  for (i in c(1:length(aids_structured$id)-1)) {
    aids_structured$tstop[i] <- aids_structured$obstime[i+1]
  }
  
  for (i in c(1:length(ids))) {
    n <- length(aids_structured$obstime[aids_structured$id == ids[i]])
    aids_structured$tstop[aids_structured$id == ids[i]][n] <- aids_structured$time[aids_structured$id == ids[i]][n]
    aids_structured$death.int[aids_structured$id == ids[i]][n] <- aids_structured$death[aids_structured$id == ids[i]][n]
  }
  
  # Fit the Cox survival model with predicted longitudinal response
  library(survival)
  survFit_structured <- coxph(Surv(time = tstart, time2 = tstop, event = death.int) ~ drug + y_hat_structured + prevOI + cluster(id), data = aids_structured, x = TRUE)
  
  # Summarise results
  outdata <- data.frame(HR = paste0(round(summary(survFit_structured)$coefficients[1,2],3),
                                    " (",
                                    round(summary(survFit_structured)$coefficients[1,4],3),
                                    ")"),
                        CI = paste0("(", round(exp(confint(survFit_structured)[1,1]), 3), 
                                    ", ", 
                                    round(exp(confint(survFit_structured)[1,2]), 3), 
                                    ")"))
  summary_table <- rbind(summary_table, outdata)
  
}

rownames(summary_table) <- c("TVC", "Two-stage", 
                             "Str Two-stage (3 visits)",
                             "Str Two-stage (5 visits)",
                             "Str Two-stage (10 visits)",
                             "Str Two-stage (100 visits)")

# Summary plots 
aids_sub2 <- aids2 %>% dplyr::filter(id == 10)
aids_structured2 <- aids2 %>% dplyr::distinct(id, .keep_all = TRUE) %>%
  dplyr::filter(id == 10) %>%
  dplyr::select(-obstime, -y_hat) %>%
  dplyr::slice(rep(1:n(), 10)) %>%
  dplyr::arrange(id) %>%
  dplyr::mutate(obstime = as.numeric(rep(seq(0,18,length.out = 10), length(unique(id))))) %>%
  dplyr::filter(obstime <= time) %>%
  dplyr::mutate(y_hat_structured = predict(lmeFit7, newdata = .))

for (i in c(1:8)) {
  aids_structured2$stop[i] <- aids_structured2$obstime[i+1]
}
aids_structured2$stop[9] <- aids_structured2$time[9]

for (i in c(1:3)) {
  aids_sub2$stop[i] <- aids_sub2$obstime[i+1]
}
aids_sub2$stop[4] <- aids_sub2$time[4]

library(ggplot2)
p5 <- ggplot(data = aids_sub2, aes(x = obstime, y = y_hat, 
                                   xend = stop, yend = y_hat)) +
  geom_segment(colour = cols[2], size = 2) +
  geom_point(aes(x = time, y = tail(y_hat, 1)),
             shape = 4, size = 6) +
  geom_segment(data = aids_structured2, aes(x = obstime, y = y_hat_structured, 
                                            xend = stop, 
                                            yend = y_hat_structured),
               colour = cols[1], size = 1) +
  geom_line(data = aids_sub2, aes(x = obstime, y = CD4, 
                                  xend = stop, yend = CD4), colour = "black", lty = 2) +
  geom_point(data = aids_structured2, aes(x = time, 
                                          y = tail(y_hat_structured, 1), 
                                          xend = stop, yend = y_hat_structured),
             shape = 4, size = 6) +
  theme_classic() +
  labs(title = "Stepwise Longitudinal Profile Assumed By Cox Model",
       x = "Follow-up time (months)", y = "CD4 Cell Count",
       caption = "Crosses (X) indicate when patient 10 was censored.")

#### Joint model approach ####
# Mixed-effects longitudinal model - use chosen model from two-stage approach
library(nlme)
lmeFit <- nlme::lme(CD4 ~ as.numeric(obstime) + drug + prevOI, data = aids2, 
                    random = ~ obstime | id,
                    control = lmeControl(maxIter = 100, msMaxIter = 100, 
                                         niterEM = 50,opt = "optim"))
summary(lmeFit)

# Cox survival model
library(survival)
survFit <- coxph(Surv(time, death) ~ drug + prevOI, data = aids.id, x = TRUE)
# include prevOI as was sign in two-stage approach
summary(survFit)

# Joint model - Maximum likelihood framework
library(JM)
set.seed(1234)
jointFit1 <- jointModel(lmeFit, survFit, timeVar = "obstime",
                        method = "Cox-PH-GH")
summary(jointFit1)

# set.seed(123)
# jointFit2 <- jointModel(lmeFit, survFit, timeVar = "obstime",
#                         method = "piecewise-PH-aGH")
# summary(jointFit2)
# 
# set.seed(123)
# jointFit3 <- jointModel(lmeFit, survFit, timeVar = "obstime",
#                         method = "Cox-PH-GH")
# summary(jointFit3)
# 
# set.seed(123)
# jointFit4 <- jointModel(lmeFit, survFit, timeVar = "obstime",
#                         method = "spline-PH-aGH")
# summary(jointFit4)

# # Joint model - Bayesian framework
# library(JMbayes2)
# jointFit2 <- jm(survFit2, lmeFit1, time_var = "obstime")
# summary(jointFit2)

#### Summary table ####
alpha_est <- jointFit1$coefficients$gammas[1]
beta_est <- jointFit1$coefficients$alpha
gamma_est <- jointFit1$coefficients$betas[3]
est <- alpha_est + beta_est*gamma_est
mean <- c(alpha_est, beta_est, gamma_est)
cov <- vcov(jointFit1)[c("T.drugddI","T.alpha","Y.drugddI"), c("T.drugddI","T.alpha","Y.drugddI")]
names(mean) <- c("alpha_est", "beta_est", "gamma_est")
colnames(cov) <- c("alpha_est", "beta_est", "gamma_est")
rownames(cov) <- c("alpha_est", "beta_est", "gamma_est")
library(msm)
ses <- deltamethod(~ x1 + x2*x3, mean = mean, cov = cov, ses = TRUE) # calculate SE for non-linear combination of parameters using delta method
cil <- est - 1.96*ses
ciu <- est + 1.96*ses

jm <- data.frame(HR = paste0(round(exp(est), 3),
                             " (",
                             round(ses, 3),
                             ")"),
                 CI = paste0("(",
                             round(exp(cil), 3),
                             ", ",
                             round(exp(ciu), 3),
                             ")"))

summary_table <- rbind(summary_table, jm)
rownames(summary_table)[7] <- c("Joint model")

# SE much larger in joint model because it is accounting for uncertainity
# in both longitudinal and survival processes.

# By increasing the number of timepoints in two-stage approach parameter estimates
# tend towards (asymptotically) the estimates of the joint model approach. But by 
# increasing the number of timepoints the SE is reduced and estimates are overly precise.

# Two-stage approach would lead to an incorrect conclusion of a statistically significant
# treatment effect (CI excludes 1). By correctly capturing the uncertainty in both processes
# using a joint model our CI includes 1. Therefore, demonstrating joint models lead to more
# reliable inferences and estimation of treatment effect.

#### Data table ####
data_table1 <- aids2 %>% dplyr::filter(id %in% c(21, 258)) %>%
  dplyr::select(id, drug, gender, prevOI, AZT, tstart, tstop, death.int, CD4) # tvc
data_table2 <- aids2 %>% dplyr::filter(id %in% c(21, 258)) %>%
  dplyr::select(id, drug, gender, prevOI, AZT, obstime, CD4) # 2-stage longitudinal & joint model longitudinal
data_table3 <- aids2 %>% dplyr::filter(id %in% c(21, 258)) %>%
  dplyr::select(id, drug, gender, prevOI, AZT, tstart, tstop, death.int, y_hat) # 2-stage survival
data_table4 <- aids.id %>% dplyr::filter(id %in% c(21, 258)) %>%
  dplyr::select(id, drug, gender, prevOI, AZT, time, death) # joint model 
i <- 10
data_table5 <- aids2 %>% dplyr::distinct(id, .keep_all = TRUE) %>%
  dplyr::select(-obstime, -y_hat, -tstart, -tstop, -death.int) %>%
  dplyr::slice(rep(1:n(), i)) %>%
  dplyr::arrange(id) %>%
  dplyr::mutate(obstime = as.numeric(rep(seq(0,18,length.out = i), length(unique(id))))) %>%
  dplyr::filter(obstime <= time) %>%
  dplyr::mutate(y_hat_structured = predict(lmeFit5, newdata = .))

# Allow for delayed entry 
ids <- unique(data_table5$id)
data_table5$tstart <- data_table5$obstime
data_table5$tstop <- NA
data_table5$death.int <- 0

for (i in c(1:length(data_table5$id)-1)) {
  data_table5$tstop[i] <- data_table5$obstime[i+1]
}

for (i in c(1:length(ids))) {
  n <- length(data_table5$obstime[data_table5$id == ids[i]])
  data_table5$tstop[data_table5$id == ids[i]][n] <- data_table5$time[data_table5$id == ids[i]][n]
  data_table5$death.int[data_table5$id == ids[i]][n] <- data_table5$death[data_table5$id == ids[i]][n]
}
data_table5 <- data_table5 %>% dplyr::filter(id %in% c(21, 258)) %>%
  dplyr::select(id, drug, gender, prevOI, AZT, tstart, tstop, death.int, y_hat_structured) # two-stage structured with 10 timepoints

#### Dynamic prediction ####
# Data for patient 454
data454 <- aids %>% dplyr::filter(id == "454")
data454$time <- data454$death <- NULL

###
sfit <- survfitJM(jointFit1, newdata = data454[1:2,])
sfit
plot(sfit, include.y = TRUE)
###

# dynamic predictions, each time an extra measurement
n <- nrow(data454)
for (i in seq_len(n)) {
  Spred <- predict(jointFit2, newdata = data454[1:i, ],
                   process = "event",
                   times = seq(0, 24, length.out = 51),
                   return_newdata = TRUE)
  Lpred <- predict(jointFit2, newdata = data454[1:i, ],
                   times = seq(0, 24, length.out = 51),
                   return_newdata = TRUE)
  
  plot(Lpred, Spred, fun_event = function(x) 1 - x,
       ylab_event = "Survival Probabilities")
}

#### Model performance ####

# ROC at 12 months with 6 month window
roc <- tvROC(jointFit1, newdata = aids, Tstart = 12, Dt = 6)
roc
plot(roc)
# AUC at 12 months with 6 month window
tvAUC(roc)

# Calibration plot at 12 months with 6 month window
calibration_plot(jointFit1, newdata = aids, Tstart = 12, Dt = 6)

# Brier score at 12 months with 6 month window
tvBrier(jointFit1, newdata = aids, Tstart = 12, Dt = 6)

