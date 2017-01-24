#/* 
# * Colin Olito. Created 13/01/2016
# * Analysis of 1st flume experiment: N-invest
# * 
# * NOTES:  This file contains all the necessary
# * 		code to read in the STAN sample files
# * 		and perform posterior predictive checks
# *         model selection, and produce regression
# *         plots for the NxRate flume data
# * 
# */

rm(list=ls())
#################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

###############
# DEPENDENCIES
source('R/functions.R')

# Import NxRate Data Set
data <- read.csv('data/NxRate_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)
# head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Rate      <-  factor(data$Rate)
data$EggPos    <-  factor(data$EggPos)
data$Lane      <-  factor(data$Lane)
data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
data$Date      <-  dmy(data$Date)
data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

# Plot Colors for X^2 Posterior Predictive Check Plots
COLS  <-  c("#009ce0",
            "#f9684a",
            "#8be36a",
            "#7a0058",
            "#9fa3ff",
            "#be0042",
            "#2d1956"
)


#####################################
# Import model results from the stan 
# sample_files for further analysis 
#####################################

csvFiles  <-  c('./output/StanFits/NxRate_m1.csv1',
                './output/StanFits/NxRate_m1.csv2',
                './output/StanFits/NxRate_m1.csv3')
m1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m2.csv1',
                './output/StanFits/NxRate_m2.csv2',
                './output/StanFits/NxRate_m2.csv3')
m2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m3.csv1',
                './output/StanFits/NxRate_m3.csv2',
                './output/StanFits/NxRate_m3.csv3')
m3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m3a.csv1',
                './output/StanFits/NxRate_m3a.csv2',
                './output/StanFits/NxRate_m3a.csv3')
m3a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m4.csv1',
                './output/StanFits/NxRate_m4.csv2',
                './output/StanFits/NxRate_m4.csv3')
m4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m4a.csv1',
                './output/StanFits/NxRate_m4a.csv2',
                './output/StanFits/NxRate_m4a.csv3')
m4a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m5.csv1',
                './output/StanFits/NxRate_m5.csv2',
                './output/StanFits/NxRate_m5.csv3')
m5        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m6.csv1',
                './output/StanFits/NxRate_m6.csv2',
                './output/StanFits/NxRate_m6.csv3')
m6        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m6a.csv1',
                './output/StanFits/NxRate_m6a.csv2',
                './output/StanFits/NxRate_m6a.csv3')
m6a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m7.csv1',
                './output/StanFits/NxRate_m7.csv2',
                './output/StanFits/NxRate_m7.csv3')
m7        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m7a.csv1',
                './output/StanFits/NxRate_m7a.csv2',
                './output/StanFits/NxRate_m7a.csv3')
m7a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m8.csv1',
                './output/StanFits/NxRate_m8.csv2',
                './output/StanFits/NxRate_m8.csv3')
m8        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m8a.csv1',
                './output/StanFits/NxRate_m8a.csv2',
                './output/StanFits/NxRate_m8a.csv3')
m8a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m9.csv1',
                './output/StanFits/NxRate_m9.csv2',
                './output/StanFits/NxRate_m9.csv3')
m9        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m9a.csv1',
                './output/StanFits/NxRate_m9a.csv2',
                './output/StanFits/NxRate_m9a.csv3')
m9a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)




############################
##  Check Convergence of
##  Fixed Effects Estimates
############################
print(m1, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m2, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m4, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m4a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m5, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));






##########################################################################
# Model: m1
##########################################################################

##############
# Diagnostics

# Model Results
#print(m1)[,9:10]
print(m1, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m1.df    <-  as.data.frame(extract(m1))  
m1.summ  <-  plyr:::adply(as.matrix(m1.df),2,MCMCsum)[-1,]

# Simple Diagnostic Plots
plot(m1, pars="beta")
pairs(m1, pars="beta")
par(mfrow=c(5,5))
rstan::traceplot(m1, pars=c("beta", "tau_run"), inc_warmup=FALSE)

# Explore Correlation structure
corrMat  <-  matrix(m1.summ[2051:3650,2], ncol=40,nrow=40)
corrplot(corrMat , method='circle', type='upper')
abline(v=20.5)
abline(h=20.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=20.5)
abline(h=20.5)

##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m1, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m1, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
plot(m1, pars="tau_run")


#########################################
# LOO Log-likelihood for model selection

m1LL     <-  extract_log_lik(m1, parameter_name = "log_lik")
m1Loo    <-  loo(m1LL)
m1WAIC   <-  waic(m1LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m1.df[1,5891:6010])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m1.df[i,5891:6010])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m1.df[,6011], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m1.df[,6012]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m1.df[,6013]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m1.df[,6014]), xlim=c(min(m1.df[,6014],sd(data$nFert)),max(m1.df[,6014],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m1, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data
X2data1   <-  as.numeric(m1.df[,6379])
X2sim1    <-  as.numeric(m1.df[,6380])


##########################################################################
# Model: m2
##########################################################################

##############
# Diagnostics

# Model Results
# print(m2)
print(m2, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m2, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m2.df    <-  as.data.frame(extract(m2))
m2.summ  <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)[-1,]

# Simple Diagnostic Plots
plot(m2, pars="beta")
pairs(m2, pars="beta")
par(mfrow=c(5,5))
rstan::traceplot(m2, c("beta"), inc_warmup=FALSE)


# Explore Correlation structure
corrMat  <-  matrix(m2.summ[1241:2140,2], ncol=30,nrow=30)
corrplot(corrMat , method='circle', type='upper')
abline(v=15.5)
abline(h=15.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=15.5)
abline(h=15.5)

##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m2, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m2, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
plot(m2, pars="tau_run")


#########################################
# LOO Log-likelihood for model selection

m2LL    <-  extract_log_lik(m2, parameter_name = "log_lik")
m2Loo   <-  loo(m2LL)
m2WAIC  <-  waic(m2LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m2.df[1,3581:3700])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m2.df[i,3581:3700])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m2.df[,3701], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m2.df[,3702]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m2.df[,3703]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m2.df[,3704]), xlim=c(min(m2.df[,3704],sd(data$nFert)),max(m2.df[,3704],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m2, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data2   <-  as.numeric(m2.df[,4069])
X2sim2    <-  as.numeric(m2.df[,4070])






##########################################################################
# Model: m3
##########################################################################

##############
# Diagnostics

# Model Results
# print(m3)
print(m3, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m3, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m3.df    <-  as.data.frame(extract(m3))
m3.summ  <-  plyr:::adply(as.matrix(m3.df),2,MCMCsum)[-1,]
m3.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3)

# Simple Diagnostic Plots
plot(m3, pars="beta")
pairs(m3, pars="beta")
rstan::traceplot(m3, c("beta"), inc_warmup=FALSE)

print(m3, pars="gamma")

# Explore Correlation structure
corrMat  <-  matrix(m3.summ[630:1029,2], ncol=20,nrow=20)
corrplot(corrMat , method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)


##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m3, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m3, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
rstan::plot(m3, pars="tau_run")
rstan::traceplot(m3,pars="tau_run", inc_warmup=FALSE)

#########################################
# LOO Log-likelihood for model selection

m3LL     <-  extract_log_lik(m3, parameter_name = "log_lik")
m3Loo    <-  loo(m3LL)
m3WAIC   <-  waic(m3LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m3.df[1,1871:1990])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m3.df[i,1871:1990])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m3.df[,1991], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1992]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1993]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1994]), xlim=c(min(m3.df[,1994],sd(data$nFert)),max(m3.df[,1994],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m3, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data3   <-  as.numeric(m3.df[,2359])
X2sim3    <-  as.numeric(m3.df[,2360])






##########################################################################
# Model: m3a
##########################################################################

##############
# Diagnostics

# Model Results
# print(m3a)
print(m3a, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m3a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m3a.df    <-  as.data.frame(extract(m3a))
m3a.summ  <-  plyr:::adply(as.matrix(m3a.df),2,MCMCsum)[-1,]
m3a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3a)

# Simple Diagnostic Plots
plot(m3a, pars="beta")
pairs(m3a, pars="gamma")
rstan::traceplot(m3a, c("beta"), inc_warmup=FALSE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m3aLL     <-  extract_log_lik(m3a, parameter_name = "log_lik")
m3aLoo    <-  loo(m3aLL)
m3aWAIC   <-  waic(m3aLL)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m3a.df[1,271:390])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m3a.df[i,271:390])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m3a.df[,391], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m3a.df[,392]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3a.df[,393]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3a.df[,394]), xlim=c(min(m3a.df[,394],sd(data$nFert)),max(m3a.df[,394],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m3a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data3a   <-  as.numeric(m3a.df[,759])
X2sim3a    <-  as.numeric(m3a.df[,760])







##########################################################################
# Model: m4
##########################################################################

##############
# Diagnostics

# Model Results
# print(m4)
print(m4, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m4, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m4.df    <-  as.data.frame(extract(m4))
m4.summ  <-  plyr:::adply(as.matrix(m4.df),2,MCMCsum)[-1,]


# Simple Diagnostic Plots
plot(m4, pars="beta")
pairs(m4, pars="beta")
rstan::traceplot(m4.mcmc, c("beta"), inc_warmup=FALSE)



# Explore Correlation structure
corrMat  <-  matrix(m4.summ[631:1030,2], ncol=16,nrow=16)
corrplot(corrMat , method='circle', type='upper')
abline(v=8.5)
abline(h=8.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=8.5)
abline(h=8.5)


##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m4, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m4, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
plot(m4, pars="tau_run")
rstan::traceplot(m4, pars="tau_run", inc_warmup=FALSE)

#########################################
# LOO Log-likelihood for model selection

m4LL  <-  extract_log_lik(m4, parameter_name = "log_lik")
m4Loo    <-  loo(m4LL)
m4WAIC   <-  waic(m4LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m4.df[1,763:880])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m4.df[i,763:880])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m4.df[,881], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m4.df[,882]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m4.df[,883]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m4.df[,884]), xlim=c(min(m4.df[,884],sd(data$nFert)),max(m4.df[,884],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m4, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data4   <-  as.numeric(m4.df[,1249])
X2sim4    <-  as.numeric(m4.df[,1250])




##########################################################################
# Model: m4a
##########################################################################

##############
# Diagnostics

# Model Results
# print(m4a)
print(m4a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m4a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m4a.df    <-  as.data.frame(extract(m4a))
m4a.summ  <-  plyr:::adply(as.matrix(m4a.df),2,MCMCsum)[-1,]


# Simple Diagnostic Plots
plot(m4a, pars="beta")
pairs(m4a, pars="beta")
rstan::traceplot(m4a, c("beta"), inc_warmup=FALSE)


#########################################
# LOO Log-likelihood for model selection

m4aLL  <-  extract_log_lik(m4a, parameter_name = "log_lik")
m4aLoo    <-  loo(m4aLL)
m4aWAIC   <-  waic(m4aLL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m4a.df[1,141:260])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m4a.df[i,141:260])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4a.summ
par(mfrow=c(2,2))
plot(density(m4a.df[,261], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m4a.df[,262]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m4a.df[,263]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m4a.df[,264]), xlim=c(min(m4a.df[,264],sd(data$nFert)),max(m4a.df[,264],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m4a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data4a   <-  as.numeric(m4a.df[,629])
X2sim4a    <-  as.numeric(m4a.df[,630])






##########################################################################
# Model: m5
##########################################################################

##############
# Diagnostics

# Model Results
# print(m5)
print(m5, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m5, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m5.df    <-  as.data.frame(extract(m5))
m5.summ  <-  plyr:::adply(as.matrix(m5.df),2,MCMCsum)[-1,]


# Simple Diagnostic Plots
plot(m5, pars="beta")
pairs(m5, pars="beta")
rstan::traceplot(m5, c("beta"), inc_warmup=FALSE)


#########################################
# LOO Log-likelihood for model selection

m5LL     <-  extract_log_lik(m5, parameter_name = "log_lik")
m5Loo    <-  loo(m5LL)
m5WAIC   <-  waic(m5LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m5.df[1,130:249])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m5.df[i,130:249])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m5.df[,250], adjust=5), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m5.df[,251]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m5.df[,252]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m5.df[,253]), xlim=c(min(m5.df[,253],sd(data$nFert)),max(m5.df[,253],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m5, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()

# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data5   <-  as.numeric(m5.df[,618])
X2sim5    <-  as.numeric(m5.df[,619])




######################################################################################################
######################################################################################################
#  Models including random effects including EggPos




##########################################################################
# Model: m6
##########################################################################

##############
# Diagnostics

# Model Results
# print(m6)
print(m6, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m6, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m6.df    <-  as.data.frame(extract(m6))
m6.summ  <-  plyr:::adply(as.matrix(m6.df),2,MCMCsum)[-1,]
m6.mcmc  <-  rstan:::as.mcmc.list.stanfit(m6)

# Simple Diagnostic Plots
plot(m6, pars="beta")
pairs(m6, pars="beta")
rstan::traceplot(m6, c("beta"), inc_warmup=FALSE)

print(m6, pars="gamma")

# Explore Correlation structure
corrMat  <-  matrix(m6.summ[630:1029,2], ncol=20,nrow=20)
corrplot(corrMat , method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)


##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m6, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m6, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
rstan::plot(m6, pars="tau_run")
rstan::traceplot(m6,pars="tau_run", inc_warmup=FALSE)

#########################################
# LOO Log-likelihood for model selection

m6LL     <-  extract_log_lik(m6, parameter_name = "log_lik")
m6Loo    <-  loo(m6LL)
m6WAIC   <-  waic(m6LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m6.df[1,1871:1990])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m6.df[i,1871:1990])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m6.df[,1991], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m6.df[,1992]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m6.df[,1993]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m6.df[,1994]), xlim=c(min(m6.df[,1994],sd(data$nFert)),max(m6.df[,1994],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m6, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data6   <-  as.numeric(m6.df[,2359])
X2sim6    <-  as.numeric(m6.df[,2360])






##########################################################################
# Model: m6a
##########################################################################

##############
# Diagnostics

# Model Results
# print(m6a)
print(m6a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m6a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m6a.df    <-  as.data.frame(extract(m6a))
m6a.summ  <-  plyr:::adply(as.matrix(m6a.df),2,MCMCsum)[-1,]
# m6a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m6a)

# Simple Diagnostic Plots
plot(m6a, pars="beta")
plot(m6a, pars="gamma")
rstan::traceplot(m6a, c("beta"), inc_warmup=FALSE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m6aLL     <-  extract_log_lik(m6a, parameter_name = "log_lik")
m6aLoo    <-  loo(m6aLL)
m6aWAIC   <-  waic(m6aLL)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m6a.df[1,291:410])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m6a.df[i,291:410])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m6a.df[,411], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m6a.df[,412]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m6a.df[,413]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m6a.df[,414]), xlim=c(min(m6a.df[,414],sd(data$nFert)),max(m6a.df[,414],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m6a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data6a   <-  as.numeric(m6a.df[,779])
X2sim6a    <-  as.numeric(m6a.df[,780])




##########################################################################
# Model: m7
##########################################################################

##############
# Diagnostics

# Model Results
# print(m7)
print(m7, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m7, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m7.df    <-  as.data.frame(extract(m7))
m7.summ  <-  plyr:::adply(as.matrix(m7.df),2,MCMCsum)[-1,]
m7.mcmc  <-  rstan:::as.mcmc.list.stanfit(m7)

# Simple Diagnostic Plots
plot(m7, pars="beta")
pairs(m7, pars="beta")
rstan::traceplot(m7, c("beta"), inc_warmup=FALSE)

print(m7, pars="gamma")

# Explore Correlation structure
corrMat  <-  matrix(m7.summ[630:1029,2], ncol=20,nrow=20)
corrplot(corrMat , method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)


##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m7, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m7, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
rstan::plot(m7, pars="tau_run")
rstan::traceplot(m7,pars="tau_run", inc_warmup=FALSE)

#########################################
# LOO Log-likelihood for model selection

m7LL     <-  extract_log_lik(m7, parameter_name = "log_lik")
m7Loo    <-  loo(m7LL)
m7WAIC   <-  waic(m7LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m7.df[1,1871:1990])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m7.df[i,1871:1990])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m7.df[,1991], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m7.df[,1992]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m7.df[,1993]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m7.df[,1994]), xlim=c(min(m7.df[,1994],sd(data$nFert)),max(m7.df[,1994],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m7, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data7   <-  as.numeric(m7.df[,2359])
X2sim7    <-  as.numeric(m7.df[,2360])






##########################################################################
# Model: m7a
##########################################################################

##############
# Diagnostics

# Model Results
# print(m7a)
print(m7a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m7a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m7a.df    <-  as.data.frame(extract(m7a))
m7a.summ  <-  plyr:::adply(as.matrix(m7a.df),2,MCMCsum)[-1,]
# m7a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m7a)

# Simple Diagnostic Plots
plot(m7a, pars="beta")
plot(m7a, pars="gamma0")
pairs(m7a, pars="gamma1")
rstan::traceplot(m7a, c("beta"), inc_warmup=FALSE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m7aLL     <-  extract_log_lik(m7a, parameter_name = "log_lik")
m7aLoo    <-  loo(m7aLL)
m7aWAIC   <-  waic(m7aLL)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m7a.df[1,281:400])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m7a.df[i,281:400])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m7a.df[,401], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m7a.df[,402]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m7a.df[,403]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m7a.df[,404]), xlim=c(min(m7a.df[,404],sd(data$nFert)),max(m7a.df[,404],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m7a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data7a   <-  as.numeric(m7a.df[,769])
X2sim7a    <-  as.numeric(m7a.df[,770])




##########################################################################
# Model: m8
##########################################################################

##############
# Diagnostics

# Model Results
# print(m8)
print(m8, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m8, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m8.df    <-  as.data.frame(extract(m8))
m8.summ  <-  plyr:::adply(as.matrix(m8.df),2,MCMCsum)[-1,]
m8.mcmc  <-  rstan:::as.mcmc.list.stanfit(m8)

# Simple Diagnostic Plots
plot(m8, pars="beta")
pairs(m8, pars="beta")
rstan::traceplot(m8, c("beta"), inc_warmup=FALSE)

print(m8, pars="gamma")

# Explore Correlation structure
corrMat  <-  matrix(m8.summ[630:1029,2], ncol=20,nrow=20)
corrplot(corrMat , method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)


##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m8, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m8, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
rstan::plot(m8, pars="tau_run")
rstan::traceplot(m8,pars="tau_run", inc_warmup=FALSE)

#########################################
# LOO Log-likelihood for model selection

m8LL     <-  extract_log_lik(m8, parameter_name = "log_lik")
m8Loo    <-  loo(m8LL)
m8WAIC   <-  waic(m8LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m8.df[1,1871:1990])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m8.df[i,1871:1990])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m8.df[,1991], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m8.df[,1992]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m8.df[,1993]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m8.df[,1994]), xlim=c(min(m8.df[,1994],sd(data$nFert)),max(m8.df[,1994],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m8, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data8   <-  as.numeric(m8.df[,2359])
X2sim8    <-  as.numeric(m8.df[,2360])






##########################################################################
# Model: m8a
##########################################################################

##############
# Diagnostics

# Model Results
# print(m8a)
print(m8a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m8a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m8a.df    <-  as.data.frame(extract(m8a))
m8a.summ  <-  plyr:::adply(as.matrix(m8a.df),2,MCMCsum)[-1,]
# m8a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m8a)

# Simple Diagnostic Plots
plot(m8a, pars="beta")
pairs(m8a, pars="gamma")
rstan::traceplot(m8a, c("beta"), inc_warmup=FALSE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m8aLL     <-  extract_log_lik(m8a, parameter_name = "log_lik")
m8aLoo    <-  loo(m8aLL)
m8aWAIC   <-  waic(m8aLL)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m8a.df[1,291:410])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m8a.df[i,291:410])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m8a.df[,411], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m8a.df[,412]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m8a.df[,413]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m8a.df[,414]), xlim=c(min(m8a.df[,414],sd(data$nFert)),max(m8a.df[,414],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m8a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data8a   <-  as.numeric(m8a.df[,779])
X2sim8a    <-  as.numeric(m8a.df[,780])




##########################################################################
# Model: m9
##########################################################################

##############
# Diagnostics

# Model Results
# print(m9)
print(m9, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m9, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m9.df    <-  as.data.frame(extract(m9))
m9.summ  <-  plyr:::adply(as.matrix(m9.df),2,MCMCsum)[-1,]
m9.mcmc  <-  rstan:::as.mcmc.list.stanfit(m9)

# Simple Diagnostic Plots
plot(m9, pars="beta")
pairs(m9, pars="beta")
rstan::traceplot(m9, c("beta"), inc_warmup=FALSE)

print(m9, pars="gamma")

# Explore Correlation structure
corrMat  <-  matrix(m9.summ[630:1029,2], ncol=20,nrow=20)
corrplot(corrMat , method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  corrMat[corrMat == 1]  = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=10.5)
abline(h=10.5)


##  Most of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m9, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m9, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
rstan::plot(m9, pars="tau_run")
rstan::traceplot(m9,pars="tau_run", inc_warmup=FALSE)

#########################################
# LOO Log-likelihood for model selection

m9LL     <-  extract_log_lik(m9, parameter_name = "log_lik")
m9Loo    <-  loo(m9LL)
m9WAIC   <-  waic(m9LL)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m9.df[1,1871:1990])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m9.df[i,1871:1990])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m9.df[,1991], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m9.df[,1992]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m9.df[,1993]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m9.df[,1994]), xlim=c(min(m9.df[,1994],sd(data$nFert)),max(m9.df[,1994],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m9, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data9   <-  as.numeric(m9.df[,2359])
X2sim9    <-  as.numeric(m9.df[,2360])






##########################################################################
# Model: m9a
##########################################################################

##############
# Diagnostics

# Model Results
# print(m9a)
print(m9a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m9a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m9a.df    <-  as.data.frame(extract(m9a))
m9a.summ  <-  plyr:::adply(as.matrix(m9a.df),2,MCMCsum)[-1,]
# m9a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m9a)

# Simple Diagnostic Plots
plot(m9a, pars="beta")
pairs(m9a, pars="gamma")
rstan::traceplot(m9a, c("beta"), inc_warmup=FALSE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m9aLL     <-  extract_log_lik(m9a, parameter_name = "log_lik")
m9aLoo    <-  loo(m9aLL)
m9aWAIC   <-  waic(m9aLL)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m9a.df[1,271:390])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m9a.df[i,271:390])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m9a.df[,391], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m9a.df[,392]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m9a.df[,393]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m9a.df[,394]), xlim=c(min(m9a.df[,394],sd(data$nFert)),max(m9a.df[,394],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m9a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data9a   <-  as.numeric(m9a.df[,759])
X2sim9a    <-  as.numeric(m9a.df[,760])




######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################




##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

##  Overall comparison of all models
looDiff   <-  compare(m1Loo, m2Loo, m3Loo, m3aLoo, m4Loo, m4aLoo, m5Loo,
#                      m6Loo, m6aLoo, m7Loo, m7aLoo, m8Loo, m8aLoo, m9Loo, m9aLoo)
                      m6aLoo, m7aLoo, m8aLoo, m9aLoo)
#waicDiff  <-  compare(m1WAIC, m2WAIC, m3WAIC, m3aWAIC, m4WAIC, m4aWAIC, m5WAIC)

print(looDiff, digits=4)
#print(waicDiff, digits=4)

# LOO Results Summary Table
LooDiff  <-  makeLooTable(looDiff)
(LooDiff)

###########################################################################
###########################################################################
## Main result from LOO model comparison:
##
##  Model m3a:  Random intercept & slopes x Run w/out covariance matrix
##               is the most parsimonious model for this analysis.
##
##  While model m2 provides a slightly better fit to the data over model
##  m3a, this does not appear to be of statistical significance (LooDiff
##  pValue = 0.676). 
##
##  Modelling the covariance structure for model m3a (m3a vs. m3) has  
##  almost no effect on the overall fit (LooDiff pValue = 0.289), and model
##  3a, in fact, fits slightly better . 
##
##  These results are reflected in the Chi-squared discrepancy posterior
##  predictive checks, where models 3a & 3 are almost perfectly overlapping. 
##  Models m1 and m2 appear to marginally improve the X^2 discrepancy 
##  compared to model 3a, but given the MUCH lower model complexity for m3a,
##  combined with the LOO results, model m3a seems the most appropriate
##  choice for the final model upon which to base our inference.
##
###########################################################################
###########################################################################




#  Plot of Chi-squared discrepancy for all models

# Plot Range for PPC plots
X2data1    <-  as.numeric(m1.df[,6379])
X2sim1     <-  as.numeric(m1.df[,6380])
plotRange  <-  c(0,max(X2data1,X2data5))

# pdf(file='./output/figs/NxRate_X2Discrepancy.pdf', width=7,height=7)
par(omi=rep(0.3, 4))
plot(X2sim5 ~ X2data5, 
    xlab=expression(paste(chi^2~discrepancy~of~observed~data)), ylab=expression(paste(chi^2~discrepancy~of~simulated~data)), 
    main=expression(paste(Posterior~predictive~check:~chi^2~"discrepancy")),
    type='n', axes=FALSE, xlim=plotRange, ylim=plotRange, xpd=NA)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2sim5 ~ X2data5, pch=21, 
        bg=transparentColor(COLS[4], 0.1),
        col=transparentColor(COLS[4], 0.3), cex=1.1)
points(X2sim4a ~ X2data4a, pch=21, 
        bg=transparentColor(COLS[6], 0.1),
        col=transparentColor(COLS[6], 0.3), cex=1.1)
points(X2sim4 ~ X2data4, pch=21, 
        bg=transparentColor(COLS[5], 0.1),
        col=transparentColor(COLS[5], 0.3), cex=1.1)
points(X2sim3 ~ X2data3, pch=21, 
        bg=transparentColor(COLS[3], 0.1),
        col=transparentColor(COLS[3], 0.3), cex=1.1)
points(X2sim3a ~ X2data3a, pch=21, 
        bg=transparentColor(COLS[7], 0.1),
        col=transparentColor(COLS[7], 0.3), cex=1.1)
points(X2sim2 ~ X2data2, pch=21, 
        bg=transparentColor(COLS[2], 0.1),
        col=transparentColor(COLS[2], 0.3), cex=1.1)
points(X2sim1 ~ X2data1, pch=21, 
        bg=transparentColor(COLS[1], 0.1),
        col=transparentColor(COLS[1], 0.3), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.17,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1)),
                      expression(paste(Model~2)),
                      expression(paste(Model~3)),
                      expression(paste(Model~"3a")),
                      expression(paste(Model~4)),
                      expression(paste(Model~"4a")),
                      expression(paste(Model~5))),
          pch     =  21,
          pt.bg   =  c(transparentColor(COLS[1],0.5), 
                       transparentColor(COLS[2],0.5),
                       transparentColor(COLS[3],0.5),
                       transparentColor(COLS[7],0.5),
                       transparentColor(COLS[5],0.5),
                       transparentColor(COLS[6],0.5),
                       transparentColor(COLS[4],0.5)),
          col     =  c(transparentColor(COLS[1], 0.7),
                       transparentColor(COLS[2], 0.7), 
                       transparentColor(COLS[3], 0.7), 
                       transparentColor(COLS[7], 0.7), 
                       transparentColor(COLS[5], 0.7), 
                       transparentColor(COLS[6], 0.7), 
                       transparentColor(COLS[4], 0.7)),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
# dev.off()




######################################################
##  Have another look at Posterior Predictive Checks
######################################################

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m3.df[,1991], adjust=2), lwd=3, col=COLS[3], main='min_y_rep (min. num. Successes)')
lines(density(m3a.df[,391], adjust=2), lwd=3, col=COLS[7])
abline(v=min(data$nFert), lwd=3, col=2)
    legend(
          x       =  8,
          y       =  0.35,
          legend  =  c(
                      expression(paste(Model~3)),
                      expression(paste(Model~"3a"))),
          pch     =  21,
          pt.bg   =  c(transparentColor(COLS[3],0.5),
                       transparentColor(COLS[7],0.5)),
          col     =  c(transparentColor(COLS[3], 0.7),
                       transparentColor(COLS[7], 0.7)),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )

plot(density(m3.df[,1992]), lwd=3, col=COLS[3], main='max_y_rep (max. num. Successes)')
lines(density(m3a.df[,392]), lwd=3, col=COLS[7])
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1993]), lwd=3, col=COLS[3], main='mean_y_rep (mean num. Successes)')
lines(density(m3a.df[,393]), lwd=3, col=COLS[7])
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1994]), xlim=c(min(m3.df[,1994],m3a.df[,394],sd(data$nFert)),max(m3.df[,1994],m3a.df[,394],sd(data$nFert))),
 lwd=3, col=COLS[3], main='sd_y_rep (sd num. Successes)')
lines(density(m3a.df[,394]), lwd=3, col=COLS[7], main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


print(m3, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));



######################################################
## Have a look at the final regression from m3 & m3a

##  Calculate Predicted Lines
m3.coef   <-  m3.summ$Mean[621:628]
m3a.coef  <-  m3a.summ$Mean[1:8]

m3Fast5   <-  inv_logit(m3.coef[1] + m3.coef[2] * data$nSperm_z)
m3Fast55  <-  inv_logit((m3.coef[1] + m3.coef[4]) + (m3.coef[2] + m3.coef[6]) * data$nSperm_z)
m3Slow5   <-  inv_logit((m3.coef[1] + m3.coef[3]) + (m3.coef[2] + m3.coef[5]) * data$nSperm_z)
m3Slow55  <-  inv_logit((m3.coef[1] + m3.coef[3] + m3.coef[7]) + (m3.coef[2] + m3.coef[5] + m3.coef[8]) * data$nSperm_z)
m3Fast    <-  inv_logit((m3.coef[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * data$nSperm_z)
m3Slow    <-  inv_logit((m3.coef[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * data$nSperm_z)

m3aFast5   <-  inv_logit(m3a.coef[1] + m3a.coef[2] * data$nSperm_z)
m3aFast55  <-  inv_logit((m3a.coef[1] + m3a.coef[4]) + (m3a.coef[2] + m3a.coef[6]) * data$nSperm_z)
m3aSlow5   <-  inv_logit((m3a.coef[1] + m3a.coef[3]) + (m3a.coef[2] + m3a.coef[5]) * data$nSperm_z)
m3aSlow55  <-  inv_logit((m3a.coef[1] + m3a.coef[3] + m3a.coef[7]) + (m3a.coef[2] + m3a.coef[5] + m3a.coef[8]) * data$nSperm_z)
m3aFast    <-  inv_logit((m3a.coef[1] + (m3a.coef[4])/2) + (m3a.coef[2] + (m3a.coef[6])/2) * data$nSperm_z)
m3aSlow    <-  inv_logit((m3a.coef[1] + m3a.coef[3] + (0.5*(m3a.coef[7]))) + (m3a.coef[2] + m3a.coef[5] + (0.5*(m3a.coef[8]))) * data$nSperm_z)

m3.low   <-  m3.summ$lower[621:628]
m3.hi   <-  m3.summ$upper[621:628]
m3Fast.low    <-  inv_logit((m3.low[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * data$nSperm_z)
m3Slow.low    <-  inv_logit((m3.low[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * data$nSperm_z)
m3Fast.hi     <-  inv_logit((m3.hi[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * data$nSperm_z)
m3Slow.hi     <-  inv_logit((m3.hi[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * data$nSperm_z)






# pdf(file='output/xRatexEggPos_m3.pdf', height=7, width=7)
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m3Fast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3Fast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3Slow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3Slow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
# dev.off()

dev.new()
# par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m3aFast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3aFast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3aSlow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3aSlow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )


#########################################
##  Compare residual plots for m3 & m3a

##  look at residuals for m3
#m3yhat        <-  inv_logit(m3.summ$Mean[1030:1149]) # mus
m3yhat        <-  (m3.summ$Mean[1999:2118])
m3.resids     <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3yhat)/sd(m3yhat)
m3.resids_z   <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3yhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m3yhat))
#m3ayhat       <-  inv_logit(m3a.summ$Mean[32:151]) # mus
m3ayhat       <-  (m3a.summ$Mean[403:522])
m3a.resids    <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3ayhat)/sd(m3ayhat)
m3a.resids_z  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3ayhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m3ayhat))

plot(m3yhat ~ m3ayhat, xlim=c(0,1), ylim=c(0,1))
abline(a = 0, b = 1)

##  Model 3 Residual Plots
par(mfrow=c(2,2))
hist(m3.resids_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(m3.resids_z ~ data$nSperm_z)
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m3.resids_z ~ seq_along(m3.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m3.resids_z)
qqline(m3.resids_z, col = 2)

dev.new()
par(mfrow=c(2,2))
hist(m3a.resids_z, breaks=40)
plot(m3a.resids_z ~ data$nSperm_z)
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m3a.resids_z ~ seq_along(m3a.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m3a.resids_z)
qqline(m3a.resids_z, col = 2)

# hist(m3.resids, breaks=40)
# plot(m3.resids ~ data$nSperm_z)
# abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
# plot(m3.resids ~ seq_along(m3.resids))
# abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
# qqnorm(m3.resids)
# qqline(m3.resids, col = 2)

# hist(m3a.resids, breaks=40)
# plot(m3a.resids ~ data$nSperm_z)
# abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
# plot(m3a.resids ~ seq_along(m3a.resids))
# abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
# qqnorm(m3a.resids)
# qqline(m3a.resids, col = 2)




##  All of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0, 
##  and yet further support for using model m3a
##  for the analysis.
print(m3, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);
m3Corrs  <-  m3.summ$Mean[grep('corrs', names(m3.df))]
plot(density(m3Corrs[m3Corrs < 0.2 & m3Corrs > -0.2]), xlim=c(-0.2,0.2))



Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:data$nSperm_z)
Znames  <-  dimnames(Z)[[2]]


m3a.allBetas   <-  as.matrix(m3a.df[2:9])
m3a.allGammas  <-  as.matrix(m3a.df[10:29])


##  Calculate Coefficients

b0Fast    <-  inv_logit((m3a.allBetas[,1] + (m3a.allBetas[,4])/2))
b0Slow    <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + (0.5*(m3a.allBetas[,7]))))
b0Fast5   <-  inv_logit(m3a.allBetas[,1])
b0Fast55  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,4]))
b0Slow5   <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3]))
b0Slow55  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + m3a.allBetas[,7]))
b1Fast    <-  inv_logit((m3a.allBetas[,2] + (m3a.allBetas[,6])/2))
b1Slow    <-  inv_logit((m3a.allBetas[,2] + m3a.allBetas[,5] + (0.5*(m3a.allBetas[,8]))))
b1Fast5   <-  inv_logit(m3a.allBetas[,2])
b1Fast55  <-  inv_logit((m3a.allBetas[,2] + m3a.allBetas[,6]))
b1Slow5   <-  inv_logit((m3a.allBetas[,2] + m3a.allBetas[,5]))
b1Slow55  <-  inv_logit((m3a.allBetas[,2] + m3a.allBetas[,5] + m3a.allBetas[,8]))

pval  <-  function(x) length(x[x < 0])/length(x)
plotContr  <-  function(Dens, name="title") {
  plot(NA, xlab=expression(paste(Delta)), type='n', axes=FALSE, ylab='Density', cex.lab=1.2, xlim=c(min(Dens$x), (max(Dens$x)+0.4*(max(Dens$x) - min(Dens$x)))), ylim=c(0, (max(Dens$y)+0.05*(max(Dens$y) - min(Dens$y)))), yaxs='i')
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  polygon(c(Dens$x), c(Dens$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
  abline(v=0, lwd=2,col=2)
  axis(1, cex.axis=0.9)
  axis(2, cex.axis=0.9, las=1)
}

c1   <-  b0Slow - b0Fast
c2   <-  b1Slow - b1Fast
c3   <-  b1Slow 
c3b  <-  b1Fast 
c4   <-  b0Fast5 - b0Fast55
c5   <-  b1Fast5 - b1Fast55
c6   <-  b0Slow5 - b0Slow55
c7   <-  b1Slow5 - b1Slow55
c8   <-  b0Fast5 - b0Slow5
c9   <-  b1Fast5 - b1Slow5
c10  <-  b0Fast55 - b0Slow55
c11  <-  b1Fast55 - b1Slow55
c12  <-  b0Fast55 - b0Slow
c13  <-  b1Fast55 - b1Slow
c14  <-  b0Fast5 - b0Slow
c15  <-  b1Fast5 - b1Slow

pval(c1)     # b0Slow - b0Fast
pval(c2)     # b1Slow - b1Fast
pval(c3)     # b1Slow - 0
pval(c3b)    # b1Fast - 0
pval(c4)     # b0Fast5 - b0Fast55
pval(c5)     # b1Fast5 - b1Fast55
pval(c6)     # b0Slow5 - b0Slow55
pval(c7)     # b1Slow5 - b1Slow55
pval(c8)     # b0Fast5 - b0Slow5
pval(c9)     # b1Fast5 - b1Slow5
pval(c10)    # b0Fast55 - b0Slow55
pval(c11)    # b1Fast55 - b1Slow55
pval(c12)    # b0Fast55 - b0Slow
pval(c13)    # b1Fast55 - b1Slow
pval(c14)    # b0Fast5 - b0Slow
pval(c15)    # b1Fast5 - b1Slow

par(mfrow=c(4,4))
plotContr(density(c1))
plotContr(density(c2))
plotContr(density(c3))
plotContr(density(c3b))
plotContr(density(c4))
plotContr(density(c5))
plotContr(density(c6))
plotContr(density(c7))
plotContr(density(c8))
plotContr(density(c9))
plotContr(density(c10))
plotContr(density(c11))

plotContr(density(c12))
plotContr(density(c13))
plotContr(density(c14))
plotContr(density(c15))


m3aFast5.neg2   <-  inv_logit(m3a.allBetas[,1] + m3a.allBetas[,2] * (-2))
m3aFast55.neg2  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,4]) + (m3a.allBetas[,2] + m3a.allBetas[,6]) * (-2))
m3aSlow5.neg2   <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3]) + (m3a.allBetas[,2] + m3a.allBetas[,5]) * (-2))
m3aSlow55.neg2  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + m3a.allBetas[,7]) + (m3a.allBetas[,2] + m3a.allBetas[,5] + m3a.allBetas[,8]) * (-2))
m3aFast.neg2    <-  inv_logit((m3a.allBetas[,1] + (m3a.allBetas[,4])/2) + (m3a.allBetas[,2] + (m3a.allBetas[,6])/2) * (-2))
m3aSlow.neg2    <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + (0.5*(m3a.allBetas[,7]))) + (m3a.allBetas[,2] + m3a.allBetas[,5] + (0.5*(m3a.allBetas[,8]))) * (-2))

m3aFast5.neg1   <-  inv_logit(m3a.allBetas[,1] + m3a.allBetas[,2] * (-1))
m3aFast55.neg1  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,4]) + (m3a.allBetas[,2] + m3a.allBetas[,6]) * (-1))
m3aSlow5.neg1   <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3]) + (m3a.allBetas[,2] + m3a.allBetas[,5]) * (-1))
m3aSlow55.neg1  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + m3a.allBetas[,7]) + (m3a.allBetas[,2] + m3a.allBetas[,5] + m3a.allBetas[,8]) * (-1))
m3aFast.neg1    <-  inv_logit((m3a.allBetas[,1] + (m3a.allBetas[,4])/2) + (m3a.allBetas[,2] + (m3a.allBetas[,6])/2) * (-1))
m3aSlow.neg1    <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + (0.5*(m3a.allBetas[,7]))) + (m3a.allBetas[,2] + m3a.allBetas[,5] + (0.5*(m3a.allBetas[,8]))) * (-1))

m3aFast5.0   <-  inv_logit(m3a.allBetas[,1] + m3a.allBetas[,2] * (0))
m3aFast55.0  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,4]) + (m3a.allBetas[,2] + m3a.allBetas[,6]) * (0))
m3aSlow5.0   <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3]) + (m3a.allBetas[,2] + m3a.allBetas[,5]) * (0))
m3aSlow55.0  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + m3a.allBetas[,7]) + (m3a.allBetas[,2] + m3a.allBetas[,5] + m3a.allBetas[,8]) * (0))
m3aFast.0    <-  inv_logit((m3a.allBetas[,1] + (m3a.allBetas[,4])/2) + (m3a.allBetas[,2] + (m3a.allBetas[,6])/2) * (0))
m3aSlow.0    <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + (0.5*(m3a.allBetas[,7]))) + (m3a.allBetas[,2] + m3a.allBetas[,5] + (0.5*(m3a.allBetas[,8]))) * (0))

m3aFast5.1   <-  inv_logit(m3a.allBetas[,1] + m3a.allBetas[,2] * (1))
m3aFast55.1  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,4]) + (m3a.allBetas[,2] + m3a.allBetas[,6]) * (1))
m3aSlow5.1   <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3]) + (m3a.allBetas[,2] + m3a.allBetas[,5]) * (1))
m3aSlow55.1  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + m3a.allBetas[,7]) + (m3a.allBetas[,2] + m3a.allBetas[,5] + m3a.allBetas[,8]) * (1))
m3aFast.1    <-  inv_logit((m3a.allBetas[,1] + (m3a.allBetas[,4])/2) + (m3a.allBetas[,2] + (m3a.allBetas[,6])/2) * (1))
m3aSlow.1    <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + (0.5*(m3a.allBetas[,7]))) + (m3a.allBetas[,2] + m3a.allBetas[,5] + (0.5*(m3a.allBetas[,8]))) * (1))

m3aFast5.2   <-  inv_logit(m3a.allBetas[,1] + m3a.allBetas[,2] * (2))
m3aFast55.2  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,4]) + (m3a.allBetas[,2] + m3a.allBetas[,6]) * (2))
m3aSlow5.2   <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3]) + (m3a.allBetas[,2] + m3a.allBetas[,5]) * (2))
m3aSlow55.2  <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + m3a.allBetas[,7]) + (m3a.allBetas[,2] + m3a.allBetas[,5] + m3a.allBetas[,8]) * (2))
m3aFast.2    <-  inv_logit((m3a.allBetas[,1] + (m3a.allBetas[,4])/2) + (m3a.allBetas[,2] + (m3a.allBetas[,6])/2) * (2))
m3aSlow.2    <-  inv_logit((m3a.allBetas[,1] + m3a.allBetas[,3] + (0.5*(m3a.allBetas[,7]))) + (m3a.allBetas[,2] + m3a.allBetas[,5] + (0.5*(m3a.allBetas[,8]))) * (2))


simpContr  <-  list(
  cSimp1   =  m3aFast5.neg2 - m3aFast55.neg2,
  cSimp2   =  m3aSlow5.neg2 - m3aSlow55.neg2,
  cSimp3   =  m3aFast5.neg2 - m3aSlow5.neg2,
  cSimp4   =  m3aFast55.neg2 - m3aSlow55.neg2,
  cSimp5   =  m3aFast.neg2 - m3aSlow.neg2,
  cSimp6   =  m3aFast5.neg1 - m3aFast55.neg1,
  cSimp7   =  m3aSlow5.neg1 - m3aSlow55.neg1,
  cSimp8   =  m3aFast5.neg1 - m3aSlow5.neg1,
  cSimp9   =  m3aFast55.neg1 - m3aSlow55.neg1,
  cSimp10  =  m3aFast.neg1 - m3aSlow.neg1,
  cSimp11  =  m3aFast5.0 - m3aFast55.0,
  cSimp12  =  m3aSlow5.0 - m3aSlow55.0,
  cSimp13  =  m3aFast5.0 - m3aSlow5.0,
  cSimp14  =  m3aFast55.0 - m3aSlow55.0,
  cSimp15  =  m3aFast.0 - m3aSlow.0,
  cSimp16  =  m3aFast5.1 - m3aFast55.1,
  cSimp17  =  m3aSlow5.1 - m3aSlow55.1,
  cSimp18  =  m3aFast5.1 - m3aSlow5.1,
  cSimp19  =  m3aFast55.1 - m3aSlow55.1,
  cSimp20  =  m3aFast.1 - m3aSlow.1,
  cSimp21  =  m3aFast5.2 - m3aFast55.2,
  cSimp22  =  m3aSlow5.2 - m3aSlow55.2,
  cSimp23  =  m3aFast5.2 - m3aSlow5.2,
  cSimp24  =  m3aFast55.2 - m3aSlow55.2,
  cSimp25  =  m3aFast.2 - m3aSlow.2
)

pdf(file="./output/contrast_histograms2.pdf", height=18, width=20)
par(mfrow=c(4,5), omi=rep(0.4,4))
# plotContr(density(simpContr[[1]]))
# proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Fast.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
# proportionalLabel(-0.25, 0.5, expression(paste('x = -2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
# plotContr(density(simpContr[[2]]))
# proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow.5 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
# plotContr(density(simpContr[[3]]))
# proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Slow.5')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
# plotContr(density(simpContr[[4]]))
# proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.55 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
# plotContr(density(simpContr[[5]]))
# proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow - Fast')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[6]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Fast.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  proportionalLabel(-0.25, 0.5, expression(paste('x = -1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[7]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow.5 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[8]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Slow.5')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[9]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.55 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[10]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast - Slow')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
for (i in 11:25) {
  plotContr(density(simpContr[[i]]))
#  if(i == 6)
#    proportionalLabel(-0.25, 0.5, expression(paste('x = -1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  if(i == 11)
    proportionalLabel(-0.25, 0.5, expression(paste('x = 0',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  if(i == 16)
    proportionalLabel(-0.25, 0.5, expression(paste('x = 1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  if(i == 21)
    proportionalLabel(-0.25, 0.5, expression(paste('x = 2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
}
dev.off()


plotContr(density(simpContr[[1]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Fast.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
proportionalLabel(-0.25, 0.5, expression(paste('x = -2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[2]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow.5 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[3]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Slow.5')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[4]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.55 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[5]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast - Slow')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)




###########################################################################
###########################################################################
## Main result from a posteriori contrasts:
##
##  There is a significant Rate x EggPos interaction for the Fast release
##  rate intercept (remember, data is Z-transformed, so intercept refers
##  to mean of nSperm_z, the center of the data). The slopes are not
##  significantly different. 
##
##  There is no corresponding significant nSperm_z x EggPos interaction for 
##  the Slow treatment. 
## 
##  There is a significant Rate effect, but this is only for the slopes,
##  not the intercepts. Accordingly, the simple contrasts between the 
##  predicted lines for Fast v. Slow treatments become significant only at
##  sigma >= 1.
##
##  Given these results, my inclination is to plot the main effect for th
##  Slow treatment (pooled regression line), but both 5cm and 55cm lines
##  for the Fast treatment... with different intercepts, but overall slope.
##  This should simplify the figure, and focus attention on the trend that 
##  is driving the overall result: The higher intercept for the 55cm Fast
##  regression line.
##
###########################################################################
###########################################################################


###########################################################################
###########################################################################
#  Plot of run-adjusted y values for the main regression plot

# Fixed effects Model Matrix
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:data$nSperm_z)
Znames  <-  dimnames(Z)[[2]]


m3a.betas    <-  m3a.summ$Mean[1:8]
m3a.gammas   <-  m3a.summ$Mean[9:28]

# Back calculations for each fixed-effect regression coefficients
# b0Fast    <-  inv_logit((m3a.betas[1] + (m3a.betas[4])/2))
# b0Slow    <-  inv_logit((m3a.betas[1] + m3a.betas[3] + (0.5*(m3a.betas[7]))))
# b0Fast5   <-  inv_logit(m3a.betas[1])
# b0Fast55  <-  inv_logit((m3a.betas[1] + m3a.betas[4]))
# b0Slow5   <-  inv_logit((m3a.betas[1] + m3a.betas[3]))
# b0Slow55  <-  inv_logit((m3a.betas[1] + m3a.betas[3] + m3a.betas[7]))
# b1Fast    <-  inv_logit((m3a.betas[2] + (m3a.betas[6])/2))
# b1Slow    <-  inv_logit((m3a.betas[2] + m3a.betas[5] + (0.5*(m3a.betas[8]))))
# b1Fast5   <-  inv_logit(m3a.betas[2])
# b1Fast55  <-  inv_logit((m3a.betas[2] + m3a.betas[6]))
# b1Slow5   <-  inv_logit((m3a.betas[2] + m3a.betas[5]))
# b1Slow55  <-  inv_logit((m3a.betas[2] + m3a.betas[5] + m3a.betas[8]))



# Predicted lines for Fast (with pooled slopes) & the overall Slow

# m3aFast5   <-  inv_logit(m3a.betas[1] + 
#                          (m3a.betas[2] * data$nSperm_z))
# m3aFast55  <-  inv_logit((m3a.betas[1] + m3a.betas[4]) +
#                          (m3a.betas[2] + (m3a.betas[6])) * data$nSperm_z)
m3aFast5   <-  inv_logit((m3a.betas[1] + 
                         (m3a.betas[2] + (m3a.betas[6])/2) * data$nSperm_z))
m3aFast55  <-  inv_logit((m3a.betas[1] + m3a.betas[4]) +
                         (m3a.betas[2] + ((m3a.betas[6])/2)) * data$nSperm_z)
m3aSlow    <-  inv_logit((m3a.betas[1] + m3a.betas[3] + (0.5*(m3a.betas[7]))) + 
                         (m3a.betas[2] + m3a.betas[5] + (0.5*(m3a.betas[8]))) * data$nSperm_z)


# Residuals for fixed-effect predicted lines
ys  <-  ((data$nFert - data$nControlFert)/data$nEggs) 
m3aFast5_Resids   <-  m3aFast5  - ys 
m3aFast55_Resids  <-  m3aFast55 - ys 
m3aSlow_Resids    <-  m3aSlow   - ys 

par(mfrow=c(2,2))
plot(density(m3aFast5_Resids))
abline(v=0,lwd=3,col=2)
plot(density(m3aFast55_Resids))
abline(v=0,lwd=3,col=2)
plot(density(m3aSlow_Resids))
abline(v=0,lwd=3,col=2)

sum(m3aFast5_Resids > 0)/length(m3aFast5_Resids)
sum(m3aFast55_Resids > 0)/length(m3aFast55_Resids)
sum(m3aSlow_Resids > 0)/length(m3aSlow_Resids)

# Adjusted y-values for each predicted line
m3aFast5_yAdj   <-  m3aFast5  + inv_logit((m3a.betas[1] + (m3a.betas[2] + (m3a.betas[6])/2) * data$nSperm_z) + Z %*% m3a.gammas) - ys 
m3aFast55_yAdj  <-  m3aFast55 + inv_logit((m3a.betas[1] + m3a.betas[4]) + (m3a.betas[2] + ((m3a.betas[6])/2)) * data$nSperm_z + Z %*% m3a.gammas) - ys
m3aSlow_yAdj    <-  m3aSlow   + inv_logit((m3a.betas[1] + m3a.betas[3] + (0.5*(m3a.betas[7]))) + 
                                          ((m3a.betas[2] + m3a.betas[5] + (0.5*(m3a.betas[8]))) * data$nSperm_z) + Z %*% m3a.gammas) - ys 



# Visually inspect, make sure they make sense
par(mfrow=c(2,2))
plot(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ m3aFast5[data$Rate == "Fast" & data$EggPos == "5"])
abline(a=0, b=1,lwd=3,col=2)
plot(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ m3aFast55[data$Rate == "Fast" & data$EggPos == "55"])
abline(a=0, b=1,lwd=3,col=2)
plot(m3aSlow_yAdj[data$Rate == "Slow"] ~ m3aSlow[data$Rate == "Slow"])
abline(a=0, b=1,lwd=3,col=2)

par(mfrow=c(2,2))
plot(density(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] - m3aFast5[data$Rate == "Fast" & data$EggPos == "5"]))
plot(density(m3aFast5_Resids))
abline(v=0,lwd=3,col=2)
plot(density(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] - m3aFast55[data$Rate == "Fast" & data$EggPos == "55"]))
abline(v=0,lwd=3,col=2)
plot(density(m3aSlow_yAdj[data$Rate == "Slow"] - m3aSlow[data$Rate == "Slow"]))
abline(v=0,lwd=3,col=2)

sum(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] - m3aFast5[data$Rate == "Fast" & data$EggPos == "5"] > 0)/length(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] - m3aFast5[data$Rate == "Fast" & data$EggPos == "5"])
sum(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] - m3aFast55[data$Rate == "Fast" & data$EggPos == "55"] > 0)/length(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] - m3aFast55[data$Rate == "Fast" & data$EggPos == "55"])
sum(m3aSlow_yAdj[data$Rate == "Slow"] - m3aSlow[data$Rate == "Slow"] > 0)/length(m3aSlow_yAdj[data$Rate == "Slow"] - m3aSlow[data$Rate == "Slow"])

# Plot 

# pdf(file='output/NxRatexEggPos_m3_yAdjusted.pdf', height=7, width=7)
par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='', ylab='', 
    type='n', axes=FALSE, ylim=c(-0.05,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m3aSlow[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3aFast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3aFast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
points(m3aSlow_yAdj[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
proportionalLabel(-0.15, 0.5, expression(paste("Adjusted Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast:~5~cm)),
                      expression(paste(Fast:~55~cm)),
                      expression(paste(Slow))),
          pch     =  c(21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
# dev.off()


#  NOW ADDING CONFIDENCE INTERVALS, TO SEE WHAT THIS LOOKS LIKE.

# Create plotting objects for each regression line
m3aFast5.plt   <-  m3aFast5.plots(m3a.betas, m3a.allBetas, data)
m3aFast55.plt  <-  m3aFast55.plots(m3a.betas, m3a.allBetas, data)
m3aSlow.plt    <-  m3aSlow.plots(m3a.betas, m3a.allBetas, data)


par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='', ylab='', 
    type='n', axes=FALSE, ylim=c(-0.05,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
polygon(x=c(m3aSlow.plt$xRaw, rev(m3aSlow.plt$xRaw)), 
        y=c(m3aSlow.plt$CIs$lower, rev(m3aSlow.plt$CIs$upper)), 
        col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
polygon(x=c(m3aFast5.plt$xRaw, rev(m3aFast5.plt$xRaw)), 
        y=c(m3aFast5.plt$CIs$lower, rev(m3aFast5.plt$CIs$upper)), 
        col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
polygon(x=c(m3aFast55.plt$xRaw, rev(m3aFast55.plt$xRaw)), 
        y=c(m3aFast55.plt$CIs$lower, rev(m3aFast55.plt$CIs$upper)), 
        col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
lines(m3aSlow.plt$y ~ m3aSlow.plt$xRaw, col='orangered1', lwd=3)
lines(m3aFast5.plt$y ~ m3aFast5.plt$xRaw, col='dodgerblue1', lwd=3)
lines(m3aFast55.plt$y ~ m3aFast55.plt$xRaw, col='dodgerblue1', lwd=3, lty=2)
points(m3aSlow_yAdj[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
proportionalLabel(-0.15, 0.5, expression(paste("Adjusted Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast:~5~cm)),
                      expression(paste(Fast:~55~cm)),
                      expression(paste(Slow))),
          pch     =  c(21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )









