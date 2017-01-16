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
source('R/dependencies.R')

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
# print(m1, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m1.df    <-  as.data.frame(extract(m1))  # Run at beginning of script to establish PPC plot ranges
# mcmc.m1  <-  as.mcmc(m1)
m1.mcmc  <-  rstan:::as.mcmc.list.stanfit(m1)
m1.summ  <-  plyr:::adply(as.matrix(m1.df),2,MCMCsum)[-1,]

# Explore Correlation structure
corrMat  <-  matrix(m1.summ[2051:3650,2], ncol=16,nrow=16)
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

# Simple Diagnostic Plots
plot(m1, pars="beta")
pairs(m1, pars="beta")
par(mfrow=c(5,5))
plot(m1.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(m1.mcmc, ask=TRUE)
#plot(density(m1.df[]))

#########################################
# LOO Log-likelihood for model selection

m1LL  <-  extract_log_lik(m1, parameter_name = "log_lik")
m1Loo    <-  loo(m1LL)
m1WAIC   <-  waic(m1LL)


########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




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
# mcmc.m2  <-  as.mcmc(m2)
m2.mcmc  <-  rstan:::as.mcmc.list.stanfit(m2)
m2.summ  <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)[-1,]

# Explore Correlation structure
corrMat  <-  matrix(m2.summ[1241:2140,2], ncol=16,nrow=16)
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

# Simple Diagnostic Plots
plot(m2, pars="beta")
pairs(m2, pars="beta")
par(mfrow=c(5,5))
plot(m2.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(m2.mcmc, c("beta"), ask=TRUE)


#########################################
# LOO Log-likelihood for model selection

m2LL  <-  extract_log_lik(m2, parameter_name = "log_lik")
m2Loo    <-  loo(m2LL)
m2WAIC   <-  waic(m2LL)


########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




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
# mcmc.m3  <-  as.mcmc(m3)
m3.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3)
m3.summ  <-  plyr:::adply(as.matrix(m3.df),2,MCMCsum)[-1,]

# Explore Correlation structure
corrMat  <-  matrix(m3.summ[631:1030,2], ncol=16,nrow=16)
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

# Simple Diagnostic Plots
plot(m3, pars="beta")
pairs(m3, pars="beta")
par(mfrow=c(5,5))
plot(m3.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(m3.mcmc, c("beta"), ask=TRUE)


#########################################
# LOO Log-likelihood for model selection

m3LL     <-  extract_log_lik(m3, parameter_name = "log_lik")
m3Loo    <-  loo(m3LL)
m3WAIC   <-  waic(m3LL)


########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




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
print(m3a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m3a, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m3a.df    <-  as.data.frame(extract(m3a))
# mcmc.m3a  <-  as.mcmc(m3a)
m3a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3a)
m3a.summ  <-  plyr:::adply(as.matrix(m3a.df),2,MCMCsum)[-1,]


# Simple Diagnostic Plots
plot(m3a, pars="gamma0")
pairs(m3a, pars="gamma1")
par(mfrow=c(5,5))
plot(m3a.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(m3a.mcmc, c("beta"), ask=TRUE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m3aLL  <-  extract_log_lik(m3a, parameter_name = "log_lik")
m3aLoo    <-  loo(m3aLL)
m3aWAIC   <-  waic(m3aLL)


########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m3a.df[1,146:265])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(m3a.df[i,146:265])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m3a.df[,266], adjust=4), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m3a.df[,267]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3a.df[,268]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3a.df[,269]), xlim=c(min(m3a.df[,269],sd(data$nFert)),max(m3a.df[,269],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m3a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data3a   <-  as.numeric(m3a.df[,634])
X2sim3a    <-  as.numeric(m3a.df[,635])







##########################################################################
# Model: m4
##########################################################################

##############
# Diagnostics

# Model Results
# print(m3a)
print(m4, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
# print(m4, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m4.df    <-  as.data.frame(extract(m4))
# mcmc.m4  <-  as.mcmc(m4)
m4.mcmc  <-  rstan:::as.mcmc.list.stanfit(m4)
m4.summ  <-  plyr:::adply(as.matrix(m4.df),2,MCMCsum)[-1,]


# Simple Diagnostic Plots
plot(m4, pars="beta")
pairs(m4, pars="beta")
par(mfrow=c(5,5))
plot(m4.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(m4.mcmc, c("beta"), ask=TRUE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m4LL  <-  extract_log_lik(m4, parameter_name = "log_lik")
m4Loo    <-  loo(m4LL)
m4WAIC   <-  waic(m4LL)


########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




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
# mcmc.m4a  <-  as.mcmc(m4a)
m4a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m4a)
m4a.summ  <-  plyr:::adply(as.matrix(m4a.df),2,MCMCsum)[-1,]


# Simple Diagnostic Plots
plot(m4a, pars="beta")
pairs(m4a, pars="beta")
par(mfrow=c(5,5))
plot(m4a.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(m4a.mcmc, c("beta"), ask=TRUE)

dev.off()

#########################################
# LOO Log-likelihood for model selection

m4aLL  <-  extract_log_lik(m4a, parameter_name = "log_lik")
m4aLoo    <-  loo(m4aLL)
m4aWAIC   <-  waic(m4aLL)


########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




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
# mcmc.m5  <-  as.mcmc(m5)
m5.mcmc  <-  rstan:::as.mcmc.list.stanfit(m5)
m5.summ  <-  plyr:::adply(as.matrix(m5.df),2,MCMCsum)[-1,]


#########################################
# LOO Log-likelihood for model selection

m5LL  <-  extract_log_lik(m5, parameter_name = "log_lik")
m5Loo    <-  loo(m5LL)
m5WAIC   <-  waic(m5LL)



########################
# Plot of main results
## !!!!!!!!!!!!!! STILL NEED TO DO THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##




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





###################################################
###################################################
#  Plot of Chi-squared discrepancy for all models

# Plot Range for PPC plots
X2data1   <-  as.numeric(m1.df[,6379])
X2sim1    <-  as.numeric(m1.df[,6380])
plotRange  <-  c(0,max(X2data1,X2data5))

#pdf(file='output/figs/NxRate_X2Discrepancy.pdf', width=7,height=7)
plot(X2sim5 ~ X2data5, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=plotRange, ylim=plotRange)
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
        bg=transparentColor(COLS[7], 0.1),
        col=transparentColor(COLS[7], 0.3), cex=1.1)
points(X2sim3a ~ X2data3a, pch=21, 
        bg=transparentColor(COLS[5], 0.1),
        col=transparentColor(COLS[5], 0.3), cex=1.1)
points(X2sim3 ~ X2data3, pch=21, 
        bg=transparentColor(COLS[3], 0.1),
        col=transparentColor(COLS[3], 0.3), cex=1.1)
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
          x       =  usr[2]*0.15,
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
                       transparentColor(COLS[5],0.5),
                       transparentColor(COLS[7],0.5),
                       transparentColor(COLS[6],0.5),
                       transparentColor(COLS[4],0.5)),
          col     =  c(transparentColor(COLS[1], 0.7),
                       transparentColor(COLS[2], 0.7), 
                       transparentColor(COLS[3], 0.7), 
                       transparentColor(COLS[5], 0.7), 
                       transparentColor(COLS[7], 0.7), 
                       transparentColor(COLS[6], 0.7), 
                       transparentColor(COLS[7], 0.7)),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
#dev.off()








##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

##  Overall comparison of all models
looDiff   <-  compare(m1Loo, m2Loo, m3Loo, m3aLoo, m4Loo, m4aLoo, m5Loo)
waicDiff  <-  compare(m1WAIC, m2WAIC, m3WAIC, m3aWAIC, m4WAIC, m4aWAIC, m5WAIC)

print(looDiff, digits=4)
print(waicDiff, digits=4)


str(looDiff)
looDiff21    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][2]
looDiff23    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][3]
looDiff23a   <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][4]
looDiff24a   <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][5]
looDiff24    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][6]
looDiff25    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][7]
looDiff13    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][3]
looDiff13a   <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][4]
looDiff14a   <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][5]
looDiff14    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][6]
looDiff15    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][7]
looDiff33a   <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][4]
looDiff34a   <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][5]
looDiff34    <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][6]
looDiff35    <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][7]
looDiff3a4a  <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][5]
looDiff3a4   <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][6]
looDiff3a5   <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][7]
looDiff4a4   <-  looDiff[,'elpd_loo'][5] - looDiff[,'elpd_loo'][6]
looDiff4a5   <-  looDiff[,'elpd_loo'][5] - looDiff[,'elpd_loo'][7]
looDiff45    <-  looDiff[,'elpd_loo'][6] - looDiff[,'elpd_loo'][7]


n  <-  length(m1Loo$pointwise[,"elpd_loo"])
selooDiff21    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"] - m1Loo$pointwise[,"elpd_loo"]))
selooDiff23    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"] - m3Loo$pointwise[,"elpd_loo"]))  
selooDiff23a   <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"] - m3aLoo$pointwise[,"elpd_loo"]))  
selooDiff24a   <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"] - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff24    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"] - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff25    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff13    <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"] - m3Loo$pointwise[,"elpd_loo"]))  
selooDiff13a   <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"] - m3aLoo$pointwise[,"elpd_loo"]))  
selooDiff14a   <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"] - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff14    <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"] - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff15    <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff33a   <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"] - m3aLoo$pointwise[,"elpd_loo"]))  
selooDiff34a   <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"] - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff34    <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"] - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff35    <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff3a4a  <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff3a4   <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff3a5   <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff4a4   <-  sqrt(n * var(m4aLoo$pointwise[,"elpd_loo"] - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff4a5   <-  sqrt(n * var(m4aLoo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff45    <-  sqrt(n * var(m4Loo$pointwise[,"elpd_loo"]  - m5Loo$pointwise[,"elpd_loo"]))  


LooDiff  <-  cbind(c(looDiff21,looDiff23,looDiff23a,looDiff24a,looDiff24,looDiff25,looDiff13,looDiff13a,looDiff14a,looDiff14,looDiff15,looDiff33a,looDiff34a,looDiff34,looDiff35,looDiff3a4a,looDiff3a4,looDiff3a5,looDiff4a4,looDiff4a5,looDiff45),
                  c(selooDiff21,selooDiff23,selooDiff23a,selooDiff24a,selooDiff24,selooDiff25,selooDiff13,selooDiff13a,selooDiff14a,selooDiff14,selooDiff15,selooDiff33a,selooDiff34a,selooDiff34,selooDiff35,selooDiff3a4a,selooDiff3a4,selooDiff3a5,selooDiff4a4,selooDiff4a5,selooDiff45))


pDiff21    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[1,1] - 0)/LooDiff[1,2])), 3))
pDiff23    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[2,1] - 0)/LooDiff[2,2])), 3))
pDiff23a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[3,1] - 0)/LooDiff[3,2])), 3))
pDiff24a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[4,1] - 0)/LooDiff[4,2])), 3))
pDiff24    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[5,1] - 0)/LooDiff[5,2])), 3))
pDiff25    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[6,1] - 0)/LooDiff[6,2])), 3))
pDiff13    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[7,1] - 0)/LooDiff[7,2])), 3))
pDiff13a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[8,1] - 0)/LooDiff[8,2])), 3))
pDiff14a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[9,1] - 0)/LooDiff[9,2])), 3))
pDiff14    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[10,1] - 0)/LooDiff[10,2])), 3))
pDiff15    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[11,1] - 0)/LooDiff[11,2])), 3))
pDiff33a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[12,1] - 0)/LooDiff[12,2])), 3))
pDiff34a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[13,1] - 0)/LooDiff[13,2])), 3))
pDiff34    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[14,1] - 0)/LooDiff[14,2])), 3))
pDiff35    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[15,1] - 0)/LooDiff[15,2])), 3))
pDiff3a4a  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[16,1] - 0)/LooDiff[16,2])), 3))
pDiff3a4   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[17,1] - 0)/LooDiff[17,2])), 3))
pDiff3a5   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[18,1] - 0)/LooDiff[18,2])), 3))
pDiff4a4   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[19,1] - 0)/LooDiff[19,2])), 3))
pDiff4a5   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[20,1] - 0)/LooDiff[20,2])), 3))
pDiff45    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[21,1] - 0)/LooDiff[21,2])), 3))

LooDiff  <-  cbind(LooDiff,c(pDiff21,pDiff23,pDiff23a,pDiff24a,pDiff24,pDiff25,pDiff13,pDiff13a,pDiff14a,pDiff14,pDiff15,pDiff33a,pDiff34a,pDiff34,pDiff35,pDiff3a4a,pDiff3a4,pDiff3a5,pDiff4a4,pDiff4a5,pDiff45))

row.names(LooDiff)  <-  c('m2 - m1',
                          'm2 - m3',
                          'm2 - m3a',
                          'm2 - m4a',
                          'm2 - m4',
                          'm2 - m5',
                          'm1 - m3',
                          'm1 - m3a',
                          'm1 - m4a',
                          'm1 - m4',
                          'm1 - m5',
                          'm3 - m3a',
                          'm3 - m4a',
                          'm3 - m4',
                          'm3 - m5',
                          'm3a - m4a',
                          'm3a - m4',
                          'm3a - m5',
                          'm4a - m4',
                          'm4a - m5',
                          'm4 - m5')
colnames(LooDiff)   <-  c("diff", "se", "p.value")
LooDiff


###########################################################################
###########################################################################
## Main result from LOO model comparison:
##
##  Model m3:  Random intercept & slopes x Run w/ COVARIANCE MATRIX model
##             is probablyh the most parsimonious model for this analysis.
##
##  Modelling the covariance structure for this model (m3 vs. m3a) only 
##  marginally improves the fit according to LOO (LooDiff pValue = 0.111). 
##  However, based on the Chi-squared discrepancy posterior predictive 
##  checks, there is a reasonably large gap between these models, with m3
##  having lower discrepancy. Thus, modelling the covariance structure here
##  seems to improve model fit pretty well, despite the marginal difference
##  in LOO.
###########################################################################
###########################################################################


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))

plot(density(m3.df[,1991], adjust=1), lwd=3, col=COLS[3], main='min_y_rep (min. num. Successes)')
lines(density(m3a.df[,266], adjust=1), lwd=3, col=COLS[5])
abline(v=min(data$nFert), lwd=3, col=2)
    legend(
          x       =  7.5,
          y       =  0.75,
          legend  =  c(
                      expression(paste(Model~3)),
                      expression(paste(Model~"3a"))),
          pch     =  21,
          pt.bg   =  c(transparentColor(COLS[3],0.5),
                       transparentColor(COLS[5],0.5)),
          col     =  c(transparentColor(COLS[3], 0.7),
                       transparentColor(COLS[5], 0.7)),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )

plot(density(m3.df[,1992]), lwd=3, col=COLS[3], main='max_y_rep (max. num. Successes)')
lines(density(m3a.df[,267]), lwd=3, col=COLS[5])
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1993]), lwd=3, col=COLS[3], main='mean_y_rep (mean num. Successes)')
lines(density(m3a.df[,268]), lwd=3, col=COLS[5])
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3.df[,1994]), xlim=c(min(m3.df[,1994],m3a.df[,269],sd(data$nFert)),max(m3.df[,1994],m3a.df[,269],sd(data$nFert))),
 lwd=3, col=COLS[3], main='sd_y_rep (sd num. Successes)')
lines(density(m3a.df[,269]), lwd=3, col=COLS[5], main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


print(m3, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3a, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));




par(mfrow=c(2,2))
plot(density(m3.df[2501:7500,1991], adjust=1), lwd=3, col=COLS[3], main='min_y_rep (min. num. Successes)')
lines(density(m3a.df[2501:7500,266], adjust=1), lwd=3, col=COLS[5])
abline(v=min(data$nFert), lwd=3, col=2)
    legend(
          x       =  7.5,
          y       =  0.75,
          legend  =  c(
                      expression(paste(Model~3)),
                      expression(paste(Model~"3a"))),
          pch     =  21,
          pt.bg   =  c(transparentColor(COLS[3],0.5),
                       transparentColor(COLS[5],0.5)),
          col     =  c(transparentColor(COLS[3], 0.7),
                       transparentColor(COLS[5], 0.7)),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )

plot(density(m3.df[2501:7500,1992]), lwd=3, col=COLS[3], main='max_y_rep (max. num. Successes)')
lines(density(m3a.df[2501:7500,267]), lwd=3, col=COLS[5])
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3.df[2501:7500,1993]), lwd=3, col=COLS[3], main='mean_y_rep (mean num. Successes)')
lines(density(m3a.df[2501:7500,268]), lwd=3, col=COLS[5])
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3.df[2501:7500,1994]), xlim=c(min(m3.df[2501:7500,1994],m3a.df[2501:7500,269],sd(data$nFert)),max(m3.df[2501:7500,1994],m3a.df[2501:7500,269],sd(data$nFert))),
 lwd=3, col=COLS[3], main='sd_y_rep (sd num. Successes)')
lines(density(m3a.df[2501:7500,269]), lwd=3, col=COLS[5], main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

