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
COLS  <-  c("#860885",
            "#006a09",
            "#028fff",
            "#ffc544",
            "#331061",
            "#ba2423"
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
print(m1)
print(m1, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m1, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m1.df    <-  as.data.frame(extract(m1))
mcmc.m1  <-  as.mcmc(m1)
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


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data1   <-  as.numeric(m1.df[,6379])
X2sim1    <-  as.numeric(m1.df[,6380])

plot(X2sim1 ~ X2data1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(X2data1)), ylim=c(0,max(X2data1)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2sim1 ~ X2data1, pch=21, 
        bg=transparentColor(COLS[1], 0.1),
        col=transparentColor(COLS[1], 0.3), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)





##########################################################################
# Model: m5
##########################################################################

##############
# Diagnostics

# Model Results
print(m5)
print(m5, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m5, c("y_rep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m5.df    <-  as.data.frame(extract(m5))
mcmc.m5  <-  as.mcmc(m5)
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


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

X2data5   <-  as.numeric(m5.df[,618])
X2sim5    <-  as.numeric(m5.df[,619])

plot(X2sim ~ X2data, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(X2data)), ylim=c(0,max(X2data)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2sim5 ~ X2data5, pch=21, 
        bg=transparentColor(COLS[6], 0.1),
        col=transparentColor(COLS[6], 0.3), cex=1.1)
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
                      expression(paste(Model~5))),
          pch     =  21,
          pt.bg   =  c(transparentColor(COLS[1],0.7), 
                       transparentColor(COLS[6],0.7)),
          col     =  c(COLS[1], 
                       COLS[6]),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
