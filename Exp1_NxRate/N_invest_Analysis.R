#/* 
# * Colin Olito. Created 13/12/2016
# * Analysis of 1st flume experiment: N-invest
# * 
# * NOTES:  This file will fit all the necessary
# * 		Stan models for the analysis of the
# * 		N-invest flume data, and write the
# * 		Stan sample files to ./output/Stanfits
# * 
# * 		These Stan sample files can then be 
# * 		imported, plotted, and analyzed using
# * 		./N_invest_Analysis
# *          
# */

rm(list=ls())
#################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

###############
# DEPENDENCIES
source('R/dependencies.R')

##############
# Import Data
data <- read.csv('data/Ninvest_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)
# head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Lane      <-  factor(data$Lane)
data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
data$Date      <-  dmy(data$Date)
data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)



#####################################
# Import model results from the stan 
# sample_files for further analysis 
#####################################

csvFiles  <-  c('./output/StanFits/N_invest_m1.csv1',
                './output/StanFits/N_invest_m1.csv2',
                './output/StanFits/N_invest_m1.csv3',
                './output/StanFits/N_invest_m1.csv4')
m1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m2.csv1',
                './output/StanFits/N_invest_m2.csv2',
                './output/StanFits/N_invest_m2.csv3',
                './output/StanFits/N_invest_m2.csv4')
m2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m2b.csv1',
                './output/StanFits/N_invest_m2b.csv2',
                './output/StanFits/N_invest_m2b.csv3',
                './output/StanFits/N_invest_m2b.csv4')
m2b        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m3.csv1',
                './output/StanFits/N_invest_m3.csv2',
                './output/StanFits/N_invest_m3.csv3',
                './output/StanFits/N_invest_m3.csv4')
m3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m4.csv1',
                './output/StanFits/N_invest_m4.csv2',
                './output/StanFits/N_invest_m4.csv3',
                './output/StanFits/N_invest_m4.csv4')
m4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)









##########################################################################
# Model: m1
##########################################################################

##############
# Diagnostics

# Model Results
print(m1)
print(m1, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m1, c("yRep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m1.df    <-  as.data.frame(extract(m1))
mcmc.m1  <-  as.mcmc(m1)
m1.mcmc  <-  rstan:::as.mcmc.list.stanfit(m1)
m1.summ  <-  plyr:::adply(as.matrix(m1.df),2,MCMCsum)[-1,]
(m1.summ)

# Simple Diagnostic Plots
plot(m1, pars="beta")
plot(m1.mcmc, ask=TRUE)
pairs(m1, pars="beta")

#########################################
# LOO Log-likelihood for model selection

m1LL  <-  extract_log_lik(m1, parameter_name = "log_lik")
m1Loo    <-  loo(m1LL)
m1WAIC   <-  waic(m1LL)



########################
# Plot of main results

#  Plot predicted line etc.
RegLine  <-  inv_logit(m1.summ$Mean[1] + m1.summ$Mean[2] * data$nSperm_z)


##  Plot showing uncertainty from Run-specific intercepts
par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ data$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(m1.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
 }, data=data, nSperm_z=data$nSperm_z)
# plot main regression line
lines(RegLine[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue4', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m1.df[1,52:99])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

#for(i in 2:(nrow(m1.df) - 1)) {
for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m1.df[i,52:99])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m1.summ
par(mfrow=c(2,2))
plot(density(m1.df[,100], adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m1.df[,101]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m1.df[,102]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m1.df[,103]), xlim=c(min(m1.df[,103],sd(data$nFert)),max(m1.df[,103],sd(data$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x1  <-  as.numeric(m1.df[,252])
y1  <-  as.numeric(m1.df[,253])

plot(y1 ~ x1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(c(x1,y1))), ylim=c(0,max(c(x1,y1))))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y1 ~ x1, pch=21, 
        bg=transparentColor('dodgerblue4', 0.2),
        col=transparentColor('dodgerblue4', 0.4), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)



##########################################################################
# Model: m2
##########################################################################

##############
# Diagnostics

# Model Results
print(m2)
print(m2, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m2, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m2.df    <-  as.data.frame(extract(m2))
mcmc.m2  <-  as.mcmc(m2)
m2.mcmc  <-  rstan:::as.mcmc.list.stanfit(m2)
m2.summ  <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)[-1,]
(m2.summ)

# Simple Diagnostic Plots
plot(m2, pars="beta")
plot(m2, pars="gamma")
plot(m2.mcmc, ask=TRUE)
pairs(m2, pars="beta")
pairs(m2, pars="gamma")


#########################################
# LOO Log-likelihood for model selection

m2LL  <-  extract_log_lik(m2, parameter_name = "log_lik")
m2Loo    <-  loo(m2LL)
m2WAIC   <-  waic(m2LL)


########################
# Plot of main results


##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[3])   + m2.summ$Mean[2] * data$nSperm_z),
               Run2  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[4])   + m2.summ$Mean[2] * data$nSperm_z),
               Run3  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[5])   + m2.summ$Mean[2] * data$nSperm_z),
               Run4  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[6])   + m2.summ$Mean[2] * data$nSperm_z),
               Run5  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[7])   + m2.summ$Mean[2] * data$nSperm_z),
               Run6  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[8])   + m2.summ$Mean[2] * data$nSperm_z),
               Run7  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[9])   + m2.summ$Mean[2] * data$nSperm_z),
               Run8  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[10])  + m2.summ$Mean[2] * data$nSperm_z)
              )

RegLine  <-  inv_logit(m2.summ$Mean[1] + m2.summ$Mean[2] * data$nSperm_z)



#pdf(file="./output/N_runVar.pdf", height=7, width=7)
#par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = data)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m2.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][data$Run == i][order(data$nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(data$nSperm_z[data$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)
#dev.off()




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m2.df[1,61:108])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m2.df[i,61:108])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m2.summ
par(mfrow=c(2,2))
plot(density(m2.df[,109], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m2.df[,110]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m2.df[,111]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m2.df[,112]), xlim=c(min(m2.df[,112],sd(data$nFert)),max(m2.df[,112],sd(data$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x2  <-  as.numeric(m2.df[,261])
y2  <-  as.numeric(m2.df[,262])

plot(y1 ~ x1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(c(x1,y1))), ylim=c(0,max(c(x1,y1))))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y1 ~ x1, pch=21, 
        bg=transparentColor('dodgerblue4', 0.1),
        col=transparentColor('dodgerblue4', 0.3), cex=1.1)
points(y2 ~ x2, pch=21, 
        bg=transparentColor('dodgerblue3', 0.1),
        col=transparentColor('dodgerblue3', 0.4), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.15,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1)),
                      expression(paste(Model~2))),
          pch     =  21,
          pt.bg   =  c(transparentColor('dodgerblue4',0.7), transparentColor('dodgerblue3',0.7)),
          col     =  c('dodgerblue4', 'dodgerblue3'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )


##########################################################################
# Model: m2b
##########################################################################

##############
# Diagnostics

# Model Results
print(m2b)
print(m2b, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m2b, c("gamma", "mu_gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m2b.df    <-  as.data.frame(extract(m2b))
mcmc.m2b  <-  as.mcmc(m2b)
m2b.mcmc  <-  rstan:::as.mcmc.list.stanfit(m2b)
m2b.summ  <-  plyr:::adply(as.matrix(m2b.df),2,MCMCsum)[-1,]
(m2b.summ)

# Simple Diagnostic Plots
plot(m2b, pars="beta")
plot(m2b, pars="gamma")
plot(m2b.mcmc, ask=TRUE)
pairs(m2b, pars="beta")
pairs(m2b, pars="gamma")


#########################################
# LOO Log-likelihood for model selection

m2bLL  <-  extract_log_lik(m2b, parameter_name = "log_lik")
m2bLoo    <-  loo(m2bLL)
m2bWAIC   <-  waic(m2bLL)


########################
# Plot of main results

##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(m2b.summ$Mean[2] + m2b.summ$Mean[1] * data$nSperm_z),
               Run2  <- inv_logit(m2b.summ$Mean[3] + m2b.summ$Mean[1] * data$nSperm_z),
               Run3  <- inv_logit(m2b.summ$Mean[4] + m2b.summ$Mean[1] * data$nSperm_z),
               Run4  <- inv_logit(m2b.summ$Mean[5] + m2b.summ$Mean[1] * data$nSperm_z),
               Run5  <- inv_logit(m2b.summ$Mean[6] + m2b.summ$Mean[1] * data$nSperm_z),
               Run6  <- inv_logit(m2b.summ$Mean[7] + m2b.summ$Mean[1] * data$nSperm_z),
               Run7  <- inv_logit(m2b.summ$Mean[8] + m2b.summ$Mean[1] * data$nSperm_z),
               Run8  <- inv_logit(m2b.summ$Mean[9] + m2b.summ$Mean[1] * data$nSperm_z)
              )

RegLine  <-  inv_logit(m2b.summ$Mean[10] + m2b.summ$Mean[1] * data$nSperm_z)



par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = data)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m2b.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma'] + x['beta'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][data$Run == i][order(data$nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(data$nSperm_z[data$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m2b.df[1,61:108])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m2b.df[i,61:108])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m2b.summ
par(mfrow=c(2,2))
plot(density(m2b.df[,109], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m2b.df[,110]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m2b.df[,111]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m2b.df[,112]), xlim=c(min(m2.df[,112],sd(data$nFert)),max(m2.df[,112],sd(data$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x2b  <-  as.numeric(m2b.df[,261])
y2b  <-  as.numeric(m2b.df[,262])

plot(y1 ~ x1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(c(x1,y1))), ylim=c(0,max(c(x1,y1))))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y1 ~ x1, pch=21, 
        bg=transparentColor('dodgerblue4', 0.1),
        col=transparentColor('dodgerblue4', 0.3), cex=1.1)
points(y2 ~ x2, pch=21, 
        bg=transparentColor('dodgerblue3', 0.1),
        col=transparentColor('dodgerblue3', 0.4), cex=1.1)
points(y2b ~ x2b, pch=21, 
        bg=transparentColor('tomato2', 0.1),
        col=transparentColor('tomato2', 0.4), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.15,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1)),
                      expression(paste(Model~2)),
                      expression(paste(Model~"2b"))),
          pch     =  21,
          pt.bg   =  c(transparentColor('dodgerblue4',0.7), 
                       transparentColor('dodgerblue3',0.7),
                       transparentColor('tomato2',0.7)),
          col     =  c('dodgerblue4', 
                       'dodgerblue3',
                       'tomato2'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##########################################################################
# Model: m3
##########################################################################

##############
# Diagnostics

# Model Results
print(m3)
print(m3, c("gamma0", "mu_gamma0", "sigma_gamma0"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma1", "mu_gamma1", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m3.df    <-  as.data.frame(extract(m3))
mcmc.m3  <-  as.mcmc(m3)
m3.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3)
m3.summ  <-  plyr:::adply(as.matrix(m3.df),2,MCMCsum)[-1,]
(m3.summ)

# Simple Diagnostic Plots
plot(m3, pars="gamma0")
plot(m3, pars="gamma1")
plot(m3.mcmc, ask=TRUE)
pairs(m3, pars="gamma0")
pairs(m3, pars="gamma1")


#########################################
# LOO Log-likelihood for model selection

m3LL  <-  extract_log_lik(m3, parameter_name = "log_lik")
m3Loo    <-  loo(m3LL)
m3WAIC   <-  waic(m3LL)

########################
# Plot of main results

##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(m3.summ$Mean[1] + m3.summ$Mean[9]  * data$nSperm_z),
               Run2  <- inv_logit(m3.summ$Mean[2] + m3.summ$Mean[10] * data$nSperm_z),
               Run3  <- inv_logit(m3.summ$Mean[3] + m3.summ$Mean[11] * data$nSperm_z),
               Run4  <- inv_logit(m3.summ$Mean[4] + m3.summ$Mean[12] * data$nSperm_z),
               Run5  <- inv_logit(m3.summ$Mean[5] + m3.summ$Mean[13] * data$nSperm_z),
               Run6  <- inv_logit(m3.summ$Mean[6] + m3.summ$Mean[14] * data$nSperm_z),
               Run7  <- inv_logit(m3.summ$Mean[7] + m3.summ$Mean[15] * data$nSperm_z),
               Run8  <- inv_logit(m3.summ$Mean[8] + m3.summ$Mean[16] * data$nSperm_z)
              )

RegLine  <-  inv_logit(m3.summ$Mean[17] + m3.summ$Mean[18] * data$nSperm_z)



par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = data)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m3.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma0'] + x['mu_gamma1'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=data$nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][data$Run == i][order(data$nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(data$nSperm_z[data$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m3.df[1,70:117])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m3.df[i,70:117])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m3.summ
par(mfrow=c(2,2))
plot(density(m3.df[,118],adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m3.df[,119]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m3.df[,120]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m3.df[,121]), xlim=c(min(m2.df[,112],sd(data$nFert)),max(m2.df[,112],sd(data$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x3   <-  as.numeric(m3.df[,270])
y3   <-  as.numeric(m3.df[,271])

plot(y1 ~ x1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(c(x1,y1))), ylim=c(0,max(c(x1,y1))))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y1 ~ x1, pch=21, 
        bg=transparentColor('dodgerblue4', 0.1),
        col=transparentColor('dodgerblue4', 0.3), cex=1.1)
points(y2 ~ x2, pch=21, 
        bg=transparentColor('dodgerblue3', 0.1),
        col=transparentColor('dodgerblue3', 0.4), cex=1.1)
points(y2b ~ x2b, pch=21, 
        bg=transparentColor('tomato2', 0.1),
        col=transparentColor('tomato2', 0.4), cex=1.1)
points(y3 ~ x3, pch=21, 
        bg=transparentColor('dodgerblue2', 0.1),
        col=transparentColor('dodgerblue2', 0.4), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.15,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1)),
                      expression(paste(Model~2)),
                      expression(paste(Model~"2b")),
                      expression(paste(Model~3))),
          pch     =  21,
          pt.bg   =  c(transparentColor('dodgerblue4',0.7), 
                       transparentColor('dodgerblue3',0.7),
                       transparentColor('tomato2',0.7),
                       transparentColor('dodgerblue2',0.7)),
          col     =  c('dodgerblue4', 
                       'dodgerblue3',
                       'tomato2',
                       'dodgerblue2'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )





##########################################################################
# Model: m4
##########################################################################

##############
# Diagnostics

# Model Results
print(m4, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
m4.df    <-  as.data.frame(extract(m4))
mcmc.m4  <-  as.mcmc(m4)
m4.mcmc  <-  rstan:::as.mcmc.list.stanfit(m4)
m4.summ  <-  plyr:::adply(as.matrix(m4.df),2,MCMCsum)[-1,]
(m4.summ)

# Explore Correlation structure
corrMat  <-  matrix(m4.summ[404:659,2], ncol=16,nrow=16)
corrplot(corrMat , method='circle', type='upper')
abline(v=8.5)
abline(h=8.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  }
}

corrplot(corrMat * 50, method='circle', type='upper')
abline(v=8.5)
abline(h=8.5)

# Simple Diagnostic Plots
plot(m4, pars="beta")
pairs(m4, pars="beta")


#########################################
# LOO Log-likelihood for model selection

m4LL  <-  extract_log_lik(m4, parameter_name = "log_lik")
m4Loo    <-  loo(m4LL)
m4WAIC   <-  waic(m4LL)




########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##

  ##  Plot predicted line etc.
  runs  <-  list(
                 Run1  <- inv_logit(m4.summ$Mean[965] + m4.summ$Mean[1029]  * data$nSperm_z),
                 Run2  <- inv_logit(m4.summ$Mean[974] + m4.summ$Mean[1038] * data$nSperm_z),
                 Run3  <- inv_logit(m4.summ$Mean[983] + m4.summ$Mean[1047] * data$nSperm_z),
                 Run4  <- inv_logit(m4.summ$Mean[992] + m4.summ$Mean[1056] * data$nSperm_z),
                 Run5  <- inv_logit(m4.summ$Mean[1001] + m4.summ$Mean[1065] * data$nSperm_z),
                 Run6  <- inv_logit(m4.summ$Mean[1010] + m4.summ$Mean[1074] * data$nSperm_z),
                 Run7  <- inv_logit(m4.summ$Mean[1019] + m4.summ$Mean[1083] * data$nSperm_z),
                 Run8  <- inv_logit(m4.summ$Mean[1028] + m4.summ$Mean[1092] * data$nSperm_z)
                )

  RegLine  <-  inv_logit(m4.summ$Mean[402] + m4.summ$Mean[403] * data$nSperm_z)



par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = data)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m3.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma0'] + x['mu_gamma1'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][data$Run == i][order(data$nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(data$nSperm_z[data$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)






##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m4.df[1,1141:1188])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m4.df[i,1141:1188])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m4.df[,1189], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m4.df[,1190]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m4.df[,1191]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m4.df[,1192]), xlim=c(min(m2.df[,112],sd(data$nFert)),max(m2.df[,112],sd(data$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x4   <-  as.numeric(m4.df[,1341])
y4   <-  as.numeric(m4.df[,1342])

plot(y1 ~ x1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(c(x1,y1))), ylim=c(0,max(c(x1,y1))))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y1 ~ x1, pch=21, 
        bg=transparentColor('dodgerblue4', 0.1),
        col=transparentColor('dodgerblue4', 0.3), cex=1.1)
points(y2 ~ x2, pch=21, 
        bg=transparentColor('dodgerblue3', 0.1),
        col=transparentColor('dodgerblue3', 0.4), cex=1.1)
points(y2b ~ x2b, pch=21, 
        bg=transparentColor('tomato2', 0.1),
        col=transparentColor('tomato2', 0.4), cex=1.1)
points(y3 ~ x3, pch=21, 
        bg=transparentColor('dodgerblue2', 0.1),
        col=transparentColor('dodgerblue2', 0.4), cex=1.1)
points(y4 ~ x4, pch=21, 
        bg=transparentColor('tomato1', 0.1),
        col=transparentColor('tomato1', 0.4), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.15,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1)),
                      expression(paste(Model~2)),
                      expression(paste(Model~"2b")),
                      expression(paste(Model~3)),
                      expression(paste(Model~4))),
          pch     =  21,
          pt.bg   =  c(transparentColor('dodgerblue4',0.7), 
                       transparentColor('dodgerblue3',0.7),
                       transparentColor('tomato2',0.7),
                       transparentColor('dodgerblue2',0.7),
                       transparentColor('tomato1',0.7)),
          col     =  c('dodgerblue4', 
                       'dodgerblue3',
                       'tomato2',
                       'dodgerblue2',
                       'tomato1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )






#########################################################
##  Final Plot for comparison of Chi-squared discrepancy 
##  posterior predictive check
##  Excludes models 2b & 4 because they provide almost
##  identical fits to models 2 & 3 respectively.

pdf(file="./output/N_invest_X2Discrepancy.pdf", height=7, width=7)
par(omi=rep(0.3, 4))
plot(y1 ~ x1, 
    xlab=expression(paste(Chi^2~discrepancy~of~observed~data)), ylab=expression(paste(Chi^2~discrepancy~of~simulated~data)), 
    type='n', axes=FALSE, xlim=c(0,max(c(x1,y1))), ylim=c(0,max(c(x1,y1))), xpd=NA)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y1 ~ x1, pch=21, 
        bg=transparentColor('dodgerblue1', 0.1),
        col=transparentColor('dodgerblue1', 0.3), cex=1.1)
points(y2 ~ x2, pch=21, 
        bg=transparentColor('dodgerblue2', 0.1),
        col=transparentColor('dodgerblue2', 0.4), cex=1.1)
points(y3 ~ x3, pch=21, 
        bg=transparentColor('dodgerblue4', 0.1),
        col=transparentColor('dodgerblue4', 0.4), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2)
axis(1)
    legend(
          x       =  usr[2]*0.572,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1:~Simple~logistic~regression)),
                      expression(paste(Model~2:~Random~intercepts)),
                      expression(paste(Model~3:~Random~intercept/slope))),
          pch     =  21,
          pt.bg   =  c(transparentColor('dodgerblue1',0.7), 
                       transparentColor('dodgerblue2',0.7),
                       transparentColor('dodgerblue4',0.7)),
          col     =  c('dodgerblue1', 
                       'dodgerblue2',
                       'dodgerblue4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
dev.off()







##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

str(m1Loo)
looDiff   <-  compare(m1Loo, m2Loo, m2bLoo, m3Loo, m4Loo)
waicDiff  <-  compare(m1WAIC, m2WAIC, m2bWAIC, m3WAIC, m4WAIC)

print(looDiff, digits=4)
print(waicDiff, digits=4)

print(compare(m1Loo, m2Loo), digits=6)
print(compare(m1Loo, m2bLoo), digits=6)
print(compare(m1Loo, m3Loo), digits=6)
print(compare(m1Loo, m4Loo), digits=6)
print(compare(m2Loo, m2bLoo), digits=6)
print(compare(m2Loo, m3Loo), digits=6)
print(compare(m2Loo, m4Loo), digits=6)
print(compare(m2bLoo, m3Loo), digits=6)
print(compare(m2bLoo, m4Loo), digits=6)
print(compare(m3Loo, m4Loo), digits=6)


looDiff34   <-  looDiff[1,3] - looDiff[2,3]
looDiff32   <-  looDiff[1,3] - looDiff[3,3]
looDiff32b  <-  looDiff[1,3] - looDiff[4,3]
looDiff31   <-  looDiff[1,3] - looDiff[5,3]
looDiff42   <-  looDiff[2,3] - looDiff[3,3]
looDiff42b  <-  looDiff[2,3] - looDiff[4,3]
looDiff41   <-  looDiff[2,3] - looDiff[5,3]
looDiff22b  <-  looDiff[3,3] - looDiff[4,3]
looDiff21   <-  looDiff[3,3] - looDiff[5,3]
looDiff2b1  <-  looDiff[4,3] - looDiff[5,3]

n  <-  length(m1Loo$pointwise[,"elpd_loo"])
selooDiff34   <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m4Loo$pointwise[,"elpd_loo"]))
selooDiff32   <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m2Loo$pointwise[,"elpd_loo"]))
selooDiff32b  <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m2bLoo$pointwise[,"elpd_loo"]))
selooDiff31   <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m1Loo$pointwise[,"elpd_loo"]))
selooDiff42   <-  sqrt(n * var(m4Loo$pointwise[,"elpd_loo"]  - m2Loo$pointwise[,"elpd_loo"]))
selooDiff42b  <-  sqrt(n * var(m4Loo$pointwise[,"elpd_loo"]  - m2bLoo$pointwise[,"elpd_loo"]))
selooDiff41   <-  sqrt(n * var(m4Loo$pointwise[,"elpd_loo"]  - m1Loo$pointwise[,"elpd_loo"]))
selooDiff22b  <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m2bLoo$pointwise[,"elpd_loo"]))
selooDiff21   <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m1Loo$pointwise[,"elpd_loo"]))
selooDiff2b1  <-  sqrt(n * var(m2bLoo$pointwise[,"elpd_loo"] - m1Loo$pointwise[,"elpd_loo"]))




LooDiff  <-  cbind(c(looDiff34,looDiff32,looDiff32b,looDiff31,looDiff42,looDiff42b,looDiff41,looDiff22b,looDiff21,looDiff2b1),
                   c(selooDiff34,selooDiff32,selooDiff32b,selooDiff31,selooDiff42,selooDiff42b,selooDiff41,selooDiff22b,selooDiff21,selooDiff2b1))

pDiff34   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[1,1] - 0)/LooDiff[1,2])), 3))
pDiff32   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[2,1] - 0)/LooDiff[2,2])), 3))
pDiff32b  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[3,1] - 0)/LooDiff[3,2])), 3))
pDiff31   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[4,1] - 0)/LooDiff[4,2])), 3))
pDiff42   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[5,1] - 0)/LooDiff[5,2])), 3))
pDiff42b  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[6,1] - 0)/LooDiff[6,2])), 3))
pDiff41   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[7,1] - 0)/LooDiff[7,2])), 3))
pDiff22b  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[8,1] - 0)/LooDiff[8,2])), 3))
pDiff21   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[9,1] - 0)/LooDiff[9,2])), 3))
pDiff2b1  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[10,1] - 0)/LooDiff[10,2])),3))

LooDiff  <-  cbind(LooDiff, c(pDiff34,pDiff32,pDiff32b,pDiff31,pDiff42,pDiff42b,pDiff41,pDiff22b,pDiff21,pDiff2b1))

row.names(LooDiff)  <-  c("m3 - m4",
                          "m3 - m2",
                          "m3 - m2b",
                          "m3 - m1",
                          "m4 - m2",
                          "m4 - m2b",
                          "m4 - m1",
                          "m2 - m2b",
                          "m2 - m1",
                          "m2b - m1")
colnames(LooDiff)   <-  c("diff", "se", "p.value")
LooDiff
LooDiff[c(1,2,4,9),]



########################################
## Main result of LOO model comparison:
#
#  Model m3 (random slopes & intercepts) is the best fitting model.
#  This is corroborated by the posterior predictive checks, especially
#  the Chi-square discrepancy graphical check. Modes m4 and m2b gave
#  nearly identical fits to models m3 and m2, and so not discussed
#  further.
#
#  Model comparison using LOO suggests that the overall fit for models
#  m3 and m2 are statistically indistinguishable (LooDiff pValue = 0.618). 
#  We therefore present the results from model m2, the random intercepts
#  model. There could be an arguement for presenting m3 based on the 
#  Chi-square discrepancy plots... but at this point in time... meh.
