#/* 
# * Colin Olito. Created 13/12/2016
# * Analysis of 1st flume experiment: N-invest
# * 
# * NOTES:  This file contains all the necessary
# * 		code to read in the STAN sample files
# * 		and perform posterior predictive checks
# *     model selection, and produce regression
# *     plots for the N_invest flume data
# * 
# */

rm(list=ls())
#################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

###########################
# DEPENDENCIES & LOAD DATA
source('R/functions.R')
source('R/loadData.R')



##########################################################################
# Model: NIm1
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm1)
print(NIm1, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm1, c("yRep"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
NIm1.df    <-  as.data.frame(extract(NIm1))
mcmc.NIm1  <-  as.mcmc(NIm1)
NIm1.mcmc  <-  rstan:::as.mcmc.list.stanfit(NIm1)
NIm1.summ  <-  plyr:::adply(as.matrix(NIm1.df),2,MCMCsum)[-1,]
(NIm1.summ)

# Simple Diagnostic Plots
plot(NIm1, pars="beta")
plot(NIm1.mcmc, ask=TRUE)
pairs(NIm1, pars="beta")

#########################################
# LOO Log-likelihood for model selection

NIm1LL  <-  extract_log_lik(NIm1, parameter_name = "log_lik")
NIm1Loo    <-  loo(NIm1LL)
NIm1WAIC   <-  waic(NIm1LL)



########################
# Plot of main results

#  Plot predicted line etc.
RegLine  <-  inv_logit(NIm1.summ$Mean[1] + NIm1.summ$Mean[2] * NinvData$nSperm_z)


##  Plot showing uncertainty from Run-specific intercepts
par(omi=rep(0.3, 4))
plot((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(NinvData$nSperm),max(NinvData$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(NIm1.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(NinvData$nSperm), max(NinvData$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
 }, data=data, nSperm_z=NinvData$nSperm_z)
# plot main regression line
lines(RegLine[order(NinvData$nSperm_z)] ~ NinvData$nSperm[order(NinvData$nSperm_z)],
                  col='black', lwd=3)
points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=21, 
        bg=transparentColor('dodgerblue4', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm1.df[1,52:99])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

#for(i in 2:(nrow(NIm1.df) - 1)) {
for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm1.df[i,52:99])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in NIm1.summ
par(mfrow=c(2,2))
plot(density(NIm1.df[,100], adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm1.df[,101]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm1.df[,102]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm1.df[,103]), xlim=c(min(NIm1.df[,103],sd(NinvData$nFert)),max(NIm1.df[,103],sd(NinvData$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x1  <-  as.numeric(NIm1.df[,252])
y1  <-  as.numeric(NIm1.df[,253])

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
# Model: NIm2
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm2)
print(NIm2, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm2, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
NIm2.df    <-  as.data.frame(extract(NIm2))
mcmc.NIm2  <-  as.mcmc(NIm2)
NIm2.mcmc  <-  rstan:::as.mcmc.list.stanfit(NIm2)
NIm2.summ  <-  plyr:::adply(as.matrix(NIm2.df),2,MCMCsum)[-1,]
(NIm2.summ)

# Simple Diagnostic Plots
plot(NIm2, pars="beta")
plot(NIm2, pars="gamma")
plot(NIm2.mcmc, ask=TRUE)
pairs(NIm2, pars="beta")
pairs(NIm2, pars="gamma")


#########################################
# LOO Log-likelihood for model selection

NIm2LL  <-  extract_log_lik(NIm2, parameter_name = "log_lik")
NIm2Loo    <-  loo(NIm2LL)
NIm2WAIC   <-  waic(NIm2LL)


########################
# Plot of main results


##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[3])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run2  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[4])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run3  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[5])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run4  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[6])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run5  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[7])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run6  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[8])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run7  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[9])   + NIm2.summ$Mean[2] * NinvData$nSperm_z),
               Run8  <- inv_logit((NIm2.summ$Mean[1] + NIm2.summ$Mean[10])  + NIm2.summ$Mean[2] * NinvData$nSperm_z)
              )

RegLine  <-  inv_logit(NIm2.summ$Mean[1] + NIm2.summ$Mean[2] * NinvData$nSperm_z)



#pdf(file="./output/N_runVar.pdf", height=7, width=7)
#par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = NinvData)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(NIm2.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(NinvData$nSperm), max(NinvData$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])] ~ NinvData$nSperm[NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(NinvData$nSperm_z)] ~ NinvData$nSperm[order(NinvData$nSperm_z)],
                  col='black', lwd=3)
points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)
#dev.off()




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm2.df[1,61:108])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm2.df[i,61:108])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in NIm2.summ
par(mfrow=c(2,2))
plot(density(NIm2.df[,109], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,110]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,111]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,112]), xlim=c(min(NIm2.df[,112],sd(NinvData$nFert)),max(NIm2.df[,112],sd(NinvData$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x2  <-  as.numeric(NIm2.df[,261])
y2  <-  as.numeric(NIm2.df[,262])

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
# Model: NIm2b
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm2b)
print(NIm2b, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm2b, c("gamma", "mu_gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
NIm2b.df    <-  as.data.frame(extract(NIm2b))
mcmc.NIm2b  <-  as.mcmc(NIm2b)
NIm2b.mcmc  <-  rstan:::as.mcmc.list.stanfit(NIm2b)
NIm2b.summ  <-  plyr:::adply(as.matrix(NIm2b.df),2,MCMCsum)[-1,]
(NIm2b.summ)

# Simple Diagnostic Plots
plot(NIm2b, pars="beta")
plot(NIm2b, pars="gamma")
plot(NIm2b.mcmc, ask=TRUE)
pairs(NIm2b, pars="beta")
pairs(NIm2b, pars="gamma")


#########################################
# LOO Log-likelihood for model selection

NIm2bLL  <-  extract_log_lik(NIm2b, parameter_name = "log_lik")
NIm2bLoo    <-  loo(NIm2bLL)
NIm2bWAIC   <-  waic(NIm2bLL)


########################
# Plot of main results

##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(NIm2b.summ$Mean[2] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run2  <- inv_logit(NIm2b.summ$Mean[3] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run3  <- inv_logit(NIm2b.summ$Mean[4] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run4  <- inv_logit(NIm2b.summ$Mean[5] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run5  <- inv_logit(NIm2b.summ$Mean[6] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run6  <- inv_logit(NIm2b.summ$Mean[7] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run7  <- inv_logit(NIm2b.summ$Mean[8] + NIm2b.summ$Mean[1] * NinvData$nSperm_z),
               Run8  <- inv_logit(NIm2b.summ$Mean[9] + NIm2b.summ$Mean[1] * NinvData$nSperm_z)
              )

RegLine  <-  inv_logit(NIm2b.summ$Mean[10] + NIm2b.summ$Mean[1] * NinvData$nSperm_z)



par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = NinvData)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(NIm2b.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(NinvData$nSperm), max(NinvData$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma'] + x['beta'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])] ~ NinvData$nSperm[NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(NinvData$nSperm_z)] ~ NinvData$nSperm[order(NinvData$nSperm_z)],
                  col='black', lwd=3)
points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm2b.df[1,61:108])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm2b.df[i,61:108])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in NIm2b.summ
par(mfrow=c(2,2))
plot(density(NIm2b.df[,109], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2b.df[,110]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2b.df[,111]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2b.df[,112]), xlim=c(min(NIm2.df[,112],sd(NinvData$nFert)),max(NIm2.df[,112],sd(NinvData$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x2b  <-  as.numeric(NIm2b.df[,261])
y2b  <-  as.numeric(NIm2b.df[,262])

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
# Model: NIm3
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm3)
print(NIm3, c("gamma0", "mu_gamma0", "sigma_gamma0"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm3, c("gamma1", "mu_gamma1", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
NIm3.df    <-  as.data.frame(extract(NIm3))
mcmc.NIm3  <-  as.mcmc(NIm3)
NIm3.mcmc  <-  rstan:::as.mcmc.list.stanfit(NIm3)
NIm3.summ  <-  plyr:::adply(as.matrix(NIm3.df),2,MCMCsum)[-1,]
(NIm3.summ)

# Simple Diagnostic Plots
plot(NIm3, pars="gamma0")
plot(NIm3, pars="gamma1")
plot(NIm3.mcmc, ask=TRUE)
pairs(NIm3, pars="gamma0")
pairs(NIm3, pars="gamma1")


#########################################
# LOO Log-likelihood for model selection

NIm3LL  <-  extract_log_lik(NIm3, parameter_name = "log_lik")
NIm3Loo    <-  loo(NIm3LL)
NIm3WAIC   <-  waic(NIm3LL)

########################
# Plot of main results

##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(NIm3.summ$Mean[1] + NIm3.summ$Mean[9]  * NinvData$nSperm_z),
               Run2  <- inv_logit(NIm3.summ$Mean[2] + NIm3.summ$Mean[10] * NinvData$nSperm_z),
               Run3  <- inv_logit(NIm3.summ$Mean[3] + NIm3.summ$Mean[11] * NinvData$nSperm_z),
               Run4  <- inv_logit(NIm3.summ$Mean[4] + NIm3.summ$Mean[12] * NinvData$nSperm_z),
               Run5  <- inv_logit(NIm3.summ$Mean[5] + NIm3.summ$Mean[13] * NinvData$nSperm_z),
               Run6  <- inv_logit(NIm3.summ$Mean[6] + NIm3.summ$Mean[14] * NinvData$nSperm_z),
               Run7  <- inv_logit(NIm3.summ$Mean[7] + NIm3.summ$Mean[15] * NinvData$nSperm_z),
               Run8  <- inv_logit(NIm3.summ$Mean[8] + NIm3.summ$Mean[16] * NinvData$nSperm_z)
              )

RegLine  <-  inv_logit(NIm3.summ$Mean[17] + NIm3.summ$Mean[18] * NinvData$nSperm_z)



par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = NinvData)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(NIm3.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(NinvData$nSperm), max(NinvData$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma0'] + x['mu_gamma1'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=NinvData$nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])] ~ NinvData$nSperm[NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(NinvData$nSperm_z)] ~ NinvData$nSperm[order(NinvData$nSperm_z)],
                  col='black', lwd=3)
points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)




##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm3.df[1,70:117])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm3.df[i,70:117])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in NIm3.summ
par(mfrow=c(2,2))
plot(density(NIm3.df[,118],adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm3.df[,119]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm3.df[,120]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm3.df[,121]), xlim=c(min(NIm2.df[,112],sd(NinvData$nFert)),max(NIm2.df[,112],sd(NinvData$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x3   <-  as.numeric(NIm3.df[,270])
y3   <-  as.numeric(NIm3.df[,271])

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
# Model: NIm4
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm4, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
NIm4.df    <-  as.data.frame(extract(NIm4))
mcmc.NIm4  <-  as.mcmc(NIm4)
NIm4.mcmc  <-  rstan:::as.mcmc.list.stanfit(NIm4)
NIm4.summ  <-  plyr:::adply(as.matrix(NIm4.df),2,MCMCsum)[-1,]
(NIm4.summ)

# Explore Correlation structure
corrMat  <-  matrix(NIm4.summ[404:659,2], ncol=16,nrow=16)
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
plot(NIm4, pars="beta")
pairs(NIm4, pars="beta")


#########################################
# LOO Log-likelihood for model selection

NIm4LL  <-  extract_log_lik(NIm4, parameter_name = "log_lik")
NIm4Loo    <-  loo(NIm4LL)
NIm4WAIC   <-  waic(NIm4LL)




########################
# Plot of main results 
## !!!!!!!!!!!!!! STILL NEED TO FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!! ##

  ##  Plot predicted line etc.
  runs  <-  list(
                 Run1  <- inv_logit(NIm4.summ$Mean[965] + NIm4.summ$Mean[1029]  * NinvData$nSperm_z),
                 Run2  <- inv_logit(NIm4.summ$Mean[974] + NIm4.summ$Mean[1038] * NinvData$nSperm_z),
                 Run3  <- inv_logit(NIm4.summ$Mean[983] + NIm4.summ$Mean[1047] * NinvData$nSperm_z),
                 Run4  <- inv_logit(NIm4.summ$Mean[992] + NIm4.summ$Mean[1056] * NinvData$nSperm_z),
                 Run5  <- inv_logit(NIm4.summ$Mean[1001] + NIm4.summ$Mean[1065] * NinvData$nSperm_z),
                 Run6  <- inv_logit(NIm4.summ$Mean[1010] + NIm4.summ$Mean[1074] * NinvData$nSperm_z),
                 Run7  <- inv_logit(NIm4.summ$Mean[1019] + NIm4.summ$Mean[1083] * NinvData$nSperm_z),
                 Run8  <- inv_logit(NIm4.summ$Mean[1028] + NIm4.summ$Mean[1092] * NinvData$nSperm_z)
                )

  RegLine  <-  inv_logit(NIm4.summ$Mean[402] + NIm4.summ$Mean[403] * NinvData$nSperm_z)



par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = NinvData)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(NIm3.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(NinvData$nSperm), max(NinvData$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma0'] + x['mu_gamma1'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])] ~ NinvData$nSperm[NinvData$Run == i][order(NinvData$nSperm_z[NinvData$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(NinvData$nSperm_z)] ~ NinvData$nSperm[order(NinvData$nSperm_z)],
                  col='black', lwd=3)
points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)






##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm4.df[1,1141:1188])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm4.df[i,1141:1188])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1) 


# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in NIm4.summ
par(mfrow=c(2,2))
plot(density(NIm4.df[,1189], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. numm Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm4.df[,1190]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm4.df[,1191]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm4.df[,1192]), xlim=c(min(NIm2.df[,112],sd(NinvData$nFert)),max(NIm2.df[,112],sd(NinvData$nFert))), lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)


# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data

x4   <-  as.numeric(NIm4.df[,1341])
y4   <-  as.numeric(NIm4.df[,1342])

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

str(NIm1Loo)
looDiff   <-  compare(NIm1Loo, NIm2Loo, NIm2bLoo, NIm3Loo, NIm4Loo)
waicDiff  <-  compare(NIm1WAIC, NIm2WAIC, NIm2bWAIC, NIm3WAIC, NIm4WAIC)

row.names(looDiff)
"NIm3Loo" "NIm4Loo" "NIm2Loo" "NIm2bLoo" "NIm1Loo"

looList  <-  lapply(row.names(looDiff),get)
names(looList)  <-  row.names(looDiff)

print(looDiff, digits=4)
print(waicDiff, digits=4)

print(compare(NIm1Loo, NIm2Loo), digits=6)
print(compare(NIm1Loo, NIm2bLoo), digits=6)
print(compare(NIm1Loo, NIm3Loo), digits=6)
print(compare(NIm1Loo, NIm4Loo), digits=6)
print(compare(NIm2Loo, NIm2bLoo), digits=6)
print(compare(NIm2Loo, NIm3Loo), digits=6)
print(compare(NIm2Loo, NIm4Loo), digits=6)
print(compare(NIm2bLoo, NIm3Loo), digits=6)
print(compare(NIm2bLoo, NIm4Loo), digits=6)
print(compare(NIm3Loo, NIm4Loo), digits=6)


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

n  <-  length(NIm1Loo$pointwise[,"elpd_loo"])
selooDiff34   <-  sqrt(n * var(NIm3Loo$pointwise[,"elpd_loo"]  - NIm4Loo$pointwise[,"elpd_loo"]))
selooDiff32   <-  sqrt(n * var(NIm3Loo$pointwise[,"elpd_loo"]  - NIm2Loo$pointwise[,"elpd_loo"]))
selooDiff32b  <-  sqrt(n * var(NIm3Loo$pointwise[,"elpd_loo"]  - NIm2bLoo$pointwise[,"elpd_loo"]))
selooDiff31   <-  sqrt(n * var(NIm3Loo$pointwise[,"elpd_loo"]  - NIm1Loo$pointwise[,"elpd_loo"]))
selooDiff42   <-  sqrt(n * var(NIm4Loo$pointwise[,"elpd_loo"]  - NIm2Loo$pointwise[,"elpd_loo"]))
selooDiff42b  <-  sqrt(n * var(NIm4Loo$pointwise[,"elpd_loo"]  - NIm2bLoo$pointwise[,"elpd_loo"]))
selooDiff41   <-  sqrt(n * var(NIm4Loo$pointwise[,"elpd_loo"]  - NIm1Loo$pointwise[,"elpd_loo"]))
selooDiff22b  <-  sqrt(n * var(NIm2Loo$pointwise[,"elpd_loo"]  - NIm2bLoo$pointwise[,"elpd_loo"]))
selooDiff21   <-  sqrt(n * var(NIm2Loo$pointwise[,"elpd_loo"]  - NIm1Loo$pointwise[,"elpd_loo"]))
selooDiff2b1  <-  sqrt(n * var(NIm2bLoo$pointwise[,"elpd_loo"] - NIm1Loo$pointwise[,"elpd_loo"]))




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

row.names(LooDiff)  <-  c("NIm3 - NIm4",
                          "NIm3 - NIm2",
                          "NIm3 - NIm2b",
                          "NIm3 - NIm1",
                          "NIm4 - NIm2",
                          "NIm4 - NIm2b",
                          "NIm4 - NIm1",
                          "NIm2 - NIm2b",
                          "NIm2 - NIm1",
                          "NIm2b - NIm1")
colnames(LooDiff)   <-  c("diff", "se", "p.value")
LooDiff
LooDiff[c(1,2,4,9),]


# or, using new makeLooTable() function:
makeLooTable(looDiff)
########################################
## Main result of LOO model comparison:
#
#  Model NIm3 (random slopes & intercepts) is the best fitting model.
#  This is corroborated by the posterior predictive checks, especially
#  the Chi-square discrepancy graphical check. Modes NIm4 and NIm2b gave
#  nearly identical fits to models NIm3 and NIm2, and so not discussed
#  further.
#
#  Model comparison using LOO suggests that the overall fit for models
#  NIm3 and NIm2 are statistically indistinguishable (LooDiff pValue = 0.618). 
#  We therefore present the results from model NIm2, the random intercepts
#  model. There could be an arguement for presenting NIm3 based on the 
#  Chi-square discrepancy plots... but at this point in time... meh.
