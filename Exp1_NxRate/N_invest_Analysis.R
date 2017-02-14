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
source('./loadData_N_invest.R')




##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

##  Model Comparisons for Binomial Error
looDiff     <-  compare(NIm1Loo, NIm2Loo, NIm3Loo)
looDiffBB   <-  compare(NIm1BBLoo, NIm2BBLoo, NIm3BBLoo)

print(looDiff, digits=4)
print(looDiffBB, digits=4)

# Pairwise LOO Comparisons Summary Table
LooDiff    <-  makeLooTable(looDiff)
LooDiffBB  <-  makeLooTable(looDiffBB)
LooDiff
LooDiffBB

# Write Loo Model Selection Results to .csv
write.csv(LooDiff, file= './output/tables/N_invest_LooDiff_Bin.csv')
write.csv(LooDiffBB, file= './output/tables/N_invest_LooDiff_BetaBin.csv')

###########################################################################
## Main result from LOO model comparison:
##
##  -- The best fitting, and most parsimonious model is NIm2 - with random
##     intercepts for Run. This result is robust to the choice of error
##     distribution.  
##
##  -- Encouragingly, the result is robust to the method of comparison for
##     the pointwise elpd differences. 
##      
##  -- Will focus analysis of posterior predictive checks on models NIm2,
##     NIm3, NIm2BB and NIm3BB... Basically between random intercepts only,
##     and random slopes and intercepts models.
##
###########################################################################




###########################################################################
###########################################################################
#  Plot of Chi-squared discrepancy for all models
###########################################################################
###########################################################################



#####################################################
# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data
#####################################################

# Model NIm3
X2data3      <-  as.numeric(NIm3.df[,316])
X2sim3       <-  as.numeric(NIm3.df[,317])

# Model NIm2
X2data2      <-  as.numeric(NIm2.df[,308])
X2sim2       <-  as.numeric(NIm2.df[,309])

# Minimal model
X2data1      <-  as.numeric(NIm1.df[,251])
X2sim1       <-  as.numeric(NIm1.df[,252])

# Model NIm3BB 
X2data3BB    <-  as.numeric(NIm3BB.df[,365])
X2sim3BB     <-  as.numeric(NIm3BB.df[,366])

# Model NIm2BB (also best fitting model)
X2data2BB    <-  as.numeric(NIm2BB.df[,357])
X2sim2BB     <-  as.numeric(NIm2BB.df[,358])

# Model NIm2BB (also best fitting model)
X2data1BB    <-  as.numeric(NIm1BB.df[,348])
X2sim1BB     <-  as.numeric(NIm1BB.df[,349])

# Plot Range for X^2 plots
plotRange    <-  c(0,max(X2data1))
plotRangeBB  <-  c(0,max(X2sim1BB))

# Plot Colors for X^2 Posterior Predictive Check Plots
COLS  <-  c("#009ce0",
            "#f9684a",
            "#8be36a",
            "#7a0058",
            "#9fa3ff",
            "#be0042",
            "#2d1956"
)


#  Plot of X^2 discrepancies for candidate models
# X^2 Discrepancy plot for Binomial error models
# pdf(file='./output/figs/N_invest_X2_Bin.pdf', width=7,height=7)
par(omi=rep(0.3, 4))
plot(X2sim1 ~ X2data1, 
    xlab=expression(paste(chi^2~discrepancy~of~observed~data)), ylab=expression(paste(chi^2~discrepancy~of~simulated~data)), 
    main=expression(paste(Posterior~predictive~check:~chi^2~"discrepancy")),
    type='n', axes=FALSE, xlim=plotRange, ylim=plotRange, xpd=NA)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2sim1 ~ X2data1, pch=21, 
        bg=transparentColor(COLS[4], 0.1),
        col=transparentColor(COLS[4], 0.3), cex=1.1)
points(X2sim2 ~ X2data2, pch=21, 
        bg=transparentColor(COLS[3], 0.1),
        col=transparentColor(COLS[3], 0.3), cex=1.1)
points(X2sim3 ~ X2data3, pch=21, 
        bg=transparentColor(COLS[7], 0.1),
        col=transparentColor(COLS[7], 0.3), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.25,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~"NIm1")),
                      expression(paste(Model~"NIm2")),
                      expression(paste(Model~"NIm3"))),
          pch     =  21,
          pt.bg   =  c(
                       transparentColor(COLS[4],0.5),
                       transparentColor(COLS[3],0.5),
                       transparentColor(COLS[7],0.5)
                       ),
          col     =  c(
                       transparentColor(COLS[4], 0.7), 
                       transparentColor(COLS[3], 0.5),
                       transparentColor(COLS[7], 0.7)
                       ),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
# dev.off()

#  Plot of X^2 discrepancies for candidate models
# X^2 Discrepancy plot for Binomial error models
# pdf(file='./output/figs/N_invest_X2BetaBin.pdf', width=7,height=7)
par(omi=rep(0.3, 4))
plot(X2sim1BB ~ X2data1BB, 
    xlab=expression(paste(chi^2~discrepancy~of~observed~data)), ylab=expression(paste(chi^2~discrepancy~of~simulated~data)), 
    main=expression(paste(Posterior~predictive~check:~chi^2~"discrepancy")),
    type='n', axes=FALSE, xlim=plotRangeBB, ylim=plotRangeBB, xpd=NA)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2sim1BB ~ X2data1BB, pch=21, 
        bg=transparentColor(COLS[4], 0.1),
        col=transparentColor(COLS[4], 0.3), cex=1.1)
points(X2sim2BB ~ X2data2BB, pch=21, 
        bg=transparentColor(COLS[3], 0.1),
        col=transparentColor(COLS[3], 0.3), cex=1.1)
points(X2sim3BB ~ X2data3BB, pch=21, 
        bg=transparentColor(COLS[7], 0.1),
        col=transparentColor(COLS[7], 0.3), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.25,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~"NIm1BB")),
                      expression(paste(Model~"NIm2BB")),
                      expression(paste(Model~"NIm3BB"))),
          pch     =  21,
          pt.bg   =  c(
                       transparentColor(COLS[4],0.5),
                       transparentColor(COLS[3],0.5),
                       transparentColor(COLS[7],0.5)
                       ),
          col     =  c(
                       transparentColor(COLS[4], 0.7), 
                       transparentColor(COLS[3], 0.5),
                       transparentColor(COLS[7], 0.7)
                       ),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
# dev.off()

# Compare Fixed effects estimates
print(NIm2BB, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm2, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));





###########################################################################
###########################################################################
#  MODEL DIAGNOSTICS AND POSTERIOR PREDICTIVE CHECKS FOR CANDIDATE MODELS
###########################################################################
###########################################################################


##########################################################################
# Model: NIm2
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm2, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm2, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(NIm2, pars="beta")
pairs(NIm2, pars="beta")
plot(NIm2, pars="gamma")
pairs(NIm2, pars="gamma")
par(mfrow=c(5,5))
rstan::traceplot(NIm2, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*10 ~ x, lwd=2, type='l')
lines(density(NIm2.df$sigma_gamma), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot((dnorm(x, sd=3)*75) ~ x, lwd=2, type='l')
lines(density(NIm2.df$beta.2), lwd=3, col='dodgerblue1')


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm2.df[1,108:155])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm2.df[i,108:155])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated NinvData, benchmarked with
#  calculated values for real NinvData
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(NIm2.df[,156], adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,157]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,158]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,159]), xlim=c(min(NIm2.df[,159],sd(NinvData$nFert)),max(NIm2.df[,159],sd(NinvData$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)

print(NIm2, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()








##########################################################################
# Model: NIm2BB
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm2BB, c("beta", "phi"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm2BB, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(NIm2BB, pars="beta")
pairs(NIm2BB, pars="beta")
pairs(NIm2BB, pars="gamma")
rstan::traceplot(NIm2BB, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*10 ~ x, lwd=2, type='l')
lines(density(NIm2BB.df$sigma_gamma), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot((dnorm(x, sd=3)*50) ~ x, lwd=2, type='l')
lines(density(NIm2BB.df$beta.2), lwd=3, col='dodgerblue1')

# Check posteriors against priors
x  <-  seq(from=0, to=100, length=500)
plot(dcauchy(x, scale=100)*10 ~ x, lwd=2, type='l')
lines(density(NIm2BB.df$phi), lwd=3, col='dodgerblue1')

##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm2BB.df[1,205:252])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm2BB.df[i,205:252])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated NinvData, benchmarked with
#  calculated values for real NinvData
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(NIm2BB.df[,253], adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2BB.df[,254]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2BB.df[,255]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2BB.df[,256]), xlim=c(min(NIm2BB.df[,256],sd(NinvData$nFert)),max(NIm2BB.df[,256],sd(NinvData$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)

print(NIm2BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()




##########################################################################
# Model: NIm3BB
##########################################################################

##############
# Diagnostics

# Model Results
print(NIm3BB, c("beta", "phi"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(NIm3BB, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(NIm3BB, pars="beta")
pairs(NIm3BB, pars="beta")
plot(NIm3BB, pars="gamma")
pairs(NIm3BB, pars="gamma")
rstan::traceplot(NIm3BB, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*15 ~ x, lwd=2, type='l')
lines(density(NIm3BB.df$sigma_gamma), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot((dnorm(x, sd=3)*17) ~ x, lwd=2, type='l')
lines(density(NIm3BB.df$beta.2), lwd=3, col='dodgerblue1')

x  <-  seq(from=0, to=100, length=500)
plot(dcauchy(x, scale=100)*7 ~ x, lwd=2, type='l')
lines(density(NIm3BB.df$phi), lwd=3, col='dodgerblue1')


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(NIm3BB.df[1,213:260])/NinvData$nEggs
x  <-  NinvData$nFert/NinvData$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(NIm3BB.df[i,213:260])/NinvData$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated NinvData, benchmarked with
#  calculated values for real NinvData
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(NIm3BB.df[,261], adjust=2), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(NinvData$nFert), lwd=3, col=2)

plot(density(NIm3BB.df[,262]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm3BB.df[,263]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm3BB.df[,264]), xlim=c(min(NIm3BB.df[,264],sd(NinvData$nFert)),max(NIm3BB.df[,264],sd(NinvData$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(NinvData$nFert), lwd=3, col=2)

print(NIm3BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()




###########################################################################
###########################################################################
##  Have a closer look at Posterior Predictive Checks for Candidate Models
###########################################################################
###########################################################################

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
# pdf(file='./output/figs/N_invest_PPDs.pdf', width=7,height=7)
par(mfrow=c(2,2))
plot(density(NIm2.df[,156],   adjust=2), ylim=c(0,0.25), xlim=c(-2,15), lwd=3, col='dodgerBlue1', main='min_y_rep (min. num. Successes)')
lines(density(NIm2BB.df[,253], adjust=2), lwd=3, col='orangered2')
lines(density(NIm3BB.df[,261], adjust=2), lwd=3, col=COLS[4])
abline(v=min(NinvData$nFert), lwd=3, col=2)
    legend(
          x       =  16,
          y       =  0.25,
          legend  =  c(
                      expression(paste(Model~"NIm2")),
                      expression(paste(Model~"NIm2BB")),
                      expression(paste(Model~"NIm3BB"))),
          lty     =  1,
          lwd     =  3,
          col     =  c(
                       'dodgerblue1',
                       'orangered2',
                       COLS[4]),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )

plot(density(NIm2.df[,157]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
lines(density(NIm2BB.df[,254]), lwd=3, col='orangered2')
lines(density(NIm3BB.df[,262]), lwd=3, col=COLS[4])
abline(v=max(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,158]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
lines(density(NIm2BB.df[,255]), lwd=3, col='orangered2')
lines(density(NIm3BB.df[,263]), lwd=3, col=COLS[4])
abline(v=mean(NinvData$nFert), lwd=3, col=2)

plot(density(NIm2.df[,159]), xlim=c(min(NIm2BB.df[,256],sd(NinvData$nFert)),max(NIm2BB.df[,256],sd(NinvData$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
lines(density(NIm2BB.df[,256]), lwd=3, col='orangered2')
lines(density(NIm3BB.df[,264]), lwd=3, col=COLS[4])
abline(v=sd(NinvData$nFert), lwd=3, col=2)

dev.off()


###########################################################################
## Main result from PPD's:
##
##  --  Not much to choose between NIm2 and NIm2BBm. NIm2 may slightly 
##      under-predict variation in data relative to NIm2BB, but does a 
##      better job of predicting the minimum, and a tighter distribution
##      around the mean. Comparable performance for max().
##     
##  --  Overall, NIm2 is the best choice. A simpler model than NIm2BB, 
##      and no issues with X^2 discrepancy.
##
###########################################################################






###########################################################################
###########################################################################
##  A priori and posteriori conrasts of interest
###########################################################################
###########################################################################

print(NIm2, c("beta", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95), digits=3);
NIm2.summ[1:2,]
NIm2.betas  <-  NIm2.df[,1:2]
NIm2.BetaSumm  <-  plyr:::adply(as.matrix(NIm2.betas),2,MCMCsum)
NIm2.BetaSumm

#  Contrasts

# Bayesian P-values for Contrasts
plotContr(density(NIm2.df[,2]))

pval(NIm2.df[,2])


NIm2.neg2   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (-2))



N_investPlot(NIm2.df, NIm2.summ, NinvData)
abline(v=c(
           mean(NinvData$nSperm) - (2*sd(NinvData$nSperm)),
           mean(NinvData$nSperm) - sd(NinvData$nSperm),
           mean(NinvData$nSperm),
           mean(NinvData$nSperm) + sd(NinvData$nSperm),
           mean(NinvData$nSperm) + (2*sd(NinvData$nSperm))
           ), col=c(1,1,2,1,1), lwd=c(1,1,3,1,1), lty=c(2,2,1,2,2))

sd(NinvData$nSperm)