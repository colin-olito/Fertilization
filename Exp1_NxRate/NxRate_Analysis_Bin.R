#/* 
# * Colin Olito. Created 13/01/2016
# * Analysis of 1st flume experiment: NxRate
# * 
# * NOTES:  This file contains all the necessary
# *     		code to read in the STAN sample files
# *     		and perform posterior predictive checks
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
source('./loadData_NxRate_Bin.R')



##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

##  Overall comparison of all models
looDiff   <-  compare(m1Loo,  m2Loo,  m3Loo,  m4Loo,  m5Loo,
                      m6Loo,  m7Loo,  m8Loo,  m9Loo,  m10Loo,
                      m11Loo, m12Loo, m13Loo, m14Loo, m15Loo,
                      m16Loo, m17Loo, m18Loo, m19Loo, m20Loo,
                      m21Loo, m22Loo, m23Loo, m24Loo, m25Loo
                     )

print(looDiff, digits=4)

# LOO Results Summary Table
LooDiff  <-  makeLooTable(looDiff)
LooDiff

# Write Loo Model Selection Results to .csv
write.csv(LooDiff, file= './output/tables/NxRate_LooDiff_Bin.csv')

###########################################################################
## Main result from LOO model comparison:
##
##  -- Model m2 provides the best overall fit to the data.
##  
##  -- 2 models stand out as being the simplest models that are statistically
##     indistinguishable from the maximal model:  Models m12 & m10.
##      
##  -- Will focus analysis of posterior predictive checks on these models.
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

# Maximal model (also best fitting model)
X2data1   <-  as.numeric(m1.df[,818])
X2sim1    <-  as.numeric(m1.df[,819])

# Simplest models that are statistically
# indistinguishable from m1
X2data12   <-  as.numeric(m12.df[,788])
X2sim12    <-  as.numeric(m12.df[,789])

X2data10   <-  as.numeric(m10.df[,798])
X2sim10    <-  as.numeric(m10.df[,799])

# Minimal model
X2data25   <-  as.numeric(m25.df[,617])
X2sim25    <-  as.numeric(m25.df[,618])



# Plot Colors for X^2 Posterior Predictive Check Plots
COLS  <-  c("#009ce0",
            "#f9684a",
            "#8be36a",
            "#7a0058",
            "#9fa3ff",
            "#be0042",
            "#2d1956"
)



# Plot Range for PPC plots
plotRange  <-  c(0,max(X2data25))

# pdf(file='./output/figs/NxRate_X2Discrepancy.pdf', width=7,height=7)
par(omi=rep(0.3, 4))
plot(X2sim25 ~ X2data25, 
    xlab=expression(paste(chi^2~discrepancy~of~observed~data)), ylab=expression(paste(chi^2~discrepancy~of~simulated~data)), 
    main=expression(paste(Posterior~predictive~check:~chi^2~"discrepancy")),
    type='n', axes=FALSE, xlim=plotRange, ylim=plotRange, xpd=NA)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2sim25 ~ X2data25, pch=21, 
        bg=transparentColor(COLS[4], 0.1),
        col=transparentColor(COLS[4], 0.3), cex=1.1)
points(X2sim10 ~ X2data10, pch=21, 
        bg=transparentColor(COLS[7], 0.1),
        col=transparentColor(COLS[7], 0.3), cex=1.1)
points(X2sim12 ~ X2data12, pch=21, 
        bg=transparentColor(COLS[2], 0.1),
        col=transparentColor(COLS[2], 0.3), cex=1.1)
points(X2sim1 ~ X2data1, pch=21, 
        bg=transparentColor('dodgerblue1', 0.1),
        col=transparentColor('dodgerblue1', 0.3), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.17,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Model~1)),
                      expression(paste(Model~"12")),
                      expression(paste(Model~"10")),
                      expression(paste(Model~"25"))),
          pch     =  21,
          pt.bg   =  c(
                       transparentColor('dodgerblue1',0.5),
                       transparentColor(COLS[2],0.5),
                       transparentColor(COLS[7],0.5),
                       transparentColor(COLS[4],0.5)
                       ),
          col     =  c(
                       transparentColor('dodgerblue1', 0.7), 
                       transparentColor(COLS[2], 0.7), 
                       transparentColor(COLS[7], 0.7), 
                       transparentColor(COLS[4], 0.7)
                       ),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
# dev.off()

###########################################################################
## Main result from X^2 discrepancy comparison:
##
##  -- Results from X^2 discrepancy are congruent with LOO results.
##  
##  -- Models m2, m12, and m10 overlay eachother almost perfectly. There isn't
##     much to choose between these models.
##
###########################################################################







###########################################################################
###########################################################################
#  MODEL DIAGNOSTICS FOR BEST CANDIDATE MODELS
###########################################################################
###########################################################################


##########################################################################
# Model: m2, the 'best-fitting model'
##########################################################################

##############
# Diagnostics

# Model Results
print(m2, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m2, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m2, pars="beta")
pairs(m2, pars="beta")
par(mfrow=c(5,5))
rstan::traceplot(m2, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*25 ~ x, lwd=2, type='l')
lines(density(m2.df$sigma_gamma), lwd=3, col='dodgerblue1')


x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m2.df[,i]), lwd=3, col=i)
}

##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m2.df[1,330:449])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m2.df[i,330:449])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m2.df[,450], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m2.df[,451]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m2.df[,452]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m2.df[,453]), xlim=c(min(m2.df[,453],sd(data$nFert)),max(m2.df[,453],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m2, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()




##########################################################################
# Model: m12
##########################################################################

##############
# Diagnostics

# Model Results
print(m12, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m12, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m12, c("beta", "gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m12, pars="beta")
pairs(m12, pars="beta")
rstan::traceplot(m12, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)
dev.off()

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m12.df$sigma_gamma), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m12.df[,i]), lwd=3, col=i)
}


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m12.df[1,300:419])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m12.df[i,300:419])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m12.df[,419], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m12.df[,420]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m12.df[,421]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m12.df[,422]), xlim=c(min(m12.df[,423],sd(data$nFert)),max(m12.df[,423],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m12, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()





##########################################################################
# Model: m10
##########################################################################

##############
# Diagnostics

# Model Results
print(m10, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m10, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m10, pars="beta")
pairs(m10, pars="beta")
rstan::traceplot(m10, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m10.df$sigma_gamma), lwd=3, col='dodgerblue1')


x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m10.df[,i]), lwd=3, col=i)
}

##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m10.df[1,310:429])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m10.df[i,310:429])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m10.summ
par(mfrow=c(2,2))
plot(density(m10.df[,430], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m10.df[,431]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m10.df[,432]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m10.df[,433]), xlim=c(min(m10.df[,433],sd(data$nFert)),max(m10.df[,433],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m10, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()



###########################################################################
###########################################################################
##  Have a closer look at Posterior Predictive Checks for Candidate Models
###########################################################################
###########################################################################

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
 pdf(file='./output/figs/NxRate_Bin_PPDs.pdf', width=7,height=7)
par(mfrow=c(2,2))
plot(density(m12.df[,420],   adjust=3), lwd=3, col='dodgerBlue1', main='min_y_rep (min. num. Successes)')
lines(density(m10.df[,430], adjust=3), lwd=3, col=COLS[7])
lines(density(m2.df[,450], adjust=3), lwd=3, col=COLS[2])
abline(v=min(data$nFert), lwd=3, col=2)
    legend(
          x       =  7.25,
          y       =  0.325,
          legend  =  c(
                      expression(paste(Mod.~12)),
                      expression(paste(Mod.~"10")),
                      expression(paste(Mod.~"2"))),
          lty     =  1,
          lwd     =  3,
          col     =  c(
                       'dodgerblue1',
                       COLS[2],
                       COLS[7]),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )

plot(density(m12.df[,421]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
lines(density(m10.df[,431]), lwd=3, col=COLS[7])
lines(density(m2.df[,451]), lwd=3, col=COLS[2])
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m12.df[,422]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
lines(density(m10.df[,432]), lwd=3, col=COLS[7])
lines(density(m2.df[,452]), lwd=3, col=COLS[2])
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m12.df[,423]), xlim=c(min(m2.df[,453],sd(data$nFert)),max(m2.df[,453],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
lines(density(m10.df[,433]), lwd=3, col=COLS[7])
lines(density(m2.df[,453]), lwd=3, col=COLS[2])
abline(v=sd(data$nFert), lwd=3, col=2)

dev.off()

print(m1, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m12, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m10, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


###########################################################################
## Main result from PPD's:
##
##  --  Not much to choose between models m2, m12, and m10. m12 may 
##      slightly under-predict max() and sd() relative to the other models,
##      but nothing of real significance.
##     
##  --  Overall, m12 is probably the best choice because it has the fewest
##      parameters, and doesn't noticeably under-perform relative to the 
##      more complex models.
##
##  --  HOWEVER, all models systematically under-predict the sd() in the 
##      data... suggesting that these data are over-dispersed for a
##      Binomial error distribution. Will re-run analyses using 
##      Beta-Binomial error to see if this improves the fit compared to 
##      model m12
##
###########################################################################




###########################################################################
###########################################################################
##  Compare Residual Plots for Candidate Models
###########################################################################
###########################################################################

##  look at residuals for m2
#m1yhat        <-  inv_logit(m1.summ$Mean[90:209]) # mus
m2yhat        <-  (m2.summ$Mean[458:577])
m2.resids     <-  (((data$nFert - data$nControlFert)/data$nEggs) - m2yhat)/sd(m2yhat)
m2.resids_z   <-  (((data$nFert - data$nControlFert)/data$nEggs) - m2yhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m2yhat))
#m12yhat       <-  inv_logit(m12.summ$Mean[60:179]) # mus
m12yhat       <-  (m12.summ$Mean[428:547])
m12.resids    <-  (((data$nFert - data$nControlFert)/data$nEggs) - m12yhat)/sd(m12yhat)
m12.resids_z  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m12yhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m12yhat))
#m10yhat        <-  inv_logit(m10.summ$Mean[70:189]) # mus
m10yhat        <-  (m10.summ$Mean[438:557])
m10.resids     <-  (((data$nFert - data$nControlFert)/data$nEggs) - m10yhat)/sd(m10yhat)
m10.resids_z   <-  (((data$nFert - data$nControlFert)/data$nEggs) - m10yhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m10yhat))

plot(m2yhat ~ m12yhat, xlim=c(0,1), ylim=c(0,1))
abline(a = 0, b = 1)
plot(m12yhat ~ m10yhat, xlim=c(0,1), ylim=c(0,1))
abline(a = 0, b = 1)

##  Model 2 Residual Plots
par(mfrow=c(2,2))
hist(m2.resids_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(m2.resids_z ~ data$nSperm_z, main="Model m1")
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m2.resids_z ~ seq_along(m2.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m2.resids_z)
qqline(m2.resids_z, col = 2)

dev.new()
par(mfrow=c(2,2))
hist(m12.resids_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(m12.resids_z ~ data$nSperm_z, main="Model m12")
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m12.resids_z ~ seq_along(m12.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m12.resids_z)
qqline(m12.resids_z, col = 2)

dev.new()
par(mfrow=c(2,2))
hist(m10.resids_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(m10.resids_z ~ data$nSperm_z, main="Model m10")
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m10.resids_z ~ seq_along(m10.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m10.resids_z)
qqline(m10.resids_z, col = 2)



###########################################################################
## Overall impression from residual plots:
##
##  -- Again, all models pretty equivocal, and the residual plots don't
##     appear to have any major problems other than slightly heavy tails.
###########################################################################












######################################################
## Have a look at the final regressions from 
##  candidate models (m1, m12, m10)

##  Calculate Predicted Lines
m2.coef   <-  m2.summ$Mean[1:8]
m10.coef  <-  m10.summ$Mean[1:8]
m12.coef  <-  m12.summ$Mean[1:8]

m2Fast5   <-  inv_logit(m2.coef[1] + m2.coef[2] * data$nSperm_z)
m2Fast55  <-  inv_logit((m2.coef[1] + m2.coef[4]) + (m2.coef[2] + m2.coef[6]) * data$nSperm_z)
m2Slow5   <-  inv_logit((m2.coef[1] + m2.coef[3]) + (m2.coef[2] + m2.coef[5]) * data$nSperm_z)
m2Slow55  <-  inv_logit((m2.coef[1] + m2.coef[3] + m2.coef[7]) + (m2.coef[2] + m2.coef[5] + m2.coef[8]) * data$nSperm_z)
m2Fast    <-  inv_logit((m2.coef[1] + (m2.coef[4])/2) + (m2.coef[2] + (m2.coef[6])/2) * data$nSperm_z)
m2Slow    <-  inv_logit((m2.coef[1] + m2.coef[3] + (0.5*(m2.coef[7]))) + (m2.coef[2] + m2.coef[5] + (0.5*(m2.coef[8]))) * data$nSperm_z)

m10Fast5   <-  inv_logit(m10.coef[1] + m10.coef[2] * data$nSperm_z)
m10Fast55  <-  inv_logit((m10.coef[1] + m10.coef[4]) + (m10.coef[2] + m10.coef[6]) * data$nSperm_z)
m10Slow5   <-  inv_logit((m10.coef[1] + m10.coef[3]) + (m10.coef[2] + m10.coef[5]) * data$nSperm_z)
m10Slow55  <-  inv_logit((m10.coef[1] + m10.coef[3] + m10.coef[7]) + (m10.coef[2] + m10.coef[5] + m10.coef[8]) * data$nSperm_z)
m10Fast    <-  inv_logit((m10.coef[1] + (m10.coef[4])/2) + (m10.coef[2] + (m10.coef[6])/2) * data$nSperm_z)
m10Slow    <-  inv_logit((m10.coef[1] + m10.coef[3] + (0.5*(m10.coef[7]))) + (m10.coef[2] + m10.coef[5] + (0.5*(m10.coef[8]))) * data$nSperm_z)

m12Fast5   <-  inv_logit(m12.coef[1] + m12.coef[2] * data$nSperm_z)
m12Fast55  <-  inv_logit((m12.coef[1] + m12.coef[4]) + (m12.coef[2] + m12.coef[6]) * data$nSperm_z)
m12Slow5   <-  inv_logit((m12.coef[1] + m12.coef[3]) + (m12.coef[2] + m12.coef[5]) * data$nSperm_z)
m12Slow55  <-  inv_logit((m12.coef[1] + m12.coef[3] + m12.coef[7]) + (m12.coef[2] + m12.coef[5] + m12.coef[8]) * data$nSperm_z)
m12Fast    <-  inv_logit((m12.coef[1] + (m12.coef[4])/2) + (m12.coef[2] + (m12.coef[6])/2) * data$nSperm_z)
m12Slow    <-  inv_logit((m12.coef[1] + m12.coef[3] + (0.5*(m12.coef[7]))) + (m12.coef[2] + m12.coef[5] + (0.5*(m12.coef[8]))) * data$nSperm_z)

# pdf(file='output/xRatexEggPos_m3.pdf', height=7, width=7)
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'),  main="Model 1",
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m2Fast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m2Fast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m2Slow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m2Slow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
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
    xlab='Sperm released', ylab=substitute('Fertilization rate'), main="Model 10",
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m10Fast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m10Fast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m10Slow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m10Slow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
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




dev.new()
# par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'),  main="Model 12",
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m12Fast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m12Fast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m12Slow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m12Slow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
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







###########################################################################
###########################################################################
##  A priori and posteriori conrasts of interest
###########################################################################
###########################################################################


m2.allBetas   <-  as.matrix(m2.df[1:8])
m2.allGammas  <-  as.matrix(m2.df[9:88])
m10.allBetas   <-  as.matrix(m10.df[1:8])
m10.allGammas  <-  as.matrix(m10.df[9:68])
m12.allBetas   <-  as.matrix(m12.df[1:8])
m12.allGammas  <-  as.matrix(m12.df[9:59])


##  Calculate Coefficients
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
(Xnames)

# for model m1
b0Fast    <-  inv_logit((m2.allBetas[,1] + (m2.allBetas[,4])/2))
b0Slow    <-  inv_logit((m2.allBetas[,1] + m2.allBetas[,3] + (0.5*(m2.allBetas[,7]))))
b0Fast5   <-  inv_logit(m2.allBetas[,1])
b0Fast55  <-  inv_logit((m2.allBetas[,1] + m2.allBetas[,4]))
b0Slow5   <-  inv_logit((m2.allBetas[,1] + m2.allBetas[,3]))
b0Slow55  <-  inv_logit((m2.allBetas[,1] + m2.allBetas[,3] + m2.allBetas[,7]))
b1Fast    <-  inv_logit((m2.allBetas[,2] + (m2.allBetas[,6])/2))
b1Slow    <-  inv_logit((m2.allBetas[,2] + m2.allBetas[,5] + (0.5*(m2.allBetas[,8]))))
b1Fast5   <-  inv_logit(m2.allBetas[,2])
b1Fast55  <-  inv_logit((m2.allBetas[,2] + m2.allBetas[,6]))
b1Slow5   <-  inv_logit((m2.allBetas[,2] + m2.allBetas[,5]))
b1Slow55  <-  inv_logit((m2.allBetas[,2] + m2.allBetas[,5] + m2.allBetas[,8]))

# for model m10
b0Fast    <-  inv_logit((m10.allBetas[,1] + (m10.allBetas[,4])/2))
b0Slow    <-  inv_logit((m10.allBetas[,1] + m10.allBetas[,3] + (0.5*(m10.allBetas[,7]))))
b0Fast5   <-  inv_logit(m10.allBetas[,1])
b0Fast55  <-  inv_logit((m10.allBetas[,1] + m10.allBetas[,4]))
b0Slow5   <-  inv_logit((m10.allBetas[,1] + m10.allBetas[,3]))
b0Slow55  <-  inv_logit((m10.allBetas[,1] + m10.allBetas[,3] + m10.allBetas[,7]))
b1Fast    <-  inv_logit((m10.allBetas[,2] + (m10.allBetas[,6])/2))
b1Slow    <-  inv_logit((m10.allBetas[,2] + m10.allBetas[,5] + (0.5*(m10.allBetas[,8]))))
b1Fast5   <-  inv_logit(m10.allBetas[,2])
b1Fast55  <-  inv_logit((m10.allBetas[,2] + m10.allBetas[,6]))
b1Slow5   <-  inv_logit((m10.allBetas[,2] + m10.allBetas[,5]))
b1Slow55  <-  inv_logit((m10.allBetas[,2] + m10.allBetas[,5] + m10.allBetas[,8]))

# for model m12
b0Fast    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2))
b0Slow    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))))
b0Fast5   <-  inv_logit(m12.allBetas[,1])
b0Fast55  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]))
b0Slow5   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]))
b0Slow55  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]))
b1Fast    <-  inv_logit((m12.allBetas[,2] + (m12.allBetas[,6])/2))
b1Slow    <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))))
b1Fast5   <-  inv_logit(m12.allBetas[,2])
b1Fast55  <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,6]))
b1Slow5   <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,5]))
b1Slow55  <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]))

b0  <-  inv_logit((m12.allBetas[,1]))
b1  <-  inv_logit((m12.allBetas[,2]))
sigma_gamma12  <-  m12.allGammas[,51]

head(sigma_gamma12)
#  Contrasts
c1   <-  b1Fast   - b1Slow
c2   <-  b1Fast5  - b1Fast55
c3   <-  b1Slow5  - b1Slow55
c4   <-  b1Fast5  - b1Slow5
c5   <-  b1Fast55 - b1Slow55
c6   <-  b1Fast55 - b1Slow
c7   <-  b1Fast5  - b1Slow

# Summarize contrasts
plyr:::adply(as.matrix(cbind(c1,c2,c3,c4,c5,c6,c7)),2,MCMCsum)


# Summarize back-transformed coefficients
plyr:::adply(as.matrix(cbind(b1Fast,b1Slow,b1Fast5,b1Fast55,b1Slow5,b1Slow55)),2,MCMCsum)
plyr:::adply(as.matrix(cbind(b0,b1,b0Fast, b0Slow, b0Fast5, b0Fast55, b0Slow5, b0Slow55, b1Fast, b1Slow, b1Fast5, b1Fast55, b1Slow5, b1Slow55, sigma_gamma12)),2,MCMCsum)




# Bayesian P-values for Contrasts
pval(c1)     # b1Slow   - b1Fast
pval(c2)     # b1Fast5  - b1Fast55
pval(c3)     # b1Slow5  - b1Slow55
pval(c4)     # b1Fast5  - b1Slow5
pval(c5)     # b1Fast55 - b1Slow55
pval(c6)     # b1Fast55 - b1Slow
pval(c7)     # b1Fast5  - b1Slow


par(mfrow=c(2,4))
plotContr(density(c1))
plotContr(density(c2))
plotContr(density(c3))
plotContr(density(c4))
plotContr(density(c5))
plotContr(density(c6))
plotContr(density(c7))



#  Simple Contrasts at -1, 0, 1, & 2 standard deviations of nSperm
m12Fast5.neg2   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (-2))
m12Fast55.neg2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (-2))
m12Slow5.neg2   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (-2))
m12Slow55.neg2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (-2))
m12Fast.neg2    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (-2))
m12Slow.neg2    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (-2))

m12Fast5.neg1   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (-1))
m12Fast55.neg1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (-1))
m12Slow5.neg1   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (-1))
m12Slow55.neg1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (-1))
m12Fast.neg1    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (-1))
m12Slow.neg1    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (-1))

m12Fast5.0   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (0))
m12Fast55.0  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (0))
m12Slow5.0   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (0))
m12Slow55.0  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (0))
m12Fast.0    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (0))
m12Slow.0    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (0))

m12Fast5.1   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (1))
m12Fast55.1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (1))
m12Slow5.1   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (1))
m12Slow55.1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (1))
m12Fast.1    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (1))
m12Slow.1    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (1))

m12Fast5.2   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (2))
m12Fast55.2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (2))
m12Slow5.2   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (2))
m12Slow55.2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (2))
m12Fast.2    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (2))
m12Slow.2    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (2))


simpContr  <-  list(
  cSimp1   =  m12Fast5.neg1 - m12Fast55.neg1,
  cSimp2   =  m12Slow5.neg1 - m12Slow55.neg1,
  cSimp3   =  m12Fast5.neg1 - m12Slow5.neg1,
  cSimp4   =  m12Fast55.neg1 - m12Slow55.neg1,
  cSimp5   =  m12Fast.neg1 - m12Slow.neg1,
  cSimp6   =  m12Fast5.0 - m12Fast55.0,
  cSimp7   =  m12Slow5.0 - m12Slow55.0,
  cSimp8   =  m12Fast5.0 - m12Slow5.0,
  cSimp9   =  m12Fast55.0 - m12Slow55.0,
  cSimp10  =  m12Fast.0 - m12Slow.0,
  cSimp11  =  m12Fast5.1 - m12Fast55.1,
  cSimp12  =  m12Slow5.1 - m12Slow55.1,
  cSimp13  =  m12Fast5.1 - m12Slow5.1,
  cSimp14  =  m12Fast55.1 - m12Slow55.1,
  cSimp15  =  m12Fast.1 - m12Slow.1,
  cSimp16  =  m12Fast5.2 - m12Fast55.2,
  cSimp17  =  m12Slow5.2 - m12Slow55.2,
  cSimp18  =  m12Fast5.2 - m12Slow5.2,
  cSimp19  =  m12Fast55.2 - m12Slow55.2,
  cSimp20  =  m12Fast.2 - m12Slow.2
)
pval(simpContr[[5]])

pdf(file="./output/figs/NxRate_SimpContrasts_Bin.pdf", height=18, width=20)
par(mfrow=c(4,5), omi=rep(0.4,4))
plotContr(density(simpContr[[1]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Fast.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  proportionalLabel(-0.25, 0.5, expression(paste('x = -1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[2]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow.5 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[3]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Slow.5')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[4]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.55 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[5]]))
  proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow - Fast')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
for (i in 6:20) {
  plotContr(density(simpContr[[i]]))
  if(i == 6)
    proportionalLabel(-0.25, 0.5, expression(paste('x = 0',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  if(i == 11)
    proportionalLabel(-0.25, 0.5, expression(paste('x = 1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
  if(i == 16)
    proportionalLabel(-0.25, 0.5, expression(paste('x = 2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
}
dev.off()



###########################################################################
###########################################################################
## Main result from a posteriori contrasts for model m12:
##
## There is a significant main effect of Rate for the slope of the 
## regression lines (p = 1.000). 
## 
## The Rate x EggPos interaction for the Fast slopes is marginal
## (p = 0.944).
## 
##  There is no corresponding interaction for the Slow treatment. 
##
##  Both the 5cm and 55cm partial regression slopes for the Fast treatment
##  are significantly steeper than the "Slow" slope (p = 0.000, 0.000) 
## 
##  The simple contrasts between the predicted lines for Fast v. Slow
##  treatments are significant at high and low values of sigma (sigma = 1, 
##  sigma > 1).
##
##
###########################################################################
###########################################################################



###########################################################################
###########################################################################
#  Plot of run-adjusted y values for the main regression plot

# Fixed effects Model Matrix
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]

# Random effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : Rate   : EggPos,
                         data = data)

# Mean coefficients
m12.betas    <-  m12.summ$Mean[1:8]
m12.gammas   <-  m12.summ$Mean[9:58]

source('R/functions-figures.R')
# Predicted lines for Fast (with pooled slopes) & the overall Slow
# Create plotting objects for each regression line
Fast5.plt   <-  Fast5.plots(m12.betas, m12.allBetas, m12.gammas, Z=Z, data=data)
Fast55.plt  <-  Fast55.plots(m12.betas, m12.allBetas, m12.gammas, Z=Z, data=data)
Slow.plt    <-  Slow.plots(m12.betas, m12.allBetas, m12.gammas, Z=Z, data=data)

# Residuals for fixed-effect predicted lines
par(mfrow=c(2,2))
plot(density(Fast5.plt$rawResids[data$Rate == "Fast" & data$EggPos == "5"]), lwd=3, col='darkolivegreen', ylim=c(0,max(density(Fast5.plt$adjResids[data$Rate == "Fast" & data$EggPos == "5"])$y)))
lines(density(Fast5.plt$adjResids[data$Rate == "Fast" & data$EggPos == "5"]), lwd=3, col='dodgerblue1')
abline(v=0,lwd=3,col=2)
plot(density(Fast55.plt$rawResids[data$Rate == "Fast" & data$EggPos == "55"]), lwd=3, col='darkolivegreen', ylim=c(0,max(density(Fast55.plt$adjResids[data$Rate == "Fast" & data$EggPos == "55"])$y)))
lines(density(Fast55.plt$adjResids[data$Rate == "Fast" & data$EggPos == "55"]), lwd=3, col='dodgerblue1')
abline(v=0,lwd=3,col=2)
plot(density(Slow.plt$rawResids[data$Rate == "Slow"]), lwd=3, col='darkolivegreen', ylim=c(0,max(density(Slow.plt$adjResids[data$Rate == "Slow"])$y)))
lines(density(Slow.plt$adjResids[data$Rate == "Slow"]), lwd=3, col='dodgerblue1')
abline(v=0,lwd=3,col=2)


sum(Fast5.plt$adjResids[data$Rate == "Fast" & data$EggPos == "5"] > 0)/length(Fast5.plt$adjResids[data$Rate == "Fast" & data$EggPos == "5"])
sum(Fast55.plt$adjResids[data$Rate == "Fast" & data$EggPos == "55"] > 0)/length(Fast55.plt$adjResids[data$Rate == "Fast" & data$EggPos == "55"])
sum(Slow.plt$adjResids[data$Rate == "Slow"] > 0)/length(Slow.plt$adjResids[data$Rate == "Slow"])


# Visually inspect, make sure they make sense
par(mfrow=c(2,2))
plot(Fast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ Fast5.plt$yHats[data$Rate == "Fast" & data$EggPos == "5"])
abline(a=0, b=1,lwd=3,col=2)
plot(Fast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ Fast55.plt$yHats[data$Rate == "Fast" & data$EggPos == "55"])
abline(a=0, b=1,lwd=3,col=2)
plot(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$yHats[data$Rate == "Slow"])
abline(a=0, b=1,lwd=3,col=2)

par(mfrow=c(2,2))
plot(Fast5.plt$rawResids[data$Rate == "Fast" & data$EggPos == "5"] ~ 
     Fast5.plt$adjResids[data$Rate == "Fast" & data$EggPos == "5"])
abline(a=0, b=1,lwd=3,col=2)
plot(Fast55.plt$rawResids[data$Rate == "Fast" & data$EggPos == "55"] ~ 
     Fast55.plt$adjResids[data$Rate == "Fast" & data$EggPos == "55"])
abline(a=0, b=1,lwd=3,col=2)
plot(Slow.plt$rawResids[data$Rate == "Slow"] ~ 
     Slow.plt$adjResids[data$Rate == "Slow"])
abline(a=0, b=1,lwd=3,col=2)


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
lines(Slow.plt$y   ~ Slow.plt$xRaw, col='orangered1', lwd=3)
lines(Fast5.plt$y  ~ Fast5.plt$xRaw, col='dodgerblue1', lwd=3)
lines(Fast55.plt$y ~ Fast55.plt$xRaw, col='dodgerblue1', lty=2, lwd=3)
points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(Fast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ Fast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(Fast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ Fast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
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
par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='', ylab='', 
    type='n', axes=FALSE, ylim=c(-0.05,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
polygon(x=c(Slow.plt$xRaw, rev(Slow.plt$xRaw)), 
        y=c(Slow.plt$CIs$lower, rev(Slow.plt$CIs$upper)), 
        col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
polygon(x=c(Fast5.plt$xRaw, rev(Fast5.plt$xRaw)), 
        y=c(Fast5.plt$CIs$lower, rev(Fast5.plt$CIs$upper)), 
        col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
polygon(x=c(Fast55.plt$xRaw, rev(Fast55.plt$xRaw)), 
        y=c(Fast55.plt$CIs$lower, rev(Fast55.plt$CIs$upper)), 
        col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
lines(Slow.plt$y   ~ Slow.plt$xRaw, col='orangered1', lwd=3)
lines(Fast5.plt$y  ~ Fast5.plt$xRaw, col='dodgerblue1', lwd=3)
lines(Fast55.plt$y ~ Fast55.plt$xRaw, col='dodgerblue1', lty=2, lwd=3)
points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(Fast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ Fast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(Fast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ Fast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
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









