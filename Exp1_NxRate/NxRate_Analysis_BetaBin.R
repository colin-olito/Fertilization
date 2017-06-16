#/* 
# * Colin Olito. Created 09/02/2017
# * Analysis of 1st flume experiment: N-invest Beta-Binomial
# * 
# * NOTES:  This file contains all the necessary
# *     		code to read in the STAN sample files
# *     		and perform posterior predictive checks
# *         model selection, and produce regression
# *         plots for the NxRate flume data fitted 
# *         using Beta-Binomial error.
# * 
# */

rm(list=ls())
#################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

###############
# DEPENDENCIES
source('./loadData_NxRate_BetaBin.R')



##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

##  Overall comparison of all models
looDiff   <-  compare(m1BBLoo,  m2BBLoo,  m3BBLoo,  m4BBLoo,  m5BBLoo,
                      m6BBLoo,  m7BBLoo,  m8BBLoo,  m9BBLoo,  m10BBLoo,
                      m11BBLoo, m12BBLoo, m13BBLoo, m14BBLoo, m15BBLoo,
                      m16BBLoo, m17BBLoo, m18BBLoo, m19BBLoo, m20BBLoo,
                      m21BBLoo, m22BBLoo, m23BBLoo, m24BBLoo, m25BBLoo
                     )
print(looDiff, digits=4)

# LOO Results Summary Table
LooDiff  <-  makeLooTable(looDiff)
LooDiff

# Write Loo Model Selection Results to .csv
write.csv(LooDiff, file= './output/tables/NxRate_LooDiff_BetaBin.csv')

###########################################################################
## Main result from LOO model comparison:
##
##  -- Model m12BB provides the best overall fit to the data. Interesting,
##     since this was the most parsimonious Binomial model.
##  
##  -- Sevaral candidate models stand out being potentially simpler models
##     than m12BB: models m21, m19BB, and m17BB.
##      
##  -- Will focus analysis of posterior predictive checks on models m12BB,
##     m19BB, and m21BB.
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

# Best fitting model
X2data12   <-  as.numeric(m12BB.df[,909])
X2sim12    <-  as.numeric(m12BB.df[,910])

# Simplest models that are statistically
# indistinguishable from m12BB
X2data21   <-  as.numeric(m21BB.df[,879])
X2sim21    <-  as.numeric(m21BB.df[,880])

X2data19   <-  as.numeric(m19BB.df[,889])
X2sim19    <-  as.numeric(m19BB.df[,890])

# Minimal model
X2data25   <-  as.numeric(m25BB.df[,858])
X2sim25    <-  as.numeric(m25BB.df[,859])



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
plotRange  <-  c(0,max(X2sim25))

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
points(X2sim19 ~ X2data19, pch=21, 
        bg=transparentColor(COLS[7], 0.1),
        col=transparentColor(COLS[7], 0.3), cex=1.1)
points(X2sim21 ~ X2data21, pch=21, 
        bg=transparentColor(COLS[2], 0.1),
        col=transparentColor(COLS[2], 0.3), cex=1.1)
points(X2sim12 ~ X2data12, pch=21, 
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
##  -- Models m12BB, m19BB overlay each other almost entirely. m21BB has 
##     slightly higher discrepancy, but the difference pretty small, in 
##     accordance with the LOO results. As in the N_invest analysis, the 
##     X^2 discrepancy for the simulated data is noticeably higher than for
##     the real data... not sure exactly why this is happening, except that
##     perhaps this is suggesting that the data is under-dispersed for a 
##     Beta-Binomial distribution.
##
###########################################################################







###########################################################################
###########################################################################
#  MODEL DIAGNOSTICS FOR BEST CANDIDATE MODELS
###########################################################################
###########################################################################

##########################################################################
# Model: m12BB, the 'best-fitting model'
##########################################################################

##############
# Diagnostics

# Model Results
print(m12BB, c("beta","phi"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m12BB, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m12BB, c("a", "b"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m12BB, pars="beta")
pairs(m12BB, pars="beta")
rstan::traceplot(m12BB, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)
dev.off()

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m12BB.df$sigma_gamma), lwd=3, col='dodgerblue1')

x  <-  seq(from=0, to=100, length=500)
plot(dcauchy(x, scale=100)*20 ~ x, lwd=2, type='l')
lines(density(m12BB.df$phi), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*22 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m12BB.df[,i]), lwd=3, col=i)
}


##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m12BB.df[1,541:660])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m12BB.df[i,541:660])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(m12BB.df[,661], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m12BB.df[,662]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m12BB.df[,663]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m12BB.df[,664]), xlim=c(min(m12BB.df[,664],sd(data$nFert)),max(m12BB.df[,664],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m12BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()





##########################################################################
# Model: m19BB
##########################################################################

##############
# Diagnostics

# Model Results
print(m19BB, c("beta", "phi"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m19BB, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m19BB, c("a", "b"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m19BB, pars="beta")
pairs(m19BB, pars="beta")
rstan::traceplot(m19BB, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m19BB.df$sigma_gamma), lwd=3, col='dodgerblue1')


x  <-  seq(from=0, to=100, length=500)
plot(dcauchy(x, scale=100)*25 ~ x, lwd=2, type='l')
lines(density(m19BB.df$phi), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m19BB.df[,i]), lwd=3, col=i)
}

##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m19BB.df[1,521:640])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m19BB.df[i,521:640])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m19BB.summ
par(mfrow=c(2,2))
plot(density(m19BB.df[,641], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m19BB.df[,642]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m19BB.df[,643]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m19BB.df[,644]), xlim=c(min(m19BB.df[,644],sd(data$nFert)),max(m19BB.df[,644],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m19BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()



##########################################################################
# Model: m21BB
##########################################################################

##############
# Diagnostics

# Model Results
print(m21BB, c("beta", "phi"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("a", "b"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m21BB, pars="beta")
pairs(m21BB, pars="beta")
rstan::traceplot(m21BB, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)

# Check posteriors against priors
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m21BB.df$sigma_gamma), lwd=3, col='dodgerblue1')


x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m21BB.df[,i]), lwd=3, col=i)
}

##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(m21BB.df[1,511:630])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m21BB.df[i,511:630])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m21BB.summ
par(mfrow=c(2,2))
plot(density(m21BB.df[,631], adjust=3), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(data$nFert), lwd=3, col=2)

plot(density(m21BB.df[,632]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m21BB.df[,633]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m21BB.df[,634]), xlim=c(min(m21BB.df[,634],sd(data$nFert)),max(m21BB.df[,634],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(data$nFert), lwd=3, col=2)

print(m21BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

dev.off()






###########################################################################
###########################################################################
##  Have a closer look at Posterior Predictive Checks for Candidate Models
###########################################################################
###########################################################################

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
 pdf(file='./output/figs/NxRate_BetaBin_PPDs.pdf', width=7,height=7)
par(mfrow=c(2,2))
plot(density(m12BB.df[,661],   adjust=3), lwd=3, col='dodgerBlue1', main='min_y_rep (min. num. Successes)')
lines(density(m19BB.df[,641], adjust=3), lwd=3, col=COLS[7])
lines(density(m21BB.df[,631], adjust=3), lwd=3, col=COLS[2])
abline(v=min(data$nFert), lwd=3, col=2)
    legend(
          x       =  7.25,
          y       =  0.325,
          legend  =  c(
                      expression(paste(Mod.~12)),
                      expression(paste(Mod.~"19")),
                      expression(paste(Mod.~"21"))),
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

plot(density(m12BB.df[,662]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
lines(density(m19BB.df[,642]), lwd=3, col=COLS[7])
lines(density(m21BB.df[,632]), lwd=3, col=COLS[2])
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m12BB.df[,663]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
lines(density(m19BB.df[,643]), lwd=3, col=COLS[7])
lines(density(m21BB.df[,633]), lwd=3, col=COLS[2])
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m12BB.df[,664]), xlim=c(min(m12BB.df[,664],sd(data$nFert)),max(m12BB.df[,664],sd(data$nFert))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
lines(density(m19BB.df[,644]), lwd=3, col=COLS[7])
lines(density(m21BB.df[,634]), lwd=3, col=COLS[2])
abline(v=sd(data$nFert), lwd=3, col=2)

dev.off()

print(m1, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m12, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


###########################################################################
## Main result from PPD's:
##
##  --  Models m12BB, m19BB, and m21BB are extremely similar. They all 
##      under-predict min(). But they also seem to do a better job of 
##      predicting sd() relative to the Binomial models.
##     
##  --  Overall, m21 is the most parsimonioius model, as it has far fewer
##      parameters, and doesn't noticeably under-perform relative to the 
##      more complex models.
##
###########################################################################




###########################################################################
###########################################################################
##  Compare Residual Plots for Candidate Models
###########################################################################
###########################################################################

##  look at residuals for m3
#m12yhat       <-  inv_logit(m12.summ$Mean[60:179]) # mus
m12BByhat       <-  (m12BB.summ$Mean[301:420])
m12BB.resids    <-  (((data$nFert - data$nControlFert)/data$nEggs) - m12BByhat)/sd(m12BByhat)
m12BB.resids_z  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m12BByhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m12BByhat))
#m21BByhat        <-  inv_logit(m21BB.summ$Mean[70:189]) # mus
m21BByhat        <-  (m21BB.summ$Mean[271:390])
m21BB.resids     <-  (((data$nFert - data$nControlFert)/data$nEggs) - m21BByhat)/sd(m21BByhat)
m21BB.resids_z   <-  (((data$nFert - data$nControlFert)/data$nEggs) - m21BByhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m21BByhat))

plot(m12BByhat ~ m21BByhat, xlim=c(0,1), ylim=c(0,1))
abline(a = 0, b = 1)

##  Model 12 Residual Plots
dev.new()
par(mfrow=c(2,2))
hist(m12BB.resids_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(m12BB.resids_z ~ data$nSperm_z, main="Model m12BB")
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m12BB.resids_z ~ seq_along(m12BB.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m12BB.resids_z)
qqline(m12BB.resids_z, col = 2)

dev.new()
par(mfrow=c(2,2))
hist(m21BB.resids_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(m21BB.resids_z ~ data$nSperm_z, main="Model m21BB")
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m21BB.resids_z ~ seq_along(m21BB.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m21BB.resids_z)
qqline(m21BB.resids_z, col = 2)



###########################################################################
## Overall impression from residual plots:
##
##  -- Again, all models pretty equivocal, and the residual plots don't
##     appear to have any major problems other than slightly heavy tails.
###########################################################################







######################################################
## Have a look at the final regressions from 
##  candidate models (m1, m12, m21BB)

##  Calculate Predicted Lines
m12BB.coef  <-  m12BB.summ$Mean[1:8]
m21BB.coef  <-  m21BB.summ$Mean[1:8]

m12BBFast5   <-  inv_logit(m12BB.coef[1] + m12BB.coef[2] * data$nSperm_z)
m12BBFast55  <-  inv_logit((m12BB.coef[1] + m12BB.coef[4]) + (m12BB.coef[2] + m12BB.coef[6]) * data$nSperm_z)
m12BBSlow5   <-  inv_logit((m12BB.coef[1] + m12BB.coef[3]) + (m12BB.coef[2] + m12BB.coef[5]) * data$nSperm_z)
m12BBSlow55  <-  inv_logit((m12BB.coef[1] + m12BB.coef[3] + m12BB.coef[7]) + (m12BB.coef[2] + m12BB.coef[5] + m12BB.coef[8]) * data$nSperm_z)
m12BBFast    <-  inv_logit((m12BB.coef[1] + (m12BB.coef[4])/2) + (m12BB.coef[2] + (m12BB.coef[6])/2) * data$nSperm_z)
m12BBSlow    <-  inv_logit((m12BB.coef[1] + m12BB.coef[3] + (0.5*(m12BB.coef[7]))) + (m12BB.coef[2] + m12BB.coef[5] + (0.5*(m12BB.coef[8]))) * data$nSperm_z)

m21BBFast5   <-  inv_logit(m21BB.coef[1] + m21BB.coef[2] * data$nSperm_z)
m21BBFast55  <-  inv_logit((m21BB.coef[1] + m21BB.coef[4]) + (m21BB.coef[2] + m21BB.coef[6]) * data$nSperm_z)
m21BBSlow5   <-  inv_logit((m21BB.coef[1] + m21BB.coef[3]) + (m21BB.coef[2] + m21BB.coef[5]) * data$nSperm_z)
m21BBSlow55  <-  inv_logit((m21BB.coef[1] + m21BB.coef[3] + m21BB.coef[7]) + (m21BB.coef[2] + m21BB.coef[5] + m21BB.coef[8]) * data$nSperm_z)
m21BBFast    <-  inv_logit((m21BB.coef[1] + (m21BB.coef[4])/2) + (m21BB.coef[2] + (m21BB.coef[6])/2) * data$nSperm_z)
m21BBSlow    <-  inv_logit((m21BB.coef[1] + m21BB.coef[3] + (0.5*(m21BB.coef[7]))) + (m21BB.coef[2] + m21BB.coef[5] + (0.5*(m21BB.coef[8]))) * data$nSperm_z)

# pdf(file='output/xRatexEggPos_m3.pdf', height=7, width=7)
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'),  main="Model 12BB",
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m12BBFast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m12BBFast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m12BBSlow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m12BBSlow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
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
    xlab='Sperm released', ylab=substitute('Fertilization rate'), main="Model 21BB",
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m21BBFast5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m21BBFast55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m21BBSlow5[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                  col='orangered1', lwd=3)
lines(m21BBSlow55[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
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


m21BB.allBetas   <-  as.matrix(m21BB.df[1:8])
m21BB.allGammas  <-  as.matrix(m21BB.df[9:28])
m12BB.allBetas   <-  as.matrix(m12BB.df[1:8])
m12BB.allGammas  <-  as.matrix(m12BB.df[9:58])


##  Calculate Coefficients
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
(Xnames)


# for model m12BB
b0Fast    <-  inv_logit((m12BB.allBetas[,1] + (m12BB.allBetas[,4])/2))
b0Slow    <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + (0.5*(m12BB.allBetas[,7]))))
b0Fast5   <-  inv_logit(m12BB.allBetas[,1])
b0Fast55  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,4]))
b0Slow5   <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3]))
b0Slow55  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + m12BB.allBetas[,7]))
b1Fast    <-  inv_logit((m12BB.allBetas[,2] + (m12BB.allBetas[,6])/2))
b1Slow    <-  inv_logit((m12BB.allBetas[,2] + m12BB.allBetas[,5] + (0.5*(m12BB.allBetas[,8]))))
b1Fast5   <-  inv_logit(m12BB.allBetas[,2])
b1Fast55  <-  inv_logit((m12BB.allBetas[,2] + m12BB.allBetas[,6]))
b1Slow5   <-  inv_logit((m12BB.allBetas[,2] + m12BB.allBetas[,5]))
b1Slow55  <-  inv_logit((m12BB.allBetas[,2] + m12BB.allBetas[,5] + m12BB.allBetas[,8]))

# for model m21BB
b0Fast    <-  inv_logit((m21BB.allBetas[,1] + (m21BB.allBetas[,4])/2))
b0Slow    <-  inv_logit((m21BB.allBetas[,1] + m21BB.allBetas[,3] + (0.5*(m21BB.allBetas[,7]))))
b0Fast5   <-  inv_logit(m21BB.allBetas[,1])
b0Fast55  <-  inv_logit((m21BB.allBetas[,1] + m21BB.allBetas[,4]))
b0Slow5   <-  inv_logit((m21BB.allBetas[,1] + m21BB.allBetas[,3]))
b0Slow55  <-  inv_logit((m21BB.allBetas[,1] + m21BB.allBetas[,3] + m21BB.allBetas[,7]))
b1Fast    <-  inv_logit((m21BB.allBetas[,2] + (m21BB.allBetas[,6])/2))
b1Slow    <-  inv_logit((m21BB.allBetas[,2] + m21BB.allBetas[,5] + (0.5*(m21BB.allBetas[,8]))))
b1Fast5   <-  inv_logit(m21BB.allBetas[,2])
b1Fast55  <-  inv_logit((m21BB.allBetas[,2] + m21BB.allBetas[,6]))
b1Slow5   <-  inv_logit((m21BB.allBetas[,2] + m21BB.allBetas[,5]))
b1Slow55  <-  inv_logit((m21BB.allBetas[,2] + m21BB.allBetas[,5] + m21BB.allBetas[,8]))

#  Contrasts
c1   <-  b1Slow   - b1Fast
c2   <-  b1Fast5  - b1Fast55
c3   <-  b1Slow5  - b1Slow55
c4   <-  b1Fast5  - b1Slow5
c5   <-  b1Fast55 - b1Slow55
c6   <-  b1Fast55 - b1Slow
c7   <-  b1Fast5  - b1Slow

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
m12BBFast5.neg2   <-  inv_logit(m12BB.allBetas[,1] + m12BB.allBetas[,2] * (-2))
m12BBFast55.neg2  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,4]) + (m12BB.allBetas[,2] + m12BB.allBetas[,6]) * (-2))
m12BBSlow5.neg2   <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5]) * (-2))
m12BBSlow55.neg2  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + m12BB.allBetas[,7]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + m12BB.allBetas[,8]) * (-2))
m12BBFast.neg2    <-  inv_logit((m12BB.allBetas[,1] + (m12BB.allBetas[,4])/2) + (m12BB.allBetas[,2] + (m12BB.allBetas[,6])/2) * (-2))
m12BBSlow.neg2    <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + (0.5*(m12BB.allBetas[,7]))) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + (0.5*(m12BB.allBetas[,8]))) * (-2))

m12BBFast5.neg1   <-  inv_logit(m12BB.allBetas[,1] + m12BB.allBetas[,2] * (-1))
m12BBFast55.neg1  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,4]) + (m12BB.allBetas[,2] + m12BB.allBetas[,6]) * (-1))
m12BBSlow5.neg1   <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5]) * (-1))
m12BBSlow55.neg1  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + m12BB.allBetas[,7]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + m12BB.allBetas[,8]) * (-1))
m12BBFast.neg1    <-  inv_logit((m12BB.allBetas[,1] + (m12BB.allBetas[,4])/2) + (m12BB.allBetas[,2] + (m12BB.allBetas[,6])/2) * (-1))
m12BBSlow.neg1    <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + (0.5*(m12BB.allBetas[,7]))) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + (0.5*(m12BB.allBetas[,8]))) * (-1))

m12BBFast5.0   <-  inv_logit(m12BB.allBetas[,1] + m12BB.allBetas[,2] * (0))
m12BBFast55.0  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,4]) + (m12BB.allBetas[,2] + m12BB.allBetas[,6]) * (0))
m12BBSlow5.0   <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5]) * (0))
m12BBSlow55.0  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + m12BB.allBetas[,7]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + m12BB.allBetas[,8]) * (0))
m12BBFast.0    <-  inv_logit((m12BB.allBetas[,1] + (m12BB.allBetas[,4])/2) + (m12BB.allBetas[,2] + (m12BB.allBetas[,6])/2) * (0))
m12BBSlow.0    <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + (0.5*(m12BB.allBetas[,7]))) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + (0.5*(m12BB.allBetas[,8]))) * (0))

m12BBFast5.1   <-  inv_logit(m12BB.allBetas[,1] + m12BB.allBetas[,2] * (1))
m12BBFast55.1  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,4]) + (m12BB.allBetas[,2] + m12BB.allBetas[,6]) * (1))
m12BBSlow5.1   <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5]) * (1))
m12BBSlow55.1  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + m12BB.allBetas[,7]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + m12BB.allBetas[,8]) * (1))
m12BBFast.1    <-  inv_logit((m12BB.allBetas[,1] + (m12BB.allBetas[,4])/2) + (m12BB.allBetas[,2] + (m12BB.allBetas[,6])/2) * (1))
m12BBSlow.1    <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + (0.5*(m12BB.allBetas[,7]))) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + (0.5*(m12BB.allBetas[,8]))) * (1))

m12BBFast5.2   <-  inv_logit(m12BB.allBetas[,1] + m12BB.allBetas[,2] * (2))
m12BBFast55.2  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,4]) + (m12BB.allBetas[,2] + m12BB.allBetas[,6]) * (2))
m12BBSlow5.2   <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5]) * (2))
m12BBSlow55.2  <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + m12BB.allBetas[,7]) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + m12BB.allBetas[,8]) * (2))
m12BBFast.2    <-  inv_logit((m12BB.allBetas[,1] + (m12BB.allBetas[,4])/2) + (m12BB.allBetas[,2] + (m12BB.allBetas[,6])/2) * (2))
m12BBSlow.2    <-  inv_logit((m12BB.allBetas[,1] + m12BB.allBetas[,3] + (0.5*(m12BB.allBetas[,7]))) + (m12BB.allBetas[,2] + m12BB.allBetas[,5] + (0.5*(m12BB.allBetas[,8]))) * (2))


simpContr  <-  list(
  cSimp1   =  m12BBFast5.neg1 - m12BBFast55.neg1,
  cSimp2   =  m12BBSlow5.neg1 - m12BBSlow55.neg1,
  cSimp3   =  m12BBFast5.neg1 - m12BBSlow5.neg1,
  cSimp4   =  m12BBFast55.neg1 - m12BBSlow55.neg1,
  cSimp5   =  m12BBFast.neg1 - m12BBSlow.neg1,
  cSimp6   =  m12BBFast5.0 - m12BBFast55.0,
  cSimp7   =  m12BBSlow5.0 - m12BBSlow55.0,
  cSimp8   =  m12BBFast5.0 - m12BBSlow5.0,
  cSimp9   =  m12BBFast55.0 - m12BBSlow55.0,
  cSimp10  =  m12BBFast.0 - m12BBSlow.0,
  cSimp11  =  m12BBFast5.1 - m12BBFast55.1,
  cSimp12  =  m12BBSlow5.1 - m12BBSlow55.1,
  cSimp13  =  m12BBFast5.1 - m12BBSlow5.1,
  cSimp14  =  m12BBFast55.1 - m12BBSlow55.1,
  cSimp15  =  m12BBFast.1 - m12BBSlow.1,
  cSimp16  =  m12BBFast5.2 - m12BBFast55.2,
  cSimp17  =  m12BBSlow5.2 - m12BBSlow55.2,
  cSimp18  =  m12BBFast5.2 - m12BBSlow5.2,
  cSimp19  =  m12BBFast55.2 - m12BBSlow55.2,
  cSimp20  =  m12BBFast.2 - m12BBSlow.2
)
pval(simpContr[[5]])

pdf(file="./output/figs/NxRate_SimpContrasts_BetaBin.pdf", height=18, width=20)
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
## Main result from a posteriori contrasts for model m12BB:
##
## There is a significant main effect of Rate for the slope of the 
## regression lines (p = 1.000). 
## 
## The Rate x EggPos interaction for the Fast slopes is no longer 
## significant (p = 0.70).
## 
##  There is still no interaction for the Slow treatment. 
##
##  Both the 5cm and 55cm partial regression slopes for the Fast treatment
##  are significantly steeper than the "Slow" slope (p = 0.023, 0.002) 
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
Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z,
                         data = data)

# Mean coefficients
m21BB.betas    <-  m21BB.summ$Mean[1:8]
m21BB.gammas   <-  m21BB.summ$Mean[9:28]

source('R/functions-figures.R')
# Predicted lines for Fast (with pooled slopes) & the overall Slow
# Create plotting objects for each regression line
Fast5.plt   <-  Fast5.plots(m21BB.betas, m21BB.allBetas, m21BB.gammas, Z=Z, data=data)
Fast55.plt  <-  Fast55.plots(m21BB.betas, m21BB.allBetas, m21BB.gammas, Z=Z, data=data)
Fast.plt    <-  Fast.plots(m21BB.betas, m21BB.allBetas, m21BB.gammas, Z=Z, data=data)
Slow.plt    <-  Slow.plots(m21BB.betas, m21BB.allBetas, m21BB.gammas, Z=Z, data=data)

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
lines(Fast.plt$y   ~ Fast.plt$xRaw, col='dodgerblue1', lwd=3)
points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(Fast.plt$yAdj[data$Rate == "Fast"] ~ Fast.plt$xReal[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
proportionalLabel(-0.15, 0.5, expression(paste("Adjusted Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue4','orangered4'),
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
polygon(x=c(Fast.plt$xRaw, rev(Fast.plt$xRaw)), 
        y=c(Fast.plt$CIs$lower, rev(Fast.plt$CIs$upper)), 
        col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
lines(Slow.plt$y   ~ Slow.plt$xRaw, col='orangered1', lwd=3)
lines(Fast.plt$y  ~ Fast.plt$xRaw, col='dodgerblue1', lwd=3)
points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(Fast.plt$yAdj[data$Rate == "Fast"] ~ Fast.plt$xReal[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
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





###########################################################################
###########################################################################
###########################################################################
##  Compare most parsimonious Binomial vs. Beta-Binomial model


# Import Binomial model m12
csvFiles  <-  c('./output/StanFits/NxRate_m12.csv1',
                './output/StanFits/NxRate_m12.csv2',
                './output/StanFits/NxRate_m12.csv3')
m12        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

m12.df    <-  as.data.frame(extract(m12))[,-1]
m12.summ  <-  plyr:::adply(as.matrix(m12.df),2,MCMCsum)




##########################################################################
# Compare models: m12 & m21BB
##########################################################################

##############
# Diagnostics

# Model Results
print(m12, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("beta", "phi"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

print(m12, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("a", "b"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(m12, pars="beta")
plot(m21BB, pars="beta")

pairs(m12, pars="beta")
pairs(m21BB, pars="beta")

rstan::traceplot(m12, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)
rstan::traceplot(m21BB, pars=c("beta", "sigma_gamma"), inc_warmup=FALSE)
dev.off()


# Check posteriors against priors
par(mfrow=c(2,3))
x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m12.df$sigma_gamma), lwd=3, col='dodgerblue1')

x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m12.df[,i]), lwd=3, col=i)
}

x  <-  seq(from=0, to=1, length=500)
plot(dcauchy(x)*20 ~ x, lwd=2, type='l')
lines(density(m21BB.df$sigma_gamma), lwd=3, col='dodgerblue1')


x  <-  seq(from=0, to=100, length=500)
plot(dcauchy(x, scale=100)*25 ~ x, lwd=2, type='l')
lines(density(m21BB.df$phi), lwd=3, col='dodgerblue1')


x  <-  seq(from=-3, to=3, length=500)
plot(dnorm(x, sd=3)*17 ~ x, lwd=2, type='l')
for(i in 1:8) {
  lines(density(m21BB.df[,i]), lwd=3, col=i)
}
dev.off()

##############################
# Posterior Predictive Checks

#  Quick self-consistency check:
#  Plot of simulated data against real data
par(mfrow=c(1,2))
y  <-  as.numeric(m12.df[1,300:419])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m12.df[i,300:419])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 


y  <-  as.numeric(m21BB.df[1,511:630])/data$nEggs
x  <-  data$nFert/data$nEggs
plot(y ~ x, xlim=c(0,1), ylim=c(0,1))

for(i in 2:1000) {
  rm(y)
  y  <-  as.numeric(m21BB.df[i,511:630])/data$nEggs
  points(y ~ jitter(x,factor=500))
}
abline(a=0,b=1, col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
 pdf(file='./output/figs/NxRate_Compare_PPDs.pdf', width=7,height=7)
par(mfrow=c(2,2))
plot(density(m12.df[,420],   adjust=3), xlim=c(-2,9), lwd=3, col='dodgerBlue1', main='min_y_rep (min. num. Successes)')
lines(density(m21BB.df[,631], adjust=3), lwd=3, col=COLS[2])
abline(v=min(data$nFert), lwd=3, col=2)
    legend(
          x       =  9.25,
          y       =  0.325,
          legend  =  c(
                      expression(paste(Mod.~12)),
                      expression(paste(Mod.~"21BB"))),
          lty     =  1,
          lwd     =  3,
          col     =  c(
                       'dodgerblue1',
                       COLS[2]),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )

plot(density(m12.df[,421]), xlim=c(65,105), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
lines(density(m21BB.df[,632]), lwd=3, col=COLS[2])
abline(v=max(data$nFert), lwd=3, col=2)

plot(density(m12.df[,422]), xlim=c(25,33), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
lines(density(m21BB.df[,633]), lwd=3, col=COLS[2])
abline(v=mean(data$nFert), lwd=3, col=2)

plot(density(m12.df[,423]), xlim=c(17,24),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
lines(density(m21BB.df[,634]), lwd=3, col=COLS[2])
abline(v=sd(data$nFert), lwd=3, col=2)

dev.off()


print(m12, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m21BB, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

###########################################################################
###########################################################################
###########################################################################
## Overall Impression Re: Using Binomial vs. Beta-Binomial
##
##  --  There really isn't much to decide between the Binomial and Beta-Binomial
##      models. While the Binomial models tend to under-predict the sd() of
##      the real data, the Beta-Binomial models tend to under-predict the min().
##      So pick your poison. Slighly under-dispersed, or zero-inflated.
##     
##  --  Making the models even more equivocal is the fact that the main
##      inference from the fixed-effects remain the same! The Rate x EggPos
##      interaction is 'more signifiant' in the Binomial models... but it's
##      still only marginal. 
##
##  --  So... my overall impression is that we ought to just stick with the
##      simpler Binomial model... but make it clear in the MS that we compared
##      our results with a Beta-Binomial analysis to see what happens when we
##      account for over-dispersion... but that the results from the models do
##      not change qualitatively.
##
###########################################################################
###########################################################################
###########################################################################
