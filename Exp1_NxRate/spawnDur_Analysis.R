#/* 
# * Colin Olito. Created 26/01/2017
# * Analysis of Spawning Duration Experiment
# * 
# * NOTES:  This file contains all the necessary
# * 		code to read in the STAN sample files
# * 		and perform posterior predictive checks
# *         model selection, and produce regression
# *         plots for the spawning duration data.
# * 
# */

rm(list=ls())
#################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

###############
# DEPENDENCIES
# source('./loadData.R')  # This is a very small analysis, so I load the data separately here,
                          # rather than running the whole loadData.R script, which also loadhistory
                          # all the data from the other analyses. However, loadData.R will also
                          # load the spawning duration data.
source('./R/functions.R')


####################################
# Import Spawning Duration Data Set
print('Importing inSitu Spawning Duration Data Set')
spwnData <- read.csv('data/inSitu-SpawnDuration.csv', header=TRUE, stringsAsFactors=FALSE)
spwnData <- data.frame(spwnData)

# Convert grouping variables to factors; Correct Dates
spwnData$Ind    <-  factor(spwnData$Ind)
spwnData$Run    <-  factor(spwnData$Run)
spwnData$Start  <-  hms(spwnData$Start)
spwnData$End    <-  hms(spwnData$End)
spwnData$Dur    <-  hms(spwnData$Dur)
# Check that duration column is equivalent to End - Start.
# as.numeric(seconds(spwnData$Dur)) == as.numeric(seconds(spwnData$End - spwnData$Start))
spwnData$start  <-  as.numeric(seconds(spwnData$Start))
spwnData$end    <-  as.numeric(seconds(spwnData$End))
spwnData$dur    <-  as.numeric(seconds(spwnData$Dur))

############################################
# Import Spawning Duration model results 
# from the stan sample_files for further analysis 
############################################
csvFiles  <-  c('./output/StanFits/spawnDur_m1.csv1',
                './output/StanFits/spawnDur_m1.csv2',
                './output/StanFits/spawnDur_m1.csv3')
SDm1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/spawnDur_m2.csv1',
                './output/StanFits/spawnDur_m2.csv2',
                './output/StanFits/spawnDur_m2.csv3')
SDm2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/spawnDur_m3.csv1',
                './output/StanFits/spawnDur_m3.csv2',
                './output/StanFits/spawnDur_m3.csv3')
SDm3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/spawnDur_m3.csv1',
                './output/StanFits/spawnDur_m3.csv2',
                './output/StanFits/spawnDur_m3.csv3')
SDm4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)




##########################################################################
# Exploratory Analysis
##########################################################################

summary(spwnData$dur)
spwnData$durCentered  <-  spwnData$dur - mean(spwnData$dur)
spwnData$dur_z        <-  (spwnData$dur - mean(spwnData$dur))/ sd(spwnData$dur)


par(mfrow=c(2,2))
plot(density(spwnData$dur, adjust=0.5), lwd=3, col='dodgerBlue1')
abline(v=mean(spwnData$dur), lwd=3, col=2)
plot(density(spwnData$durCentered), lwd=3, col='dodgerBlue1')
abline(v=mean(spwnData$durCentered), lwd=3, col=2)
plot(density(spwnData$dur_z), lwd=3, col='dodgerBlue1')
abline(v=mean(spwnData$dur_z), lwd=3, col=2)
hist(spwnData$dur_z, breaks=50, xlim=c(-3,3))
abline(v=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))

##  Residual Plots
par(mfrow=c(2,2))
hist(spwnData$durCentered, breaks=40)
abline(v=c(-2,2), lty=2)
plot(spwnData$durCentered ~ as.numeric(spwnData$Ind), ylim=c(-300,300))
abline(h=c((-2*sd(spwnData$durCentered)),0,(2*sd(spwnData$durCentered))), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(spwnData$durCentered ~ seq_along(spwnData$durCentered), ylim=c(-300,300))
abline(h=c((-2*sd(spwnData$durCentered)),0,(2*sd(spwnData$durCentered))), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(spwnData$durCentered)
qqline(spwnData$durCentered, col = 2)


##  Residual Plots
par(mfrow=c(2,2))
hist(spwnData$dur_z, breaks=40)
abline(v=c(-2,2), lty=2)
plot(spwnData$dur_z ~ as.numeric(spwnData$Ind), ylim=c(-3,3))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(spwnData$dur_z ~ seq_along(spwnData$dur_z), ylim=c(-3,3))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(spwnData$dur_z)
qqline(spwnData$dur_z, col = 2)

##  Have a look at among-run variation
COLS  <-  c("#860885",
            "#006a09",
            "#028fff",
            "#ffc544",
            "#331061",
            "#ba2423"
            )

plot(density(spwnData$dur), lwd=3, ylim=c(0,0.01), col='dodgerBlue1')
for (i in 1:5) {
	lines(density(spwnData$dur[spwnData$Run == i]), lwd=2, col=COLS[(i)])
}
abline(v=mean(spwnData$dur), lwd=3, col=2)

########################################################
##  Normal error distribution will probably work
##  just fine for this analysis. Most of the runs
##  look like they have homogenous means/variances,
##  except for run #1... Will run a few mixed models
##  to account for among-run variation, and also see
##  what an exponential error distribution looks like.
########################################################


##########################################################################
# Model: SDm1
##########################################################################

##############
# Diagnostics

# Model Results
SDm1.df     <-  as.data.frame(extract(SDm1))  
SDm1.summ   <-  plyr:::adply(as.matrix(SDm1.df),2,MCMCsum)[-1,]

print(SDm1, c("beta", "gamma", "sigma_gamma", "sigma_y"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(SDm1, pars=c("beta", "gamma"))
pairs(SDm1, pars="gamma")
rstan::traceplot(SDm1, c("beta", "gamma", "sigma_gamma", "sigma_y"), inc_warmup=FALSE)

dev.off()

##############################
# Posterior Predictive Checks
##############################

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(SDm1.df[1, 62:87])
x  <-  spwnData$dur
plot(y ~ x, xlim = c((min(spwnData$dur)*0.95), (max(spwnData$dur)*1.05)),
	        ylim = c(min(SDm1.df[,62:87]), max(SDm1.df[,62:87])))
for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(SDm1.df[i, 62:87])
  points(y ~ jitter(x,factor=500))
}
abline(h=mean(spwnData$dur), col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(SDm1.df[,88], adjust=1), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(spwnData$dur), lwd=3, col=2)

plot(density(SDm1.df[,89]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(spwnData$dur), lwd=3, col=2)

plot(density(SDm1.df[,90]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(spwnData$dur), lwd=3, col=2)

plot(density(SDm1.df[,91]), xlim=c(min(SDm1.df[,91],sd(spwnData$dur)),max(SDm1.df[,91],sd(spwnData$dur))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(spwnData$dur), lwd=3, col=2)

print(SDm1, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


#########################################
# LOO Log-likelihood for model selection
#########################################
SDm1LL     <-  extract_log_lik(SDm1, parameter_name = "log_lik")
SDm1Loo    <-  loo(SDm1LL)
SDm1WAIC   <-  waic(SDm1LL)

#####################################################
# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data
#####################################################
X2dataSD1   <-  as.numeric(SDm1.df[,174])
X2simSD1    <-  as.numeric(SDm1.df[,175])









##########################################################################
# Model: SDm2
##########################################################################

##############
# Diagnostics

# Model Results
SDm2.df     <-  as.data.frame(extract(SDm2))  
SDm2.summ   <-  plyr:::adply(as.matrix(SDm2.df),2,MCMCsum)[-1,]

print(SDm2, c("beta", "sigma_y"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(SDm2, pars=c("beta", "sigma_y"))
rstan::traceplot(SDm2, c("beta", "sigma_y"), inc_warmup=FALSE)

dev.off()

##############################
# Posterior Predictive Checks
##############################

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(SDm2.df[1, 56:81])
x  <-  spwnData$dur
plot(y ~ x, xlim = c((min(spwnData$dur)*0.95), (max(spwnData$dur)*1.05)),
	        ylim = c(min(SDm2.df[,56:81]), max(SDm2.df[,56:81])))
for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(SDm2.df[i, 56:81])
  points(y ~ jitter(x,factor=500))
}
abline(h=mean(spwnData$dur), col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(SDm2.df[,82], adjust=1), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(spwnData$dur), lwd=3, col=2)

plot(density(SDm2.df[,83]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(spwnData$dur), lwd=3, col=2)

plot(density(SDm2.df[,84]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(spwnData$dur), lwd=3, col=2)

plot(density(SDm2.df[,85]), xlim=c(min(SDm2.df[,85],sd(spwnData$dur)),max(SDm2.df[,85],sd(spwnData$dur))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(spwnData$dur), lwd=3, col=2)

print(SDm2, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


#########################################
# LOO Log-likelihood for model selection
#########################################
SDm2LL     <-  extract_log_lik(SDm2, parameter_name = "log_lik")
SDm2Loo    <-  loo(SDm2LL)
SDm2WAIC   <-  waic(SDm2LL)

#####################################################
# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data
#####################################################
X2dataSD2   <-  as.numeric(SDm2.df[,168])
X2simSD2    <-  as.numeric(SDm2.df[,169])













##########################################################################
# Model: SDm3
##########################################################################

##############
# Diagnostics

# Model Results
SDm3.df     <-  as.data.frame(extract(SDm3))  
SDm3.summ   <-  plyr:::adply(as.matrix(SDm3.df),2,MCMCsum)[-1,]

print(SDm3, c("beta", "gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

# Simple Diagnostic Plots
plot(SDm3, pars=c("beta", "gamma"))
rstan::traceplot(SDm3, c("beta", "gamma"), inc_warmup=FALSE)

dev.off()

##############################
# Posterior Predictive Checks
##############################

#  Quick self-consistency check:
#  Plot of simulated data against real data

y  <-  as.numeric(1/SDm3.df[1, 61:86])
x  <-  spwnData$dur
plot(y ~ x, xlim = c((1/min(spwnData$dur)*0.95), (1/max(spwnData$dur)*1.05)),
          ylim = c(1/min(SDm3.df[,61:86]), 1/max(SDm3.df[,61:86])))
for(i in 2:500) {
  rm(y)
  y  <-  as.numeric(SDm3.df[i, 61:86])
  points(y ~ jitter(x,factor=500))
}
abline(h=mean(spwnData$dur), col=2, lwd=3) 

# Density plots of min, max, mean, sd
#  of replicated data, benchmarked with
#  calculated values for real data
#  Find associated p-values in m4.summ
par(mfrow=c(2,2))
plot(density(SDm3.df[,87], adjust=1), lwd=3, col='dodgerBlue3', main='min_y_rep (min. num. Successes)')
abline(v=min(spwnData$dur), lwd=3, col=2)

plot(density(SDm3.df[,88]), lwd=3, col='dodgerBlue3', main='max_y_rep (max. num. Successes)')
abline(v=max(spwnData$dur), lwd=3, col=2)

plot(density(SDm3.df[,89]), lwd=3, col='dodgerBlue3', main='mean_y_rep (mean num. Successes)')
abline(v=mean(spwnData$dur), lwd=3, col=2)

plot(density(SDm3.df[,90]), xlim=c(min(SDm3.df[,90],sd(spwnData$dur)),max(SDm3.df[,90],sd(spwnData$dur))),
 lwd=3, col='dodgerBlue3', main='sd_y_rep (sd num. Successes)')
abline(v=sd(spwnData$dur), lwd=3, col=2)

print(SDm3, c("p_min","p_max","p_mean","p_sd"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


#########################################
# LOO Log-likelihood for model selection
#########################################
SDm3LL     <-  extract_log_lik(SDm3, parameter_name = "log_lik")
SDm3Loo    <-  loo(SDm3LL)
SDm3WAIC   <-  waic(SDm3LL)

#####################################################
# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data
#####################################################
X2dataSD3   <-  as.numeric(SDm3.df[,173])
X2simSD3    <-  as.numeric(SDm3.df[,174])

















##########################################################################################
##########################################################################################
# Model selection using LOO cross-validation
##########################################################################################
##########################################################################################

##  Overall comparison of all models
looDiff   <-  compare(SDm1Loo, SDm2Loo)
#waicDiff  <-  compare(m1WAIC, m2WAIC, m3WAIC, m3aWAIC, m4WAIC, m4aWAIC, m5WAIC)

print(looDiff, digits=4)
#print(waicDiff, digits=4)

# LOO Results Summary Table
# Check sample sizes
    if(length(SDm1Loo$pointwise[,"elpd_loo"]) != length(SDm1Loo$pointwise[,"elpd_loo"]))
         stop("sample sizes (chain lengths) differ between models")
    n  <-  length(SDm1Loo$pointwise[,"elpd_loo"])
# Pairwise elpd_loo differences for all models
    elpd_pair  <-  SDm1Loo$elpd_loo - SDm2Loo$elpd_loo
# Standard Error of pairwise elpd_loo differences
    selooDiff  <-  sqrt(n * var(SDm1Loo$pointwise[,"elpd_loo"]  - SDm3Loo$pointwise[,"elpd_loo"]))
# P-values for pairwise elpd_loo differences
    pDiff  <-  as.numeric(rounded(2*pnorm(-abs((elpd_pair - 0)/selooDiff)), 3))
    LooDiff  <-  cbind(elpd_pair, selooDiff, pDiff)
(LooDiff)


########################################################
#  MAIN CONCLUSION
#
#  Perhaps not surprisingly, accounting for among-run
#  variation in spawning time did not significantly 
#  improve model fit. Thus, there is really no good
#  reason to provide a Bayesian estimation of the mean
#  and s.e. for this analysis.  Can stick with basic
#  summary(spwnData$dur).
########################################################




# Plot Range for PPC plots
plotRange  <-  c(0,max(c(X2dataSD2, X2simSD2)))

# pdf(file='./output/figs/NxRate_X2Discrepancy.pdf', width=7,height=7)
par(omi=rep(0.3, 4))
plot(X2simSD2 ~ X2dataSD2, 
    xlab=expression(paste(chi^2~discrepancy~of~observed~data)), ylab=expression(paste(chi^2~discrepancy~of~simulated~data)), 
    main=expression(paste(Posterior~predictive~check:~chi^2~"discrepancy")),
    type='n', axes=FALSE, xlim=plotRange, ylim=plotRange, xpd=NA)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(X2simSD2 ~ X2dataSD2, pch=21, 
        bg=transparentColor('dodgerBlue3', 0.1),
        col=transparentColor('dodgerBlue3', 0.3), cex=1.1)
points(X2simSD1 ~ X2dataSD1, pch=21, 
        bg=transparentColor('orangered1', 0.1),
        col=transparentColor('orangered1', 0.3), cex=1.1)
abline(a=0,b=1, lwd=2) 
axis(2, las=1)
axis(1)
#    legend(
#          x       =  usr[2]*0.17,
#          y       =  usr[4],
#          legend  =  c(
#                      expression(paste(Model~1)),
#                      expression(paste(Model~2)),
#                      expression(paste(Model~3)),
#                      expression(paste(Model~"3a")),
#                      expression(paste(Model~4)),
#                      expression(paste(Model~"4a")),
#                      expression(paste(Model~5))),
#          pch     =  21,
#          pt.bg   =  c(transparentColor(COLS[1],0.5), 
#                       transparentColor(COLS[2],0.5),
#                       transparentColor(COLS[3],0.5),
#                       transparentColor(COLS[7],0.5),
#                       transparentColor(COLS[5],0.5),
#                       transparentColor(COLS[6],0.5),
#                       transparentColor(COLS[4],0.5)),
#          col     =  c(transparentColor(COLS[1], 0.7),
#                       transparentColor(COLS[2], 0.7), 
#                       transparentColor(COLS[3], 0.7), 
#                       transparentColor(COLS[7], 0.7), 
#                       transparentColor(COLS[5], 0.7), 
#                       transparentColor(COLS[6], 0.7), 
#                       transparentColor(COLS[4], 0.7)),
#          cex     =  1,
#          xjust   =  1,
#          yjust   =  1,
#          bty     =  'n',
#          border  =  NA
#    )
# dev.off()
