#/* 
# * Colin Olito. Created 28/05/2016.
# * 
# * NOTES: 2nd Flume Experiment
# *         crossing N x Rate; with 2 egg patches
# *          Fitting using Beta-Binomial error dist.
# *          
# */

## Beta-Binomial Mixed Effects Regression Analysis

rm(list=ls())
################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

#******************
# DEPENDENCIES
source('R/dependencies.R')

#*******************
# Import Data
data <- read.csv('data/NxRate_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)
head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Rate      <-  factor(data$Rate)
data$EggPos    <-  factor(data$EggPos)
data$Lane      <-  factor(data$Lane)
data$Date      <-  dmy(data$Date)

# Centered and rescaled nsperm variable for easier estimation
nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)






##################################################################
##################################################################
##  Fit all models
##################################################################
##################################################################

##  Options for all analyses
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)



#  Fixed Effects Model Matrix (Same for all models)
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

################################################
## Model #1: MAXIMAL MODEL W/ COVARIANCE MATRIX
################################################




################################################
## Model #2: Maximal model W/o COVARIANCE MATRIX
################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z) #+
#                                data$Run:data$Rate) 
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

##  Assemble data for stan
data.list  <-  list(N    =  nrow(data),
                    P    =  ncol(X),
                    K    =  ncol(Z),
                    nT   =  data$nEggs - data$nControlFert,
                    nS   =  data$nFert,
                    X    =  X,
                    Z    =  Z
                   )

## Call to STAN
m2 <- stan(data     =  data.list,
             file     =  './Stan/mat-BetaBin-1Z.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )

# Model Results
m2.df    <-  as.data.frame(extract(m2))
mcmc.m2  <-  as.mcmc(m2)
m2.mcmc  <-  rstan:::as.mcmc.list.stanfit(m2)
m2.summ  <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)

##  LOO Log-likelihood for model selection
m2LL     <-  extract_log_lik(m2, parameter_name = "log_lik")
m2Loo    <-  loo(m2LL)
m2WAIC   <-  waic(m2LL)

## Garbage Collection
rm(Z)
rm(data.list)

##  Notification
system('notify-send "Sampling for m2 complete"')
print("Sampling for m2 complete")




print(m2, c("beta","sigma_gamma", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
plot(m2.mcmc, ask=TRUE)

m2.summ$X1
m2yhat  <-  inv_logit(m2.summ$Mean[271:390])
m2.resids  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m2yhat)/sd(m2yhat)
m2.resids_z  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m2yhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m2yhat))


m2.coef   <-  m2.summ$Mean[1:8]

m2Fast5   <-  inv_logit(m2.coef[1] + m2.coef[2] * nSperm_z)
m2Fast55  <-  inv_logit((m2.coef[1] + m2.coef[4]) + (m2.coef[2] + m2.coef[6]) * nSperm_z)
m2Slow5   <-  inv_logit((m2.coef[1] + m2.coef[3]) + (m2.coef[2] + m2.coef[5]) * nSperm_z)
m2Slow55  <-  inv_logit((m2.coef[1] + m2.coef[3] + m2.coef[7]) + (m2.coef[2] + m2.coef[5] + m2.coef[8]) * nSperm_z)
m2Fast    <-  inv_logit((m2.coef[1] + (m2.coef[4])/2) + (m2.coef[2] + (m2.coef[6])/2) * nSperm_z)
m2Slow    <-  inv_logit((m2.coef[1] + m2.coef[3] + (0.5*(m2.coef[7]))) + (m2.coef[2] + m2.coef[5] + (0.5*(m2.coef[8]))) * nSperm_z)

m2.low   <-  m2.summ$lower[1:8]
m2.hi   <-  m2.summ$upper[1:8]
m2Fast.low    <-  inv_logit((m2.low[1] + (m2.coef[4])/2) + (m2.coef[2] + (m2.coef[6])/2) * nSperm_z)
m2Slow.low    <-  inv_logit((m2.low[1] + m2.coef[3] + (0.5*(m2.coef[7]))) + (m2.coef[2] + m2.coef[5] + (0.5*(m2.coef[8]))) * nSperm_z)
m2Fast.hi     <-  inv_logit((m2.hi[1] + (m2.coef[4])/2) + (m2.coef[2] + (m2.coef[6])/2) * nSperm_z)
m2Slow.hi     <-  inv_logit((m2.hi[1] + m2.coef[3] + (0.5*(m2.coef[7]))) + (m2.coef[2] + m2.coef[5] + (0.5*(m2.coef[8]))) * nSperm_z)





par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m2Fast5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m2Fast55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m2Slow5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(m2Slow55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
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






##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m2.df, 1, function(x, data, nSperm_z){
#     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue1',0.01))
#     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered1',0.01))
# }, data=data, nSperm_z=nSperm_z)

polygon(c(data$nSperm[order(nSperm_z)], rev(data$nSperm[order(nSperm_z)])), 
        c(m2Slow.low[order(nSperm_z)], rev(m2Slow.hi[order(nSperm_z)])), col=transparentColor('orangered2', 0.5), border='orangered2')
polygon(c(data$nSperm[order(nSperm_z)], rev(data$nSperm[order(nSperm_z)])), 
        c(m2Fast.low[order(nSperm_z)], rev(m2Fast.hi[order(nSperm_z)])), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
lines(m2Slow[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
lines(m2Fast[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)

points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )

