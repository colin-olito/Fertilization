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
                    J    =  max(as.numeric(as.factor(data$Run))),
                    K    =  ncol(Z),
                    grp  =  as.numeric(as.factor(data$Run)),
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
