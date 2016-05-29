#/* 
# * Colin Olito. Created 26/04/2016.
# * 
# * NOTES: 2nd Flume Experiment
# *         crossing N x Rate; with 2 egg patches
# *          
# *          
# */

## Logistic Mixed Effects Regression Analysis

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

#  Random Effects Model Matrix
#  Not enough observations to include EggPos ixns! Results in only 
#  3 observations per run.
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z +
                                data$Run:data$Rate +
                                data$Run:nSperm_z:data$Rate)
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
m1 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  3,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )

# Model Results
m1.df    <-  as.data.frame(extract(m1))
mcmc.m1  <-  as.mcmc(m1)
m1.mcmc  <-  rstan:::as.mcmc.list.stanfit(m1)
m1.summ  <-  plyr:::adply(as.matrix(m1.df),2,MCMCsum)

#  LOO Log-likelihood for model selection
m1LL  <-  extract_log_lik(m1, parameter_name = "log_lik")
m1Loo    <-  loo(m1LL)
m1WAIC   <-  waic(m1LL)

## Garbage Collection
rm(Z)
rm(data.list)

##  Notification
system('notify-send "Sampling for m1 complete"')
print("Sampling for m1 complete")


################################################
## Model #2: Remove Run x Rate x nSperm ixn 
##			 W/ COVARIANCE MATRIX
################################################

##  Fixed Effects Model Matrix Remains the same
# X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
# Xnames  <-  dimnames(X)[[2]]
# X       <-  unname(X)
# attr(X,"assign") <- NULL
# str(X)
# head(X)

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z +
                                data$Run:data$Rate) 
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
             file     =  './Stan/mat-logistic-1Z-cov.stan',
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

################################################
## Model #3: Random intercept & slopes x Run
##			 W/ COVARIANCE MATRIX
################################################

##  Fixed Effects Model Matrix Remains the same
# X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
# Xnames  <-  dimnames(X)[[2]]
# X       <-  unname(X)
# attr(X,"assign") <- NULL
# str(X)
# head(X)

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z)
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
m3 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )

## Model Results
m3.df    <-  as.data.frame(extract(m3))
mcmc.m3  <-  as.mcmc(m3)
m3.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3)
m3.summ  <-  plyr:::adply(as.matrix(m3.df),2,MCMCsum)

##  LOO Log-likelihood for model selection
m3LL     <-  extract_log_lik(m3, parameter_name = "log_lik")
m3Loo    <-  loo(m3LL)
m3WAIC   <-  waic(m3LL)

## Garbage Collection
rm(Z)
rm(data.list)

##  Notification
system('notify-send "Sampling for m3 complete"')
print("Sampling for m3 complete")


################################################
## Model #3a: Random intercept & slopes x Run
##			 W/ NO COVARIANCE MATRIX
################################################

##  Fixed Effects Model Matrix Remains the same
# X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
# Xnames  <-  dimnames(X)[[2]]
# X       <-  unname(X)
# attr(X,"assign") <- NULL
# str(X)
# head(X)

##  Random Slopes Model Matrix
Z0       <-  model.matrix(~ -1 + data$Run, data=data)
Z0names  <-  dimnames(Z0)[[2]]
Z0       <-  unname(Z0)
attr(Z0,"assign") <- NULL

##  Random Intercepts Model Matrix
Z1  <-  model.matrix(~ -1 +  data$Run:nSperm_z , data=data)
Z1       <-  unname(Z1)
attr(Z1,"assign") <- NULL

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z0),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    Z0  =  Z0,
                    Z1  =  Z1
                   )

## Call to STAN
m3a <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-2Z.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
             )

# Model Results
m3a.df    <-  as.data.frame(extract(m3a))
mcmc.m3a  <-  as.mcmc(m3a)
m3a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3a)
m3a.summ  <-  plyr:::adply(as.matrix(m3a.df),2,MCMCsum)

#  LOO Log-likelihood for model selection
m3aLL  <-  extract_log_lik(m3a, parameter_name = "log_lik")
m3aLoo    <-  loo(m3aLL)
m3aWAIC   <-  waic(m3aLL)

## Garbage Collection
rm(Z0)
rm(Z1)
rm(data.list)

##  Notification
system('notify-send "Sampling for m3a complete"')
print("Sampling for m3a complete")


################################################
## Model #4: Random intercepts
##			 W/ COVARIANCE MATRIX
################################################

##  Fixed Effects Model Matrix Remains the same
# X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
# Xnames  <-  dimnames(X)[[2]]
# X       <-  unname(X)
# attr(X,"assign") <- NULL
# str(X)
# head(X)

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run)
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
m4 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )

## Model Results
m4.df    <-  as.data.frame(extract(m4))
mcmc.m4  <-  as.mcmc(m4)
m4.mcmc  <-  rstan:::as.mcmc.list.stanfit(m4)
m4.summ  <-  plyr:::adply(as.matrix(m4.df),2,MCMCsum)

##  LOO Log-likelihood for model selection
m4LL     <-  extract_log_lik(m4, parameter_name = "log_lik")
m4Loo    <-  loo(m4LL)
m4WAIC   <-  waic(m4LL)

## Garbage Collection
rm(Z)
rm(data.list)

##  Notification
system('notify-send "Sampling for m4 complete"')
print("Sampling for m4 complete")


################################################
## Model #4a: Random intercepts
##			 W/ NO COVARIANCE MATRIX
################################################

##  Fixed Effects Model Matrix Remains the same
# X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
# Xnames  <-  dimnames(X)[[2]]
# X       <-  unname(X)
# attr(X,"assign") <- NULL
# str(X)
# head(X)

Z       <-  model.matrix(~ -1 + data$Run, data=data)
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m4a <- stan(data     =  data.list,
            file     =  './Stan/mat-logistic-1Z.stan',
            chains   =  nChains,
            iter     =  numSavedSteps,
            thin     =  thinSteps,
            save_dso =  TRUE
           )

# Model Results
m4a.df    <-  as.data.frame(extract(m4a))
mcmc.m4a  <-  as.mcmc(m4a)
m4a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m4a)
m4a.summ  <-  plyr:::adply(as.matrix(m4a.df),2,MCMCsum)

#  LOO Log-likelihood for model selection
m4aLL  <-  extract_log_lik(m4a, parameter_name = "log_lik")
m4aLoo    <-  loo(m4aLL)
m4aWAIC   <-  waic(m4aLL)

## Garbage Collection
rm(Z)
rm(data.list)

##  Notification
system('notify-send "Sampling for m4a complete"')
print("Sampling for m4a complete")


################################################
## Model #5: Simple Logistic Regression
################################################

##  Fixed Effects Model Matrix Remains the same
# X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
# Xnames  <-  dimnames(X)[[2]]
# X       <-  unname(X)
# attr(X,"assign") <- NULL
# str(X)
# head(X)

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X
                   )

## Call to STAN
m5 <- stan(data     =  data.list,
                 file     =  './Stan/mat-logistic-reg.stan',
                 chains   =  nChains,
                 iter     =  numSavedSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                )

# Model Results
m5.df    <-  as.data.frame(extract(m5))
mcmc.m5  <-  as.mcmc(m5)
m5.mcmc  <-  rstan:::as.mcmc.list.stanfit(m5)
m5.summ  <-  plyr:::adply(as.matrix(m5.df),2,MCMCsum)

#  LOO Log-likelihood for model selection
m5LL     <-  extract_log_lik(m5, parameter_name = "log_lik")
m5Loo    <-  loo(m5LL)
m5WAIC   <-  waic(m5LL)

## Garbage Collection
rm(data.list)

##  Notification
system('notify-send "Sampling for m5 complete"')
print("Sampling for m5 complete")
