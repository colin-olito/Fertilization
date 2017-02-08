#/* 
# * Colin Olito. Created 08/02/2017
# * Analysis of 1st flume experiment: N-invest
# * using Beta-Binomial error distribution
# * 
# * NOTES:  This file will fit all the necessary
# * 		Stan models for the analysis of the
# * 		N-invest flume data, and write the
# * 		Stan sample files to ./output/Stanfits
# * 
# * 		These Stan sample files can then be 
# * 		imported, plotted, and analyzed using
# * 		./N_invest_analyze
# *          
# */

rm(list=ls())
################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

#******************
# DEPENDENCIES
source('R/functions.R')

# str(data)


##################
# Import Data Sets

#********************
#  N_invest Data Set
print('Importing N_invest Data Set')
NinvData <- read.csv('data/Ninvest_master.csv', header=TRUE, stringsAsFactors=FALSE)
NinvData <- data.frame(NinvData)

# Convert grouping variables to factors; Correct Dates
NinvData$Run       <-  factor(NinvData$Run)
NinvData$Colony    <-  factor(NinvData$Colony)
NinvData$N         <-  factor(NinvData$N)
NinvData$Lane      <-  factor(NinvData$Lane)
NinvData$nSperm_c  <-  NinvData$nSperm - mean(NinvData$nSperm)
NinvData$Date      <-  dmy(NinvData$Date)
NinvData$nSperm_z  <-  (NinvData$nSperm - mean(NinvData$nSperm))/sd(NinvData$nSperm)


########################################################
########################################################
## FIT THE MODELS
########################################################
########################################################


#  Fixed Effects Model Matrix (Same for all models)
X       <-  model.matrix(~ 1 + nSperm_z, data=NinvData)
X       <-  unname(X)
attr(X,"assign") <- NULL


# Options for the STAN analyses
nChains       = 3
thinSteps     = 1
nIter         = 2000 #for each chain
burnInSteps   = nIter / 2
nSavedSteps  = (nIter/2)*nChains
(nSavedSteps)

########################################################
#  m1BB: Simple Logistic regression 
#      --  FertRate ~ nSperm_z
#	   --  Complete pooling of observations
########################################################

# create data.list
data.list  <-  list(N   =  nrow(NinvData),
                    P   =  ncol(X), 
                    nT  = NinvData$nEggs,
                    nS  =  NinvData$nFert,
                    X   =  X
                   )

# Call to STAN
m1BB <- stan(data         =  data.list,
           seed         =  123456789,
           file         =  './Stan/mat-BetaBin.stan',
           sample_file  =  './output/StanFits/N_invest_m1BB.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m1BB")
system("notify-send \"STAN has finished fitting model m1BB\"")

# garbage collection
rm(data.list)
rm(m1BB)

########################################################
#  m2BB: Logistic regression w/ Random Intercept ~ Run 
#      --  FertRate ~ nSperm_z + (1 | Run)
########################################################

# 'random effects'
Z  <-  unname(model.matrix(~ NinvData$Run -1, data=NinvData))
attr(Z,"assign") <- NULL
# str(Z)
# head(Z)

# create data.list
data.list  <-  list(N   =  nrow(NinvData),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  NinvData$nEggs,
                    nS  =  NinvData$nFert,
                    X   =  X,
                    Z   =  Z
                   )

# Call to STAN
m2BB <- stan(data         =  data.list,
           seed         =  234567891,
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/N_invest_m2BB.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2BB")
system("notify-send \"N-invest analysis model fitting:\n
					  STAN has finished fitting model m2BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m2BB)


########################################################
#  m3BB: Logistic Mixed Effects Regression
#       --  FertRate ~ nSperm_z + ...
## -- random intercept for RUN
## -- random slopes for Run x nSperm
## -- Estimate covariance matrix
########################################################

Z  <-  unname(model.matrix(~ NinvData$Run *nSperm_z -1, data=NinvData))
attr(Z,"assign") <- NULL
# str(Z)
# Z[1:20,]

# create data.list
data.list  <-  list(N    =  nrow(NinvData),
                    P    =  ncol(X),
                    J    =  max(as.numeric(NinvData$Run)),
                    K    =  ncol(Z),
                    grp  =  as.numeric(NinvData$Run),
                    nT   =  NinvData$nEggs,
                    nS   =  NinvData$nFert,
                    X    =  X,
                    Z    =  Z
                   )

# inits  <-  list(
#   list(L_run    =  matrix(runif(16^2,-0.5,0.5), nrow=16,ncol=16)),
#   list(tau_run  =  runif(16,0.1,2)),
#   list(u        =  matrix(runif((8*16),-1,1),nrow=8,ncol=16)),
#   list(beta     =  runif(2,-1,1)),
#   list(sigma_y  =  runif(1,0.1,2))
#   )

# Call to STAN
m3BB <- stan(data         =  data.list,
           seed         =  567891234,
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/N_invest_m3BB.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE,
           control      =  list(adapt_delta = 0.95) # default adapt_delta of 0.8 threw many divergent transitions. 
          )

# message
message("STAN has finished fitting model m3BB")
system("notify-send \"STAN has finished fitting model m3BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m3BB)
