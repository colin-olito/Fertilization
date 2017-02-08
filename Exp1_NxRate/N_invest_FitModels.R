#/* 
# * Colin Olito. Created 12/12/2016
# * Analysis of 1st flume experiment: N-invest
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


##################
# Import Data

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
#  m1: Simple Logistic regression 
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
m1 <- stan(data         =  data.list,
           seed         =  123456789,
           file         =  './Stan/mat-logistic-reg.stan',
           sample_file  =  './output/StanFits/N_invest_m1.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m1")
system("notify-send \"STAN has finished fitting model m1\"")

# garbage collection
rm(data.list)
rm(m1)

########################################################
#  m2: Logistic regression w/ Random Intercept ~ Run 
#      --  FertRate ~ nSperm_z + (1 | Run)
########################################################

# 'random effects'
Z  <-  unname(model.matrix(~ -1 + Run, data=NinvData))
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
m2 <- stan(data         =  data.list,
           seed         =  234567891,
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/N_invest_m2.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2")
system("notify-send \"N-invest analysis model fitting:\n
					  STAN has finished fitting model m2\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m2)


########################################################
#  m3: Logistic Mixed Effects Regression
#       --  FertRate ~ nSperm_z + ...
#       --  random intercept for RUN
#       --  random slopes for Run x nSperm
########################################################

Z  <-  unname(model.matrix(~ -1 + Run +
                                  Run:nSperm_z , data=NinvData))
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
	m3 <- stan(data         =  data.list,
	           seed         =  345678912,
	           file         =  './Stan/mat-logistic-1Z.stan',
	           sample_file  =  './output/StanFits/N_invest_m3.csv',
	           chains       =  nChains,
	           warmup       =  burnInSteps,
	           iter         =  nIter,
	           thin         =  thinSteps,
	           save_dso     =  TRUE,
	           control      =  list(adapt_delta = 0.99) # default adapt_delta value threw ~60 divergent transitions.
	          )

# message
message("STAN has finished fitting model m3")
system("notify-send \"STAN has finished fitting model m3\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m3)



########################################################
#  m1BB: Simple Logistic regression 
#      --  FertRate ~ nSperm_z
#    --  Complete pooling of observations
########################################################

# create data.list
data.list  <-  list(N   =  nrow(NinvData),
                    P   =  ncol(X), 
                    nT  = NinvData$nEggs,
                    nS  =  NinvData$nFert,
                    X   =  X
                   )

# Call to STAN
m1BB <- stan(data       =  data.list,
           seed         =  56789123,
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
Z  <-  unname(model.matrix(~ -1 + Run, data=NinvData))
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
           seed         =  678912345,
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

Z  <-  unname(model.matrix(~ -1 + Run +
                                  Run:nSperm_z, data=NinvData))
attr(Z,"assign") <- NULL

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

# Call to STAN
m3BB <- stan(data         =  data.list,
             seed         =  789123456,
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
