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
source('R/dependencies.R')

#*******************
# Import Data
data <- read.csv('data/Ninvest_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)
# head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Lane      <-  factor(data$Lane)
data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
data$Date      <-  dmy(data$Date)
data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

# str(data)



########################################################
########################################################
## FIT THE MODELS
########################################################
########################################################

########################################################
#  m1: Simple Logistic regression 
#      --  FertRate ~ nSperm_z
#	   --  Complete pooling of observations
########################################################


# model matrix
X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
# str(X)
# head(X)

# create data.list
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    X   =  X
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 1
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
#nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# Call to STAN
m1 <- stan(data         =  data.list,
           seed         =  123456789,
           file         =  './Stan/mat-logistic-reg.stan',
           sample_file  =  './output/StanFits/N_invest_m1.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  numSavedSteps,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m1")
system("notify-send \"STAN has finished fitting model m1\"")

# garbage collection
rm(X)
rm(data.list)
rm(m1)

########################################################
#  m2: Logistic regression w/ Random Intercept ~ Run 
#      --  FertRate ~ nSperm_z + (1 | Run)
########################################################

# model matrices
# 'fixed effects'
X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
# str(X)
# head(X)

# 'random effects'
Z  <-  unname(model.matrix(~ data$Run -1, data=data))
attr(Z,"assign") <- NULL
# str(Z)
# head(Z)

# create data.list
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    X   =  X,
                    Z   =  Z
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 1
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# Call to STAN
m2 <- stan(data         =  data.list,
           seed         =  234567891,
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/N_invest_m2.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  numSavedSteps,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2")
system("notify-send \"N-invest analysis model fitting:\n
					  STAN has finished fitting model m2\"")

# garbage collection
rm(X)
rm(Z)
rm(data.list)
rm(m2)


########################################################
#  m2b: Logistic Mixed Effects Regression w/ random intercept for RUN
#       --  FertRate ~ nSperm_z + (1 | Run)
#       --  Alternative cell-mean model specification
########################################################

# model matrices
X  <-  unname(model.matrix(~ nSperm_z -1, data=data))
attr(X,"assign") <- NULL
# str(X)
# head(X)

Z  <-  unname(model.matrix(~ data$Run -1, data=data))
attr(Z,"assign") <- NULL
# str(Z)
# head(Z)

# create data.list
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    X   =  X,
                    Z   =  Z
                   )


#  Options for the analysis
nChains        = 4
thinSteps      = 1
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
# nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# Call to STAN
m2b <- stan(data         =  data.list,
           seed         =  345678912,
           file         =  './Stan/mat-logistic-1Z-cellmean.stan',
           sample_file  =  './output/StanFits/N_invest_m2b.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  numSavedSteps,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2b")
system("notify-send \"STAN has finished fitting model m2b\"")

# garbage collection
rm(X)
rm(Z)
rm(data.list)
rm(m2b)



########################################################
#  m3: Logistic Mixed Effects Regression
#       --  FertRate ~ nSperm_z + ...
#       --  random intercept for RUN
#       --  random slopes for Run x nSperm
########################################################

# model matrices
X  <-  unname(model.matrix(~ nSperm_z -1, data=data))
attr(X,"assign") <- NULL
# str(X)
# head(X)

Z0  <-  unname(model.matrix(~ data$Run -1, data=data))[,-c(9:16)]
attr(Z0,"assign") <- NULL
# str(Z0)
# head(Z0)

Z1  <-  unname(model.matrix(~ data$Run * nSperm_z , data=data))[,-c(1:8)]
attr(Z1,"assign") <- NULL
# str(Z1)
Z1[7:nrow(Z1),1]  <-  0
# Z1[1:20,]

# create data.list
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K0  =  ncol(Z0),
                    K1  =  ncol(Z1),
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    Z0  =  Z0,
                    Z1  =  Z1
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
# nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)


	# Call to STAN
	m3 <- stan(data         =  data.list,
	           seed         =  456789123,
	           file         =  './Stan/mat-logistic-allZ.stan',
	#           sample_file  =  './output/StanFits/N_invest_m3.csv',
	           chains       =  nChains,
	           warmup       =  burnInSteps,
	           iter         =  numSavedSteps,
	           thin         =  thinSteps,
	           save_dso     =  TRUE,
	           control      =  list(adapt_delta = 0.99)
	          )

# message
message("STAN has finished fitting model m3")
system("notify-send \"STAN has finished fitting model m3\"")

# garbage collection
rm(Z0)
rm(Z1)
rm(data.list)
rm(m3)


########################################################
#  m4: Logistic Mixed Effects Regression
#       --  FertRate ~ nSperm_z + ...
## -- random intercept for RUN
## -- random slopes for Run x nSperm
## -- Estimate covariance matrix
########################################################


# model matrices
X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
# str(X)
# head(X)

Z  <-  unname(model.matrix(~ data$Run *nSperm_z -1, data=data))
attr(Z,"assign") <- NULL
# str(Z)
# Z[1:20,]

# create data.list
data.list  <-  list(N    =  nrow(data),
                    P    =  ncol(X),
                    J    =  max(as.numeric(data$Run)),
                    K    =  ncol(Z),
                    grp  =  as.numeric(data$Run),
                    nT   =  data$nEggs,
                    nS   =  data$nFert,
                    X    =  X,
                    Z    =  Z
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 1
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
# nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# inits  <-  list(
#   list(L_run    =  matrix(runif(16^2,-0.5,0.5), nrow=16,ncol=16)),
#   list(tau_run  =  runif(16,0.1,2)),
#   list(u        =  matrix(runif((8*16),-1,1),nrow=8,ncol=16)),
#   list(beta     =  runif(2,-1,1)),
#   list(sigma_y  =  runif(1,0.1,2))
#   )

# Call to STAN
m4 <- stan(data         =  data.list,
           seed         =  567891234,
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/N_invest_m4.csv',
           chains       =  nChains,
           warmup       =  burnInSteps,
           iter         =  numSavedSteps,
           thin         =  thinSteps,
           save_dso     =  TRUE,
           control      =  list(adapt_delta = 0.95) # default adapt_delta of 0.8 threw many divergent transitions. 
          )

# message
message("STAN has finished fitting model m4")
system("notify-send \"STAN has finished fitting model m4\"")

# garbage collection
rm(X)
rm(Z)
rm(data.list)
rm(m4)
