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
head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Lane      <-  factor(data$Lane)
data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
data$Date      <-  dmy(data$Date)
data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

str(data)



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
str(X)
head(X)

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
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# Call to STAN
m1 <- stan(data         =  data.list,
           seed         =  123456789,
           file         =  './Stan/mat-logistic-reg.stan',
           sample_file  =  './output/StanFits/N_invest_m1.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )



########################################################
#  m2: Logistic regression w/ Random Intercept ~ Run 
#      --  FertRate ~ nSperm_z + (1 | Run)
########################################################

# model matrices
# 'fixed effects'
X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
str(X)
head(X)

# 'random effects'
Z  <-  unname(model.matrix(~ data$Run -1, data=data))
attr(Z,"assign") <- NULL
str(Z)
head(Z)

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
thinSteps      = 5
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# Call to STAN
m2 <- stan(data         =  data.list,
           seed         =  234567891,
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/N_invest_m2.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )



########################################################
#  m2b: Logistic Mixed Effects Regression w/ random intercept for RUN
#       --  FertRate ~ nSperm_z + (1 | Run)
#       --  Alternative cell-mean model specification
########################################################

# model matrices
X  <-  unname(model.matrix(~ nSperm_z -1, data=data))
attr(X,"assign") <- NULL
str(X)
head(X)

Z  <-  unname(model.matrix(~ data$Run -1, data=data))
attr(Z,"assign") <- NULL
str(Z)
head(Z)

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
thinSteps      = 5
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# Call to STAN
m2 <- stan(data         =  data.list,
           seed         =  345678912,
           file         =  './Stan/mat-logistic-1Z-cellmean.stan',
           sample_file  =  './output/StanFits/N_invest_m2b.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )
