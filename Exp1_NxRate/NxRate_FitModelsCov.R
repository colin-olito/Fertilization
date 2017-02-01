#/* 
# * Colin Olito. Created 12/01/2017
# * Analysis of 2nd flume experiment: N x Rate
# * 
# * NOTES:  This file will fit all the necessary
# * 		Stan models for the analysis of the
# * 		NxRate flume data, and write the
# * 		Stan sample files to ./output/Stanfits
# * 
# * 		These Stan sample files can then be 
# * 		imported, plotted, and analyzed using
# * 		./NxRate_analyze
# *          
# */

rm(list=ls())
################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

#******************
# DEPENDENCIES
source('R/functions.R')

#*******************
# Import Data
data <- read.csv('data/NxRate_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)
# head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Rate      <-  factor(data$Rate)
data$EggPos    <-  factor(data$EggPos)
data$Lane      <-  factor(data$Lane)
data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
data$Date      <-  dmy(data$Date)
data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)






########################################################
########################################################
## FIT THE MODELS
########################################################
########################################################

#  Fixed Effects Model Matrix (Same for all models)
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL


# Options for the STAN analyses
nChains       = 3
thinSteps     = 1
nIter         = 2000 * thinSteps #for each chain
burnInSteps   = nIter / 2
nSavedSteps  = ((nIter/thinSteps)/2)*nChains
(nSavedSteps)

# Create pseudo-random seeds for Stan 
set.seed(234567891)
randos  <- as.integer(runif(n=49) * 1e8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################
##  NESTED MODEL SET
##   Modeling Covariance Structure Using 
##   Spherical Random Effects
########################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

########################################################
#  Model m1a: MAXIMAL MODEL
#      --  Random Effects: (1 + nSperm_z | Run : Rate : EggPos)
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run                     +
                                Run : nSperm_z          +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z +
                                Run : Rate   : EggPos   +
                                Run : Rate   : EggPos  : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

#  Create data.list for stan
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
m1a <- stan(data         =  data.list,
           seed         =  randos[1],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m1a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9), # increase adapt_delta above 0.8 if many divergent transitions. 
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m1a")
system("notify-send \"STAN has finished fitting model m1a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m1a)

########################################################
#  Model m2a: 
#      --  Maximal model, remove Run : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run                     +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z +
                                Run : Rate   : EggPos   +
                                Run : Rate   : EggPos  : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m2a <- stan(data         =  data.list,
           seed         =  randos[2],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m2a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2a")
system("notify-send \"STAN has finished fitting model m2a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m2a)

########################################################
#  Model m3a: 
#      --  remove 4-way interaction term
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z          +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m3a <- stan(data         =  data.list,
           seed         =  randos[3],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m3a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m3a")
system("notify-send \"STAN has finished fitting model m3a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m3a)

########################################################
#  Model m4a: 
#      --  remove 4-way interaction term
#      --  remove Run : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m4a <- stan(data         =  data.list,
            seed         =  randos[4],
            file         =  './Stan/mat-logistic-1Z-cov.stan',
            sample_file  =  './output/StanFits/NxRate_m4a.csv',
            chains       =  nChains,
            iter         =  nIter,
            thin         =  thinSteps,
#            control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
            save_dso     =  TRUE
           )

# message
message("STAN has finished fitting model m4a")
system("notify-send \"STAN has finished fitting model m4a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m4a)

########################################################
#  Model m5a: 
#      --  #  Remove ONE 3-way term: Run : Rate : EggPos
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z          +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m5a <- stan(data         =  data.list,
           seed         =  randos[5],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m5a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m5a")
system("notify-send \"STAN has finished fitting model m5a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m5a)

########################################################
#  Model m6a: 
#      --  #  Remove ONE 3-way term: Run : EggPos : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z          +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m6a <- stan(data         =  data.list,
           seed         =  randos[6],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m6a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m6a")
system("notify-send \"STAN has finished fitting model m6a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m6a)

########################################################
#  Model m7a: 
#      --  #  Remove ONE 3-way term: Run : Rate   : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run                     +
                                Run : nSperm_z          +
                                Run : Rate              +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m7a <- stan(data         =  data.list,
           seed         =  randos[7],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m7a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m7a")
system("notify-send \"STAN has finished fitting model m7a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m7a)

########################################################
#  Model m8a: 
#      --  #  Remove ONE 3-way term: Run : Rate : EggPos
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m8a <- stan(data         =  data.list,
           seed         =  randos[8],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m8a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m8a")
system("notify-send \"STAN has finished fitting model m8a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m8a)

########################################################
#  Model m9a: 
#      --  #  Remove ONE 3-way term: Run : EggPos : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos            +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m9a <- stan(data         =  data.list,
           seed         =  randos[9],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m9a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m9a")
system("notify-send \"STAN has finished fitting model m9a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m9a)

########################################################
#  Model m10a: 
#      --  #  Remove ONE 3-way term: Run : Rate   : nSperm_z
#      --  #  Remove Run : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run                     +
                                Run : Rate              +
                                Run : EggPos            +
                                Run : EggPos : nSperm_z +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m10a <- stan(data         =  data.list,
           seed         =  randos[10],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m10a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m10a")
system("notify-send \"STAN has finished fitting model m10a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m10a)

########################################################
#  Model m11a: 
#      --  #  Remove TWO 3-way terms :
#                Run : Rate : EggPos
#                Run : EggPos : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run                     +
                                Run : nSperm_z          +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m11a <- stan(data         =  data.list,
           seed         =  randos[11],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m11a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m11a")
system("notify-send \"STAN has finished fitting model m11a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m11a)

########################################################
#  Model m12a: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : nSperm_z
#                Run : EggPos : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m12a <- stan(data         =  data.list,
           seed         =  randos[12],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m12a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m12a")
system("notify-send \"STAN has finished fitting model m12a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m12a)

########################################################
#  Model m13a: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : EggPos
#                Run : Rate : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : EggPos : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m13a <- stan(data         =  data.list,
           seed         =  randos[13],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m13a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m13a")
system("notify-send \"STAN has finished fitting model m13a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m13a)

########################################################
#  Model m14a: 
#      --  #  Remove TWO 3-way terms :
#                Run : Rate : EggPos
#                Run : EggPos : nSperm_z
#      --  #  Remove Run : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run                     +
                                Run : Rate              +
                                Run : Rate   : nSperm_z +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m14a <- stan(data         =  data.list,
           seed         =  randos[14],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m14a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m14a")
system("notify-send \"STAN has finished fitting model m14a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m14a)

########################################################
#  Model m15a: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : nSperm_z
#                Run : EggPos : nSperm_z
#      --  #  Remove Run : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : Rate   : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m15a <- stan(data         =  data.list,
           seed         =  randos[15],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m15a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m15a")
system("notify-send \"STAN has finished fitting model m15a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m15a)

########################################################
#  Model m16a: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : EggPos
#                Run : Rate : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : EggPos : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m16a <- stan(data         =  data.list,
           seed         =  randos[16],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m16a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m16a")
system("notify-send \"STAN has finished fitting model m16a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m16a)

########################################################
#  Model m17a: 
#      --  #  Remove ALL 3-way terms
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate     +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m17a <- stan(data         =  data.list,
           seed         =  randos[17],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m17a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m17a")
system("notify-send \"STAN has finished fitting model m17a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m17a)

########################################################
#  Model m18a: 
#      --  #  Remove ONE 2-way term : Run: EggPos
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m18a <- stan(data         =  data.list,
           seed         =  randos[18],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m18a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m18a")
system("notify-send \"STAN has finished fitting model m18a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m18a)

########################################################
#  Model m19a: 
#      --  #  Remove ONE 2-way term : Run: Rate
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m19a <- stan(data         =  data.list,
           seed         =  randos[19],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m19a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m19a")
system("notify-send \"STAN has finished fitting model m19a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m19a)

########################################################
#  Model m20a: 
#      --  #  Remove ONE 2-way term : Run: nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run        +
                                Run : Rate +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m20a <- stan(data         =  data.list,
           seed         =  randos[20],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m20a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m20a")
system("notify-send \"STAN has finished fitting model m20a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m20a)

########################################################
#  Model m21a: 
#      --  #  Remove Two 2-way terms: 
#               Run : Rate
#               Run : EggPos
#      --  Random intercept & slope ~ Run
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m21a <- stan(data         =  data.list,
           seed         =  randos[21],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m21a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m21a")
system("notify-send \"STAN has finished fitting model m21a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m21a)

########################################################
#  Model m22a: 
#      --  #  Remove Two 2-way terms: 
#               Run : Rate
#               Run : nSperm_z
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m22a <- stan(data         =  data.list,
           seed         =  randos[22],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m22a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m22a")
system("notify-send \"STAN has finished fitting model m22a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m22a)

########################################################
#  Model m23a: 
#      --  #  Remove Two 2-way terms: 
#               Run : nSperm_z
#               Run : EggPos
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run + 
                                Run : Rate,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m23a <- stan(data         =  data.list,
           seed         =  randos[23],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m23a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m23a")
system("notify-send \"STAN has finished fitting model m23a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m23a)

########################################################
#  Model m24a: 
#      --  #  Remove All 2-way terms: 
#               Random intercept ~ Run
#      --  w/ ESTIMATED COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run, data = data)
Znames  <-  dimnames(Z)[[2]]
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
m24a <- stan(data         =  data.list,
           seed         =  randos[24],
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m24a.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m24a")
system("notify-send \"STAN has finished fitting model m24a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m24a)