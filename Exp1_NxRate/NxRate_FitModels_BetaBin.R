#/* 
# * Colin Olito. Created 03/02/2017
# * Analysis of 2nd flume experiment: N x Rate
# * 
# * NOTES:  This file will fit all the necessary
# * 		Stan models for the analysis of the
# * 		NxRate flume data, using Beta-Binomial 
# *     regression, and write the Stan sample 
# *     files to ./output/Stanfits.
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
nChains       =  3
thinSteps     =  1
nIter         =  2000 * thinSteps #for each chain
burnInSteps   =  nIter / 2
nSavedSteps   =  ((nIter/thinSteps)/2)*nChains
(nSavedSteps)

# Create pseudo-random seeds for Stan 
set.seed(12345678)
randos  <- as.integer(runif(n=49) * 1e8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################
##  NESTED MODEL SET, WITHOUT MODELING COVARIANCE
##  STRUCTURE
########################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


########################################################
#  Model m1BB: MAXIMAL MODEL
#      --  Random Effects: (1 + nSperm_z | Run : Rate : EggPos)
#      --  NO COVARIANCE MATRIX
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
m1BB <- stan(data         =  data.list,
           seed         =  randos[1],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m1BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m1BB")
system("notify-send \"STAN has finished fitting model m1BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m1BB)

########################################################
#  Model m2BB: 
#      --  Maximal model, remove Run : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m2BB <- stan(data         =  data.list,
           seed         =  randos[2],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m2BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2BB")
system("notify-send \"STAN has finished fitting model m2BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m2BB)

########################################################
#  Model m3BB: 
#      --  remove 4-way interaction term
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m3BB <- stan(data         =  data.list,
           seed         =  randos[3],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m3BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m3BB")
system("notify-send \"STAN has finished fitting model m3BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m3BB)

########################################################
#  Model m4BB: 
#      --  remove 4-way interaction term
#      --  remove Run : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m4BB <- stan(data         =  data.list,
           seed         =  randos[4],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m4BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m4BB")
system("notify-send \"STAN has finished fitting model m4BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m4BB)

########################################################
#  Model m5BB: 
#      --  #  Remove ONE 3-way term: Run : Rate : EggPos
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m5BB <- stan(data         =  data.list,
           seed         =  randos[5],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m5BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m5BB")
system("notify-send \"STAN has finished fitting model m5BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m5BB)

########################################################
#  Model m6BB: 
#      --  #  Remove ONE 3-way term: Run : EggPos : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m6BB <- stan(data         =  data.list,
           seed         =  randos[6],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m6BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m6BB")
system("notify-send \"STAN has finished fitting model m6BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m6BB)

########################################################
#  Model m7BB: 
#      --  #  Remove ONE 3-way term: Run : Rate   : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m7BB <- stan(data         =  data.list,
           seed         =  randos[7],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m7BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m7BB")
system("notify-send \"STAN has finished fitting model m7BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m7BB)

########################################################
#  Model m8BB: 
#      --  #  Remove ONE 3-way term: Run : Rate : EggPos
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m8BB <- stan(data         =  data.list,
           seed         =  randos[8],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m8BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m8BB")
system("notify-send \"STAN has finished fitting model m8BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m8BB)

########################################################
#  Model m9BB: 
#      --  #  Remove ONE 3-way term: Run : EggPos : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m9BB <- stan(data         =  data.list,
           seed         =  randos[9],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m9BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m9BB")
system("notify-send \"STAN has finished fitting model m9BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m9BB)

########################################################
#  Model m10BB: 
#      --  #  Remove ONE 3-way term: Run : Rate   : nSperm_z
#      --  #  Remove Run : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m10BB <- stan(data         =  data.list,
           seed         =  randos[10],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m10BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m10BB")
system("notify-send \"STAN has finished fitting model m10BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m10BB)

########################################################
#  Model m11BB: 
#      --  #  Remove TWO 3-way terms :
#                Run : Rate : EggPos
#                Run : EggPos : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m11BB <- stan(data         =  data.list,
           seed         =  randos[11],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m11BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m11BB")
system("notify-send \"STAN has finished fitting model m11BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m11BB)

########################################################
#  Model m12BB: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : nSperm_z
#                Run : EggPos : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m12BB <- stan(data         =  data.list,
           seed         =  randos[12],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m12BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m12BB")
system("notify-send \"STAN has finished fitting model m12BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m12BB)

########################################################
#  Model m13BB: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : EggPos
#                Run : Rate : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m13BB <- stan(data         =  data.list,
           seed         =  randos[13],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m13BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m13BB")
system("notify-send \"STAN has finished fitting model m13BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m13BB)

########################################################
#  Model m14BB: 
#      --  #  Remove TWO 3-way terms :
#                Run : Rate : EggPos
#                Run : EggPos : nSperm_z
#      --  #  Remove Run : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m14BB <- stan(data         =  data.list,
           seed         =  randos[14],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m14BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m14BB")
system("notify-send \"STAN has finished fitting model m14BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m14BB)

########################################################
#  Model m15BB: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : nSperm_z
#                Run : EggPos : nSperm_z
#      --  #  Remove Run : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m15BB <- stan(data         =  data.list,
           seed         =  randos[15],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m15BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m15BB")
system("notify-send \"STAN has finished fitting model m15BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m15BB)

########################################################
#  Model m16BB: 
#      --  #  Remove TWO 3-way terms
#                Run : Rate : EggPos
#                Run : Rate : nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m16BB <- stan(data         =  data.list,
           seed         =  randos[16],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m16BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m16BB")
system("notify-send \"STAN has finished fitting model m16BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m16BB)

########################################################
#  Model m17BB: 
#      --  #  Remove ALL 3-way terms
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m17BB <- stan(data         =  data.list,
           seed         =  randos[17],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m17BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m17BB")
system("notify-send \"STAN has finished fitting model m17BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m17BB)

########################################################
#  Model m18BB: 
#      --  #  Remove ONE 2-way term : Run: EggPos
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m18BB <- stan(data         =  data.list,
           seed         =  randos[18],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m18BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m18BB")
system("notify-send \"STAN has finished fitting model m18BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m18BB)

########################################################
#  Model m19BB: 
#      --  #  Remove ONE 2-way term : Run: Rate
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m19BB <- stan(data         =  data.list,
           seed         =  randos[19],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m19BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m19BB")
system("notify-send \"STAN has finished fitting model m19BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m19BB)

########################################################
#  Model m20BB: 
#      --  #  Remove ONE 2-way term : Run: nSperm_z
#      --  NO COVARIANCE MATRIX
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
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m20BB <- stan(data         =  data.list,
           seed         =  randos[20],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m20BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m20BB")
system("notify-send \"STAN has finished fitting model m20BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m20BB)

########################################################
#  Model m21BB: 
#      --  #  Remove Two 2-way terms: 
#               Run : Rate
#               Run : EggPos
#      --  Random intercept & slope ~ Run
#      --  NO COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m21BB <- stan(data         =  data.list,
           seed         =  randos[21],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m21BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m21BB")
system("notify-send \"STAN has finished fitting model m21BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m21BB)

########################################################
#  Model m22BB: 
#      --  #  Remove Two 2-way terms: 
#               Run : Rate
#               Run : nSperm_z
#      --  NO COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run +
                                Run : EggPos,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m22BB <- stan(data         =  data.list,
           seed         =  randos[22],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m22BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m22BB")
system("notify-send \"STAN has finished fitting model m22BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m22BB)

########################################################
#  Model m23BB: 
#      --  #  Remove Two 2-way terms: 
#               Run : nSperm_z
#               Run : EggPos
#      --  NO COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run + 
                                Run : Rate,
                         data = data)
Znames  <-  dimnames(Z)[[2]]
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
m23BB <- stan(data         =  data.list,
           seed         =  randos[23],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m23BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m23BB")
system("notify-send \"STAN has finished fitting model m23BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m23BB)

########################################################
#  Model m24BB: 
#      --  #  Remove All 2-way terms: 
#               Random intercept ~ Run
#      --  NO COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + Run, data = data)
Znames  <-  dimnames(Z)[[2]]
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
m24BB <- stan(data         =  data.list,
           seed         =  randos[24],
           file         =  './Stan/mat-BetaBin-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m24BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # increase adapt_delta above 0.8 if many divergent transitions.
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m24BB")
system("notify-send \"STAN has finished fitting model m24BB\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m24BB)

########################################################
#  Model m25BB: Simple Beta-Binomial Logistic Regression 
#       --  NO Random Effects
#       --  NO COVARIANCE MATRIX
########################################################


#  Assemble data.list for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X
                   )

# Call to STAN
m25BB <- stan(data         =  data.list,
           seed         =  randos[25],
           file         =  './Stan/mat-BetaBin.stan',
           sample_file  =  './output/StanFits/NxRate_m25BB.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m25BB")
system("notify-send \"STAN has finished fitting model m25BB\"")

# garbage collection
rm(data.list)
rm(m25BB)
