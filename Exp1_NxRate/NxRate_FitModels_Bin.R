#/* 
# * Author: XXXX XXXX. Created 12/01/2017
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
nSavedSteps   = ((nIter/thinSteps)/2)*nChains
print('nSavedSteps')
nSavedSteps

# Create pseudo-random seeds for Stan 
set.seed(123456789)
randos  <- as.integer(runif(n=49) * 1e8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################
##  NESTED MODEL SET, WITHOUT MODELING COVARIANCE
##  STRUCTURE
########################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


########################################################
#  Model m1: MAXIMAL MODEL
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
m1 <- stan(data         =  data.list,
           seed         =  randos[1],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m1.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m1")
system("notify-send \"STAN has finished fitting model m1\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m1)

########################################################
#  Model m2: 
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
m2 <- stan(data         =  data.list,
           seed         =  randos[2],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m2.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m2")
system("notify-send \"STAN has finished fitting model m2\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m2)

########################################################
#  Model m3: 
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
m3 <- stan(data         =  data.list,
           seed         =  randos[3],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m3.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m3")
system("notify-send \"STAN has finished fitting model m3\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m3)

########################################################
#  Model m4: 
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
m4 <- stan(data         =  data.list,
           seed         =  randos[4],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m4.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m4")
system("notify-send \"STAN has finished fitting model m4\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m4)

########################################################
#  Model m5: 
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
m5 <- stan(data         =  data.list,
           seed         =  randos[5],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m5.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m5")
system("notify-send \"STAN has finished fitting model m5\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m5)

########################################################
#  Model m6: 
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
m6 <- stan(data         =  data.list,
           seed         =  randos[6],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m6.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m6")
system("notify-send \"STAN has finished fitting model m6\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m6)

########################################################
#  Model m7: 
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
m7 <- stan(data         =  data.list,
           seed         =  randos[7],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m7.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m7")
system("notify-send \"STAN has finished fitting model m7\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m7)

########################################################
#  Model m8: 
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
m8 <- stan(data         =  data.list,
           seed         =  randos[8],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m8.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m8")
system("notify-send \"STAN has finished fitting model m8\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m8)

########################################################
#  Model m9: 
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
m9 <- stan(data         =  data.list,
           seed         =  randos[9],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m9.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m9")
system("notify-send \"STAN has finished fitting model m9\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m9)

########################################################
#  Model m10: 
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
m10 <- stan(data         =  data.list,
           seed         =  randos[10],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m10.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m10")
system("notify-send \"STAN has finished fitting model m10\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m10)

########################################################
#  Model m11: 
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
m11 <- stan(data         =  data.list,
           seed         =  randos[11],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m11.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m11")
system("notify-send \"STAN has finished fitting model m11\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m11)

########################################################
#  Model m12: 
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
m12 <- stan(data         =  data.list,
           seed         =  randos[12],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m12.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m12")
system("notify-send \"STAN has finished fitting model m12\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m12)

########################################################
#  Model m13: 
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
m13 <- stan(data         =  data.list,
           seed         =  randos[13],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m13.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m13")
system("notify-send \"STAN has finished fitting model m13\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m13)

########################################################
#  Model m14: 
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
m14 <- stan(data         =  data.list,
           seed         =  randos[14],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m14.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m14")
system("notify-send \"STAN has finished fitting model m14\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m14)

########################################################
#  Model m15: 
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
m15 <- stan(data         =  data.list,
           seed         =  randos[15],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m15.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m15")
system("notify-send \"STAN has finished fitting model m15\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m15)

########################################################
#  Model m16: 
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
m16 <- stan(data         =  data.list,
           seed         =  randos[16],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m16.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m16")
system("notify-send \"STAN has finished fitting model m16\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m16)

########################################################
#  Model m17: 
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
m17 <- stan(data         =  data.list,
           seed         =  randos[17],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m17.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m17")
system("notify-send \"STAN has finished fitting model m17\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m17)

########################################################
#  Model m18: 
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
m18 <- stan(data         =  data.list,
           seed         =  randos[18],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m18.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m18")
system("notify-send \"STAN has finished fitting model m18\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m18)

########################################################
#  Model m19: 
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
m19 <- stan(data         =  data.list,
           seed         =  randos[19],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m19.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m19")
system("notify-send \"STAN has finished fitting model m19\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m19)

########################################################
#  Model m20: 
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
m20 <- stan(data         =  data.list,
           seed         =  randos[20],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m20.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m20")
system("notify-send \"STAN has finished fitting model m20\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m20)

########################################################
#  Model m21: 
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
m21 <- stan(data         =  data.list,
           seed         =  randos[21],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m21.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m21")
system("notify-send \"STAN has finished fitting model m21\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m21)

########################################################
#  Model m22: 
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
m22 <- stan(data         =  data.list,
           seed         =  randos[22],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m22.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m22")
system("notify-send \"STAN has finished fitting model m22\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m22)

########################################################
#  Model m23: 
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
m23 <- stan(data         =  data.list,
           seed         =  randos[23],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m23.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m23")
system("notify-send \"STAN has finished fitting model m23\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m23)

########################################################
#  Model m24: 
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
m24 <- stan(data         =  data.list,
           seed         =  randos[24],
           file         =  './Stan/mat-logistic-1Z.stan',
           sample_file  =  './output/StanFits/NxRate_m24.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m24")
system("notify-send \"STAN has finished fitting model m24\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m24)

########################################################
#  Model m25: Simple Logistic Regression 
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
m25 <- stan(data         =  data.list,
           seed         =  randos[25],
           file         =  './Stan/mat-logistic-reg.stan',
           sample_file  =  './output/StanFits/NxRate_m25.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model m25")
system("notify-send \"STAN has finished fitting model m25\"")

# garbage collection
rm(data.list)
rm(m25)
