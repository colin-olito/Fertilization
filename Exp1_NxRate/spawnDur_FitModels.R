#/* 
# * Colin Olito. Created 26/01/2017
# * Analysis of spawning duration data
# * 
# * NOTES:  This file will fit all the necessary
# * 		Stan models for the analysis of the
# * 		spawning duration data, and write the
# * 		Stan sample files to ./output/Stanfits
# * 
# * 		These Stan sample files can then be 
# * 		imported, plotted, and analyzed using
# * 		./spawnDur_Analysis.R
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
spwnData <- read.csv('data/inSitu-SpawnDuration.csv', header=TRUE, stringsAsFactors=FALSE)
spwnData <- data.frame(spwnData)

# Convert grouping variables to factors; Correct Dates
spwnData$Ind    <-  factor(spwnData$Ind)
spwnData$Run    <-  factor(spwnData$Run)
spwnData$Start  <-  hms(spwnData$Start)
spwnData$End    <-  hms(spwnData$End)
spwnData$Dur    <-  hms(spwnData$Dur)
# Check that duration column is equivalent to End - Start.
# as.numeric(seconds(spwnData$Dur)) == as.numeric(seconds(spwnData$End - spwnData$Start))
spwnData$start  <-  as.numeric(seconds(spwnData$Start))
spwnData$end    <-  as.numeric(seconds(spwnData$End))
spwnData$dur    <-  as.numeric(seconds(spwnData$Dur))






########################################################
########################################################
## FIT THE MODELS
########################################################
########################################################

#  Fixed Effects Model Matrix (Same for all models)
X       <-  model.matrix(~ 1 , data=spwnData)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL


# Options for the STAN analyses
nChains       = 3
thinSteps     = 5
nIter         = 2000 * thinSteps #for each chain
burnInSteps   = nIter / 2
nSavedSteps  = ((nIter/thinSteps)/2)*nChains
(nSavedSteps)




########################################################
#  Model SDm1: Gaussian Error
#			 --  Random Effects: Run
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + spwnData$Run)
Znames  <-  dimnames(Z)[[2]]
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

##  Assemble data for stan
data.list  <-  list(N   =  nrow(spwnData),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    Y   =  spwnData$dur,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
SDm1 <- stan(data       =  data.list,
           seed         =  123456,
           file         =  './Stan/mat-gauss-1Z.stan',
           sample_file  =  './output/StanFits/spawnDur_m1.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model SDm1")
system("notify-send \"STAN has finished fitting model SDm1 \"")

# garbage collection
rm(Z)
rm(data.list)
rm(SD.m1)




########################################################
#  Model SDm2: Gaussian Error
#			 --  NO Random Effects
########################################################

##  Assemble data for stan
data.list  <-  list(N   =  nrow(spwnData),
                    P   =  ncol(X), 
                    Y   =  spwnData$dur,
                    X   =  X
                   )

## Call to STAN
SDm2 <- stan(data       =  data.list,
           seed         =  234561,
           file         =  './Stan/mat-gauss.stan',
           sample_file  =  './output/StanFits/spawnDur_m2.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model SDm2")
system("notify-send \"STAN has finished fitting model SDm2 \"")

# garbage collection
rm(data.list)
rm(SDm2)




########################################################
#  Model SDm3: Exponential Error
#			 --  Random Effects: Run
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + spwnData$Run)
Znames  <-  dimnames(Z)[[2]]
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

##  Assemble data for stan
data.list  <-  list(N   =  nrow(spwnData),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    Y   =  spwnData$dur,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
SDm3 <- stan(data       =  data.list,
           seed         =  345612,
           file         =  './Stan/mat-exponential-1Z.stan',
           sample_file  =  './output/StanFits/spawnDur_m3.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

# message
message("STAN has finished fitting model SDm3")
system("notify-send \"STAN has finished fitting model SDm3 \"")

# garbage collection
rm(Z)
rm(data.list)
rm(SDm3)



########################################################
#  Model SDm4: Exponential Error
#			 --  NO Random Effects
########################################################

##  Assemble data for stan
data.list  <-  list(N   =  nrow(spwnData),
                    P   =  ncol(X), 
                    Y   =  spwnData$dur,
                    X   =  X
                   )

## Call to STAN
SDm4 <- stan(data       =  data.list,
           seed         =  456123,
           file         =  './Stan/mat-exponential.stan',
           sample_file  =  './output/StanFits/spawnDur_m4.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )


# message
message("STAN has finished fitting model SDm4")
system("notify-send \"STAN has finished fitting model SDm4 \"")

# garbage collection
rm(data.list)
rm(SDm4)