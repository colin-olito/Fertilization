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
thinSteps     = 5
nIter         = 2000 * thinSteps #for each chain
burnInSteps   = nIter / 2
nSavedSteps  = ((nIter/thinSteps)/2)*nChains
(nSavedSteps)

########################################################
#  Model m1: MAXIMAL MODEL 
#			 --  Random Effects: Run*nSperm_z*Rate
#			 --  w/ COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
#  Not enough observations to include EggPos ixns! Results in only 
#  3 observations per run.
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:data$nSperm_z +
                                data$Run:data$Rate +
                                data$Run:data$nSperm_z:data$Rate)
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


# Call to STAN
m1 <- stan(data         =  data.list,
           seed         =  123456789,
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m1.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
           save_dso     =  TRUE
          )

#  Notification
system('notify-send "Sampling for model m1 complete"')
print("Sampling for model m1 complete")

# garbage collection
rm(data.list)
rm(Z)
rm(m1)





########################################################
#  Model m2: Remove Run x Rate x nSperm ixn 
#			 --  Random Effects: Run + Run*nSperm_z + Run*Rate
#			 --  w/ COVARIANCE MATRIX
########################################################


#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:data$nSperm_z +
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
m2 <- stan(data         =  data.list,
           seed         =  234567891,
           file         =  './Stan/mat-logistic-1Z-cov.stan',
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
#  Model m3: Random intercept & slopes x Run 
#			 --  Random Effects: Run + Run*nSperm_z
#			 --  w/ COVARIANCE MATRIX
########################################################

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:data$nSperm_z)
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

#  Create data.list for STAN
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

# Call to STAN
m3 <- stan(data         =  data.list,
           seed         =  345678912,
           file         =  './Stan/mat-logistic-1Z-cov.stan',
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
#  Model m3a: Random intercept & slopes x Run 
#			 --  Random Effects: Run + Run*nSperm_z
#			 --  NO COVARIANCE MATRIX
########################################################

# NOTE: STAN threw a single error: 
#         Exception thrown at line 47: normal_log: Scale parameter is 0 but must be > 0

#  Random Slopes Model Matrix
Z0       <-  model.matrix(~ -1 + data$Run, data=data)
Z0names  <-  dimnames(Z0)[[2]]
Z0       <-  unname(Z0)
attr(Z0,"assign") <- NULL

#  Random Intercepts Model Matrix
Z1  <-  model.matrix(~ -1 +  data$Run:nSperm_z , data=data)
Z1       <-  unname(Z1)
attr(Z1,"assign") <- NULL

#  Assemble data.list for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z0),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    Z0  =  Z0,
                    Z1  =  Z1
                   )

# Call to STAN
m3a <- stan(data         =  data.list,
            seed         =  456789123,
            file         =  './Stan/mat-logistic-2Z.stan',
            sample_file  =  './output/StanFits/NxRate_m3a.csv',
            chains       =  nChains,
            iter         =  nIter,
            thin         =  thinSteps,
            save_dso     =  TRUE
           )

# message
message("STAN has finished fitting model m3a")
system("notify-send \"STAN has finished fitting model m3a\"")

# garbage collection
rm(Z0)
rm(Z1)
rm(data.list)
rm(m3a)



########################################################
#  Model m4: Random intercept x Run 
#			 --  Random Effects: Run
#			 --  w/ COVARIANCE MATRIX
########################################################

# NOTE: STAN threw a few numerical problems for this analysis... 
#        maybe play with the control settings?

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
m4 <- stan(data         =  data.list,
           seed         =  567891234,
           file         =  './Stan/mat-logistic-1Z-cov.stan',
           sample_file  =  './output/StanFits/NxRate_m4.csv',
           chains       =  nChains,
           iter         =  nIter,
           thin         =  thinSteps,
#           control      =  list(adapt_delta = 0.9) # default adapt_delta of 0.8 threw many divergent transitions. 
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
#  Model m4a: Random intercept x Run 
#			  --  Random Effects: Run
#			  --  NO COVARIANCE MATRIX
########################################################


#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run, data=data)
Z       <-  unname(Z)
attr(Z,"assign") <- NULL

#  Assemble data.list for STAN
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X,
                    Z   =  Z
                   )

## Call to STAN
m4a <- stan(data        =  data.list,
            seed         =  678912345,
            file        =  './Stan/mat-logistic-1Z.stan',
            sample_file  =  './output/StanFits/NxRate_m4a.csv',
            chains      =  nChains,
            iter        =  nIter,
            thin        =  thinSteps,
            save_dso    =  TRUE
           )


# message
message("STAN has finished fitting model m4a")
system("notify-send \"STAN has finished fitting model m4a\"")

# garbage collection
rm(Z)
rm(data.list)
rm(m4a)




########################################################
#  Model m5: Simple Logistic Regression 
#			  --  NO Random Effects
#			  --  NO COVARIANCE MATRIX
########################################################


#  Assemble data.list for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
                    X   =  X
                   )

# Call to STAN
m5 <- stan(data         =  data.list,
           seed         =  678912345,
           file         =  './Stan/mat-logistic-reg.stan',
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
rm(data.list)
rm(m5)



