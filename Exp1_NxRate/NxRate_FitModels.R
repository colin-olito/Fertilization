#/* 
# * Colin Olito. Created 12/041/2016.
# * 
# * NOTES: 2nd Flume Experiment
# *         crossing N x Rate; with 2 egg patches
# *          
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

#  Options for all analyses
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)







