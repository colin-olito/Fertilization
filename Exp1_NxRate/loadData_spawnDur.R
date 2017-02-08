#/* 
# * Colin Olito. Created 09/02/2017
# * 
# * Load all Data Sets
# * 
# */

rm(list=ls())
#################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

###############
# DEPENDENCIES
source('R/functions.R')
source('R/functions-figures.R')

##################
# Import Data Set

#*****************************
#  Spawning Duration Data Set
print('Importing inSitu Spawning Duration Data Set')
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


############################################
# Import Spawning Duration model results 
# from the stan sample_files for further analysis 
############################################
print('Loading Spawning Duration stanfits')
pb <- txtProgressBar(min=0,max=3, style=3)
csvFiles  <-  c('./output/StanFits/spawnDur_m1.csv1',
                './output/StanFits/spawnDur_m1.csv2',
                './output/StanFits/spawnDur_m1.csv3')
SDm1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 1)

csvFiles  <-  c('./output/StanFits/spawnDur_m2.csv1',
                './output/StanFits/spawnDur_m2.csv2',
                './output/StanFits/spawnDur_m2.csv3')
SDm2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 2)

# csvFiles  <-  c('./output/StanFits/spawnDur_m2b.csv1',
#                 './output/StanFits/spawnDur_m2b.csv2',
#                 './output/StanFits/spawnDur_m2b.csv3')
# SDm3        <-  read_stan_csv(csvFiles, col_major = TRUE)
# rm(csvFiles)

csvFiles  <-  c('./output/StanFits/spawnDur_m3.csv1',
                './output/StanFits/spawnDur_m3.csv2',
                './output/StanFits/spawnDur_m3.csv3')
SDm4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 3)
close(pb)
