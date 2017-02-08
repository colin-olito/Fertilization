#/* 
# * Colin Olito. Created 08/02/2017
# * 
# * Load all Data and Models for Analysis of 
# * N-invest experiments
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
# Import Data Sets

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



###############################################
# Import N_invest model results from
#  the stan  sample_files for further analysis 
###############################################
print('Loading N_invest stanfits')
pb <- txtProgressBar(min=0, max=6, style=3)
csvFiles  <-  c('./output/StanFits/N_invest_m1.csv1',
                './output/StanFits/N_invest_m1.csv2',
                './output/StanFits/N_invest_m1.csv3')
NIm1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 1)

csvFiles  <-  c('./output/StanFits/N_invest_m2.csv1',
                './output/StanFits/N_invest_m2.csv2',
                './output/StanFits/N_invest_m2.csv3')
NIm2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 2)

csvFiles  <-  c('./output/StanFits/N_invest_m3.csv1',
                './output/StanFits/N_invest_m3.csv2',
                './output/StanFits/N_invest_m3.csv3')
NIm3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 3)

csvFiles  <-  c('./output/StanFits/N_invest_m1BB.csv1',
                './output/StanFits/N_invest_m1BB.csv2',
                './output/StanFits/N_invest_m1BB.csv3')
NIm1BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 4)

csvFiles  <-  c('./output/StanFits/N_invest_m2BB.csv1',
                './output/StanFits/N_invest_m2BB.csv2',
                './output/StanFits/N_invest_m2BB.csv3')
NIm2BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 5)

csvFiles  <-  c('./output/StanFits/N_invest_m3BB.csv1',
                './output/StanFits/N_invest_m3BB.csv2',
                './output/StanFits/N_invest_m3BB.csv3')
NIm3BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 6)
close(pb)


################################################
##  Extract Samples into data frames, summarize
################################################
print('Extracting stanfits into data frames, etc.')
pb <- txtProgressBar(min=0, max=6, style=3)

NIm1.df    <-  as.data.frame(extract(NIm1))[,-1]
NIm1.summ  <-  plyr:::adply(as.matrix(NIm1.df),2,MCMCsum)
setTxtProgressBar(pb, 1)

NIm2.df    <-  as.data.frame(extract(NIm2))[,-1]
NIm2.summ  <-  plyr:::adply(as.matrix(NIm2.df),2,MCMCsum)
setTxtProgressBar(pb, 2)

NIm3.df    <-  as.data.frame(extract(NIm3))[,-1]
NIm3.summ  <-  plyr:::adply(as.matrix(NIm3.df),2,MCMCsum)
setTxtProgressBar(pb, 3)

NIm1BB.df    <-  as.data.frame(extract(NIm1BB))[,-1]
NIm1BB.summ  <-  plyr:::adply(as.matrix(NIm1BB.df),2,MCMCsum)
setTxtProgressBar(pb, 4)

NIm2BB.df    <-  as.data.frame(extract(NIm2BB))[,-1]
NIm2BB.summ  <-  plyr:::adply(as.matrix(NIm2BB.df),2,MCMCsum)
setTxtProgressBar(pb, 5)

NIm3BB.df    <-  as.data.frame(extract(NIm3BB))[,-1]
NIm3BB.summ  <-  plyr:::adply(as.matrix(NIm3BB.df),2,MCMCsum)
setTxtProgressBar(pb, 6)
close(pb)


#########################################
# LOO Log-likelihood for model selection
#########################################
print('Calculating LOO for model comparison')
pb <- txtProgressBar(min=0, max=6, style=3)
NIm1LL  <-  extract_log_lik(NIm1, parameter_name = "log_lik")
NIm1Loo    <-  loo(NIm1LL)
NIm1WAIC   <-  waic(NIm1LL)
setTxtProgressBar(pb, 1)

NIm2LL  <-  extract_log_lik(NIm2, parameter_name = "log_lik")
NIm2Loo    <-  loo(NIm2LL)
NIm2WAIC   <-  waic(NIm2LL)
setTxtProgressBar(pb, 2)

NIm3LL  <-  extract_log_lik(NIm3, parameter_name = "log_lik")
NIm3Loo    <-  loo(NIm3LL)
NIm3WAIC   <-  waic(NIm3LL)
setTxtProgressBar(pb, 3)

NIm1BBLL  <-  extract_log_lik(NIm1BB, parameter_name = "log_lik")
NIm1BBLoo    <-  loo(NIm1BBLL)
NIm1BBWAIC   <-  waic(NIm1BBLL)
setTxtProgressBar(pb, 4)

NIm2BBLL  <-  extract_log_lik(NIm2BB, parameter_name = "log_lik")
NIm2BBLoo    <-  loo(NIm2BBLL)
NIm2BBWAIC   <-  waic(NIm2BBLL)
setTxtProgressBar(pb, 5)

NIm3BBLL  <-  extract_log_lik(NIm3BB, parameter_name = "log_lik")
NIm3BBLoo    <-  loo(NIm3BBLL)
NIm3BBWAIC   <-  waic(NIm3BBLL)
setTxtProgressBar(pb, 6)
close(pb)
