#/* 
# * Colin Olito. Created 25/01/2017
# * 
# * Load all Data and Models for Analysis of 
# * 1st flume experiment: N-invest
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

#******************
#  NxRate Data Set
print('Importing NxRate Data Set')
data <- read.csv('data/NxRate_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)

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

############################################
# Import N_invest model results from the stan 
# sample_files for further analysis 
############################################
print('Loading N_invest stanfits')
pb <- txtProgressBar(min=0,max=5, style=3)
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

csvFiles  <-  c('./output/StanFits/N_invest_m2b.csv1',
                './output/StanFits/N_invest_m2b.csv2',
                './output/StanFits/N_invest_m2b.csv3')
NIm2b        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 3)

csvFiles  <-  c('./output/StanFits/N_invest_m3.csv1',
                './output/StanFits/N_invest_m3.csv2',
                './output/StanFits/N_invest_m3.csv3')
NIm3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 4)

csvFiles  <-  c('./output/StanFits/N_invest_m4.csv1',
                './output/StanFits/N_invest_m4.csv2',
                './output/StanFits/N_invest_m4.csv3')
NIm4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 5)
close(pb)

############################################
# Import NxRate model results from the stan 
# sample_files for further analysis 
############################################
print('Loading NxRate stanfits')
pb <- txtProgressBar(min=0,max=25, style=3)
setTxtProgressBar(pb, 0)
csvFiles  <-  c('./output/StanFits/NxRate_m1.csv1',
                './output/StanFits/NxRate_m1.csv2',
                './output/StanFits/NxRate_m1.csv3')
m1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 1)

csvFiles  <-  c('./output/StanFits/NxRate_m2.csv1',
                './output/StanFits/NxRate_m2.csv2',
                './output/StanFits/NxRate_m2.csv3')
m2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 2)

csvFiles  <-  c('./output/StanFits/NxRate_m3.csv1',
                './output/StanFits/NxRate_m3.csv2',
                './output/StanFits/NxRate_m3.csv3')
m3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 3)

csvFiles  <-  c('./output/StanFits/NxRate_m4.csv1',
                './output/StanFits/NxRate_m4.csv2',
                './output/StanFits/NxRate_m4.csv3')
m4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 4)

csvFiles  <-  c('./output/StanFits/NxRate_m5.csv1',
                './output/StanFits/NxRate_m5.csv2',
                './output/StanFits/NxRate_m5.csv3')
m5        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 5)

csvFiles  <-  c('./output/StanFits/NxRate_m6.csv1',
                './output/StanFits/NxRate_m6.csv2',
                './output/StanFits/NxRate_m6.csv3')
m6        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 6)

csvFiles  <-  c('./output/StanFits/NxRate_m7.csv1',
                './output/StanFits/NxRate_m7.csv2',
                './output/StanFits/NxRate_m7.csv3')
m7        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 7)

csvFiles  <-  c('./output/StanFits/NxRate_m8.csv1',
                './output/StanFits/NxRate_m8.csv2',
                './output/StanFits/NxRate_m8.csv3')
m8        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 8)

csvFiles  <-  c('./output/StanFits/NxRate_m9.csv1',
                './output/StanFits/NxRate_m9.csv2',
                './output/StanFits/NxRate_m9.csv3')
m9        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 9)

csvFiles  <-  c('./output/StanFits/NxRate_m10.csv1',
                './output/StanFits/NxRate_m10.csv2',
                './output/StanFits/NxRate_m10.csv3')
m10        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 10)

csvFiles  <-  c('./output/StanFits/NxRate_m11.csv1',
                './output/StanFits/NxRate_m11.csv2',
                './output/StanFits/NxRate_m11.csv3')
m11        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 11)

csvFiles  <-  c('./output/StanFits/NxRate_m12.csv1',
                './output/StanFits/NxRate_m12.csv2',
                './output/StanFits/NxRate_m12.csv3')
m12        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 12)

csvFiles  <-  c('./output/StanFits/NxRate_m13.csv1',
                './output/StanFits/NxRate_m13.csv2',
                './output/StanFits/NxRate_m13.csv3')
m13        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 13)

csvFiles  <-  c('./output/StanFits/NxRate_m14.csv1',
                './output/StanFits/NxRate_m14.csv2',
                './output/StanFits/NxRate_m14.csv3')
m14        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 14)

csvFiles  <-  c('./output/StanFits/NxRate_m15.csv1',
                './output/StanFits/NxRate_m15.csv2',
                './output/StanFits/NxRate_m15.csv3')
m15        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 15)

csvFiles  <-  c('./output/StanFits/NxRate_m16.csv1',
                './output/StanFits/NxRate_m16.csv2',
                './output/StanFits/NxRate_m16.csv3')
m16        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 16)

csvFiles  <-  c('./output/StanFits/NxRate_m17.csv1',
                './output/StanFits/NxRate_m17.csv2',
                './output/StanFits/NxRate_m17.csv3')
m17        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 17)

csvFiles  <-  c('./output/StanFits/NxRate_m18.csv1',
                './output/StanFits/NxRate_m18.csv2',
                './output/StanFits/NxRate_m18.csv3')
m18        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 18)

csvFiles  <-  c('./output/StanFits/NxRate_m19.csv1',
                './output/StanFits/NxRate_m19.csv2',
                './output/StanFits/NxRate_m19.csv3')
m19        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 19)

csvFiles  <-  c('./output/StanFits/NxRate_m20.csv1',
                './output/StanFits/NxRate_m20.csv2',
                './output/StanFits/NxRate_m20.csv3')
m20        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 20)

csvFiles  <-  c('./output/StanFits/NxRate_m21.csv1',
                './output/StanFits/NxRate_m21.csv2',
                './output/StanFits/NxRate_m21.csv3')
m21        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 21)

csvFiles  <-  c('./output/StanFits/NxRate_m22.csv1',
                './output/StanFits/NxRate_m22.csv2',
                './output/StanFits/NxRate_m22.csv3')
m22        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 22)

csvFiles  <-  c('./output/StanFits/NxRate_m23.csv1',
                './output/StanFits/NxRate_m23.csv2',
                './output/StanFits/NxRate_m23.csv3')
m23        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 23)

csvFiles  <-  c('./output/StanFits/NxRate_m24.csv1',
                './output/StanFits/NxRate_m24.csv2',
                './output/StanFits/NxRate_m24.csv3')
m24        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 24)

csvFiles  <-  c('./output/StanFits/NxRate_m25.csv1',
                './output/StanFits/NxRate_m25.csv2',
                './output/StanFits/NxRate_m25.csv3')
m25        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 25)
close(pb)






################################################
##  Extract Samples into data frames, summarize
################################################
print('Extracting stanfits into data frames, etc.')
pb <- txtProgressBar(min=0,max=25, style=3)
m1.df     <-  as.data.frame(extract(m1))  [,-1]
m1.summ   <-  plyr:::adply(as.matrix(m1.df),2,MCMCsum)
setTxtProgressBar(pb, 1)
m2.df     <-  as.data.frame(extract(m2))[,-1]
m2.summ   <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)
setTxtProgressBar(pb, 2)
m3.df     <-  as.data.frame(extract(m3))[,-1]
m3.summ   <-  plyr:::adply(as.matrix(m3.df),2,MCMCsum)
setTxtProgressBar(pb, 3)
m4.df     <-  as.data.frame(extract(m4))[,-1]
m4.summ   <-  plyr:::adply(as.matrix(m4.df),2,MCMCsum)
setTxtProgressBar(pb, 4)
m5.df     <-  as.data.frame(extract(m5))[,-1]
m5.summ   <-  plyr:::adply(as.matrix(m5.df),2,MCMCsum)
setTxtProgressBar(pb, 5)
m6.df     <-  as.data.frame(extract(m6))[,-1]
m6.summ   <-  plyr:::adply(as.matrix(m6.df),2,MCMCsum)
setTxtProgressBar(pb, 6)
m7.df     <-  as.data.frame(extract(m7))[,-1]
m7.summ   <-  plyr:::adply(as.matrix(m7.df),2,MCMCsum)
setTxtProgressBar(pb, 7)
m8.df     <-  as.data.frame(extract(m8))[,-1]
m8.summ   <-  plyr:::adply(as.matrix(m8.df),2,MCMCsum)
setTxtProgressBar(pb, 8)
m9.df     <-  as.data.frame(extract(m9))[,-1]
m9.summ   <-  plyr:::adply(as.matrix(m9.df),2,MCMCsum)
setTxtProgressBar(pb, 9)
m10.df    <-  as.data.frame(extract(m10))[,-1]
m10.summ  <-  plyr:::adply(as.matrix(m10.df),2,MCMCsum)
setTxtProgressBar(pb, 10)
m11.df    <-  as.data.frame(extract(m11))[,-1]
m11.summ  <-  plyr:::adply(as.matrix(m11.df),2,MCMCsum)
setTxtProgressBar(pb, 11)
m12.df    <-  as.data.frame(extract(m12))[,-1]
m12.summ  <-  plyr:::adply(as.matrix(m12.df),2,MCMCsum)
setTxtProgressBar(pb, 12)
m13.df    <-  as.data.frame(extract(m13))[,-1]
m13.summ  <-  plyr:::adply(as.matrix(m13.df),2,MCMCsum)
setTxtProgressBar(pb, 13)
m14.df    <-  as.data.frame(extract(m14))[,-1]
m14.summ  <-  plyr:::adply(as.matrix(m14.df),2,MCMCsum)
setTxtProgressBar(pb, 14)
m15.df    <-  as.data.frame(extract(m15))[,-1]
m15.summ  <-  plyr:::adply(as.matrix(m15.df),2,MCMCsum)
setTxtProgressBar(pb, 15)
m16.df    <-  as.data.frame(extract(m16))[,-1]
m16.summ  <-  plyr:::adply(as.matrix(m16.df),2,MCMCsum)
setTxtProgressBar(pb, 16)
m17.df    <-  as.data.frame(extract(m17))[,-1]
m17.summ  <-  plyr:::adply(as.matrix(m17.df),2,MCMCsum)
setTxtProgressBar(pb, 17)
m18.df    <-  as.data.frame(extract(m18))[,-1]
m18.summ  <-  plyr:::adply(as.matrix(m18.df),2,MCMCsum)
setTxtProgressBar(pb, 18)
m19.df    <-  as.data.frame(extract(m19))[,-1]
m19.summ  <-  plyr:::adply(as.matrix(m19.df),2,MCMCsum)
setTxtProgressBar(pb, 19)
m20.df    <-  as.data.frame(extract(m20))[,-1]
m20.summ  <-  plyr:::adply(as.matrix(m20.df),2,MCMCsum)
setTxtProgressBar(pb, 20)
m21.df    <-  as.data.frame(extract(m21))[,-1]
m21.summ  <-  plyr:::adply(as.matrix(m21.df),2,MCMCsum)
setTxtProgressBar(pb, 21)
m22.df    <-  as.data.frame(extract(m22))[,-1]
m22.summ  <-  plyr:::adply(as.matrix(m22.df),2,MCMCsum)
setTxtProgressBar(pb, 22)
m23.df    <-  as.data.frame(extract(m23))[,-1]
m23.summ  <-  plyr:::adply(as.matrix(m23.df),2,MCMCsum)
setTxtProgressBar(pb, 23)
m24.df    <-  as.data.frame(extract(m24))[,-1]
m24.summ  <-  plyr:::adply(as.matrix(m24.df),2,MCMCsum)
setTxtProgressBar(pb, 24)
m25.df    <-  as.data.frame(extract(m25))[,-1]
m25.summ  <-  plyr:::adply(as.matrix(m25.df),2,MCMCsum)
setTxtProgressBar(pb, 25)
close(pb)




#########################################
# LOO Log-likelihood for model selection
#########################################
print('Calculating LOO for model comparison')
pb <- txtProgressBar(min=0,max=25, style=3)
m1LL     <-  extract_log_lik(m1, parameter_name = "log_lik")
m1Loo    <-  loo(m1LL)
m1WAIC   <-  waic(m1LL)
setTxtProgressBar(pb, 1)
m2LL     <-  extract_log_lik(m2, parameter_name = "log_lik")
m2Loo    <-  loo(m2LL)
m2WAIC   <-  waic(m2LL)
setTxtProgressBar(pb, 2)
m3LL     <-  extract_log_lik(m3, parameter_name = "log_lik")
m3Loo    <-  loo(m3LL)
m3WAIC   <-  waic(m3LL)
setTxtProgressBar(pb, 3)
m4LL     <-  extract_log_lik(m4, parameter_name = "log_lik")
m4Loo    <-  loo(m4LL)
m4WAIC   <-  waic(m4LL)
setTxtProgressBar(pb, 4)
m5LL     <-  extract_log_lik(m5, parameter_name = "log_lik")
m5Loo    <-  loo(m5LL)
m5WAIC   <-  waic(m5LL)
setTxtProgressBar(pb, 5)
m6LL     <-  extract_log_lik(m6, parameter_name = "log_lik")
m6Loo    <-  loo(m6LL)
m6WAIC   <-  waic(m6LL)
setTxtProgressBar(pb, 6)
m7LL     <-  extract_log_lik(m7, parameter_name = "log_lik")
m7Loo    <-  loo(m7LL)
m7WAIC   <-  waic(m7LL)
setTxtProgressBar(pb, 7)
m8LL     <-  extract_log_lik(m8, parameter_name = "log_lik")
m8Loo    <-  loo(m8LL)
m8WAIC   <-  waic(m8LL)
setTxtProgressBar(pb, 8)
m9LL     <-  extract_log_lik(m9, parameter_name = "log_lik")
m9Loo    <-  loo(m9LL)
m9WAIC   <-  waic(m9LL)
setTxtProgressBar(pb, 9)
m10LL    <-  extract_log_lik(m10, parameter_name = "log_lik")
m10Loo   <-  loo(m10LL)
m10WAIC  <-  waic(m10LL)
setTxtProgressBar(pb, 10)
m11LL    <-  extract_log_lik(m11, parameter_name = "log_lik")
m11Loo   <-  loo(m11LL)
m11WAIC  <-  waic(m11LL)
setTxtProgressBar(pb, 11)
m12LL    <-  extract_log_lik(m12, parameter_name = "log_lik")
m12Loo   <-  loo(m12LL)
m12WAIC  <-  waic(m12LL)
setTxtProgressBar(pb, 12)
m13LL    <-  extract_log_lik(m13, parameter_name = "log_lik")
m13Loo   <-  loo(m13LL)
m13WAIC  <-  waic(m13LL)
setTxtProgressBar(pb, 13)
m14LL    <-  extract_log_lik(m14, parameter_name = "log_lik")
m14Loo   <-  loo(m14LL)
m14WAIC  <-  waic(m14LL)
setTxtProgressBar(pb, 14)
m15LL    <-  extract_log_lik(m15, parameter_name = "log_lik")
m15Loo   <-  loo(m15LL)
m15WAIC  <-  waic(m15LL)
setTxtProgressBar(pb, 15)
m16LL    <-  extract_log_lik(m16, parameter_name = "log_lik")
m16Loo   <-  loo(m16LL)
m16WAIC  <-  waic(m16LL)
setTxtProgressBar(pb, 16)
m17LL    <-  extract_log_lik(m17, parameter_name = "log_lik")
m17Loo   <-  loo(m17LL)
m17WAIC  <-  waic(m17LL)
setTxtProgressBar(pb, 17)
m18LL    <-  extract_log_lik(m18, parameter_name = "log_lik")
m18Loo   <-  loo(m18LL)
m18WAIC  <-  waic(m18LL)
setTxtProgressBar(pb, 18)
m19LL    <-  extract_log_lik(m19, parameter_name = "log_lik")
m19Loo   <-  loo(m19LL)
m19WAIC  <-  waic(m19LL)
setTxtProgressBar(pb, 19)
m20LL    <-  extract_log_lik(m20, parameter_name = "log_lik")
m20Loo   <-  loo(m20LL)
m20WAIC  <-  waic(m20LL)
setTxtProgressBar(pb, 20)
m21LL    <-  extract_log_lik(m21, parameter_name = "log_lik")
m21Loo   <-  loo(m21LL)
m21WAIC  <-  waic(m21LL)
setTxtProgressBar(pb, 21)
m22LL    <-  extract_log_lik(m22, parameter_name = "log_lik")
m22Loo   <-  loo(m22LL)
m22WAIC  <-  waic(m22LL)
setTxtProgressBar(pb, 22)
m23LL    <-  extract_log_lik(m23, parameter_name = "log_lik")
m23Loo   <-  loo(m23LL)
m23WAIC  <-  waic(m23LL)
setTxtProgressBar(pb, 23)
m24LL    <-  extract_log_lik(m24, parameter_name = "log_lik")
m24Loo   <-  loo(m24LL)
m24WAIC  <-  waic(m24LL)
setTxtProgressBar(pb, 24)
m25LL    <-  extract_log_lik(m25, parameter_name = "log_lik")
m25Loo   <-  loo(m25LL)
m25WAIC  <-  waic(m25LL)
setTxtProgressBar(pb, 25)
close(pb)


m12Z       <-  model.matrix(~ -1 + Run            +
                                   Run : nSperm_z +
                                   Run : Rate     +
                                   Run : EggPos   +
                                   Run : Rate   : EggPos,
                            data = data)
