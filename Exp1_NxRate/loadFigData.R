#/* 
# * Colin Olito. Created 08/02/2017
# * 
# * Load all Data Sets and Models necessary
# * for making figures.
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






############################################
# Import NxRate model results from the stan 
# sample_files for further analysis 
############################################
print('Loading NxRate Binomial stanfits')
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




############################################
# Import NxRate Beta-Binomial model results 
# from the stan sample_files for further 
# analysis 
############################################
print('Loading NxRate Beta-Binomial stanfits')
pb <- txtProgressBar(min=0,max=25, style=3)
setTxtProgressBar(pb, 0)
csvFiles  <-  c('./output/StanFits/NxRate_m1BB.csv1',
                './output/StanFits/NxRate_m1BB.csv2',
                './output/StanFits/NxRate_m1BB.csv3')
m1BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 1)

csvFiles  <-  c('./output/StanFits/NxRate_m2BB.csv1',
                './output/StanFits/NxRate_m2BB.csv2',
                './output/StanFits/NxRate_m2BB.csv3')
m2BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 2)

csvFiles  <-  c('./output/StanFits/NxRate_m3BB.csv1',
                './output/StanFits/NxRate_m3BB.csv2',
                './output/StanFits/NxRate_m3BB.csv3')
m3BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 3)

csvFiles  <-  c('./output/StanFits/NxRate_m4BB.csv1',
                './output/StanFits/NxRate_m4BB.csv2',
                './output/StanFits/NxRate_m4BB.csv3')
m4BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 4)

csvFiles  <-  c('./output/StanFits/NxRate_m5BB.csv1',
                './output/StanFits/NxRate_m5BB.csv2',
                './output/StanFits/NxRate_m5BB.csv3')
m5BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 5)

csvFiles  <-  c('./output/StanFits/NxRate_m6BB.csv1',
                './output/StanFits/NxRate_m6BB.csv2',
                './output/StanFits/NxRate_m6BB.csv3')
m6BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 6)

csvFiles  <-  c('./output/StanFits/NxRate_m7BB.csv1',
                './output/StanFits/NxRate_m7BB.csv2',
                './output/StanFits/NxRate_m7BB.csv3')
m7BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 7)

csvFiles  <-  c('./output/StanFits/NxRate_m8BB.csv1',
                './output/StanFits/NxRate_m8BB.csv2',
                './output/StanFits/NxRate_m8BB.csv3')
m8BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 8)

csvFiles  <-  c('./output/StanFits/NxRate_m9BB.csv1',
                './output/StanFits/NxRate_m9BB.csv2',
                './output/StanFits/NxRate_m9BB.csv3')
m9BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 9)

csvFiles  <-  c('./output/StanFits/NxRate_m10BB.csv1',
                './output/StanFits/NxRate_m10BB.csv2',
                './output/StanFits/NxRate_m10BB.csv3')
m10BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 10)

csvFiles  <-  c('./output/StanFits/NxRate_m11BB.csv1',
                './output/StanFits/NxRate_m11BB.csv2',
                './output/StanFits/NxRate_m11BB.csv3')
m11BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 11)

csvFiles  <-  c('./output/StanFits/NxRate_m12BB.csv1',
                './output/StanFits/NxRate_m12BB.csv2',
                './output/StanFits/NxRate_m12BB.csv3')
m12BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 12)

csvFiles  <-  c('./output/StanFits/NxRate_m13BB.csv1',
                './output/StanFits/NxRate_m13BB.csv2',
                './output/StanFits/NxRate_m13BB.csv3')
m13BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 13)

csvFiles  <-  c('./output/StanFits/NxRate_m14BB.csv1',
                './output/StanFits/NxRate_m14BB.csv2',
                './output/StanFits/NxRate_m14BB.csv3')
m14BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 14)

csvFiles  <-  c('./output/StanFits/NxRate_m15BB.csv1',
                './output/StanFits/NxRate_m15BB.csv2',
                './output/StanFits/NxRate_m15BB.csv3')
m15BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 15)

csvFiles  <-  c('./output/StanFits/NxRate_m16BB.csv1',
                './output/StanFits/NxRate_m16BB.csv2',
                './output/StanFits/NxRate_m16BB.csv3')
m16BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 16)

csvFiles  <-  c('./output/StanFits/NxRate_m17BB.csv1',
                './output/StanFits/NxRate_m17BB.csv2',
                './output/StanFits/NxRate_m17BB.csv3')
m17BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 17)

csvFiles  <-  c('./output/StanFits/NxRate_m18BB.csv1',
                './output/StanFits/NxRate_m18BB.csv2',
                './output/StanFits/NxRate_m18BB.csv3')
m18BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 18)

csvFiles  <-  c('./output/StanFits/NxRate_m19BB.csv1',
                './output/StanFits/NxRate_m19BB.csv2',
                './output/StanFits/NxRate_m19BB.csv3')
m19BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 19)

csvFiles  <-  c('./output/StanFits/NxRate_m20BB.csv1',
                './output/StanFits/NxRate_m20BB.csv2',
                './output/StanFits/NxRate_m20BB.csv3')
m20BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 20)

csvFiles  <-  c('./output/StanFits/NxRate_m21BB.csv1',
                './output/StanFits/NxRate_m21BB.csv2',
                './output/StanFits/NxRate_m21BB.csv3')
m21BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 21)

csvFiles  <-  c('./output/StanFits/NxRate_m22BB.csv1',
                './output/StanFits/NxRate_m22BB.csv2',
                './output/StanFits/NxRate_m22BB.csv3')
m22BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 22)

csvFiles  <-  c('./output/StanFits/NxRate_m23BB.csv1',
                './output/StanFits/NxRate_m23BB.csv2',
                './output/StanFits/NxRate_m23BB.csv3')
m23BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 23)

csvFiles  <-  c('./output/StanFits/NxRate_m24BB.csv1',
                './output/StanFits/NxRate_m24BB.csv2',
                './output/StanFits/NxRate_m24BB.csv3')
m24BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 24)

csvFiles  <-  c('./output/StanFits/NxRate_m25BB.csv1',
                './output/StanFits/NxRate_m25BB.csv2',
                './output/StanFits/NxRate_m25BB.csv3')
m25BB        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)
setTxtProgressBar(pb, 25)
close(pb)




################################################
##  Extract Samples into data frames, summarize
################################################
print('Extracting stanfits into data frames, etc.')
pb <- txtProgressBar(min=0,max=25, style=3)
m1BB.df     <-  as.data.frame(extract(m1BB))  [,-1]
m1BB.summ   <-  plyr:::adply(as.matrix(m1BB.df),2,MCMCsum)
setTxtProgressBar(pb, 1)
m2BB.df     <-  as.data.frame(extract(m2BB))[,-1]
m2BB.summ   <-  plyr:::adply(as.matrix(m2BB.df),2,MCMCsum)
setTxtProgressBar(pb, 2)
m3BB.df     <-  as.data.frame(extract(m3BB))[,-1]
m3BB.summ   <-  plyr:::adply(as.matrix(m3BB.df),2,MCMCsum)
setTxtProgressBar(pb, 3)
m4BB.df     <-  as.data.frame(extract(m4BB))[,-1]
m4BB.summ   <-  plyr:::adply(as.matrix(m4BB.df),2,MCMCsum)
setTxtProgressBar(pb, 4)
m5BB.df     <-  as.data.frame(extract(m5BB))[,-1]
m5BB.summ   <-  plyr:::adply(as.matrix(m5BB.df),2,MCMCsum)
setTxtProgressBar(pb, 5)
m6BB.df     <-  as.data.frame(extract(m6BB))[,-1]
m6BB.summ   <-  plyr:::adply(as.matrix(m6BB.df),2,MCMCsum)
setTxtProgressBar(pb, 6)
m7BB.df     <-  as.data.frame(extract(m7BB))[,-1]
m7BB.summ   <-  plyr:::adply(as.matrix(m7BB.df),2,MCMCsum)
setTxtProgressBar(pb, 7)
m8BB.df     <-  as.data.frame(extract(m8BB))[,-1]
m8BB.summ   <-  plyr:::adply(as.matrix(m8BB.df),2,MCMCsum)
setTxtProgressBar(pb, 8)
m9BB.df     <-  as.data.frame(extract(m9BB))[,-1]
m9BB.summ   <-  plyr:::adply(as.matrix(m9BB.df),2,MCMCsum)
setTxtProgressBar(pb, 9)
m10BB.df    <-  as.data.frame(extract(m10BB))[,-1]
m10BB.summ  <-  plyr:::adply(as.matrix(m10BB.df),2,MCMCsum)
setTxtProgressBar(pb, 10)
m11BB.df    <-  as.data.frame(extract(m11BB))[,-1]
m11BB.summ  <-  plyr:::adply(as.matrix(m11BB.df),2,MCMCsum)
setTxtProgressBar(pb, 11)
m12BB.df    <-  as.data.frame(extract(m12BB))[,-1]
m12BB.summ  <-  plyr:::adply(as.matrix(m12BB.df),2,MCMCsum)
setTxtProgressBar(pb, 12)
m13BB.df    <-  as.data.frame(extract(m13BB))[,-1]
m13BB.summ  <-  plyr:::adply(as.matrix(m13BB.df),2,MCMCsum)
setTxtProgressBar(pb, 13)
m14BB.df    <-  as.data.frame(extract(m14BB))[,-1]
m14BB.summ  <-  plyr:::adply(as.matrix(m14BB.df),2,MCMCsum)
setTxtProgressBar(pb, 14)
m15BB.df    <-  as.data.frame(extract(m15BB))[,-1]
m15BB.summ  <-  plyr:::adply(as.matrix(m15BB.df),2,MCMCsum)
setTxtProgressBar(pb, 15)
m16BB.df    <-  as.data.frame(extract(m16BB))[,-1]
m16BB.summ  <-  plyr:::adply(as.matrix(m16BB.df),2,MCMCsum)
setTxtProgressBar(pb, 16)
m17BB.df    <-  as.data.frame(extract(m17BB))[,-1]
m17BB.summ  <-  plyr:::adply(as.matrix(m17BB.df),2,MCMCsum)
setTxtProgressBar(pb, 17)
m18BB.df    <-  as.data.frame(extract(m18BB))[,-1]
m18BB.summ  <-  plyr:::adply(as.matrix(m18BB.df),2,MCMCsum)
setTxtProgressBar(pb, 18)
m19BB.df    <-  as.data.frame(extract(m19BB))[,-1]
m19BB.summ  <-  plyr:::adply(as.matrix(m19BB.df),2,MCMCsum)
setTxtProgressBar(pb, 19)
m20BB.df    <-  as.data.frame(extract(m20BB))[,-1]
m20BB.summ  <-  plyr:::adply(as.matrix(m20BB.df),2,MCMCsum)
setTxtProgressBar(pb, 20)
m21BB.df    <-  as.data.frame(extract(m21BB))[,-1]
m21BB.summ  <-  plyr:::adply(as.matrix(m21BB.df),2,MCMCsum)
setTxtProgressBar(pb, 21)
m22BB.df    <-  as.data.frame(extract(m22BB))[,-1]
m22BB.summ  <-  plyr:::adply(as.matrix(m22BB.df),2,MCMCsum)
setTxtProgressBar(pb, 22)
m23BB.df    <-  as.data.frame(extract(m23BB))[,-1]
m23BB.summ  <-  plyr:::adply(as.matrix(m23BB.df),2,MCMCsum)
setTxtProgressBar(pb, 23)
m24BB.df    <-  as.data.frame(extract(m24BB))[,-1]
m24BB.summ  <-  plyr:::adply(as.matrix(m24BB.df),2,MCMCsum)
setTxtProgressBar(pb, 24)
m25BB.df    <-  as.data.frame(extract(m25BB))[,-1]
m25BB.summ  <-  plyr:::adply(as.matrix(m25BB.df),2,MCMCsum)
setTxtProgressBar(pb, 25)
close(pb)
