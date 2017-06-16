#/* 
# * Colin Olito. Created 03/02/2017
# * 
# * Load all Data and Models for Beta-Binomial Analysis of 
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



#########################################
# LOO Log-likelihood for model selection
#########################################
print('Calculating LOO for model comparison')
pb <- txtProgressBar(min=0,max=25, style=3)
m1BBLL     <-  extract_log_lik(m1BB, parameter_name = "log_lik")
m1BBLoo    <-  loo(m1BBLL)
m1BBWAIC   <-  waic(m1BBLL)
setTxtProgressBar(pb, 1)
m2BBLL     <-  extract_log_lik(m2BB, parameter_name = "log_lik")
m2BBLoo    <-  loo(m2BBLL)
m2BBWAIC   <-  waic(m2BBLL)
setTxtProgressBar(pb, 2)
m3BBLL     <-  extract_log_lik(m3BB, parameter_name = "log_lik")
m3BBLoo    <-  loo(m3BBLL)
m3BBWAIC   <-  waic(m3BBLL)
setTxtProgressBar(pb, 3)
m4BBLL     <-  extract_log_lik(m4BB, parameter_name = "log_lik")
m4BBLoo    <-  loo(m4BBLL)
m4BBWAIC   <-  waic(m4BBLL)
setTxtProgressBar(pb, 4)
m5BBLL     <-  extract_log_lik(m5BB, parameter_name = "log_lik")
m5BBLoo    <-  loo(m5BBLL)
m5BBWAIC   <-  waic(m5BBLL)
setTxtProgressBar(pb, 5)
m6BBLL     <-  extract_log_lik(m6BB, parameter_name = "log_lik")
m6BBLoo    <-  loo(m6BBLL)
m6BBWAIC   <-  waic(m6BBLL)
setTxtProgressBar(pb, 6)
m7BBLL     <-  extract_log_lik(m7BB, parameter_name = "log_lik")
m7BBLoo    <-  loo(m7BBLL)
m7BBWAIC   <-  waic(m7BBLL)
setTxtProgressBar(pb, 7)
m8BBLL     <-  extract_log_lik(m8BB, parameter_name = "log_lik")
m8BBLoo    <-  loo(m8BBLL)
m8BBWAIC   <-  waic(m8BBLL)
setTxtProgressBar(pb, 8)
m9BBLL     <-  extract_log_lik(m9BB, parameter_name = "log_lik")
m9BBLoo    <-  loo(m9BBLL)
m9BBWAIC   <-  waic(m9BBLL)
setTxtProgressBar(pb, 9)
m10BBLL    <-  extract_log_lik(m10BB, parameter_name = "log_lik")
m10BBLoo   <-  loo(m10BBLL)
m10BBWAIC  <-  waic(m10BBLL)
setTxtProgressBar(pb, 10)
m11BBLL    <-  extract_log_lik(m11BB, parameter_name = "log_lik")
m11BBLoo   <-  loo(m11BBLL)
m11BBWAIC  <-  waic(m11BBLL)
setTxtProgressBar(pb, 11)
m12BBLL    <-  extract_log_lik(m12BB, parameter_name = "log_lik")
m12BBLoo   <-  loo(m12BBLL)
m12BBWAIC  <-  waic(m12BBLL)
setTxtProgressBar(pb, 12)
m13BBLL    <-  extract_log_lik(m13BB, parameter_name = "log_lik")
m13BBLoo   <-  loo(m13BBLL)
m13BBWAIC  <-  waic(m13BBLL)
setTxtProgressBar(pb, 13)
m14BBLL    <-  extract_log_lik(m14BB, parameter_name = "log_lik")
m14BBLoo   <-  loo(m14BBLL)
m14BBWAIC  <-  waic(m14BBLL)
setTxtProgressBar(pb, 14)
m15BBLL    <-  extract_log_lik(m15BB, parameter_name = "log_lik")
m15BBLoo   <-  loo(m15BBLL)
m15BBWAIC  <-  waic(m15BBLL)
setTxtProgressBar(pb, 15)
m16BBLL    <-  extract_log_lik(m16BB, parameter_name = "log_lik")
m16BBLoo   <-  loo(m16BBLL)
m16BBWAIC  <-  waic(m16BBLL)
setTxtProgressBar(pb, 16)
m17BBLL    <-  extract_log_lik(m17BB, parameter_name = "log_lik")
m17BBLoo   <-  loo(m17BBLL)
m17BBWAIC  <-  waic(m17BBLL)
setTxtProgressBar(pb, 17)
m18BBLL    <-  extract_log_lik(m18BB, parameter_name = "log_lik")
m18BBLoo   <-  loo(m18BBLL)
m18BBWAIC  <-  waic(m18BBLL)
setTxtProgressBar(pb, 18)
m19BBLL    <-  extract_log_lik(m19BB, parameter_name = "log_lik")
m19BBLoo   <-  loo(m19BBLL)
m19BBWAIC  <-  waic(m19BBLL)
setTxtProgressBar(pb, 19)
m20BBLL    <-  extract_log_lik(m20BB, parameter_name = "log_lik")
m20BBLoo   <-  loo(m20BBLL)
m20BBWAIC  <-  waic(m20BBLL)
setTxtProgressBar(pb, 20)
m21BBLL    <-  extract_log_lik(m21BB, parameter_name = "log_lik")
m21BBLoo   <-  loo(m21BBLL)
m21BBWAIC  <-  waic(m21BBLL)
setTxtProgressBar(pb, 21)
m22BBLL    <-  extract_log_lik(m22BB, parameter_name = "log_lik")
m22BBLoo   <-  loo(m22BBLL)
m22BBWAIC  <-  waic(m22BBLL)
setTxtProgressBar(pb, 22)
m23BBLL    <-  extract_log_lik(m23BB, parameter_name = "log_lik")
m23BBLoo   <-  loo(m23BBLL)
m23BBWAIC  <-  waic(m23BBLL)
setTxtProgressBar(pb, 23)
m24BBLL    <-  extract_log_lik(m24BB, parameter_name = "log_lik")
m24BBLoo   <-  loo(m24BBLL)
m24BBWAIC  <-  waic(m24BBLL)
setTxtProgressBar(pb, 24)
m25BBLL    <-  extract_log_lik(m25BB, parameter_name = "log_lik")
m25BBLoo   <-  loo(m25BBLL)
m25BBWAIC  <-  waic(m25BBLL)
setTxtProgressBar(pb, 25)
close(pb)
