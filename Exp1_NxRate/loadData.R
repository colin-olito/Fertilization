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

# Plot Colors for X^2 Posterior Predictive Check Plots
COLS  <-  c("#009ce0",
            "#f9684a",
            "#8be36a",
            "#7a0058",
            "#9fa3ff",
            "#be0042",
            "#2d1956"
)


############################################
# Import N_invest model results from the stan 
# sample_files for further analysis 
############################################

csvFiles  <-  c('./output/StanFits/N_invest_m1.csv1',
                './output/StanFits/N_invest_m1.csv2',
                './output/StanFits/N_invest_m1.csv3')
NIm1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m2.csv1',
                './output/StanFits/N_invest_m2.csv2',
                './output/StanFits/N_invest_m2.csv3')
NIm2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m2b.csv1',
                './output/StanFits/N_invest_m2b.csv2',
                './output/StanFits/N_invest_m2b.csv3')
NIm2b        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m3.csv1',
                './output/StanFits/N_invest_m3.csv2',
                './output/StanFits/N_invest_m3.csv3')
NIm3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/N_invest_m4.csv1',
                './output/StanFits/N_invest_m4.csv2',
                './output/StanFits/N_invest_m4.csv3')
NIm4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

############################################
# Import NxRate model results from the stan 
# sample_files for further analysis 
############################################

csvFiles  <-  c('./output/StanFits/NxRate_m1.csv1',
                './output/StanFits/NxRate_m1.csv2',
                './output/StanFits/NxRate_m1.csv3')
m1        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m2.csv1',
                './output/StanFits/NxRate_m2.csv2',
                './output/StanFits/NxRate_m2.csv3')
m2        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m3.csv1',
                './output/StanFits/NxRate_m3.csv2',
                './output/StanFits/NxRate_m3.csv3')
m3        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m3a.csv1',
                './output/StanFits/NxRate_m3a.csv2',
                './output/StanFits/NxRate_m3a.csv3')
m3a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m4.csv1',
                './output/StanFits/NxRate_m4.csv2',
                './output/StanFits/NxRate_m4.csv3')
m4        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m4a.csv1',
                './output/StanFits/NxRate_m4a.csv2',
                './output/StanFits/NxRate_m4a.csv3')
m4a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m5.csv1',
                './output/StanFits/NxRate_m5.csv2',
                './output/StanFits/NxRate_m5.csv3')
m5        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

# csvFiles  <-  c('./output/StanFits/NxRate_m6.csv1',
#                 './output/StanFits/NxRate_m6.csv2',
#                 './output/StanFits/NxRate_m6.csv3')
# m6        <-  read_stan_csv(csvFiles, col_major = TRUE)
# rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m6a.csv1',
                './output/StanFits/NxRate_m6a.csv2',
                './output/StanFits/NxRate_m6a.csv3')
m6a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m7.csv1',
                './output/StanFits/NxRate_m7.csv2',
                './output/StanFits/NxRate_m7.csv3')
m7        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m7a.csv1',
                './output/StanFits/NxRate_m7a.csv2',
                './output/StanFits/NxRate_m7a.csv3')
m7a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

# csvFiles  <-  c('./output/StanFits/NxRate_m8.csv1',
#                 './output/StanFits/NxRate_m8.csv2',
#                 './output/StanFits/NxRate_m8.csv3')
# m8        <-  read_stan_csv(csvFiles, col_major = TRUE)
# rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m8a.csv1',
                './output/StanFits/NxRate_m8a.csv2',
                './output/StanFits/NxRate_m8a.csv3')
m8a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)

# csvFiles  <-  c('./output/StanFits/NxRate_m9.csv1',
#                 './output/StanFits/NxRate_m9.csv2',
#                 './output/StanFits/NxRate_m9.csv3')
# m9        <-  read_stan_csv(csvFiles, col_major = TRUE)
# rm(csvFiles)

csvFiles  <-  c('./output/StanFits/NxRate_m9a.csv1',
                './output/StanFits/NxRate_m9a.csv2',
                './output/StanFits/NxRate_m9a.csv3')
m9a        <-  read_stan_csv(csvFiles, col_major = TRUE)
rm(csvFiles)




################################################
##  Extract Samples into data frames, summarize
################################################
m1.df     <-  as.data.frame(extract(m1))  
m1.summ   <-  plyr:::adply(as.matrix(m1.df),2,MCMCsum)[-1,]
m2.df     <-  as.data.frame(extract(m2))
m2.summ   <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)[-1,]
m3.df     <-  as.data.frame(extract(m3))
m3.summ   <-  plyr:::adply(as.matrix(m3.df),2,MCMCsum)[-1,]
# m3.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3)
m3a.df    <-  as.data.frame(extract(m3a))
m3a.summ  <-  plyr:::adply(as.matrix(m3a.df),2,MCMCsum)[-1,]
# m3a.mcmc  <-  rstan:::as.mcmc.list.stanfit(m3a)
m4.df     <-  as.data.frame(extract(m4))
m4.summ   <-  plyr:::adply(as.matrix(m4.df),2,MCMCsum)[-1,]
m4a.df    <-  as.data.frame(extract(m4a))
m4a.summ  <-  plyr:::adply(as.matrix(m4a.df),2,MCMCsum)[-1,]
m5.df     <-  as.data.frame(extract(m5))
m5.summ   <-  plyr:::adply(as.matrix(m5.df),2,MCMCsum)[-1,]
# m6.df    <-  as.data.frame(extract(m6))
# m6.summ  <-  plyr:::adply(as.matrix(m6.df),2,MCMCsum)[-1,]
m6a.df    <-  as.data.frame(extract(m6a))
m6a.summ  <-  plyr:::adply(as.matrix(m6a.df),2,MCMCsum)[-1,]
m7.df     <-  as.data.frame(extract(m7))
m7.summ   <-  plyr:::adply(as.matrix(m7.df),2,MCMCsum)[-1,]
m7a.df    <-  as.data.frame(extract(m7a))
m7a.summ  <-  plyr:::adply(as.matrix(m7a.df),2,MCMCsum)[-1,]
# m8.df    <-  as.data.frame(extract(m8))
# m8.summ  <-  plyr:::adply(as.matrix(m8.df),2,MCMCsum)[-1,]
m8a.df    <-  as.data.frame(extract(m8a))
m8a.summ  <-  plyr:::adply(as.matrix(m8a.df),2,MCMCsum)[-1,]
# m9.df     <-  as.data.frame(extract(m9))
# m9.summ   <-  plyr:::adply(as.matrix(m9.df),2,MCMCsum)[-1,]
m9a.df    <-  as.data.frame(extract(m9a))
m9a.summ  <-  plyr:::adply(as.matrix(m9a.df),2,MCMCsum)[-1,]



#########################################
# LOO Log-likelihood for model selection
#########################################

m1LL     <-  extract_log_lik(m1, parameter_name = "log_lik")
m1Loo    <-  loo(m1LL)
m1WAIC   <-  waic(m1LL)
m2LL     <-  extract_log_lik(m2, parameter_name = "log_lik")
m2Loo    <-  loo(m2LL)
m2WAIC   <-  waic(m2LL)
m3LL     <-  extract_log_lik(m3, parameter_name = "log_lik")
m3Loo    <-  loo(m3LL)
m3WAIC   <-  waic(m3LL)
m3aLL    <-  extract_log_lik(m3a, parameter_name = "log_lik")
m3aLoo   <-  loo(m3aLL)
m3aWAIC  <-  waic(m3aLL)
m4LL     <-  extract_log_lik(m4, parameter_name = "log_lik")
m4Loo    <-  loo(m4LL)
m4WAIC   <-  waic(m4LL)
m4aLL    <-  extract_log_lik(m4a, parameter_name = "log_lik")
m4aLoo   <-  loo(m4aLL)
m4aWAIC  <-  waic(m4aLL)
m5LL     <-  extract_log_lik(m5, parameter_name = "log_lik")
m5Loo    <-  loo(m5LL)
m5WAIC   <-  waic(m5LL)
# m6LL     <-  extract_log_lik(m6, parameter_name = "log_lik")
# m6Loo    <-  loo(m6LL)
# m6WAIC   <-  waic(m6LL)
m6aLL    <-  extract_log_lik(m6a, parameter_name = "log_lik")
m6aLoo   <-  loo(m6aLL)
m6aWAIC  <-  waic(m6aLL)
m7LL     <-  extract_log_lik(m7, parameter_name = "log_lik")
m7Loo    <-  loo(m7LL)
m7WAIC   <-  waic(m7LL)
m7aLL    <-  extract_log_lik(m7a, parameter_name = "log_lik")
m7aLoo   <-  loo(m7aLL)
m7aWAIC  <-  waic(m7aLL)
# m8LL     <-  extract_log_lik(m8, parameter_name = "log_lik")
# m8Loo    <-  loo(m8LL)
# m8WAIC   <-  waic(m8LL)
m8aLL    <-  extract_log_lik(m8a, parameter_name = "log_lik")
m8aLoo   <-  loo(m8aLL)
m8aWAIC  <-  waic(m8aLL)
# m9LL     <-  extract_log_lik(m9, parameter_name = "log_lik")
# m9Loo    <-  loo(m9LL)
# m9WAIC   <-  waic(m9LL)
m9aLL    <-  extract_log_lik(m9a, parameter_name = "log_lik")
m9aLoo   <-  loo(m9aLL)
m9aWAIC  <-  waic(m9aLL)





#####################################################
# Chi-squared goodness of fit measure of discrepancy
# for simulated data ~ real data
#####################################################
X2data1   <-  as.numeric(m1.df[,6379])
X2sim1    <-  as.numeric(m1.df[,6380])
X2data2   <-  as.numeric(m2.df[,4069])
X2sim2    <-  as.numeric(m2.df[,4070])
X2data3   <-  as.numeric(m3.df[,2359])
X2sim3    <-  as.numeric(m3.df[,2360])
X2data3a  <-  as.numeric(m3a.df[,759])
X2sim3a   <-  as.numeric(m3a.df[,760])
X2data4   <-  as.numeric(m4.df[,1249])
X2sim4    <-  as.numeric(m4.df[,1250])
X2data4a  <-  as.numeric(m4a.df[,629])
X2sim4a   <-  as.numeric(m4a.df[,630])
X2data5    <-  as.numeric(m5.df[,618])
X2sim5    <-  as.numeric(m5.df[,619])
# X2data6   <-  as.numeric(m6.df[,2359])
# X2sim6    <-  as.numeric(m6.df[,2360])
X2data6a  <-  as.numeric(m6a.df[,779])
X2sim6a   <-  as.numeric(m6a.df[,780])
X2data7   <-  as.numeric(m7.df[,4069])
X2sim7    <-  as.numeric(m7.df[,4070])
X2data7a  <-  as.numeric(m7a.df[,769])
X2sim7a   <-  as.numeric(m7a.df[,770])
# X2data8   <-  as.numeric(m8.df[,2359])
# X2sim8    <-  as.numeric(m8.df[,2360])
X2data8a  <-  as.numeric(m8a.df[,779])
X2sim8a   <-  as.numeric(m8a.df[,780])
# X2data9   <-  as.numeric(m9.df[,2359])
# X2sim9    <-  as.numeric(m9.df[,2360])
X2data9a  <-  as.numeric(m9a.df[,759])
X2sim9a   <-  as.numeric(m9a.df[,760])
