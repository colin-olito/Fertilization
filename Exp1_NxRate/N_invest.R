#/* 
# * Colin Olito. Created 21/10/2015.
# * BODY MASS X SPERM MORPHOLOGY CORRELATION
# * NOTES: 
# *          
# */

rm(list=ls())
################
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

#******************
# DEPENDENCIES
source('R/dependencies_N_invest.R')

#*******************
# Import Data
data <- read.csv('data/Ninvest_master.csv', header=TRUE, stringsAsFactors=FALSE)
data <- data.frame(data)
head(data)

# Convert grouping variables to factors; Correct Dates
data$Run       <-  factor(data$Run)
data$Colony    <-  factor(data$Colony)
data$N         <-  factor(data$N)
data$Lane      <-  factor(data$Lane)
data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
data$Date      <-  dmy(data$Date)

str(data)


#########################################
# A few exploratory plots
#########################################

blue1  <-  adjustcolor("dodgerblue3",  alpha.f=0.5)
blue2  <-  adjustcolor("dodgerblue6", alpha.f=1)
red1   <-  adjustcolor("orangered",  alpha.f=0.5)
red2   <-  adjustcolor("orangered3", alpha.f=1)

Runcols   <-  c(blue1, 
                red1,
                adjustcolor("darkolivegreen", alpha.f=1),
                adjustcolor("salmon1", alpha.f=1),
                adjustcolor("purple2", alpha.f=1),
                adjustcolor("magenta", alpha.f=1),
                adjustcolor("grey70", alpha.f=1),
                adjustcolor("orangered3", alpha.f=1)
            ) 

Lanecols   <-  c(adjustcolor("dodgerblue", alpha.f=1),
                 adjustcolor("darkolivegreen", alpha.f=1),
                 adjustcolor("salmon1", alpha.f=1),
                 adjustcolor("purple2", alpha.f=1),
                 adjustcolor("dodgerblue4", alpha.f=1),
                 adjustcolor("red", alpha.f=1)
            ) 

# plot of fertilization rate ~ sperm
par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm, data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(1)
axis(2, las=1)



# plot of fertilization rate ~ sperm, grouped by run
par(omi=rep(0.3, 4))
plot((nFert[Run == 1]/nEggs[Run == 1]) ~ nSperm[Run == 1], data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points((data$nFert[data$Run == 1]/data$nEggs[data$Run == 1]) ~ data$nSperm[data$Run == 1], pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
for (i in 2:max(as.numeric(data$Run))){
  points((data$nFert[data$Run == i]/data$nEggs[data$Run == i]) ~ data$nSperm[data$Run == i], pch=21, 
          bg=Runcols[i],
          col=Runcols[i], cex=1.1)
}
axis(1)
axis(2, las=1)


# plot of fertilization rate ~ sperm, grouped by lane
par(omi=rep(0.3, 4))
plot((nFert[Lane == 1]/nEggs[Lane == 1]) ~ nSperm[Lane == 1], data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points((data$nFert[data$Lane == 1]/data$nEggs[data$Lane == 1]) ~ data$nSperm[data$Lane == 1], pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
for (i in 2:max(as.numeric(data$Lane))){
  points((data$nFert[data$Lane == i]/data$nEggs[data$Lane == i]) ~ data$nSperm[data$Lane == i], pch=21, 
          bg=Lanecols[i],
          col=Lanecols[i], cex=1.1)
}
axis(1)
axis(2, las=1)


# plot of fertilization rate ~ sperm, grouped by colony
par(omi=rep(0.3, 4))
plot((nFert[Colony == 1]/nEggs[Colony == 1]) ~ nSperm[Colony == 1], data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points((data$nFert[data$Colony == 1]/data$nEggs[data$Colony == 1]) ~ data$nSperm[data$Colony == 1], pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
for (i in 2:max(as.numeric(data$Colony))){
  points((data$nFert[data$Colony == i]/data$nEggs[data$Colony == i]) ~ data$nSperm[data$Colony == i], pch=21, 
          bg=Runcols[i],
          col=Runcols[i], cex=1.1)
}
axis(1)
axis(2, las=1)



# plot of fertilization rate ~ Vol
par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ Vol, data=data, 
    xlab=substitute('Volume sperm (mL)'), ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$Vol),max(data$Vol)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points((data$nFert/data$nEggs) ~ data$Vol, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(1)
axis(2, las=1)






########################################################
#  Simple Logistic regression ~ nSperm. Complete Pooling 
#  Call to STAN:
########################################################
nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

data.list  <-  list(N       =  nrow(data),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    nSperm  =  nSperm_z)

#  Options for the analysis
nChains        = 4
burnInSteps    = 0
thinSteps      = 5
numSavedSteps  = 10000 #across all chains
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

Logistic1 <- stan(data    =  data.list,
                 file     =  './Stan/logistic-reg-pool.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(Logistic1)
print(Logistic1, c("theta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
Logistic1.df    <-  as.data.frame(extract(Logistic1))
mcmc.Logistic1  <-  as.mcmc(Logistic1)
Logistic1.mcmc  <-  rstan:::as.mcmc.list.stanfit(Logistic1)
Logistic1.summary <- plyr:::adply(as.matrix(Logistic1.df),2,MCMCsum)
(Logistic1.summary)

# Simple Diagnostic Plots
plot(Logistic1, pars="theta")
par(mfrow=c(2,2))
plot(Logistic1.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(Logistic1.mcmc, ask=TRUE)
pairs(Logistic1, pars="theta")




##  Plot predicted line etc.
pred     <-  Logistic1.summary[-c(1,2,51),]
RegLine  <-  inv_logit(Logistic1.summary$Mean[1] + Logistic1.summary$Mean[2] * nSperm_z)
RegLow   <-  inv_logit(Logistic1.summary$lower[1] + Logistic1.summary$lower[2] * nSperm_z)
RegHi    <-  inv_logit(Logistic1.summary$upper[1] + Logistic1.summary$upper[2] * nSperm_z)

par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
lines((pred$Mean[order(nSperm_z)] / data$nEggs[order(nSperm_z)]) ~ 
        data$nSperm[order(nSperm_z)], type='l', lwd=3)
lines((pred$upper[order(nSperm_z)] / data$nEggs[order(nSperm_z)]) ~ 
        data$nSperm[order(nSperm_z)], type='l', lty=2)
lines((pred$lower[order(nSperm_z)] / data$nEggs[order(nSperm_z)]) ~ 
        data$nSperm[order(nSperm_z)], type='l', lty=2)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)

# RegLine matches pred !!!  YAY!!!
# lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
#                   col='black', lwd=3)
# lines(RegHi[order(nSperm_z)]   ~ data$nSperm[order(nSperm_z)], type='l', lty=2)
# lines(RegLow[order(nSperm_z)]  ~ data$nSperm[order(nSperm_z)], type='l', lty=2)
# cbind(RegHi[order(nSperm_z)], 
#       RegLine[order(nSperm_z)],
#       RegLow[order(nSperm_z)])














########################################################
#  Logistic regression ~ nSperm. Random Intercept ~ Run. 
#  Call to STAN:
########################################################
nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

data.list  <-  list(N       =  nrow(data),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    Run     =  as.numeric(data$Run),
                    nSperm  =  nSperm_z)

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 10000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

mlLogistic1 <- stan(data    =  data.list,
                 file     =  './Stan/logistic-reg-Run.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 warmup   =  burnInSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic1)
print(mlLogistic1, c("theta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mlLogistic1, c("a"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mlLogistic1.df    <-  as.data.frame(extract(mlLogistic1))
mcmc.mlLogistic1  <-  as.mcmc(mlLogistic1)
mlLogistic1.mcmc  <-  rstan:::as.mcmc.list.stanfit(mlLogistic1)
mlLogistic1.summary <- plyr:::adply(as.matrix(mlLogistic1.df),2,MCMCsum)
(mlLogistic1.summary)
head(mlLogistic1.df)
dim(mlLogistic1.df)

# Simple Diagnostic Plots
plot(mlLogistic1, pars="theta")
par(mfrow=c(2,2))
plot(mlLogistic1.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic1.mcmc, ask=TRUE)
#pairs(mlLogistic1, pars="theta")
pairs(mlLogistic1, pars="a")





##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(mlLogistic1.summary$Mean[1] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run2  <- inv_logit(mlLogistic1.summary$Mean[2] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run3  <- inv_logit(mlLogistic1.summary$Mean[3] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run4  <- inv_logit(mlLogistic1.summary$Mean[4] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run5  <- inv_logit(mlLogistic1.summary$Mean[5] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run6  <- inv_logit(mlLogistic1.summary$Mean[6] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run7  <- inv_logit(mlLogistic1.summary$Mean[7] + mlLogistic1.summary$Mean[9] * nSperm_z),
               Run8  <- inv_logit(mlLogistic1.summary$Mean[8] + mlLogistic1.summary$Mean[9] * nSperm_z)
              )

mlLogistic1.df


mlLogistic1.summary[1:11,]
RegLine  <-  inv_logit(mlLogistic1.summary$Mean[11] + mlLogistic1.summary$Mean[9] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mlLogistic1.df, 1, function(x, data, nSperm_z){
     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['mu_a'] + x['theta'] * xrange2), col=transparentColor('grey90',0.01))
 }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
# for(i in 1:8) {
#   lines(runs[[i]][data$Run == i][order(nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(nSperm_z[data$Run == i])],
#                   col='grey75', lwd=3)
# }
# plot main regression line
lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)








########################################################
#  Logistic regression ~ nSperm. Random Intercept ~ Run. 
#  
#  ALTERNATIVE PARAMETERIZATION... NOW, MODELED WITH 
#  AN OVERALL INERCEPT, AND RUN-SPECIFIC DEVIATIONS FROM 
#  THE OVERALL INTERCEPT. GIVES AN IDENTICAL FIT, BUT IS
#  A POORER SPECIFICATION OF THE MODEL BECAUSE IT IMPOSES
#  STRONG CORRELATION BETWEEN THE OVERALL INTERCEPT AND THE 
#  RUN-SPECIFIC DEVIATIONS.  TAKES LONGER TO FIT.
#  Call to STAN:
########################################################
nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

data.list  <-  list(N       =  nrow(data),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    Run     =  as.numeric(data$Run),
                    nSperm  =  nSperm_z)

#  Options for the analysis
nChains        = 3
thinSteps      = 5
numSavedSteps  = 1000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

mlLogistic1b <- stan(data    =  data.list,
                 file     =  './Stan/logistic-reg-Run-Int.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 warmup   =  burnInSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic1b)
print(mlLogistic1b, c("theta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95), digits_summary=4);
print(mlLogistic1b, c("a"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95), digits_summary=4);
mlLogistic1b.df    <-  as.data.frame(extract(mlLogistic1b))
mcmc.mlLogistic1b  <-  as.mcmc(mlLogistic1b)
mlLogistic1b.mcmc  <-  rstan:::as.mcmc.list.stanfit(mlLogistic1b)
mlLogistic1b.summary <- plyr:::adply(as.matrix(mlLogistic1b.df),2,MCMCsum)
head(mlLogistic1b.summary)

plot(mlLogistic1b.df$a.8 ~ mlLogistic1b.df$theta.1)
plot(mlLogistic1b.df$a.1 ~ mlLogistic1b.df$a.2)

x<-y<-numeric(8)
for(i in 1:8) {
x[i]  <-  mean(mlLogistic1b.df[[paste0('a.',i)]] + mlLogistic1b.df$theta.1)
y[i]  <-  mean(mlLogistic1.df[[paste0('a.',i)]])

}
plot(x~y)
abline(0,1)

head(mlLogistic1.df)

is(print(mlLogistic1b))
# Simple Diagnostic Plots
plot(mlLogistic1b, pars="theta")
par(mfrow=c(2,2))
plot(mlLogistic1b.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic1b.mcmc, ask=TRUE)
#pairs(mlLogistic1b, pars="theta")
pairs(mlLogistic1b, pars="a")
pairs(mlLogistic1b, pars="theta")





##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[1] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run2  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[2] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run3  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[3] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run4  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[4] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run5  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[5] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run6  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[6] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run7  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[7] + mlLogistic1b.summary$Mean[10] * nSperm_z),
               Run8  <- inv_logit(mlLogistic1b.summary$Mean[9] + mlLogistic1b.summary$Mean[8] + mlLogistic1b.summary$Mean[10] * nSperm_z)
              )

mlLogistic1b.summary[1:10,]
RegLine  <-  inv_logit(mlLogistic1b.summary$Mean[9]  + mlLogistic1b.summary$Mean[10]  * nSperm_z)


par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mlLogistic1b.df, 1, function(x, data, nSperm_z){
     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['theta.1'] + x['theta.2'] * xrange2), 
           col=transparentColor('black',0.01))
 }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
for(i in 1:8) {
  lines(runs[[i]][data$Run == i][order(nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(nSperm_z[data$Run == i])],
                  col='grey80', lwd=3)
}
# plot overall regression line
lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)















########################################################
#  Logistic regression ~ nSperm. Random Intercept ~ Lane 
#  Call to STAN:
########################################################
nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

data.list  <-  list(N       =  nrow(data),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    Lane    =  as.numeric(data$Lane),
                    nSperm  =  nSperm_z)

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

mlLogistic2 <- stan(data    =  data.list,
                 file     =  './Stan/logistic-reg-Lane.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic2)
print(mlLogistic2, c("theta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mlLogistic2, c("a"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mlLogistic2.df    <-  as.data.frame(extract(mlLogistic2))
mcmc.mlLogistic2  <-  as.mcmc(mlLogistic2)
mlLogistic2.mcmc  <-  rstan:::as.mcmc.list.stanfit(mlLogistic2)
mlLogistic2.summary <- plyr:::adply(as.matrix(mlLogistic2.df),2,MCMCsum)
(mlLogistic2.summary)
head(mlLogistic2.df)
dim(mlLogistic2.df)

# Simple Diagnostic Plots
plot(mlLogistic2, pars="theta")
par(mfrow=c(2,2))
plot(mlLogistic2.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic2.mcmc, ask=TRUE)
#pairs(mlLogistic2, pars="theta")
pairs(mlLogistic2, pars="a")
pairs(mlLogistic2, pars="a")



mlLogistic2.summary[1:10,]
##  Plot predicted line etc.
lanes  <-  list(
               Lane1  <- inv_logit(mlLogistic2.summary$Mean[1] + mlLogistic2.summary$Mean[7] * nSperm_z),
               Lane2  <- inv_logit(mlLogistic2.summary$Mean[2] + mlLogistic2.summary$Mean[7] * nSperm_z),
               Lane3  <- inv_logit(mlLogistic2.summary$Mean[3] + mlLogistic2.summary$Mean[7] * nSperm_z),
               Lane4  <- inv_logit(mlLogistic2.summary$Mean[4] + mlLogistic2.summary$Mean[7] * nSperm_z),
               Lane5  <- inv_logit(mlLogistic2.summary$Mean[5] + mlLogistic2.summary$Mean[7] * nSperm_z),
               Lane6  <- inv_logit(mlLogistic2.summary$Mean[6] + mlLogistic2.summary$Mean[7] * nSperm_z)
              )

mlLogistic2.summary[1:10,]
RegLine  <-  inv_logit(mlLogistic2.summary$Mean[9] + mlLogistic2.summary$Mean[7] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# Plot all regressions from MCMC chains
# apply(mlLogistic2.df, 1, function(x, data, nSperm_z){
#      xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#      xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#      lines(xrange, inv_logit(x['mu_a'] + x['theta'] * xrange2), 
#            col=transparentColor('grey90',0.01))
#  }, data=data, nSperm_z=nSperm_z)
# Plot Lane-specific regression lines
  for(i in 1:8) {
   lines(lanes[[i]][data$Lane == i][order(nSperm_z[data$Lane == i])] ~ data$nSperm[data$Lane == i][order(nSperm_z[data$Lane == i])],
                   col='grey80', lwd=3)
 }
lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)












########################################################
#  Logistic regression ~ nSperm. Random Intercept ~ Colony 
#  Call to STAN:
########################################################
nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

data.list  <-  list(N       =  nrow(data),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    Colony    =  as.numeric(data$Colony),
                    nSperm  =  nSperm_z)

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

mlLogistic3 <- stan(data    =  data.list,
                 file     =  './Stan/logistic-reg-Colony.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic3)
print(mlLogistic3, c("theta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mlLogistic3, c("a"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mlLogistic3.df    <-  as.data.frame(extract(mlLogistic3))
mcmc.mlLogistic3  <-  as.mcmc(mlLogistic3)
mlLogistic3.mcmc  <-  rstan:::as.mcmc.list.stanfit(mlLogistic3)
mlLogistic3.summary <- plyr:::adply(as.matrix(mlLogistic3.df),2,MCMCsum)
(mlLogistic3.summary)
head(mlLogistic3.df)
dim(mlLogistic3.df)

# Simple Diagnostic Plots
plot(mlLogistic3, pars="theta")
par(mfrow=c(2,2))
plot(mlLogistic3.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic3.mcmc, ask=TRUE)
#pairs(mlLogistic3, pars="theta")
pairs(mlLogistic3, pars="a")



mlLogistic3.summary[1:10,]
##  Plot predicted line etc.
Colonys  <-  list(
               Colony1  <- inv_logit(mlLogistic3.summary$Mean[1] + mlLogistic3.summary$Mean[4] * nSperm_z),
               Colony2  <- inv_logit(mlLogistic3.summary$Mean[2] + mlLogistic3.summary$Mean[4] * nSperm_z),
               Colony3  <- inv_logit(mlLogistic3.summary$Mean[3] + mlLogistic3.summary$Mean[4] * nSperm_z)
              )

mlLogistic3.summary[1:10,]
RegLine  <-  inv_logit(mlLogistic3.summary$Mean[6] + mlLogistic3.summary$Mean[4] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# Plot all regressions from MCMC chains
# apply(mlLogistic3.df, 1, function(x, data, nSperm_z){
#      xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#      xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#      lines(xrange, inv_logit(x['mu_a'] + x['theta'] * xrange2), 
#            col=transparentColor('grey60',0.01))
#  }, data=data, nSperm_z=nSperm_z)
# Plot Colony-specific regression lines
  for(i in 1:8) {
   lines(Colonys[[i]][data$Colony == i][order(nSperm_z[data$Colony == i])] ~ data$nSperm[data$Colony == i][order(nSperm_z[data$Colony == i])],
                   col='grey80', lwd=3)
 }
lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)









#########################################
#  TRAINER MODELS!!! 
#########################################
# Start with a simple logistic regression
# models
#########################################

head(data)
glmfit1  <-  glm(cbind(nFert, (nEggs-nFert)) ~ log(nSperm), data=data, 
                 family='binomial')
summary(glmfit1)

# Input variables for STAN


data.list  <-  list(N       =  nrow(data),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    nSperm  =  data$nSperm)

# OPTIONS FOR THE ANALYSES
nChains        = 3
burnInSteps    = 500
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)




#  Grand Mean - complete pooling. 
#  Call to STAN:
mlLogistic1 <- stan(data    =  data.list,
                 file     =  './Stan/hier-logistic-pool.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 warmup   =  burnInSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic1)
mlLogistic1.df <-as.data.frame(extract(mlLogistic1))
mcmc.mlLogistic1 <- as.mcmc(mlLogistic1)
mlLogistic1.mcmc<-rstan:::as.mcmc.list.stanfit(mlLogistic1)

# Simple Diagnostic Plots
plot(mlLogistic1)
par(mfrow=c(2,2))
plot(mlLogistic1.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic1.mcmc, ask=TRUE)
mlLogistic1.summary <- plyr:::adply(as.matrix(mlLogistic1.df),2,MCMCsum)
(mlLogistic1.summary)





#  subject-specific means. NO POOLING. 
#  Call to STAN:
mlLogistic2 <- stan(data    =  data.list,
                 file     =  './Stan/hier-logistic-nopool.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 warmup   =  burnInSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic2)
mlLogistic2.df <-as.data.frame(extract(mlLogistic2))
mcmc.mlLogistic2 <- as.mcmc(mlLogistic2)
mlLogistic2.mcmc<-rstan:::as.mcmc.list.stanfit(mlLogistic2)

# Simple Diagnostic Plots
plot(mlLogistic2)
par(mfrow=c(2,2))
plot(mlLogistic2.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic2.mcmc, ask=TRUE)
mlLogistic2.summary <- plyr:::adply(as.matrix(mlLogistic2.df),2,MCMCsum)
(mlLogistic2.summary)






#  subject-specific means. Partial POOLING. 
#  Using Beta-Binomial parameterization
#  Call to STAN:
mlLogistic3 <- stan(data    =  data.list,
                 file     =  './Stan/hier-logistic.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 warmup   =  burnInSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )

# Model Results
print(mlLogistic3)
mlLogistic3.df <-as.data.frame(extract(mlLogistic3))
mcmc.mlLogistic3 <- as.mcmc(mlLogistic3)
mlLogistic3.mcmc<-rstan:::as.mcmc.list.stanfit(mlLogistic3)
#print(mlLogistic3, c("theta", "kappa", "phi"), probs=c(0.1, 0.5, 0.9));
ss_hier  <-  extract(mlLogistic3);

# Simple Diagnostic Plots
plot(mlLogistic3)
par(mfrow=c(2,2))
plot(mlLogistic3.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(mlLogistic3.mcmc, ask=TRUE)
mlLogistic3.summary <- plyr:::adply(as.matrix(mlLogistic3.df),2,MCMCsum)
(mlLogistic3.summary)


##  Figure from Ben Carpenters vignette showing the fitted
##  values for ϕ and κ on the unconstrained scale, which is
##  the space over which Stan is sampling.
df_bda3_fig_5_3 <- with(ss_hier,
                        data.frame(x = log(phi / (1 - phi)),
                                   y = log(kappa)));
phi_sim <- ss_hier$phi;
kappa_sim <- ss_hier$kappa;
df_bda3_fig_5_3 <- data.frame(x = log(phi_sim / (1 - phi_sim)),
                              y = log(kappa_sim));

plot(y ~ x, data=df_bda3_fig_5_3, 
    xlab='logit(phi) = log(alpha / beta)', 
    ylab='log(kappa) = log(alpha + beta)', 
    type='n', axes=FALSE)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(y ~ x, data=df_bda3_fig_5_3,
       pch=21, 
       bg=transparentColor('dodgerblue3', 0.3),
       col=transparentColor('dodgerblue1', 0.3), cex=1.1)
axis(1)
axis(2, las=1)










# Call to STAN
PoiReg1 <- stan(data    =  data.list,
                 file     =  './Stan/poisson_offset.stan',
                 chains   =  nChains,
                 iter     =  nIter,
#                 warmup   =  burnInSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                 )



# Model Results
print(PoiReg1)
PoiReg1.df <-as.data.frame(extract(PoiReg1))
mcmc.PoiReg1 <- as.mcmc(PoiReg1)
PoiReg1.mcmc<-rstan:::as.mcmc.list.stanfit(PoiReg1)


# Simple Diagnostic Plots
plot(PoiReg1)
par(mfrow=c(2,2))
plot(PoiReg1.mcmc, ask=TRUE)
par(mfrow=c(3,2))
traceplot(PoiReg1.mcmc, ask=TRUE)
PoiReg1.summary <- plyr:::adply(as.matrix(PoiReg1.df),2,MCMCsum)
(PoiReg1.summary)



##  make data for Diego
test = cbind(nFert=data$nFert, nEggs=data$nEggs, nSperm=data$nSperm)
write.csv(data.frame(test), file="dat.csv")

N          <-  nrow(data)
data.list  <-  list(N       =  N,
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    nSperm  =  data$nSperm)

# OPTIONS FOR THE ANALYSIS
nChains = 3
burnInSteps = 500
thinSteps = 5
numSavedSteps = 5000 #across all chains
nIter = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)








#########################################
# THE PRIMARY QUESTION OF INTEREST IS WHETHER
# THERE IS A SIGNIFICANT DIFFERENCE  IN DELTA-VO2
# MAX BETWEEN SPAWNING & NON-SPAWNING INDIVIDUALS
# AFTER RECEIVING KCl INJECTIONS
#########################################

#### Start with a naive ANOVA model for reference ####
fm1 <- lm(delta.vo2max ~ sex*spawn, data=data)

par(mfrow=c(2,2))
plot(fm1)
summary(fm1)


# Planned contrast of interest
# A) Spawn vs NoSpawn for each sex
lsex=levels(data$sex)
cc <- contrast(fm1,
         a = list(spawn = "N", sex = lsex),
         b = list(spawn = "Y", sex = lsex)
         )
print(cc,X=TRUE)
cc$X
summary(glht(fm1, linfct=cc$X))


# PLOTTING
newdata = expand.grid(SEX=levels(data$sex), SPAWN=levels(data$spawn))
coefs = coef(fm1)
Xmat <- model.matrix(~SEX*SPAWN , data=newdata)
fit = t(coefs %*% t(Xmat))
fit
SE = sqrt(diag(Xmat %*% vcov(fm1) %*% t(Xmat)))
CI = qt(0.975, length(data$sex) - 2)*SE
CI
newdata = cbind(newdata, data.frame(fit=fit, se=SE, lower=fit-CI, upper=fit+CI))
newdata

p_naive_anova <-  ggplot(data=newdata, aes(y=fit, x=as.numeric(SPAWN)+c(-0.05,0.05))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05) +
#  geom_linerange(aes(ymin=lower, ymax=upper)) +
  geom_line(aes(linetype=SEX)) +
  geom_point(aes(shape=SEX, fill=SEX), size=4) +
  scale_fill_manual('Sex', values=c('grey','black')) +
  scale_shape_manual('Sex', values=c(21,21,21)) +
  scale_linetype_discrete('Sex') +
  scale_y_continuous(expression(Delta~VO[2]~max)) +
  scale_x_continuous('Spawn', breaks=c(1,2), labels=c('No','Yes'), limits=c(0.5,2.5)) +
  theme_bw() +
  theme(
    text=element_text(size=10),
    axis.title.y=element_text(vjust=2, size=rel(2)),
    axis.title.x=element_text(vjust=-1, size=rel(2)),
    legend.position=c(1,0), legend.justification=c(1,0))
p_naive_anova










########################################
#### Mixed Effects Models ####
# Exploring interactions between sex,spawn, met.date/run
# use lme4, remember to detatch("package:nlme") because of conflicts


# All random effects
fmm1 <- lmer(delta.vo2max ~ sex*spawn +
                  (1|sex:m.date) +
                  (1|spawn:m.date) +
                  (1|sex:spawn:m.date) +
                  (1|m.date/run),
                      data=data)
plot(fmm1)

fmm1.prof <- profile(fmm1)
xyplot(fmm1.prof, aspect=1.3)
confint(fmm1.prof)
splom(fmm1.prof)

summary(fmm1)

# No sex*spawn*met.date interaction
fmm2 <- lmer(delta.vo2max ~ sex*spawn +
                  (1|sex:met.date) +
                  (1|spawn:met.date) +
                  (1|met.date/run),
                      data=data)
plot(fmm2)

fmm2.prof <- profile(fmm2)
xyplot(fmm2.prof, aspect=1.3)
confint(fmm2.prof)
splom(fmm2.prof)

summary(fmm2)

# Sex*met.date random effect only
fmm3 <- lmer(delta.vo2max ~ sex*spawn +
                  (1|sex:m.date) +
                  (1|m.date/run),
                      data=data)
plot(fmm3)

fmm3.prof <- profile(fmm3)
xyplot(fmm3.prof, aspect=1.3)
confint(fmm3.prof)
splom(fmm3.prof)

summary(fmm3)

# Spawn*met.date random effect only
fmm4 <- lmer(delta.vo2max ~ sex*spawn +
                  (1|spawn:m.date) +
                  (1|m.date/run),
                      data=data)
plot(fmm4)

fmm4.prof <- profile(fmm4)
xyplot(fmm4.prof, aspect=1.3)
confint(fmm4.prof)
splom(fmm4.prof)

summary(fmm4)


# met.date/run random effect only
fmm5 <- lmer(delta.vo2max ~ sex*spawn +
                  (1|m.date/run),
                      data=data)
plot(fmm5)

fmm5.prof <- profile(fmm5)
xyplot(fmm5.prof, aspect=1.3)
confint(fmm5.prof)
splom(fmm5.prof)

summary(fmm5)

####
#### Model Comparison to determine usefulness of random effects
#### Using AICc
#library(AICcmodavg, lib.loc="C:/R_home/R-3.1.2/library")

aictab(list(fmm1,fmm2,fmm3,fmm4,fmm5),
       modnames=c("fmm1","fmm2","fmm3","fmm4","fmm5"),
       sort=TRUE
       )

anova(fmm1,fmm2,fmm3,fmm4,fmm5)

# Using parametric bootstrapping
PBmodcomp(fmm1,fmm2)
PBmodcomp(fmm2,fmm3)
PBmodcomp(fmm2,fmm4)
PBmodcomp(fmm3,fmm5)
PBmodcomp(fmm4,fmm5)





# met.date/run random effect only
fmm6 <- lmer(delta.vo2max ~ sex*spawn +
                 (1|m.date) +
                     (1|m.date:run),
                      data=data)




#### A look at the Best Fitting Model... re-fit using REML ####

fmm.Best <- fmm5

# Diagnostic plots
plot(fmm.Best)
fmm.Best.prof <- profile(fmm.Best)
xyplot(fmm.Best.prof, aspect=1.3)
confint(fmm.Best.prof, method="boot")
splom(fmm.Best.prof)
dotplot(ranef(fmm.Best, condVar=TRUE))
qqmath(ranef(fmm.Best, condVar=TRUE))

# Summary
summary(fmm.Best)

# Proportion of variation attributable to between vs. within run effects
(VC <- VarCorr(fmm.Best))
v1 <- as.numeric(attr(VC$'m.date','stddev'))
v2 <- as.numeric(attr(VC$'run:m.date','stddev'))
r <- as.numeric(attr(VC,'sc'))
varattr <- c(v1,v2,r)/ sum(c(v1,v2,r))
names(varattr) <- c("m.date","among run/m.date", "within run")
varattr

# Planned contrast of interest
# A) Spawn vs NoSpawn for each sex
lm1 <- lm(delta.vo2max ~ sex*spawn, data=data)
head(model.matrix(lm1))

# Had trouble getting glht to work with covariate, so used resulting
# contrast matrix to construct my own for the comparison of interest
# (sex*spawn simple effect)
lsex=levels(data$sex)
cc <- contrast(lm1,
         a = list(spawn = "N", sex = lsex),
         b = list(spawn = "Y", sex = lsex)
         )
print(cc,X=TRUE)
cc$X
summary(glht(fmm.Best, linfct=cc$X))



# PLOTTING
newdata = expand.grid(SEX=levels(data$sex), SPAWN=levels(data$spawn))#, coll.date=unique(data$coll.date))
coefs = fixef(fmm.Best)
Xmat <- model.matrix(~SEX*SPAWN , data=newdata)
fit = t(coefs %*% t(Xmat))
fit
SE = sqrt(diag(Xmat %*% vcov(fmm.Best) %*% t(Xmat)))
CI = qt(0.975, length(data$sex) - 2)*SE
CI
newdata = cbind(newdata, data.frame(fit=fit, se=SE, lower=fit-CI, upper=fit+CI))
newdata

p_fmm.Best<-  ggplot(data=newdata, aes(y=fit, x=as.numeric(SPAWN)+c(-0.05,0.05))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05) +
#  geom_linerange(aes(ymin=lower, ymax=upper)) +
  geom_hline(aes(intercept=0),size=0.75,lty=2) +
  geom_line(aes(linetype=SEX)) +
  geom_point(aes(shape=SEX, fill=SEX), size=4) +
  scale_fill_manual('Sex', values=c('grey','black')) +
  scale_shape_manual('Sex', values=c(21,21,21)) +
  scale_linetype_discrete('Sex') +
  scale_y_continuous(expression(Delta~VO[2]~max)) +
  scale_x_continuous('Spawn', breaks=c(1,2), labels=c('No','Yes'), limits=c(0.5,2.5)) +
  theme_bw() +
  theme(
    text=element_text(size=10),
    axis.title.y=element_text(vjust=1.5, size=rel(2)),
    axis.title.x=element_text(vjust=0, size=rel(2)),
    legend.position=c(1,0), legend.justification=c(1,0))
p_fmm.Best

pdf(file="lmer.spawn.pdf",width=7,height=5)
p_fmm.Best
graphics.off()










#### A NOTE ON TWO ALTERNATIVE MODELING APPROACHES ####
#
# i) Instead of focusing on the variance structure
#    arising from m.date, instead mean-adjust
#    sex*spawn interactions by including m.date as
#    a fixed effect, but keep m.date:run variance
#    structure
#
# ii) Same as above, but include m.date as a
#     continuous covariate instead of as a factor.
#
# These alternative modeling approaches did not qualitatively
# alter the outcome of the models. Consequently, I have opted
# for the simpler model above primarily because it requires far
# fewer paramters to be estimated.
#
# See preliminary analyses for full exploration of these models




#### Including GWM as a covariate??? ####


# met.date/run random effect only
fmm10 <- lmer(delta.vo2max ~ sex*spawn*twm +
                  (1|m.date/run),
                      data=data)
plot(fmm10)

fmm10.prof <- profile(fmm10)
xyplot(fmm10.prof, aspect=1.3)
confint(fmm10.prof)
splom(fmm10.prof)

summary(fmm10)


# Planned contrast of interest
# A) Spawn vs NoSpawn for each sex
lm1 <- lm(delta.vo2max ~ sex*spawn*twm, data=data)
head(model.matrix(lm1))

# Had trouble getting glht to work with covariate, so used resulting
# contrast matrix to construct my own for the comparison of interest
# (sex*spawn simple effect)
lsex=levels(data$sex)
ltwm=mean(data$twm)
cc <- contrast(lm1,
         a = list(spawn = "N", sex = lsex, twm=ltwm),
         b = list(spawn = "Y", sex = lsex, twm=ltwm)
         )
print(cc,X=TRUE)
cc$X
summary(glht(fmm10, linfct=cc$X))




vcov(fmm10)

K=rbind("F: Y v N cat"=c(0,0,1,0,0,0,0,0),
         "M: Y v N cat"=c(0,1,0,0,1,0,0,0),
         "F: Y v N"=c(0,0,0,0,0,0,1,0),
         "M: Y v N"=c(0,0,0,0,0,1,0,1),
         "F: twm=0"=c(0,0,0,1,0,0,0.5,0),
         "M: twm=0"=c(0,0,0,1,0,1,0.5,0.5),
         "Overall twm=0"=c(0,0,0,1,0,0.5,0.5,0.5))
summary(glht(fmm10,linfct=K))



lm1 <- lm(delta.vo2max ~ twm, data=data)
lm1 <- lmer(delta.vo2max ~ twm + (1|m.date/run), data=data)
summary(lm1)


cc <- cbind("F: Y v N"=c(0,0,0,0,0,0,1,0),
             "M: Y v N"=c(0,0,0,0,0,1,0,1),
             "Overall twm=0"=c(1,0,0,1,0,0.5,0.5,0.5))

cbind(0,t(cc) %*% contr.treatment(8))






# PLOTTING
newdata = expand.grid(SEX=levels(data$sex), SPAWN=levels(data$spawn))#, coll.date=unique(data$coll.date))
coefs = fixef(fmm10)[c(1,2,3,5)]
Xmat <- model.matrix(~SEX*SPAWN , data=newdata)
fit = t(coefs %*% t(Xmat))
fit
SE = sqrt(diag(Xmat %*% vcov(fmm10)[c(1:3,5),c(1:3,5)] %*% t(Xmat)))
CI = qt(0.975, length(data$sex) - 2)*SE
CI
newdata = cbind(newdata, data.frame(fit=fit, se=SE, lower=fit-CI, upper=fit+CI))
newdata

p_fmm.10 <-  ggplot(data=newdata, aes(y=fit, x=as.numeric(SPAWN)+c(-0.05,0.05))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05) +
#  geom_linerange(aes(ymin=lower, ymax=upper)) +
  geom_line(aes(linetype=SEX)) +
  geom_point(aes(shape=SEX, fill=SEX), size=4) +
  scale_fill_manual('Sex', values=c('grey','black')) +
  scale_shape_manual('Sex', values=c(21,21,21)) +
  scale_linetype_discrete('Sex') +
  scale_y_continuous(expression(Delta~VO[2]~max)) +
  scale_x_continuous('Spawn', breaks=c(1,2), labels=c('No','Yes'), limits=c(0.5,2.5)) +
  theme_bw() +
  theme(
    text=element_text(size=10),
    axis.title.y=element_text(vjust=2, size=rel(2)),
    axis.title.x=element_text(vjust=-1, size=rel(2)),
    legend.position=c(1,0), legend.justification=c(1,0))
p_fmm.10









#*********************************************
#*********************************************
#*********************************************
#*********************************************
#  PER-GAMETE DELTA VO2 ANALYSES

names(data)
head(data)
subdata = subset(data, data$no.gametes != 'NA')
head(subdata)

subdata$pergamete <- subdata$delta.vo2max / subdata$no.gametes
mdat <- subset(subdata, subdata$sex == 'M')
# Remove a couple outlier males
mdat <- mdat[mdat$ind != 54,]
mdat <- mdat[mdat$ind != 59,]

fdat <- subset(subdata, subdata$sex == 'F')



# Some exploratory plots

# delta.VO2 x no.gamete
p.no_dv02<- ggplot(data=subdata, aes(y = delta.vo2max, x=no.gametes)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.no_dv02



p.Mno_dv02<- ggplot(data=mdat, aes(y = delta.vo2max, x=no.gametes)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw()

#pdf(file="maledv02.ng.pdf",width=7,height=5)
p.Mno_dv02
#graphics.off()

p.Fno_dv02<- ggplot(data=fdat, aes(y = delta.vo2max, x=no.gametes)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw()

#pdf(file="femaledv02.ng.pdf",width=7,height=5)
p.Fno_dv02
#graphics.off()


# Per-Gamete cost x Total Wet Mass

p.Mpg_twm<- ggplot(data=mdat, aes(y = pergamete, x=twm)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.Mpg_twm

p.Fpg_twm <- ggplot(data=fdat, aes(y = pergamete, x=twm)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.Fpg_twm

# Per-Gamete cost x Gonad Wet Mass
p.Mpg_gwm<- ggplot(data=mdat, aes(y = pergamete, x=log(gwm))) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.Mpg_gwm

p.Fpg_gwm <- ggplot(data=fdat, aes(y = pergamete, x=log(gwm))) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.Fpg_gwm



# Per-Gamete cost x Test Diameter
p.Mpg_daim<- ggplot(data=mdat, aes(y = pergamete, x=diam.c)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.Mpg_daim

p.Fpg_daim <- ggplot(data=fdat, aes(y = pergamete, x=diam.c)) +
              geom_point() +
              stat_smooth(method='lm') +
              theme_bw() +
              facet_grid(~sex)
p.Fpg_daim


# a few preliminary analyses
lm.mtwm <- lm(pergamete ~ twm, data=mdat)
summary(lm.mtwm)

lm.ftwm <- lm(pergamete ~ twm, data=fdat)
summary(lm.ftwm)

lm.mgwm <- lm(pergamete ~ gwm, data=mdat)
summary(lm.mgwm)

lm.fgwm <- lm(pergamete ~ gwm, data=fdat)
summary(lm.fgwm)

lm.mdiam <- lm(pergamete ~ diam.c, data=mdat)
summary(lm.mdiam)

lm.fdiam <- lm(pergamete ~ diam.c, data=fdat)
summary(lm.fdiam)








# A COMPARABLE MIXED EFFECTS MODEL FOR THESE DATA
# met.date/run random effect only


# SEX x TOTAL WET MASS
fmm.pg1 <- lmer(pergamete ~ sex * twm +
                  (1|m.date/run),
                      data=subdata)
plot(fmm.pg1)

# Diagnostic Plots
fmm.pg1.prof <- profile(fmm.pg1)
xyplot(fmm5.prof, aspect=1.3)
confint(fmm.pg1.prof)
splom(fmm.pg1.prof)
dotplot(ranef(fmm.pg1, condVar=TRUE))
qqmath(ranef(fmm.pg1, condVar=TRUE))

# Summary
summary(fmm.pg1)


# Proportion of variation attributable to between vs. within run effects
(VC <- VarCorr(fmm.pg1))
v1 <- as.numeric(attr(VC$'m.date','stddev'))
v2 <- as.numeric(attr(VC$'run:m.date','stddev'))
r <- as.numeric(attr(VC,'sc'))
varattr <- c(v1,v2,r)/ sum(c(v1,v2,r))
names(varattr) <- c("m.date","among run/m.date", "within run")
varattr

# Planned contrast of interest
# A) twm Slope = 0 for each sex
summary(fmm.pg1)

K <- rbind("F: twm coefficient" = c(0,0,1,0),
            "M: twm coefficient" = c(0,0,0,1))
summary(glht(fmm.pg1, linfct=K))




# SEX x GWM WET MASS
fmm.pg2 <- lmer(pergamete ~ sex * log(gwm) +
                  (1|m.date/run),
                      data=subdata)
plot(fmm.pg2)

# Diagnostic Plots
fmm.pg2.prof <- profile(fmm.pg2)
xyplot(fmm.pg2.prof, aspect=1.3)
confint(fmm.pg2.prof)
splom(fmm.pg2.prof)
dotplot(ranef(fmm.pg2, condVar=TRUE))
qqmath(ranef(fmm.pg2, condVar=TRUE))

# Summary
summary(fmm.pg2)


# Proportion of variation attributable to between vs. within run effects
(VC <- VarCorr(fmm.pg2))
v1 <- as.numeric(attr(VC$'m.date','stddev'))
v2 <- as.numeric(attr(VC$'run:m.date','stddev'))
r <- as.numeric(attr(VC,'sc'))
varattr <- c(v1,v2,r)/ sum(c(v1,v2,r))
names(varattr) <- c("m.date","among run/m.date", "within run")
varattr

# Planned contrast of interest
# A) twm Slope = 0 for each sex
summary(fmm.pg2)

K <- rbind("F: twm coefficient" = c(0,0,1,0),
            "M: twm coefficient" = c(0,0,0,1))
summary(glht(fmm.pg2, linfct=K))










# ????????????????????????????????????????????????
# PLOTTING
newdata = expand.grid(SEX=levels(data$sex))#, coll.date=unique(data$coll.date))
coefs = fixef(fmm.pg1)
Xmat <- model.matrix(~SEX*TWM , data=newdata)
fit = t(coefs %*% t(Xmat))
fit
SE = sqrt(diag(Xmat %*% vcov(fmm.Best) %*% t(Xmat)))
CI = qt(0.975, length(data$sex) - 2)*SE
CI
newdata = cbind(newdata, data.frame(fit=fit, se=SE, lower=fit-CI, upper=fit+CI))
newdata

p_fmm.Best<-  ggplot(data=newdata, aes(y=fit, x=as.numeric(SPAWN)+c(-0.05,0.05))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05) +
#  geom_linerange(aes(ymin=lower, ymax=upper)) +
  geom_hline(aes(intercept=0),size=0.75,lty=2) +
  geom_line(aes(linetype=SEX)) +
  geom_point(aes(shape=SEX, fill=SEX), size=4) +
  scale_fill_manual('Sex', values=c('grey','black')) +
  scale_shape_manual('Sex', values=c(21,21,21)) +
  scale_linetype_discrete('Sex') +
  scale_y_continuous(expression(Delta~VO[2]~max)) +
  scale_x_continuous('Spawn', breaks=c(1,2), labels=c('No','Yes'), limits=c(0.5,2.5)) +
  theme_bw() +
  theme(
    text=element_text(size=10),
    axis.title.y=element_text(vjust=1.5, size=rel(2)),
    axis.title.x=element_text(vjust=0, size=rel(2)),
    legend.position=c(1,0), legend.justification=c(1,0))
p_fmm.Best






head(subdata)
table(subdata$sex)

# DELTA.VO2MAX X NO.GAMETES
fmm.Mng <- lmer(delta.vo2max ~ no.gametes +
                  (1|m.date/run),
                      data=mdat)
plot(fmm.Mng)

# Diagnostic Plots
fmm.ng.prof <- profile(fmm.ng)
xyplot(fmm.ng.prof, aspect=1.3)
confint(fmm.ng.prof)
splom(fmm.ng.prof)
dotplot(ranef(fmm.ng, condVar=TRUE))
qqmath(ranef(fmm.ng, condVar=TRUE))

# Summary
summary(fmm.Mng)


# Proportion of variation attributable to between vs. within run effects
(VC <- VarCorr(fmm.ng))
v1 <- as.numeric(attr(VC$'m.date','stddev'))
v2 <- as.numeric(attr(VC$'run:m.date','stddev'))
r <- as.numeric(attr(VC,'sc'))
varattr <- c(v1,v2,r)/ sum(c(v1,v2,r))
names(varattr) <- c("m.date","among run/m.date", "within run")
varattr

# Planned contrast of interest
# A) twm Slope = 0 for each sex
summary(fmm.ng)

K <- rbind("F: twm coefficient" = c(0,0,1,0),
            "M: twm coefficient" = c(0,0,0,1))
summary(glht(fmm.ng, linfct=K))



k <- rbind("coef"=c(0,1))
summary(glht(fmm.Mng, linfct=k))
