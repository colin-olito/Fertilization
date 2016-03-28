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
source('R/dependencies.R')

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

blue1  <-  adjustcolor("dodgerblue4",  alpha.f=0.5)
blue2  <-  adjustcolor("dodgerblue4", alpha.f=0.8)
red1   <-  adjustcolor("orangered3",  alpha.f=0.5)
red2   <-  adjustcolor("orangered3", alpha.f=0.8)

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
axis(1, las=1)
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
numSavedSteps  = 5000 #across all chains
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
     lines(xrange, inv_logit(x['mu_a'] + x['theta'] * xrange2), col=transparentColor('grey68',0.01))
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
                    Colony  =  as.numeric(data$Colony),
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
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

mlLogistic4 <- stan(data    =  data.list,
                    file     =  './Stan/logistic-reg-Run-int-slope.stan',
                    chains   =  nChains,
                    iter     =  nIter,
                    thin     =  thinSteps,
                    save_dso =  TRUE
                 )

# Model Results
print(mlLogistic4)
print(mlLogistic4, c("b", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95), digits_summary=3);
print(mlLogistic4, c("a"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95), digits_summary=3);
mlLogistic4.df    <-  as.data.frame(extract(mlLogistic4))
mcmc.mlLogistic4  <-  as.mcmc(mlLogistic4)
mlLogistic4.mcmc  <-  rstan:::as.mcmc.list.stanfit(mlLogistic4)
mlLogistic4.summary <- plyr:::adply(as.matrix(mlLogistic4.df),2,MCMCsum)
(mlLogistic4.summary)
head(mlLogistic4.df)
dim(mlLogistic4.df)

# Simple Diagnostic Plots
plot(mlLogistic4, pars="b")
par(mfrow=c(2,2))
plot(mlLogistic4.mcmc, ask=TRUE)
pairs(mlLogistic4, pars="b")
pairs(mlLogistic4, pars="a")




##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(mlLogistic4.summary$Mean[1] + mlLogistic4.summary$Mean[9]  * nSperm_z),
               Run2  <- inv_logit(mlLogistic4.summary$Mean[2] + mlLogistic4.summary$Mean[10] * nSperm_z),
               Run3  <- inv_logit(mlLogistic4.summary$Mean[3] + mlLogistic4.summary$Mean[11] * nSperm_z),
               Run4  <- inv_logit(mlLogistic4.summary$Mean[4] + mlLogistic4.summary$Mean[12] * nSperm_z),
               Run5  <- inv_logit(mlLogistic4.summary$Mean[5] + mlLogistic4.summary$Mean[13] * nSperm_z),
               Run6  <- inv_logit(mlLogistic4.summary$Mean[6] + mlLogistic4.summary$Mean[14] * nSperm_z),
               Run7  <- inv_logit(mlLogistic4.summary$Mean[7] + mlLogistic4.summary$Mean[15] * nSperm_z),
               Run8  <- inv_logit(mlLogistic4.summary$Mean[8] + mlLogistic4.summary$Mean[16] * nSperm_z)
              )

RegLine  <-  inv_logit(mlLogistic4.summary$Mean[17] + mlLogistic4.summary$Mean[19] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mlLogistic4.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['mu_a'] + x['mu_b'] * xrange2), col=transparentColor('grey50',0.01))
 }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][data$Run == i][order(nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(nSperm_z[data$Run == i])],
                   col='grey75', lwd=3)
 }
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
#  -- Estimate covariance matrix for a x theta
#  Call to STAN:
########################################################










########################################################
#  Logistic regression ~ nSperm. 
#  -- Nested random effects:  Run/Colony
#  -- Random intercept only.
#  Call to STAN:
########################################################

colIndex <- unique(data[c("Run","Colony")])[,"Colony"]
runIndex <- unique(data[c("Run","Colony")])[,"Run"]

model.matrix(~Run/Colony, data=data)

nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

data.list  <-  list(N       =  nrow(data),
                    Nr      =  max(unique(as.numeric(data$Run))),
                    Nc      =  max(unique(as.numeric(data$Colony))),
                    nFert   =  data$nFert, 
                    nEggs   =  data$nEggs,
                    Run     =  as.numeric(data$Run),
                    Colony  =  as.numeric(colIndex),
                    nSperm  =  nSperm_z)

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

nestInt1 <- stan(data    =  data.list,
                 file     =  './Stan/logistic-reg-nest-Run-Colony-int.stan',
                 chains   =  nChains,
                 iter     =  nIter,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                )

# Model Results
print(nestInt1)
print(nestInt1, c("theta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(nestInt1, c("a"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(nestInt1, c("b"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
nestInt1.df    <-  as.data.frame(extract(nestInt1))
mcmc.nestInt1  <-  as.mcmc(nestInt1)
nestInt1.mcmc  <-  rstan:::as.mcmc.list.stanfit(nestInt1)
nestInt1.summ <- plyr:::adply(as.matrix(nestInt1.df),2,MCMCsum)
(nestInt1.summ)
head(nestInt1.df)
dim(nestInt1.df)

# Simple Diagnostic Plots
plot(nestInt1, pars="theta")
par(mfrow=c(2,2))
plot(nestInt1.mcmc, ask=TRUE)
 par(mfrow=c(3,2))
traceplot(nestInt1.mcmc, ask=TRUE)
#pairs(nestInt1, pars="theta")
pairs(nestInt1, pars="a")



##  MAYBE MAKE A PLOT WITH COLONY-SPECIFIC P

##  Plot predicted line etc.
Colonys  <-  list(
               Colony1  <- inv_logit(nestInt1.summ$Mean[9]  + nestInt1.summ$Mean[12] * nSperm_z),
               Colony2  <- inv_logit(nestInt1.summ$Mean[10] + nestInt1.summ$Mean[12] * nSperm_z),
               Colony3  <- inv_logit(nestInt1.summ$Mean[11] + nestInt1.summ$Mean[12] * nSperm_z)
              )
RegLine  <-  inv_logit(nestInt1.summ$Mean[13] + nestInt1.summ$Mean[12] * nSperm_z)




 + nestInt1.summ$Mean[9]
 + nestInt1.summ$Mean[10]
 + nestInt1.summ$Mean[10]
 + nestInt1.summ$Mean[10]
 + nestInt1.summ$Mean[10]
 + nestInt1.summ$Mean[10]
 + nestInt1.summ$Mean[11]
 + nestInt1.summ$Mean[11]

 + nestInt1.summ$Mean[15]

##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(nestInt1.summ$Mean[1] + nestInt1.summ$Mean[12] * nSperm_z),
               Run2  <- inv_logit(nestInt1.summ$Mean[2] + nestInt1.summ$Mean[12] * nSperm_z),
               Run3  <- inv_logit(nestInt1.summ$Mean[3] + nestInt1.summ$Mean[12] * nSperm_z),
               Run4  <- inv_logit(nestInt1.summ$Mean[4] + nestInt1.summ$Mean[12] * nSperm_z),
               Run5  <- inv_logit(nestInt1.summ$Mean[5] + nestInt1.summ$Mean[12] * nSperm_z),
               Run6  <- inv_logit(nestInt1.summ$Mean[6] + nestInt1.summ$Mean[12] * nSperm_z),
               Run7  <- inv_logit(nestInt1.summ$Mean[7] + nestInt1.summ$Mean[12] * nSperm_z),
               Run8  <- inv_logit(nestInt1.summ$Mean[8] + nestInt1.summ$Mean[12] * nSperm_z)
              )

RegLine  <-  inv_logit(nestInt1.summ$Mean[13] + nestInt1.summ$Mean[12] * nSperm_z)


##  Plot showing uncertainty from Run-specific intercepts
par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(nestInt1.df, 1, function(x, data, nSperm_z){
     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['mu_run'] + x['theta'] * xrange2), col=transparentColor('grey68',0.01))
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

##  Plot showing uncertainty from Colony-specific intercepts
par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(nestInt1.df, 1, function(x, data, nSperm_z){
     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['mu_col'] + x['theta'] * xrange2), col=transparentColor('grey68',0.03))
 }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
# for(i in 1:3) {
#   lines(Colonys[[i]][data$Colony == i][order(nSperm_z[data$Colony == i])] ~ data$nSperm[data$Colony == i][order(nSperm_z[data$Colony == i])],
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



















##################################################################
##################################################################
##  SAME ANALYSES, USING MATRIX NOTATION MODELS.
##  ALSO INCLUDING THE MAXIMAL MODEL ESTIMATING COVARIANCE MATRIX. 
##################################################################
##################################################################

# Centered and rescaled nsperm variable for easier estimation

nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)


#############################
## Simple Logistic Regression
#############################

head(data)

X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
str(X)
head(X)

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    X   =  X
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)


## Call to STAN
mat1 <- stan(data     =  data.list,
                 file     =  './Stan/mat-logistic-reg.stan',
                 chains   =  nChains,
                 iter     =  numSavedSteps,
                 thin     =  thinSteps,
                 save_dso =  TRUE
                )


# Model Results
print(mat1)
print(mat1, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat1.df    <-  as.data.frame(extract(mat1))
mcmc.mat1  <-  as.mcmc(mat1)
mat1.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat1)
mat1.summ  <-  plyr:::adply(as.matrix(mat1.df),2,MCMCsum)
(mat1.summ)

# Simple Diagnostic Plots
plot(mat1, pars="beta")
plot(mat1.mcmc, ask=TRUE)
pairs(mat1, pars="beta")


#  LOO Log-likelihood for model selection
mat1LL  <-  extract_log_lik(mat1, parameter_name = "log_lik")
mat1Loo    <-  loo(mat1LL)
mat1WAIC   <-  waic(mat1LL)


##  Plot predicted line etc.
RegLine  <-  inv_logit(mat1.summ$Mean[1] + mat1.summ$Mean[2] * nSperm_z)


##  Plot showing uncertainty from Run-specific intercepts
par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat1.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
 }, data=data, nSperm_z=nSperm_z)
# plot main regression line
lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue4', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)









####################################
## Logistic Mixed Effects Regression
## -- random intercept for RUN
####################################

head(data)

X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
str(X)
head(X)

Z  <-  unname(model.matrix(~ data$Run -1, data=data))
attr(Z,"assign") <- NULL
str(Z)
head(Z)

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    X   =  X,
                    Z   =  Z
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)


## Call to STAN
mat2 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
                )


# Model Results
print(mat2)
print(mat2, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat2, c("gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat2.df    <-  as.data.frame(extract(mat2))
mcmc.mat2  <-  as.mcmc(mat2)
mat2.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat2)
mat2.summ  <-  plyr:::adply(as.matrix(mat2.df),2,MCMCsum)
(mat2.summ)

# Simple Diagnostic Plots
plot(mat2, pars="beta")
plot(mat2, pars="gamma")
plot(mat2.mcmc, ask=TRUE)
pairs(mat2, pars="beta")
pairs(mat2, pars="gamma")


#  LOO Log-likelihood for model selection
mat2LL  <-  extract_log_lik(mat2, parameter_name = "log_lik")
mat2Loo    <-  loo(mat2LL)
mat2WAIC   <-  waic(mat2LL)


##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(mat2.summ$Mean[3]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run2  <- inv_logit(mat2.summ$Mean[4]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run3  <- inv_logit(mat2.summ$Mean[5]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run4  <- inv_logit(mat2.summ$Mean[6]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run5  <- inv_logit(mat2.summ$Mean[7]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run6  <- inv_logit(mat2.summ$Mean[8]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run7  <- inv_logit(mat2.summ$Mean[9]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
               Run8  <- inv_logit(mat2.summ$Mean[10]  + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z)
              )

RegLine  <-  inv_logit(mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat2.df, 1, function(x, data, nSperm_z){
     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
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





####################################
## Logistic Mixed Effects Regression
## -- alternative cell mean specification
## -- random intercept for RUN
####################################

head(data)

X  <-  unname(model.matrix(~ nSperm_z -1, data=data))
attr(X,"assign") <- NULL
str(X)
head(X)

Z  <-  unname(model.matrix(~ data$Run -1, data=data))
attr(Z,"assign") <- NULL
str(Z)
head(Z)

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    X   =  X,
                    Z   =  Z
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)


## Call to STAN
mat2b <- stan(data     =  data.list,
              file     =  './Stan/mat-logistic-1Z-cellmean.stan',
              chains   =  nChains,
              iter     =  numSavedSteps,
              thin     =  thinSteps,
              save_dso =  TRUE
              )


# Model Results
print(mat2b)
print(mat2b, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat2b, c("gamma", "mu_gamma", "sigma_gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat2b.df    <-  as.data.frame(extract(mat2b))
mcmc.mat2b  <-  as.mcmc(mat2b)
mat2b.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat2b)
mat2b.summ  <-  plyr:::adply(as.matrix(mat2b.df),2,MCMCsum)
(mat2b.summ)

# Simple Diagnostic Plots
plot(mat2b, pars="beta")
plot(mat2b, pars="gamma")
plot(mat2b.mcmc, ask=TRUE)
pairs(mat2b, pars="beta")
pairs(mat2b, pars="gamma")


#  LOO Log-likelihood for model selection
mat2bLL  <-  extract_log_lik(mat2b, parameter_name = "log_lik")
mat2bLoo    <-  loo(mat2bLL)
mat2bWAIC   <-  waic(mat2bLL)


##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(mat2b.summ$Mean[2] + mat2b.summ$Mean[1] * nSperm_z),
               Run2  <- inv_logit(mat2b.summ$Mean[3] + mat2b.summ$Mean[1] * nSperm_z),
               Run3  <- inv_logit(mat2b.summ$Mean[4] + mat2b.summ$Mean[1] * nSperm_z),
               Run4  <- inv_logit(mat2b.summ$Mean[5] + mat2b.summ$Mean[1] * nSperm_z),
               Run5  <- inv_logit(mat2b.summ$Mean[6] + mat2b.summ$Mean[1] * nSperm_z),
               Run6  <- inv_logit(mat2b.summ$Mean[7] + mat2b.summ$Mean[1] * nSperm_z),
               Run7  <- inv_logit(mat2b.summ$Mean[8] + mat2b.summ$Mean[1] * nSperm_z),
               Run8  <- inv_logit(mat2b.summ$Mean[9] + mat2b.summ$Mean[1] * nSperm_z)
              )

RegLine  <-  inv_logit(mat2b.summ$Mean[10] + mat2b.summ$Mean[1] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat2b.df, 1, function(x, data, nSperm_z){
     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit(x['mu_gamma'] + x['beta'] * xrange2), col=transparentColor('grey68',0.1))
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













####################################
## Logistic Mixed Effects Regression
## -- random intercept for RUN
## -- random slopes for Run x nSperm
####################################

head(data)

Z0  <-  unname(model.matrix(~ data$Run -1, data=data))[,-c(9:16)]
attr(Z0,"assign") <- NULL
str(Z0)
head(Z0)

Z1  <-  unname(model.matrix(~ data$Run * nSperm_z , data=data))[,-c(1:8)]
attr(Z1,"assign") <- NULL
str(Z1)
Z1[7:nrow(Z1),1]  <-  0
Z1[1:20,]

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K0  =  ncol(Z0),
                    K1  =  ncol(Z1),
                    nT  =  data$nEggs,
                    nS  =  data$nFert,
                    Z0  =  Z0,
                    Z1  =  Z1
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 5000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)


## Call to STAN
mat3 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-allZ.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
             )


# Model Results
print(mat3)
print(mat3, c("gamma0", "mu_gamma0", "sigma_gamma0"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3, c("gamma1", "mu_gamma1", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat3.df    <-  as.data.frame(extract(mat3))
mcmc.mat3  <-  as.mcmc(mat3)
mat3.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat3)
mat3.summ  <-  plyr:::adply(as.matrix(mat3.df),2,MCMCsum)
(mat3.summ)

# Simple Diagnostic Plots
plot(mat3, pars="gamma0")
plot(mat3, pars="gamma1")
plot(mat3.mcmc, ask=TRUE)
pairs(mat3, pars="gamma0")
pairs(mat3, pars="gamma1")


#  LOO Log-likelihood for model selection
mat3LL  <-  extract_log_lik(mat3, parameter_name = "log_lik")
mat3Loo    <-  loo(mat3LL)
mat3WAIC   <-  waic(mat3LL)


##  Plot predicted line etc.
runs  <-  list(
               Run1  <- inv_logit(mat3.summ$Mean[1] + mat3.summ$Mean[9] * nSperm_z),
               Run2  <- inv_logit(mat3.summ$Mean[2] + mat3.summ$Mean[10] * nSperm_z),
               Run3  <- inv_logit(mat3.summ$Mean[3] + mat3.summ$Mean[11] * nSperm_z),
               Run4  <- inv_logit(mat3.summ$Mean[4] + mat3.summ$Mean[12] * nSperm_z),
               Run5  <- inv_logit(mat3.summ$Mean[5] + mat3.summ$Mean[13] * nSperm_z),
               Run6  <- inv_logit(mat3.summ$Mean[6] + mat3.summ$Mean[14] * nSperm_z),
               Run7  <- inv_logit(mat3.summ$Mean[7] + mat3.summ$Mean[15] * nSperm_z),
               Run8  <- inv_logit(mat3.summ$Mean[8] + mat3.summ$Mean[16] * nSperm_z)
              )

RegLine  <-  inv_logit(mat3.summ$Mean[17] + mat3.summ$Mean[18] * nSperm_z)



par(omi=rep(0.3, 4))
plot((data$nFert/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(mat3.df, 1, function(x, data, nSperm_z){
#     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit(x['mu_gamma0'] + x['mu_gamma1'] * xrange2), col=transparentColor('grey68',0.1))
# }, data=data, nSperm_z=nSperm_z)
# plot run-specific regression lines
 for(i in 1:8) {
   lines(runs[[i]][data$Run == i][order(nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(nSperm_z[data$Run == i])],
                   col='grey75', lwd=3)
 }
# plot main regression line
lines(RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='black', lwd=3)
points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
        bg=transparentColor('dodgerblue3', 0.7),
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
axis(2, las=1)
axis(1)













####################################
## Logistic Mixed Effects Regression
## -- random intercept for RUN
## -- random slopes for Run x nSperm
## -- Estimate covariance matrix
####################################


head(data)

X  <-  unname(model.matrix(~ 1 + nSperm_z, data=data))
attr(X,"assign") <- NULL
str(X)
head(X)

Z  <-  unname(model.matrix(~ data$Run *nSperm_z -1, data=data))
attr(Z,"assign") <- NULL
str(Z)
Z[1:20,]

##  Assemble data for stan
data.list  <-  list(N    =  nrow(data),
                    P    =  ncol(X),
                    J    =  max(as.numeric(data$Run)),
                    K    =  ncol(Z),
                    grp  =  as.numeric(data$Run),
                    nT   =  data$nEggs,
                    nS   =  data$nFert,
                    X    =  X,
                    Z    =  Z
                   )

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 2000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

# inits  <-  list(
#   list(L_run    =  matrix(runif(16^2,-0.5,0.5), nrow=16,ncol=16)),
#   list(tau_run  =  runif(16,0.1,2)),
#   list(u        =  matrix(runif((8*16),-1,1),nrow=8,ncol=16)),
#   list(beta     =  runif(2,-1,1)),
#   list(sigma_y  =  runif(1,0.1,2))
#   )

## Call to STAN
mat4b <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
             )


# Model Results
#print(mat4)
print(mat4, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat4.df    <-  as.data.frame(extract(mat4))
mcmc.mat4  <-  as.mcmc(mat4)
mat4.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat4)
mat4.summ  <-  plyr:::adply(as.matrix(mat4.df),2,MCMCsum)
(mat4.summ)

# Explore Correlation structure
corrMat  <-  matrix(mat4.summ[404:659,2], ncol=16,nrow=16)
corrplot(corrMat , method='circle', type='upper')
abline(v=8.5)
abline(h=8.5)

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  }
}

corrplot(corrMat * 15, method='circle', type='upper')
abline(v=8.5)
abline(h=8.5)

# Simple Diagnostic Plots
plot(mat4, pars="beta")
pairs(mat4, pars="beta")


#  LOO Log-likelihood for model selection
mat4LL  <-  extract_log_lik(mat4, parameter_name = "log_lik")
mat4Loo    <-  loo(mat4LL)
mat4WAIC   <-  waic(mat4LL)










#######################
##  MODEL COMPARISONS
#######################

str(mat1Loo)
looDiff   <-  compare(mat1Loo, mat2Loo, mat2bLoo, mat3Loo, mat4Loo)
waicDiff  <-  compare(mat1WAIC, mat2WAIC, mat2bWAIC, mat3WAIC, mat4WAIC)

print(looDiff, digits=4)
print(waicDiff, digits=4)

print(compare(mat1Loo, mat2Loo), digits=6)
print(compare(mat1Loo, mat2bLoo), digits=6)
print(compare(mat1Loo, mat3Loo), digits=6)
print(compare(mat1Loo, mat4Loo), digits=6)
print(compare(mat2Loo, mat2bLoo), digits=6)
print(compare(mat2Loo, mat3Loo), digits=6)
print(compare(mat2Loo, mat4Loo), digits=6)
print(compare(mat2bLoo, mat3Loo), digits=6)
print(compare(mat2bLoo, mat4Loo), digits=6)
print(compare(mat3Loo, mat4Loo), digits=6)

looDiff34   <-  looDiff[1,3] - looDiff[2,3]
looDiff32   <-  looDiff[1,3] - looDiff[3,3]
looDiff32b  <-  looDiff[1,3] - looDiff[4,3]
looDiff31   <-  looDiff[1,3] - looDiff[5,3]
looDiff42   <-  looDiff[2,3] - looDiff[3,3]
looDiff42b  <-  looDiff[2,3] - looDiff[4,3]
looDiff41   <-  looDiff[2,3] - looDiff[5,3]
looDiff22b  <-  looDiff[3,3] - looDiff[4,3]
looDiff21   <-  looDiff[3,3] - looDiff[5,3]
looDiff2b1  <-  looDiff[4,3] - looDiff[5,3]

n  <-  length(mat1Loo$pointwise[,"elpd_loo"])
selooDiff34   <-  sqrt(n * var(mat3Loo$pointwise[,"elpd_loo"]  - mat4Loo$pointwise[,"elpd_loo"]))
selooDiff32   <-  sqrt(n * var(mat3Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"]))
selooDiff32b  <-  sqrt(n * var(mat3Loo$pointwise[,"elpd_loo"]  - mat2bLoo$pointwise[,"elpd_loo"]))
selooDiff31   <-  sqrt(n * var(mat3Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff42   <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"]))
selooDiff42b  <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat2bLoo$pointwise[,"elpd_loo"]))
selooDiff41   <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff22b  <-  sqrt(n * var(mat2Loo$pointwise[,"elpd_loo"]  - mat2bLoo$pointwise[,"elpd_loo"]))
selooDiff21   <-  sqrt(n * var(mat2Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff2b1  <-  sqrt(n * var(mat2bLoo$pointwise[,"elpd_loo"] - mat1Loo$pointwise[,"elpd_loo"]))







LooDiff  <-  cbind(c(looDiff34,looDiff32,looDiff32b,looDiff31,looDiff42,looDiff42b,looDiff41,looDiff22b,looDiff21,looDiff2b1),
                   c(selooDiff34,selooDiff32,selooDiff32b,selooDiff31,selooDiff42,selooDiff42b,selooDiff41,selooDiff22b,selooDiff21,selooDiff2b1))

pDiff34   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[1,1] - 0)/LooDiff[1,2])), 3))
pDiff32   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[2,1] - 0)/LooDiff[2,2])), 3))
pDiff32b  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[3,1] - 0)/LooDiff[3,2])), 3))
pDiff31   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[4,1] - 0)/LooDiff[4,2])), 3))
pDiff42   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[5,1] - 0)/LooDiff[5,2])), 3))
pDiff42b  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[6,1] - 0)/LooDiff[6,2])), 3))
pDiff41   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[7,1] - 0)/LooDiff[7,2])), 3))
pDiff22b  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[8,1] - 0)/LooDiff[8,2])), 3))
pDiff21   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[9,1] - 0)/LooDiff[9,2])), 3))
pDiff2b1  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[10,1] - 0)/LooDiff[10,2])),3))

LooDiff  <-  cbind(LooDiff, c(pDiff34,pDiff32,pDiff32b,pDiff31,pDiff42,pDiff42b,pDiff41,pDiff22b,pDiff21,pDiff2b1))

row.names(LooDiff)  <-  c("mat3 - mat4",
                          "mat3 - mat2",
                          "mat3 - mat2b",
                          "mat3 - mat1",
                          "mat4 - mat2",
                          "mat4 - mat2b",
                          "mat4 - mat1",
                          "mat2 - mat2b",
                          "mat2 - mat1",
                          "mat2b - mat1")
colnames(LooDiff)   <-  c("diff", "se", "p.value")
LooDiff


str(mat3Loo)

##  Plot differences 
m3m4   <-  density(mat3Loo$pointwise[,"elpd_loo"] - mat4Loo$pointwise[,"elpd_loo"])
m3m2   <-  density(mat3Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"])
m3m2b  <-  density(mat3Loo$pointwise[,"elpd_loo"]  - mat2bLoo$pointwise[,"elpd_loo"])
m3m1   <-  density(mat3Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"])

allx   <-  c(m3m4$x,m3m2$x,m3m2b$x,m3m1$x)
ally   <-  c(m3m4$y,m3m2$y,m3m2b$y,m3m1$y)



plot(NA, xlab=expression(paste(Delta[LOOic])), type='n', axes=FALSE, ylab='Density', cex.lab=1.2, 
     xlim=c(min(allx), (max(allx)+0.4*(max(allx) - min(allx, allx)))), 
     ylim=c(0, (max(ally, ally)+0.05*(max(ally, ally) - min(ally, ally)))), yaxs='i')
proportionalLabel(0.5, 1.1, expression(paste('Distribution of', Delta[elpd-loo])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.2)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
polygon(c(m3m4$x), c(m3m4$y),   col=transparentColor('dodgerblue1', 0.5), border='dodgerblue1')

polygon(c(m3m2$x), c(m3m2$y),   col=transparentColor('dodgerblue1', 0.5), border='dodgerblue1')
polygon(c(m3m2b$x), c(m3m2b$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue1')
polygon(c(m3m1$x), c(m3m1$y),   col=transparentColor('dodgerblue3', 0.5), border='dodgerblue1')
abline(v=0, col=transparentColor('red', 0.7), lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9, las=1)







m3  <-  density(mat3Loo$pointwise[,"looic"])
m4  <-  density(mat4Loo$pointwise[,"looic"])
m2  <-  density(mat2Loo$pointwise[,"looic"])
m2b  <-  density(mat2bLoo$pointwise[,"looic"])
m1  <-  density(mat1Loo$pointwise[,"looic"])

allx   <-  c(m3$x,m4$x,m2$x,m2b$x,m1$x)
ally   <-  c(m3$y,m4$y,m2$y,m2b$y,m1$y)

plot(NA, xlab=expression(paste(Delta[LOOic])), type='n', axes=FALSE, ylab='Density', cex.lab=1.2, 
     xlim=c(min(allx), (max(allx)+0.4*(max(allx) - min(allx, allx)))), 
     ylim=c(0, (max(ally, ally)+0.05*(max(ally, ally) - min(ally, ally)))), yaxs='i')
proportionalLabel(0.5, 1.1, expression(paste('Distribution of', Delta[elpd-loo])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.2)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
polygon(c(m3$x), c(m3$y),   col=transparentColor('dodgerblue1', 0.5), border='dodgerblue1')
polygon(c(m4$x), c(m4$y),   col=transparentColor('dodgerblue4', 0.5), border='dodgerblue1')
polygon(c(m2$x), c(m2$y),   col=transparentColor('dodgerblue4', 0.5), border='dodgerblue1')
polygon(c(m2b$x), c(m2b$y),   col=transparentColor('dodgerblue4', 0.5), border='dodgerblue1')
polygon(c(m1$x), c(m1$y),   col=transparentColor('dodgerblue4', 0.5), border='dodgerblue1')
abline(v=0, col=transparentColor('red', 0.7), lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9, las=1)









