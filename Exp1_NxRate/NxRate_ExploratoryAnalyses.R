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
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.7), cex=1.1)
axis(1, las=1)
axis(2, las=1)



# plot of fertilization rate ~ sperm X Rate
par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm, data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == 'Fast'] ~ data$nSperm[data$Rate == 'Fast'], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.7), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == 'Slow'] ~ data$nSperm[data$Rate == 'Slow'], pch=21, 
        bg=transparentColor('orangeRed1', 0.7),
        col=transparentColor('orangeRed3', 0.7), cex=1.1)
axis(1, las=1)
axis(2, las=1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  21,
          pt.bg   =  c(transparentColor('dodgerblue1',0.7), transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue1', 'orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



# plot of fertilization rate ~ sperm X Rate x EggPos
par(omi=rep(0.3, 4))
plot((nFert/nEggs) ~ nSperm, data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == 'Fast' & data$EggPos == 5] ~ data$nSperm[data$Rate == 'Fast' & data$EggPos == 5], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue3', 0.7), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == 'Slow' & data$EggPos == 5] ~ data$nSperm[data$Rate == 'Slow' & data$EggPos == 5], pch=21, 
        bg=transparentColor('OrangeRed1', 0.7),
        col=transparentColor('OrangeRed3', 0.7), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == 'Fast' & data$EggPos == 55] ~ data$nSperm[data$Rate == 'Fast' & data$EggPos == 55], pch=21, 
        col=transparentColor('dodgerblue1', 0.7), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == 'Slow' & data$EggPos == 55] ~ data$nSperm[data$Rate == 'Slow' & data$EggPos == 55], pch=21, 
        col=transparentColor('OrangeRed1', 0.7), cex=1.1)
axis(1, las=1)
axis(2, las=1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),NA,transparentColor('orangered1',0.7),NA),
          col     =  c('dodgerblue3','dodgerblue3','orangered3','orangered3'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )


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





##################################################################
##################################################################
##  MATRIX NOTATION MODELS.
##################################################################
##################################################################

# Centered and rescaled nsperm variable for easier estimation

nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)


################################
## Simple Logistic Regression
##  -- nSperm_z x Rate x eggPos
##  -- nSperm_z is continuous
################################

head(data)

X  <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X  <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
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


# Compare with frequentest GLM
mer1  <-  glm(cbind(nFert,nEggs) ~ nSperm_z*Rate*EggPos, 
                family='binomial', data=data)
summary(mer1)$coefficients
print(mat1, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


#  LOO Log-likelihood for model selection
mat1LL  <-  extract_log_lik(mat1, parameter_name = "log_lik")
mat1Loo    <-  loo(mat1LL)
mat1WAIC   <-  waic(mat1LL)


##  Calculate Predicted Lines
Xnames[1]  <-  "RateFast(Intercept)"
Xnames
Coeff  <-  mat1.summ$Mean
Fast5_RegLine   <-  inv_logit(Coeff[1] + Coeff[2] * nSperm_z)
Fast55_RegLine  <-  inv_logit((Coeff[1] + Coeff[4]) + (Coeff[2] + Coeff[6]) * nSperm_z)
Slow5_RegLine   <-  inv_logit((Coeff[1] + Coeff[3]) + (Coeff[2] + Coeff[5]) * nSperm_z)
Slow55_RegLine  <-  inv_logit((Coeff[1] + Coeff[3] + Coeff[7]) + (Coeff[2] + Coeff[5] + Coeff[8]) * nSperm_z)
Fast_RegLine    <-  inv_logit((Coeff[1] + (Coeff[4])/2) + (Coeff[2] + (Coeff[6])/2) * nSperm_z)
Slow_RegLine    <-  inv_logit((Coeff[1] + Coeff[3] + (0.5*(Coeff[7]))) + (Coeff[2] + Coeff[5] + (0.5*(Coeff[8]))) * nSperm_z)


##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )




##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
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
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue3',0.01))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered3',0.01))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(16,16),
          col     =  c('dodgerblue1','orangered1'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )







####################################
## Logistic Mixed Effects Regression
##  -- nSperm_z x Rate x eggPos
##  -- nSperm_z is continuous
## -- random intercept for RUN
####################################

head(data)

X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

Z       <-  model.matrix(~ -1 + data$Run, data=data)
Znames  <-  dimnames(Z)[[2]]
Z       <-  unname(Z)
attr(Z,"assign") <- NULL
str(Z)
head(Z)

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
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
mat2.summ[1:16,1:2]
mat1.summ[1:16,1:2]
Xnames


##  Have to think about what 'run effect' I want to visualize here... 
##  in the meantime... will leave out and just plot overall regressions
#
# runs  <-  list(
#                Run1   <- inv_logit(mat2.summ$Mean[3]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run2   <- inv_logit(mat2.summ$Mean[4]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run3   <- inv_logit(mat2.summ$Mean[5]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run4   <- inv_logit(mat2.summ$Mean[6]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run5   <- inv_logit(mat2.summ$Mean[7]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run6   <- inv_logit(mat2.summ$Mean[8]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run7   <- inv_logit(mat2.summ$Mean[9]   + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run8   <- inv_logit(mat2.summ$Mean[10]  + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run9   <- inv_logit(mat2.summ$Mean[10]  + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z),
#                Run10  <- inv_logit(mat2.summ$Mean[10]  + mat2.summ$Mean[1] + mat2.summ$Mean[2] * nSperm_z)
#               )

##  Calculate Predicted Lines
Xnames[1]  <-  "RateFast(Intercept)"
Xnames
Coeff2  <-  mat2.summ$Mean
Fast5_RegLine   <-  inv_logit(Coeff2[1] + Coeff2[2] * nSperm_z)
Fast55_RegLine  <-  inv_logit((Coeff2[1] + Coeff2[4]) + (Coeff2[2] + Coeff2[6]) * nSperm_z)
Slow5_RegLine   <-  inv_logit((Coeff2[1] + Coeff2[3]) + (Coeff2[2] + Coeff2[5]) * nSperm_z)
Slow55_RegLine  <-  inv_logit((Coeff2[1] + Coeff2[3] + Coeff2[7]) + (Coeff2[2] + Coeff2[5] + Coeff2[8]) * nSperm_z)
Fast_RegLine    <-  inv_logit((Coeff2[1] + (Coeff2[4])/2) + (Coeff2[2] + (Coeff2[6])/2) * nSperm_z)
Slow_RegLine    <-  inv_logit((Coeff2[1] + Coeff2[3] + (0.5*(Coeff2[7]))) + (Coeff2[2] + Coeff2[5] + (0.5*(Coeff2[8]))) * nSperm_z)

Fast5_RegLine   <-  inv_logit(Coeff2[1] + Coeff2[2] * nSperm_z)
plot(density(mat2.df[,2]))
Fast55_RegLine  <-  inv_logit((Coeff2[1] + Coeff2[4]) + (Coeff2[2] + Coeff2[6]) * nSperm_z)
plot(density(mat2.df[,6]), lwd=3)
abline(v=0, lwd=3,col=2)
sum(mat2.df[,6] - 0)/length(mat2.df[,6])

##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.1),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.1),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat2.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue1',0.01))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered1',0.01))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )


# Compare with frequentest GLM
mer2  <-  glmer(cbind(nFert,nEggs) ~ nSperm_z*Rate*EggPos + (1|Run), 
                family='binomial', data=data)
summary(mer2)$coefficients
print(mat2, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));



####################################
## Logistic Mixed Effects Regression
## -- random intercept for RUN
## -- random slopes for Run x nSperm
####################################

head(data)

## Fixed Effects Model Matrix
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

##  Random Slopes Model Matrix
Z0       <-  model.matrix(~ -1 + data$Run, data=data)
Z0names  <-  dimnames(Z0)[[2]]
Z0       <-  unname(Z0)
attr(Z0,"assign") <- NULL
str(Z0)
head(Z0)

##  Random Intercepts Model Matrix
Z1  <-  model.matrix(~ -1 +  data$Run * nSperm_z , data=data)[,-c(1:10)]
head(Z1)
Z1names  <-  dimnames(Z1)[[2]]
Z1       <-  unname(Z1)
attr(Z1,"assign") <- NULL
str(Z1)
Z1[13:nrow(Z1),1]  <-  0
Z1[1:20,]

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z0),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
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
             file     =  './Stan/mat-logistic-2Z.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
             )


# Model Results
print(mat3)
print(mat3, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3, c("gamma0", "sigma_gamma0"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3, c("gamma1", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat3.df    <-  as.data.frame(extract(mat3))
mcmc.mat3  <-  as.mcmc(mat3)
mat3.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat3)
mat3.summ  <-  plyr:::adply(as.matrix(mat3.df),2,MCMCsum)
(mat3.summ)


# Simple Diagnostic Plots
plot(mat3, pars="gamma0")
plot(mat3, pars="gamma1")
plot(mat3.mcmc, ask=TRUE)
pairs(mat3, pars="beta")
pairs(mat3, pars="gamma0")
pairs(mat3, pars="gamma1")


#  LOO Log-likelihood for model selection
mat3LL  <-  extract_log_lik(mat3, parameter_name = "log_lik")
mat3Loo    <-  loo(mat3LL)
mat3WAIC   <-  waic(mat3LL)


##  Plot predicted line etc.
# runs  <-  list(
#                Run1   <- inv_logit(mat3.summ$Mean[1] + mat3.summ$Mean[7]  * nSperm_z),
#                Run2   <- inv_logit(mat3.summ$Mean[2] + mat3.summ$Mean[8]  * nSperm_z),
#                Run3   <- inv_logit(mat3.summ$Mean[3] + mat3.summ$Mean[9]  * nSperm_z),
#                Run4   <- inv_logit(mat3.summ$Mean[4] + mat3.summ$Mean[10] * nSperm_z),
#                Run5   <- inv_logit(mat3.summ$Mean[5] + mat3.summ$Mean[11] * nSperm_z),
#                Run6   <- inv_logit(mat3.summ$Mean[6] + mat3.summ$Mean[12] * nSperm_z),
#                Run7   <- inv_logit(mat3.summ$Mean[7] + mat3.summ$Mean[15] * nSperm_z),
#                Run8   <- inv_logit(mat3.summ$Mean[8] + mat3.summ$Mean[16] * nSperm_z),
#                Run9   <- inv_logit(mat3.summ$Mean[8] + mat3.summ$Mean[16] * nSperm_z),
#                Run10  <- inv_logit(mat3.summ$Mean[8] + mat3.summ$Mean[16] * nSperm_z),
#               )


##  Calculate Predicted Lines
Xnames[1]  <-  "RateFast(Intercept)"
Xnames
Coeff3  <-  mat3.summ$Mean

Fast5_RegLine   <-  inv_logit(Coeff3[1] + Coeff3[2] * nSperm_z)
Fast55_RegLine  <-  inv_logit((Coeff3[1] + Coeff3[4]) + (Coeff3[2] + Coeff3[6]) * nSperm_z)
Slow5_RegLine   <-  inv_logit((Coeff3[1] + Coeff3[3]) + (Coeff3[2] + Coeff3[5]) * nSperm_z)
Slow55_RegLine  <-  inv_logit((Coeff3[1] + Coeff3[3] + Coeff3[7]) + (Coeff3[2] + Coeff3[5] + Coeff3[8]) * nSperm_z)
Fast_RegLine    <-  inv_logit((Coeff3[1] + (Coeff3[4])/2) + (Coeff3[2] + (Coeff3[6])/2) * nSperm_z)
Slow_RegLine    <-  inv_logit((Coeff3[1] + Coeff3[3] + (0.5*(Coeff3[7]))) + (Coeff3[2] + Coeff3[5] + (0.5*(Coeff3[8]))) * nSperm_z)




##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat3.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue1',0.01))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered1',0.01))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )


# Compare with frequentest GLM
mer3  <-  glmer(cbind(nFert,nEggs) ~ nSperm_z*Rate*EggPos + (nSperm_z|Run), 
                family='binomial', data=data)
summary(mer3)$coefficients
print(mat3, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
coefficients(mer3)
print(mat3, c("gamma0", "sigma_gamma0"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3, c("gamma1", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

Coeff3_freq  <- summary(mer3)$coefficients[,1]
Fast5_freq   <-  inv_logit(Coeff3_freq[1] + Coeff3_freq[2] * nSperm_z)
Fast55_freq  <-  inv_logit((Coeff3_freq[1] + Coeff3_freq[4]) + (Coeff3_freq[2] + Coeff3_freq[6]) * nSperm_z)
Slow5_freq   <-  inv_logit((Coeff3_freq[1] + Coeff3_freq[3]) + (Coeff3_freq[2] + Coeff3_freq[5]) * nSperm_z)
Slow55_freq  <-  inv_logit((Coeff3_freq[1] + Coeff3_freq[3] + Coeff3_freq[7]) + (Coeff3_freq[2] + Coeff3_freq[5] + Coeff3_freq[8]) * nSperm_z)
Fast_freq    <-  inv_logit((Coeff3_freq[1] + (Coeff3_freq[4])/2) + (Coeff3_freq[2] + (Coeff3_freq[6])/2) * nSperm_z)
Slow_freq    <-  inv_logit((Coeff3_freq[1] + Coeff3_freq[3] + (0.5*(Coeff3_freq[7]))) + (Coeff3_freq[2] + Coeff3_freq[5] + (0.5*(Coeff3_freq[8]))) * nSperm_z)


##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), main='Frequentist Version',
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_freq[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_freq[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_freq[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_freq[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



####################################
## Logistic Mixed Effects Regression
##    * alternative parameterization
##    * using 1 Z matrix, with both 
##    * intercepts and slopes together;
##    * Check for equivalence with the 
##    * 2 Z matrix version (mat3)
## -- random intercept for RUN
## -- random slopes for Run x nSperm
####################################


head(data)

X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z, data=data)
Znames  <-  dimnames(Z)[[2]]
Z       <-  unname(Z)
attr(Z,"assign") <- NULL
str(Z)
head(Z)
Z[13:nrow(Z1),11]  <-  0
Z[1:20,]

##  Assemble data for stan
data.list  <-  list(N   =  nrow(data),
                    P   =  ncol(X), 
                    K   =  ncol(Z),
                    nT  =  data$nEggs,
                    nS  =  data$nFert - data$nControlFert,
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
mat3b <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-int-slope.stan',
             chains   =  nChains,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
                )


# Model Results
print(mat3b)
print(mat3b, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3b, c("gamma", "sigma_gamma0", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat3b.df    <-  as.data.frame(extract(mat3b))
mcmc.mat3b  <-  as.mcmc(mat3b)
mat3b.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat3b)
mat3b.summ  <-  plyr:::adply(as.matrix(mat3b.df),2,MCMCsum)
(mat3b.summ)

# Simple Diagnostic Plots
plot(mat3b, pars="beta")
plot(mat3b, pars="gamma")
plot(mat3b.mcmc, ask=TRUE)
pairs(mat3b, pars="beta")
pairs(mat3b, pars="gamma")


#  LOO Log-likelihood for model selection
mat3bLL  <-  extract_log_lik(mat3b, parameter_name = "log_lik")
mat3bLoo    <-  loo(mat3bLL)
mat3bWAIC   <-  waic(mat3bLL)


# Compare with previous parameterization
print(mat3, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3b, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3, c("gamma0", "sigma_gamma0", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat3b, c("gamma", "sigma_gamma0", "sigma_gamma1"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));



##  Plot predicted line etc.
mat3b.summ[1:30,1:2]
Xnames




##  Calculate Predicted Lines
Xnames[1]  <-  "RateFast(Intercept)"
Xnames
Coeff3b  <-  mat3b.summ$Mean
Fast5_RegLineb   <-  inv_logit(Coeff3b[1] + Coeff3b[2] * nSperm_z)
Fast55_RegLineb  <-  inv_logit((Coeff3b[1] + Coeff3b[4]) + (Coeff3b[2] + Coeff3b[6]) * nSperm_z)
Slow5_RegLineb   <-  inv_logit((Coeff3b[1] + Coeff3b[3]) + (Coeff3b[2] + Coeff3b[5]) * nSperm_z)
Slow55_RegLineb  <-  inv_logit((Coeff3b[1] + Coeff3b[3] + Coeff3b[7]) + (Coeff3b[2] + Coeff3b[5] + Coeff3b[8]) * nSperm_z)
Fast_RegLineb    <-  inv_logit((Coeff3b[1] + (Coeff3b[4])/2) + (Coeff3b[2] + (Coeff3b[6])/2) * nSperm_z)
Slow_RegLineb    <-  inv_logit((Coeff3b[1] + Coeff3b[3] + (0.5*(Coeff3b[7]))) + (Coeff3b[2] + Coeff3b[5] + (0.5*(Coeff3b[8]))) * nSperm_z)




##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLineb[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLineb[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLineb[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLineb[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat3.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue1',0.01))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered1',0.01))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )







##############################################################
##############################################################
##  MODELING COVARIANCE STRUCTURE


#####################################
## Logistic Mixed Effects Regression
## MAXIMAL MODEL W/ COVARIANCE MATRIX
## -- random intercept for RUN
## -- random slopes for Run x nSperm
## -- Estimate covariance matrix
#####################################

head(data)

#  Fixed Effects Model Matrix
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

#  Random Effects Model Matrix
#Za       <-  model.matrix(~ -1 + data$Run*nSperm_z*data$Rate*data$EggPos)
#Zanames  <-  dimnames(Z)[[2]]

##  Not enough observations to include EggPos ixns! Results in only 
##  3 observations per run.
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z +
                                data$Run:data$Rate +
#                                data$Run:data$EggPos +
                                data$Run:nSperm_z:data$Rate) #+
#                                data$Run:nSperm_z:data$EggPos +
#                                data$Run:data$Rate:data$EggPos +
#                                data$Run:nSperm_z:data$Rate:data$EggPos)
Znames  <-  dimnames(Z)[[2]]
any(Znames != Zanames)
any(Z != Za)
dim(Z)

Z       <-  unname(Z)
attr(Z,"assign") <- NULL
str(Z)
head(Z)
Z[,11]#  <-  0
Z[1:20,]

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
mat4 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  1,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )


# Model Results
#print(mat4)
print(mat4, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat4, c("tau_run","corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat4.df    <-  as.data.frame(extract(mat4))
mcmc.mat4  <-  as.mcmc(mat4)
mat4.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat4)
mat4.summ  <-  plyr:::adply(as.matrix(mat4.df),2,MCMCsum)
#(mat4.summ)
mat4.summ$X1

# Explore Correlation structure
hist(mat4.summ[2050:3649,2][mat4.summ[2050:3649,2] < 0.2],breaks=100)
corrMat  <-  matrix(mat4.summ[2050:3649,2], ncol=40,nrow=40)
corrplot(corrMat , method='circle', type='upper')
abline(v=c(10.5,20.5,30.5))
abline(h=c(10.5,20.5,30.5))

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  }
}

Znames
corrplot(corrMat * 15, method='circle', type='upper')
abline(v=c(10.5,20.5,30.5))
abline(h=c(10.5,20.5,30.5))

# Simple Diagnostic Plots
plot(mat4, pars="beta")
pairs(mat4, pars="beta")


#  LOO Log-likelihood for model selection
mat4LL  <-  extract_log_lik(mat4, parameter_name = "log_lik")
mat4Loo    <-  loo(mat4LL)
mat4WAIC   <-  waic(mat4LL)








print(mat3, c("beta", "lp__"), probs=c( 0.5));
print(mat4, c("beta","sigma_y", "lp__"), probs=c(0.5));


##  Calculate Predicted Lines
Xnames
head(mat4.summ)
Coeff4  <-  mat4.summ$Mean[2041:2048]
Fast5_RegLine   <-  inv_logit(Coeff4[1] + Coeff4[2] * nSperm_z)
Fast55_RegLine  <-  inv_logit((Coeff4[1] + Coeff4[4]) + (Coeff4[2] + Coeff4[6]) * nSperm_z)
Slow5_RegLine   <-  inv_logit((Coeff4[1] + Coeff4[3]) + (Coeff4[2] + Coeff4[5]) * nSperm_z)
Slow55_RegLine  <-  inv_logit((Coeff4[1] + Coeff4[3] + Coeff4[7]) + (Coeff4[2] + Coeff4[5] + Coeff4[8]) * nSperm_z)
Fast_RegLine    <-  inv_logit((Coeff4[1] + (Coeff4[4])/2) + (Coeff4[2] + (Coeff4[6])/2) * nSperm_z)
Slow_RegLine    <-  inv_logit((Coeff4[1] + Coeff4[3] + (0.5*(Coeff4[7]))) + (Coeff4[2] + Coeff4[5] + (0.5*(Coeff4[8]))) * nSperm_z)




##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat4.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), 
            col=transparentColor('dodgerblue1',0.1))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), 
            col=transparentColor('orangered1',0.1))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )












#####################################
## Logistic Mixed Effects Regression
##  2ndary ixns MODEL W/ COVARIANCE MATRIX
## -- random intercept for RUN
## -- random slopes for Run x nSperm
## -- Estimate covariance matrix
#####################################

head(data)

#  Fixed Effects Model Matrix
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

#  Random Effects Model Matrix
#Za       <-  model.matrix(~ -1 + data$Run*nSperm_z*data$Rate*data$EggPos)
#Zanames  <-  dimnames(Z)[[2]]
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z +
                                data$Run:data$Rate) #+
#                                data$Run:data$EggPos +
#                                data$Run:nSperm_z:data$Rate) +
#                                data$Run:nSperm_z:data$EggPos +
#                                data$Run:data$Rate:data$EggPos +
#                                data$Run:nSperm_z:data$Rate:data$EggPos)
Znames  <-  dimnames(Z)[[2]]
any(Znames != Zanames)
any(Z != Za)
dim(Z)

Z       <-  unname(Z)
attr(Z,"assign") <- NULL
str(Z)
head(Z)
Z[,11]#  <-  0
Z[1:20,]

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

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 2000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

## Call to STAN
mat5 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  1,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )




# Model Results
#print(mat5)
print(mat5, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat5, c("tau_run","corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat5.df    <-  as.data.frame(extract(mat5))
mcmc.mat5  <-  as.mcmc(mat5)
mat5.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat5)
mat5.summ  <-  plyr:::adply(as.matrix(mat5.df),2,MCMCsum)
#(mat5.summ)
mat5.summ$X1

# Explore Correlation structure
hist(mat5.summ[1231:1238,2][mat5.summ[1231:1238,2] < 0.2],breaks=100)
corrMat  <-  matrix(mat5.summ[1240:2139,2], ncol=30,nrow=30)
corrplot(corrMat , method='circle', type='upper')
abline(v=c(10.5,20.5,30.5))
abline(h=c(10.5,20.5,30.5))

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  }
}

Znames
corrplot(corrMat * 15, method='circle', type='upper')
abline(v=c(10.5,20.5,30.5))
abline(h=c(10.5,20.5,30.5))

# Simple Diagnostic Plots
plot(mat5, pars="beta")
pairs(mat5, pars="beta")


#  LOO Log-likelihood for model selection
mat5LL  <-  extract_log_lik(mat5, parameter_name = "log_lik")
mat5Loo    <-  loo(mat5LL)
mat5WAIC   <-  waic(mat5LL)


##  Calculate Predicted Lines
mat5.summ$X1
Coeff5  <-  mat5.summ$Mean[1231:1238]
Fast5_RegLine   <-  inv_logit(Coeff5[1] + Coeff5[2] * nSperm_z)
Fast55_RegLine  <-  inv_logit((Coeff5[1] + Coeff5[4]) + (Coeff5[2] + Coeff5[6]) * nSperm_z)
Slow5_RegLine   <-  inv_logit((Coeff5[1] + Coeff5[3]) + (Coeff5[2] + Coeff5[5]) * nSperm_z)
Slow55_RegLine  <-  inv_logit((Coeff5[1] + Coeff5[3] + Coeff5[7]) + (Coeff5[2] + Coeff5[5] + Coeff5[8]) * nSperm_z)
Fast_RegLine    <-  inv_logit((Coeff5[1] + (Coeff5[4])/2) + (Coeff5[2] + (Coeff5[6])/2) * nSperm_z)
Slow_RegLine    <-  inv_logit((Coeff5[1] + Coeff5[3] + (0.5*(Coeff5[7]))) + (Coeff5[2] + Coeff5[5] + (0.5*(Coeff5[8]))) * nSperm_z)




##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat4.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), 
            col=transparentColor('dodgerblue1',0.1))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), 
            col=transparentColor('orangered1',0.1))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )














#####################################
## Logistic Mixed Effects Regression
## -- random intercept for RUN
## -- random slopes for Run x nSperm
## -- Estimate covariance matrix
#####################################

head(data)

#  Fixed Effects Model Matrix
X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
X       <-  unname(X)
attr(X,"assign") <- NULL
str(X)
head(X)

#  Random Effects Model Matrix
Z       <-  model.matrix(~ -1 + data$Run +
                                data$Run:nSperm_z) # +
#                                data$Run:data$Rate) +
#                                data$Run:data$EggPos +
#                                data$Run:nSperm_z:data$Rate) +
#                                data$Run:nSperm_z:data$EggPos +
#                                data$Run:data$Rate:data$EggPos +
#                                data$Run:nSperm_z:data$Rate:data$EggPos)
Znames  <-  dimnames(Z)[[2]]
any(Znames != Zanames)
any(Z != Za)
dim(Z)

Z       <-  unname(Z)
attr(Z,"assign") <- NULL
str(Z)
head(Z)
Z[,11]#  <-  0
Z[1:20,]

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

#  Options for the analysis
nChains        = 4
thinSteps      = 5
numSavedSteps  = 2000 #across all chains
burnInSteps    = numSavedSteps / 2
nIter          = ceiling(burnInSteps+(numSavedSteps * thinSteps)/nChains)

## Call to STAN
mat6 <- stan(data     =  data.list,
             file     =  './Stan/mat-logistic-1Z-cov.stan',
             chains   =  1,
             iter     =  numSavedSteps,
             thin     =  thinSteps,
             save_dso =  TRUE
            )




# Model Results
#print(mat6)
print(mat6, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(mat6, c("tau_run","corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
mat6.df    <-  as.data.frame(extract(mat6))
mcmc.mat6  <-  as.mcmc(mat6)
mat6.mcmc  <-  rstan:::as.mcmc.list.stanfit(mat6)
mat6.summ  <-  plyr:::adply(as.matrix(mat6.df),2,MCMCsum)
#(mat6.summ)
mat6.summ$X1

# Explore Correlation structure
hist(mat6.summ[630:1029,2][mat6.summ[630:1029,2] < 0.2],breaks=100)
corrMat  <-  matrix(mat6.summ[630:1029,2], ncol=20,nrow=20)
corrplot(corrMat , method='circle', type='upper')
abline(v=c(10.5,20.5))
abline(h=c(10.5,20.5))

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  }
}

Znames
corrplot(corrMat * 15, method='circle', type='upper')
abline(v=c(10.5,20.5))
abline(h=c(10.5,20.5))

# Simple Diagnostic Plots
plot(mat6, pars="beta")
pairs(mat6, pars="beta")


#  LOO Log-likelihood for model selection
mat6LL  <-  extract_log_lik(mat6, parameter_name = "log_lik")
mat6Loo    <-  loo(mat6LL)
mat6WAIC   <-  waic(mat6LL)


##  Calculate Predicted Lines
mat6.summ$X1
Coeff6  <-  mat6.summ$Mean[621:628]
Fast5_RegLine   <-  inv_logit(Coeff6[1] + Coeff6[2] * nSperm_z)
Fast55_RegLine  <-  inv_logit((Coeff6[1] + Coeff6[4]) + (Coeff6[2] + Coeff6[6]) * nSperm_z)
Slow5_RegLine   <-  inv_logit((Coeff6[1] + Coeff6[3]) + (Coeff6[2] + Coeff6[5]) * nSperm_z)
Slow55_RegLine  <-  inv_logit((Coeff6[1] + Coeff6[3] + Coeff6[7]) + (Coeff6[2] + Coeff6[5] + Coeff6[8]) * nSperm_z)
Fast_RegLine    <-  inv_logit((Coeff6[1] + (Coeff6[4])/2) + (Coeff6[2] + (Coeff6[6])/2) * nSperm_z)
Slow_RegLine    <-  inv_logit((Coeff6[1] + Coeff6[3] + (0.5*(Coeff6[7]))) + (Coeff6[2] + Coeff6[5] + (0.5*(Coeff6[8]))) * nSperm_z)




##  Plot of all 4 regression lines for Rate x EggPos
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(Fast5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(Fast55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(Slow5_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(Slow55_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.3,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(16,21,16,21),
          col     =  c('dodgerblue1','dodgerblue1','orangered1','orangered1'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )



##  Plot of Rate x nSperm effect.
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
 apply(mat4.df, 1, function(x, data, nSperm_z){
     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), 
            col=transparentColor('dodgerblue1',0.1))
     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), 
            col=transparentColor('orangered1',0.1))
 }, data=data, nSperm_z=nSperm_z)
lines(Fast_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(Slow_RegLine[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
    legend(
          x       =  usr[2]*0.2,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(Fast)),
                      expression(paste(Slow))),
          pch     =  c(21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('orangered1',0.7)),
          col     =  c('dodgerblue3','orangered3'),
          cex     =  1,
          xjust   =  1,
          bty     =  'n',
          border  =  NA
    )
























#######################
##  MODEL COMPARISONS
#######################

str(mat1Loo)
looDiff   <-  compare(mat1Loo, mat2Loo, mat3Loo, mat4Loo, mat5Loo, mat6Loo)
waicDiff  <-  compare(mat1WAIC, mat2WAIC, mat3WAIC, mat4WAIC, mat5WAIC, mat6WAIC)

print(looDiff, digits=4)
print(waicDiff, digits=4)

print(compare(mat1Loo, mat2Loo), digits=6)
print(compare(mat1Loo, mat3Loo), digits=6)
print(compare(mat1Loo, mat4Loo), digits=6)
print(compare(mat2Loo, mat3Loo), digits=6)
print(compare(mat2Loo, mat4Loo), digits=6)
print(compare(mat3Loo, mat4Loo), digits=6)


str(looDiff)
looDiff54  <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][2]
looDiff56  <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][3]
looDiff53  <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][4]
looDiff52  <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][5]
looDiff51  <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][6]
looDiff46  <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][3]
looDiff43  <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][4]
looDiff42  <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][5]
looDiff41  <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][6]
looDiff63  <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][4]
looDiff62  <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][5]
looDiff61  <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][6]
looDiff32  <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][5]
looDiff31  <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][6]
looDiff21  <-  looDiff[,'elpd_loo'][5] - looDiff[,'elpd_loo'][6]

n  <-  length(mat1Loo$pointwise[,"elpd_loo"])
selooDiff54  <-  sqrt(n * var(mat5Loo$pointwise[,"elpd_loo"]  - mat4Loo$pointwise[,"elpd_loo"]))
selooDiff56  <-  sqrt(n * var(mat5Loo$pointwise[,"elpd_loo"]  - mat6Loo$pointwise[,"elpd_loo"]))
selooDiff53  <-  sqrt(n * var(mat5Loo$pointwise[,"elpd_loo"]  - mat3Loo$pointwise[,"elpd_loo"]))
selooDiff52  <-  sqrt(n * var(mat5Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"]))
selooDiff51  <-  sqrt(n * var(mat5Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff46  <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat6Loo$pointwise[,"elpd_loo"]))
selooDiff43  <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat3Loo$pointwise[,"elpd_loo"]))
selooDiff42  <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"]))
selooDiff41  <-  sqrt(n * var(mat4Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff63  <-  sqrt(n * var(mat6Loo$pointwise[,"elpd_loo"]  - mat3Loo$pointwise[,"elpd_loo"]))
selooDiff62  <-  sqrt(n * var(mat6Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"]))
selooDiff61  <-  sqrt(n * var(mat6Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff32  <-  sqrt(n * var(mat3Loo$pointwise[,"elpd_loo"]  - mat2Loo$pointwise[,"elpd_loo"]))
selooDiff31  <-  sqrt(n * var(mat3Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))
selooDiff21  <-  sqrt(n * var(mat2Loo$pointwise[,"elpd_loo"]  - mat1Loo$pointwise[,"elpd_loo"]))


LooDiff  <-  cbind(c(looDiff54,looDiff56,looDiff53,looDiff52,looDiff51,looDiff46,looDiff43,looDiff42,looDiff41,looDiff63,looDiff62,looDiff61,looDiff32,looDiff31,looDiff21),
                   c(selooDiff54,selooDiff56,selooDiff53,selooDiff52,selooDiff51,selooDiff46,selooDiff43,selooDiff42,selooDiff41,selooDiff63,selooDiff62,selooDiff61,selooDiff32,selooDiff31,selooDiff21))

pDiff43   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[1,1] - 0)/LooDiff[1,2])), 3))
pDiff42   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[2,1] - 0)/LooDiff[2,2])), 3))
pDiff41   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[3,1] - 0)/LooDiff[3,2])), 3))
pDiff32   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[4,1] - 0)/LooDiff[4,2])), 3))
pDiff31   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[5,1] - 0)/LooDiff[5,2])), 3))
pDiff21   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[6,1] - 0)/LooDiff[6,2])), 3))

pDiff54  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[1,1] - 0)/LooDiff[1,2])), 3))
pDiff56  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[2,1] - 0)/LooDiff[2,2])), 3))
pDiff53  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[3,1] - 0)/LooDiff[3,2])), 3))
pDiff52  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[4,1] - 0)/LooDiff[4,2])), 3))
pDiff51  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[5,1] - 0)/LooDiff[5,2])), 3))
pDiff46  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[6,1] - 0)/LooDiff[6,2])), 3))
pDiff43  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[7,1] - 0)/LooDiff[7,2])), 3))
pDiff42  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[8,1] - 0)/LooDiff[8,2])), 3))
pDiff41  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[9,1] - 0)/LooDiff[9,2])), 3))
pDiff63  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[10,1] - 0)/LooDiff[10,2])), 3))
pDiff62  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[11,1] - 0)/LooDiff[11,2])), 3))
pDiff61  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[12,1] - 0)/LooDiff[12,2])), 3))
pDiff32  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[13,1] - 0)/LooDiff[13,2])), 3))
pDiff31  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[14,1] - 0)/LooDiff[14,2])), 3))
pDiff21  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[15,1] - 0)/LooDiff[15,2])), 3))



LooDiff  <-  cbind(LooDiff, c(pDiff54,pDiff56,pDiff53,pDiff52,pDiff51,pDiff46,pDiff43,pDiff42,pDiff41,pDiff63,pDiff62,pDiff61,pDiff32,pDiff31,pDiff21))


row.names(LooDiff)  <-  c('mat5 - mat4',
                          'mat5 - mat6',
                          'mat5 - mat3',
                          'mat5 - mat2',
                          'mat5 - mat1',
                          'mat4 - mat6',
                          'mat4 - mat3',
                          'mat4 - mat2',
                          'mat4 - mat1',
                          'mat6 - mat3',
                          'mat6 - mat2',
                          'mat6 - mat1',
                          'mat3 - mat2',
                          'mat3 - mat1',
                          'mat2 - mat1')
colnames(LooDiff)   <-  c("diff", "se", "p.value")
LooDiff


str(mat3Loo)

##  Plot differences 
m3m4   <-  density(mat3Loo$pointwise[,"elpd_loo"] - mat4Loo$pointwise[,"elpd_loo"])
m3m2   <-  density(mat3Loo$pointwise[,"elpd_loo"] - mat2Loo$pointwise[,"elpd_loo"])
m3m2b  <-  density(mat3Loo$pointwise[,"elpd_loo"] - mat2bLoo$pointwise[,"elpd_loo"])
m3m1   <-  density(mat3Loo$pointwise[,"elpd_loo"] - mat1Loo$pointwise[,"elpd_loo"])

allx   <-  c(m3m4$x,m3m2$x,m3m2b$x,m3m1$x)
ally   <-  c(m3m4$y,m3m2$y,m3m2b$y,m3m1$y)



##  Not really sure what's going on here... I don't think these 
##  Density plots should look the way they do... should be more 
##  spread out... right?
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









