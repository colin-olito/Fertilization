#/* 
# * Colin Olito. Created 26/04/2016.
# * 
# * NOTES: Final Analyses and model selection 
# *         for the 2nd Flume Experiment
# *         crossing N x Rate; with 2 egg patches
# *          
# *          
# */

## Logistic Mixed Effects Regression Analysis
rm(list=ls())

###################
##  GLOBAL OPTIONS
options("menu.graphics"=FALSE)

#################
##  DEPENDENCIES
source('R/dependencies.R')


####################
##  Fit All Models
####################
source('./NxRate_FitAllModels.R')


######################
##  Check Convergence 
##  for Fixed Effects
######################
print(m1, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m2, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m4, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m4a, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m5, c("beta", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));


hist(m1.summ[2050:3649,2][m1.summ[2050:3649,2] < 0.2],breaks=100)
corrMat  <-  matrix(m1.summ[2050:3649,2], ncol=40,nrow=40)
corrplot(corrMat , method='ellipse', type='upper')
abline(v=c(10.5,20.5,30.5))
abline(h=c(10.5,20.5,30.5))

for (i in 1:nrow(corrMat)) {
  for (j in 1:ncol(corrMat)) {
    if(i == j)
      corrMat[i,j] = 0
  }
}

Znames
corrplot(corrMat * 15, method='ellipse', type='upper')
abline(v=c(10.5,20.5,30.5))
abline(h=c(10.5,20.5,30.5))

#######################
##  MODEL COMPARISONS
#######################

##  Overall comparison of all models
looDiff   <-  compare(m1Loo, m2Loo, m3Loo, m3aLoo, m4Loo, m4aLoo, m5Loo)
waicDiff  <-  compare(m1WAIC, m2WAIC, m3WAIC, m3aWAIC, m4WAIC, m4aWAIC, m5WAIC)

print(looDiff, digits=4)
print(waicDiff, digits=4)

str(looDiff)
looDiff21    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][2]
looDiff23a   <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][3]
looDiff23    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][4]
looDiff24    <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][5]
looDiff24a   <-  looDiff[,'elpd_loo'][1] - looDiff[,'elpd_loo'][6]
looDiff25    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][7]
looDiff13a   <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][3]
looDiff13    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][4]
looDiff14    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][5]
looDiff14a   <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][6]
looDiff15    <-  looDiff[,'elpd_loo'][2] - looDiff[,'elpd_loo'][7]
looDiff3a3   <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][4]
looDiff3a4   <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][5]
looDiff3a4a  <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][6]
looDiff3a5   <-  looDiff[,'elpd_loo'][3] - looDiff[,'elpd_loo'][7]
looDiff34    <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][5]
looDiff34a   <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][6]
looDiff35    <-  looDiff[,'elpd_loo'][4] - looDiff[,'elpd_loo'][7]
looDiff44a   <-  looDiff[,'elpd_loo'][5] - looDiff[,'elpd_loo'][6]
looDiff45    <-  looDiff[,'elpd_loo'][5] - looDiff[,'elpd_loo'][7]
looDiff4a5   <-  looDiff[,'elpd_loo'][6] - looDiff[,'elpd_loo'][7]


n  <-  length(m1Loo$pointwise[,"elpd_loo"])
selooDiff21    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m1Loo$pointwise[,"elpd_loo"]))
selooDiff23a   <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m3aLoo$pointwise[,"elpd_loo"]))  
selooDiff23    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m3Loo$pointwise[,"elpd_loo"]))  
selooDiff24    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff24a   <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff25    <-  sqrt(n * var(m2Loo$pointwise[,"elpd_loo"]  - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff13a   <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"]  - m3aLoo$pointwise[,"elpd_loo"]))  
selooDiff13    <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"]  - m3Loo$pointwise[,"elpd_loo"]))  
selooDiff14    <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"]  - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff14a   <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"]  - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff15    <-  sqrt(n * var(m1Loo$pointwise[,"elpd_loo"]  - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff3a3   <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m3Loo$pointwise[,"elpd_loo"]))  
selooDiff3a4   <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff3a4a  <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff3a5   <-  sqrt(n * var(m3aLoo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff34    <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff34a   <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m4aLoo$pointwise[,"elpd_loo"]))  
selooDiff35    <-  sqrt(n * var(m3Loo$pointwise[,"elpd_loo"]  - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff44a   <-  sqrt(n * var(m4Loo$pointwise[,"elpd_loo"]  - m4Loo$pointwise[,"elpd_loo"]))  
selooDiff45    <-  sqrt(n * var(m4Loo$pointwise[,"elpd_loo"]  - m5Loo$pointwise[,"elpd_loo"]))  
selooDiff4a5   <-  sqrt(n * var(m4aLoo$pointwise[,"elpd_loo"] - m5Loo$pointwise[,"elpd_loo"]))  


LooDiff  <-  cbind(c(looDiff21,looDiff23a,looDiff23,looDiff24,looDiff24a,looDiff25,looDiff13a,looDiff13,looDiff14,looDiff14a,looDiff15,looDiff3a3,looDiff3a4,looDiff3a4a,looDiff3a5,looDiff34,looDiff34a,looDiff35,looDiff44a,looDiff45,looDiff4a5),
                   c(selooDiff21,selooDiff23a,selooDiff23,selooDiff24,selooDiff24a,selooDiff25,selooDiff13a,selooDiff13,selooDiff14,selooDiff14a,selooDiff15,selooDiff3a3,selooDiff3a4,selooDiff3a4a,selooDiff3a5,selooDiff34,selooDiff34a,selooDiff35,selooDiff44a,selooDiff45,selooDiff4a5))


pDiff21    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[1,1] - 0)/LooDiff[1,2])), 3))
pDiff23a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[2,1] - 0)/LooDiff[2,2])), 3))
pDiff23    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[3,1] - 0)/LooDiff[3,2])), 3))
pDiff24    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[4,1] - 0)/LooDiff[4,2])), 3))
pDiff24a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[5,1] - 0)/LooDiff[5,2])), 3))
pDiff25    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[6,1] - 0)/LooDiff[6,2])), 3))
pDiff13a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[7,1] - 0)/LooDiff[7,2])), 3))
pDiff13    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[8,1] - 0)/LooDiff[8,2])), 3))
pDiff14    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[9,1] - 0)/LooDiff[9,2])), 3))
pDiff14a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[10,1] - 0)/LooDiff[10,2])), 3))
pDiff15    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[11,1] - 0)/LooDiff[11,2])), 3))
pDiff3a3   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[12,1] - 0)/LooDiff[12,2])), 3))
pDiff3a4   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[13,1] - 0)/LooDiff[13,2])), 3))
pDiff3a4a  <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[14,1] - 0)/LooDiff[14,2])), 3))
pDiff3a5   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[15,1] - 0)/LooDiff[15,2])), 3))
pDiff34    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[16,1] - 0)/LooDiff[16,2])), 3))
pDiff34a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[17,1] - 0)/LooDiff[17,2])), 3))
pDiff35    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[18,1] - 0)/LooDiff[18,2])), 3))
pDiff44a   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[19,1] - 0)/LooDiff[19,2])), 3))
pDiff45    <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[20,1] - 0)/LooDiff[20,2])), 3))
pDiff4a5   <-  as.numeric(rounded(2*pnorm(-abs((LooDiff[21,1] - 0)/LooDiff[21,2])), 3))

LooDiff  <-  cbind(LooDiff, c(pDiff21,pDiff23a,pDiff23,pDiff24,pDiff24a,pDiff25,pDiff13a,pDiff13,pDiff14,pDiff14a,pDiff15,pDiff3a3,pDiff3a4,pDiff3a4a,pDiff3a5,pDiff34,pDiff34a,pDiff35,pDiff44a,pDiff45,pDiff4a5))
row.names(LooDiff)  <-  c('m2 - m1',
						  'm2 - m3a',
						  'm2 - m3',
						  'm2 - m4',
						  'm2 - m4a',
						  'm2 - m5',
						  'm1 - m3a',
						  'm1 - m3',
						  'm1 - m4',
						  'm1 - m4a',
						  'm1 - m5',
						  'm3a - m3',
						  'm3a - m4',
						  'm3a - m4a',
						  'm3a - m5',
						  'm3 - m4',
						  'm3 - m4a',
						  'm3 - m5',
						  'm4 - m4a',
						  'm4 - m5',
						  'm4a - m5')
colnames(LooDiff)   <-  c("diff", "se", "p.value")
LooDiff


#########################################
##  m3 & m3a are best compromise between 
##  model fit and number of terms
##
##  Explore both, see whether including 
##  covariance structure helps or not.
#########################################

print(m3, c("beta","sigma_y", "lp__"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3a, c("beta"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3a, c("gamma0", "gamma1", "sigma_gamma0", "sigma_gamma1"), 
	  probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

##  Correlation structure for m3
hist(m3.summ[630:1029,2][m3.summ[630:1029,2] < 1],breaks=100)
corrMat  <-  matrix(m3.summ[630:1029,2], ncol=20,nrow=20)
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



##  All of the estimated covariances lie between 
##  -0.015 and 0.015... with standard deviations
##  in the neighborhood of 0.21... providing strong
##  evidence that these correlations could be 0
print(m3, c("corrs"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95),digits=3);

##  But the standard deviations of unconditional 
##  random effect distributions appear to be 
##  different from 0
print(m3, c("tau_run"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
plot(m3, pars="tau_run")



##  Have a look at Fixed Effect Coefficients
plot(m3, pars="beta")
plot(m3a, pars="beta")
pairs(m3, pars="beta")
pairs(m3a, pars="beta")


X       <-  model.matrix(~ 1 + nSperm_z*Rate*EggPos, data=data)
Xnames  <-  dimnames(X)[[2]]
Xnames
m3.summ

##  Calculate Predicted Lines
m3.coef   <-  m3.summ$Mean[621:628]
m3a.coef  <-  m3a.summ$Mean[1:8]

m3Fast5   <-  inv_logit(m3.coef[1] + m3.coef[2] * nSperm_z)
m3Fast55  <-  inv_logit((m3.coef[1] + m3.coef[4]) + (m3.coef[2] + m3.coef[6]) * nSperm_z)
m3Slow5   <-  inv_logit((m3.coef[1] + m3.coef[3]) + (m3.coef[2] + m3.coef[5]) * nSperm_z)
m3Slow55  <-  inv_logit((m3.coef[1] + m3.coef[3] + m3.coef[7]) + (m3.coef[2] + m3.coef[5] + m3.coef[8]) * nSperm_z)
m3Fast    <-  inv_logit((m3.coef[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * nSperm_z)
m3Slow    <-  inv_logit((m3.coef[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * nSperm_z)

m3aFast5   <-  inv_logit(m3a.coef[1] + m3a.coef[2] * nSperm_z)
m3aFast55  <-  inv_logit((m3a.coef[1] + m3a.coef[4]) + (m3a.coef[2] + m3a.coef[6]) * nSperm_z)
m3aSlow5   <-  inv_logit((m3a.coef[1] + m3a.coef[3]) + (m3a.coef[2] + m3a.coef[5]) * nSperm_z)
m3aSlow55  <-  inv_logit((m3a.coef[1] + m3a.coef[3] + m3a.coef[7]) + (m3a.coef[2] + m3a.coef[5] + m3a.coef[8]) * nSperm_z)
m3aFast    <-  inv_logit((m3a.coef[1] + (m3a.coef[4])/2) + (m3a.coef[2] + (m3a.coef[6])/2) * nSperm_z)
m3aSlow    <-  inv_logit((m3a.coef[1] + m3a.coef[3] + (0.5*(m3a.coef[7]))) + (m3a.coef[2] + m3a.coef[5] + (0.5*(m3a.coef[8]))) * nSperm_z)

m3.low   <-  m3.summ$lower[621:628]
m3.hi   <-  m3.summ$upper[621:628]
m3Fast.low    <-  inv_logit((m3.low[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * nSperm_z)
m3Slow.low    <-  inv_logit((m3.low[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * nSperm_z)
m3Fast.hi     <-  inv_logit((m3.hi[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * nSperm_z)
m3Slow.hi     <-  inv_logit((m3.hi[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * nSperm_z)


##  MODEL 3
##  Plot of all 4 regression lines for Rate x EggPos

pdf(file='./output/NxRatexEggPos_prelim_m3.pdf', height=7, width=7)
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m3Fast5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3Fast55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3Slow5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3Slow55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
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
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
dev.off()
embed_fonts("./output/NxRatexEggPos_prelim_m3.pdf", outfile="./output/NxRatexEggPos_prelim_m3_embedded.pdf")


##  MODEL 3a
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
lines(m3aFast5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3aFast55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3aSlow5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3aSlow55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
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
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )


##  look at residuals for m3
m3yhat  <-  inv_logit(m3.summ$Mean[1030:1149])
m3.resids  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3yhat)/sd(m3yhat)
m3.resids_z  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3yhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m3yhat))

##  Model 3 Residual Plots
par(mfrow=c(2,2))
hist(m3.resids, breaks=40)
plot(m3.resids ~ nSperm_z)
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m3.resids ~ seq_along(m3.resids))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m3.resids)
qqline(m3.resids, col = 2)

hist(m3.resids, breaks=40)
plot(m3.resids_z ~ nSperm_z)
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m3.resids_z ~ seq_along(m3.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m3.resids_z)
qqline(m3.resids_z, col = 2)


#################################
##  Need yhat values for Model 3a

# ##  look at residuals for m3a
m3ayhat  <-  inv_logit(m3a.summ$Mean[31:150])
m3a.resids  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3ayhat)/sqrt(m3ayhat)
m3a.resids_z  <-  (((data$nFert - data$nControlFert)/data$nEggs) - m3ayhat)/sd((((data$nFert - data$nControlFert)/data$nEggs) - m3ayhat))

##  Model 3 Residual Plots
par(mfrow=c(2,2))
hist(m3a.resids, breaks=40)
plot(m3a.resids ~ nSperm_z)
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m3a.resids ~ seq_along(m3a.resids))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m3a.resids)
qqline(m3a.resids, col = 2)

hist(m3a.resids, breaks=40)
plot(m3a.resids_z ~ nSperm_z)
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
plot(m3a.resids_z ~ seq_along(m3a.resids_z))
abline(h=c(-2,0,2), lwd=c(1,3,1), col=c(1,2,1), lty=c(2,1,2))
qqnorm(m3a.resids_z)
qqline(m3a.resids_z, col = 2)



###################################################################
###################################################################
##  Plot of Rate x nSperm effect.
pdf(file="./output/NxRate_prelim_m3.pdf", height=7, width=7)
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m3.df, 1, function(x, data, nSperm_z){
#     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), 
#     	   col=transparentColor('dodgerblue1',0.05))
#     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), 
#     	   col=transparentColor('orangered1',0.05))
# }, data=data, nSperm_z=nSperm_z)
lines(m3Fast[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(m3Slow[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
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
dev.off()


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
# apply(m3a.df, 1, function(x, data, nSperm_z){
#     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue1',0.01))
#     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered1',0.01))
# }, data=data, nSperm_z=nSperm_z)
lines(m3aFast[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(m3aSlow[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
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



###################################################################
###################################################################
##  Plot of PER-GAMETE Rate x nSperm effect.

newy  <-  ((data$nFert - data$nControlFert)/data$nEggs)/data$nSperm

pdf(file="./output/NxRatexEggPos_perCapita.pdf", height=7, width=7)
par(omi=c(0.3, 0.5, 0.3, 0.3))
plot(newy ~ nSperm_z, 
    xlab='Sperm released', ylab='',
    type='n', axes=FALSE, ylim=c(min(newy)*1.1,max(newy)*1.2), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot regression lines
lines(m3Fast5[order(nSperm_z)]/data$nSperm[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3Fast55[order(nSperm_z)]/data$nSperm[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3Slow5[order(nSperm_z)]/data$nSperm[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3Slow55[order(nSperm_z)]/data$nSperm[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points(newy[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(newy[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points(newy[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points(newy[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.2),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
proportionalLabel(-0.175,0.5,'per-Sperm Fertilization rate', xpd=NA, srt=90, adj=0.5)
    legend(
          x       =  usr[2]*1,
          y       =  usr[4],
          legend  =  c(
                      expression(paste(5~cm:~Fast)),
                      expression(paste(55~cm:~Fast)),
                      expression(paste(5~cm:~Slow)),
                      expression(paste(55~cm:~Slow))),
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
 dev.off()


##  Plot of Rate x nSperm effect.
pdf(file="./output/NxRate_perCapita.pdf", height=7, width=7)
par(omi=c(0.3,0.4,0.3,0.3))
plot(newy ~ nSperm_z, 
    xlab='Sperm released', ylab='', 
    type='n', axes=FALSE, ylim=c(min(newy)*1.1,max(newy)*1.1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
# plot all regression lines from MCMC chains
# apply(m3a.df, 1, function(x, data, nSperm_z){
#     xrange   <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
#     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
#     lines(xrange, inv_logit((x['beta.1'] + (x['beta.4'])/2) + (x['beta.2'] + (x['beta.6'])/2) * xrange2), col=transparentColor('dodgerblue1',0.01))
#     lines(xrange, inv_logit((x['beta.1'] + x['beta.3'] + (0.5*(x['beta.7']))) + (x['beta.2'] + x['beta.5'] + (0.5*(x['beta.8']))) * xrange2), col=transparentColor('orangered1',0.01))
# }, data=data, nSperm_z=nSperm_z)
lines(m3aFast[order(nSperm_z)]/data$nSperm[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(m3aSlow[order(nSperm_z)]/data$nSperm[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
points((newy)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((newy)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
axis(2, las=1)
axis(1)
proportionalLabel(-0.175,0.5,'per-Sperm Fertilization rate', xpd=NA, srt=90, adj=0.5)
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
dev.off()

##############################################
##  CONTRASTS OF INTEREST
##
##  1: Fast vs. Slow across eggpos
##  2: 5 vs. 55 for each level of Rate
##  3: Fast vs. slow for each level of eggPos
##############################################

m3.betas   <-  m3.df[,621:628]
m3a.betas  <-  m3a.df[,1:8]
Xnames

b0Fast    <-  inv_logit((m3.betas[,1] + (m3.betas[,4])/2))
b0Slow    <-  inv_logit((m3.betas[,1] + m3.betas[,3] + (0.5*(m3.betas[,7]))))
b0Fast5   <-  inv_logit(m3.betas[,1])
b0Fast55  <-  inv_logit((m3.betas[,1] + m3.betas[,4]))
b0Slow5   <-  inv_logit((m3.betas[,1] + m3.betas[,3]))
b0Slow55  <-  inv_logit((m3.betas[,1] + m3.betas[,3] + m3.betas[,7]))
b1Fast    <-  inv_logit((m3.betas[,2] + (m3.betas[,6])/2))
b1Slow    <-  inv_logit((m3.betas[,2] + m3.betas[,5] + (0.5*(m3.betas[,8]))))
b1Fast5   <-  inv_logit(m3.betas[,2])
b1Fast55  <-  inv_logit((m3.betas[,2] + m3.betas[,6]))
b1Slow5   <-  inv_logit((m3.betas[,2] + m3.betas[,5]))
b1Slow55  <-  inv_logit((m3.betas[,2] + m3.betas[,5] + m3.betas[,8]))

b0Fast    <-  inv_logit((m3a.betas[,1] + (m3a.betas[,4])/2))
b0Slow    <-  inv_logit((m3a.betas[,1] + m3a.betas[,3] + (0.5*(m3a.betas[,7]))))
b0Fast5   <-  inv_logit(m3a.betas[,1])
b0Fast55  <-  inv_logit((m3a.betas[,1] + m3a.betas[,4]))
b0Slow5   <-  inv_logit((m3a.betas[,1] + m3a.betas[,3]))
b0Slow55  <-  inv_logit((m3a.betas[,1] + m3a.betas[,3] + m3a.betas[,7]))
b1Fast    <-  inv_logit((m3a.betas[,2] + (m3a.betas[,6])/2))
b1Slow    <-  inv_logit((m3a.betas[,2] + m3a.betas[,5] + (0.5*(m3a.betas[,8]))))
b1Fast5   <-  inv_logit(m3a.betas[,2])
b1Fast55  <-  inv_logit((m3a.betas[,2] + m3a.betas[,6]))
b1Slow5   <-  inv_logit((m3a.betas[,2] + m3a.betas[,5]))
b1Slow55  <-  inv_logit((m3a.betas[,2] + m3a.betas[,5] + m3a.betas[,8]))

pval  <-  function(x) length(x[x < 0])/length(x)
plotContr  <-  function(Dens, name="title") {
	plot(NA, xlab=expression(paste(Delta)), type='n', axes=FALSE, ylab='Density', cex.lab=1.2, xlim=c(min(Dens$x), (max(Dens$x)+0.4*(max(Dens$x) - min(Dens$x)))), ylim=c(0, (max(Dens$y)+0.05*(max(Dens$y) - min(Dens$y)))), yaxs='i')
	usr  <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
	box()
	polygon(c(Dens$x), c(Dens$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
	abline(v=0, lwd=2,col=2)
	axis(1, cex.axis=0.9)
	axis(2, cex.axis=0.9, las=1)
}


c1   <-  b0Slow - b0Fast
c2   <-  b1Slow - b1Fast
c3   <-  b1Slow - 0
c3b  <-  b1Fast
c4   <-  b0Fast5 - b0Fast55
c5   <-  b1Fast5 - b1Fast55
c6   <-  b0Slow5 - b0Slow55
c7   <-  b1Slow5 - b1Slow55
c8   <-  b0Fast5 - b0Slow5
c9   <-  b1Fast5 - b1Slow5
c10  <-  b0Fast55 - b0Slow55
c11  <-  b1Fast55 - b1Slow55


pval(c1)
pval(c2)
pval(c3b)
pval(c4)
pval(c5)
pval(c6)
pval(c7)
pval(c8)
pval(c9)
pval(c10)
pval(c11)

plotContr(density(c3))



#################################################
##  It may also be of INTEREST to make simple 
##  comparisons of the various regression lines 
##  at different values of nSperm
#################################################


m3Fast5.neg2   <-  inv_logit(m3.betas[1] + m3.betas[2] * (-2))
m3Fast55.neg2  <-  inv_logit((m3.betas[1] + m3.betas[4]) + (m3.betas[2] + m3.betas[6]) * (-2))
m3Slow5.neg2   <-  inv_logit((m3.betas[1] + m3.betas[3]) + (m3.betas[2] + m3.betas[5]) * (-2))
m3Slow55.neg2  <-  inv_logit((m3.betas[1] + m3.betas[3] + m3.betas[7]) + (m3.betas[2] + m3.betas[5] + m3.betas[8]) * (-2))
m3Fast.neg2    <-  inv_logit((m3.betas[1] + (m3.betas[4])/2) + (m3.betas[2] + (m3.betas[6])/2) * (-2))
m3Slow.neg2    <-  inv_logit((m3.betas[1] + m3.betas[3] + (0.5*(m3.betas[7]))) + (m3.betas[2] + m3.betas[5] + (0.5*(m3.betas[8]))) * (-2))

m3Fast5.neg1   <-  inv_logit(m3.betas[1] + m3.betas[2] * (-1))
m3Fast55.neg1  <-  inv_logit((m3.betas[1] + m3.betas[4]) + (m3.betas[2] + m3.betas[6]) * (-1))
m3Slow5.neg1   <-  inv_logit((m3.betas[1] + m3.betas[3]) + (m3.betas[2] + m3.betas[5]) * (-1))
m3Slow55.neg1  <-  inv_logit((m3.betas[1] + m3.betas[3] + m3.betas[7]) + (m3.betas[2] + m3.betas[5] + m3.betas[8]) * (-1))
m3Fast.neg1    <-  inv_logit((m3.betas[1] + (m3.betas[4])/2) + (m3.betas[2] + (m3.betas[6])/2) * (-1))
m3Slow.neg1    <-  inv_logit((m3.betas[1] + m3.betas[3] + (0.5*(m3.betas[7]))) + (m3.betas[2] + m3.betas[5] + (0.5*(m3.betas[8]))) * (-1))

m3Fast5.0   <-  inv_logit(m3.betas[1] + m3.betas[2] * (0))
m3Fast55.0  <-  inv_logit((m3.betas[1] + m3.betas[4]) + (m3.betas[2] + m3.betas[6]) * (0))
m3Slow5.0   <-  inv_logit((m3.betas[1] + m3.betas[3]) + (m3.betas[2] + m3.betas[5]) * (0))
m3Slow55.0  <-  inv_logit((m3.betas[1] + m3.betas[3] + m3.betas[7]) + (m3.betas[2] + m3.betas[5] + m3.betas[8]) * (0))
m3Fast.0    <-  inv_logit((m3.betas[1] + (m3.betas[4])/2) + (m3.betas[2] + (m3.betas[6])/2) * (0))
m3Slow.0    <-  inv_logit((m3.betas[1] + m3.betas[3] + (0.5*(m3.betas[7]))) + (m3.betas[2] + m3.betas[5] + (0.5*(m3.betas[8]))) * (0))

m3Fast5.1   <-  inv_logit(m3.betas[1] + m3.betas[2] * (1))
m3Fast55.1  <-  inv_logit((m3.betas[1] + m3.betas[4]) + (m3.betas[2] + m3.betas[6]) * (1))
m3Slow5.1   <-  inv_logit((m3.betas[1] + m3.betas[3]) + (m3.betas[2] + m3.betas[5]) * (1))
m3Slow55.1  <-  inv_logit((m3.betas[1] + m3.betas[3] + m3.betas[7]) + (m3.betas[2] + m3.betas[5] + m3.betas[8]) * (1))
m3Fast.1    <-  inv_logit((m3.betas[1] + (m3.betas[4])/2) + (m3.betas[2] + (m3.betas[6])/2) * (1))
m3Slow.1    <-  inv_logit((m3.betas[1] + m3.betas[3] + (0.5*(m3.betas[7]))) + (m3.betas[2] + m3.betas[5] + (0.5*(m3.betas[8]))) * (1))

m3Fast5.2   <-  inv_logit(m3.betas[1] + m3.betas[2] * (2))
m3Fast55.2  <-  inv_logit((m3.betas[1] + m3.betas[4]) + (m3.betas[2] + m3.betas[6]) * (2))
m3Slow5.2   <-  inv_logit((m3.betas[1] + m3.betas[3]) + (m3.betas[2] + m3.betas[5]) * (2))
m3Slow55.2  <-  inv_logit((m3.betas[1] + m3.betas[3] + m3.betas[7]) + (m3.betas[2] + m3.betas[5] + m3.betas[8]) * (2))
m3Fast.2    <-  inv_logit((m3.betas[1] + (m3.betas[4])/2) + (m3.betas[2] + (m3.betas[6])/2) * (2))
m3Slow.2    <-  inv_logit((m3.betas[1] + m3.betas[3] + (0.5*(m3.betas[7]))) + (m3.betas[2] + m3.betas[5] + (0.5*(m3.betas[8]))) * (2))


simpContr  <-  list(
	cSimp1   =  m3Fast5.neg2[,1] - m3Fast55.neg2[,1],
	cSimp2   =  m3Slow5.neg2[,1] - m3Slow55.neg2[,1],
	cSimp3   =  m3Fast5.neg2[,1] - m3Slow5.neg2[,1],
	cSimp4   =  m3Fast55.neg2[,1] - m3Slow55.neg2[,1],
	cSimp5   =  m3Fast.neg2[,1] - m3Slow.neg2[,1],
	cSimp6   =  m3Fast5.neg1[,1] - m3Fast55.neg1[,1],
	cSimp7   =  m3Slow5.neg1[,1] - m3Slow55.neg1[,1],
	cSimp8   =  m3Fast5.neg1[,1] - m3Slow5.neg1[,1],
	cSimp9   =  m3Fast55.neg1[,1] - m3Slow55.neg1[,1],
	cSimp10  =  m3Fast.neg1[,1] - m3Slow.neg1[,1],
	cSimp11  =  m3Fast5.0[,1] - m3Fast55.0[,1],
	cSimp12  =  m3Slow5.0[,1] - m3Slow55.0[,1],
	cSimp13  =  m3Fast5.0[,1] - m3Slow5.0[,1],
	cSimp14  =  m3Fast55.0[,1] - m3Slow55.0[,1],
	cSimp15  =  m3Fast.0[,1] - m3Slow.0[,1],
	cSimp16  =  m3Fast5.1[,1] - m3Fast55.1[,1],
	cSimp17  =  m3Slow5.1[,1] - m3Slow55.1[,1],
	cSimp18  =  m3Fast5.1[,1] - m3Slow5.1[,1],
	cSimp19  =  m3Fast55.1[,1] - m3Slow55.1[,1],
	cSimp20  =  m3Fast.1[,1] - m3Slow.1[,1],
	cSimp21  =  m3Fast5.2[,1] - m3Fast55.2[,1],
	cSimp22  =  m3Slow5.2[,1] - m3Slow55.2[,1],
	cSimp23  =  m3Fast5.2[,1] - m3Slow5.2[,1],
	cSimp24  =  m3Fast55.2[,1] - m3Slow55.2[,1],
	cSimp25  =  m3Fast.2[,1] - m3Slow.2[,1]
)



##  Histograms of simple comparisons between regression lines accross the 
##  nSperm gradient


pdf(file="./output/contrast_histograms.pdf", height=18, width=20)
par(mfrow=c(5,5), omi=rep(0.4,4))
plotContr(density(simpContr[[1]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Fast.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
proportionalLabel(-0.25, 0.5, expression(paste('x = -2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[2]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow.5 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[3]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Slow.5')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[4]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.55 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
plotContr(density(simpContr[[5]]))
proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast - Fast')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
for (i in 6:25) {
	plotContr(density(simpContr[[i]]))
	if(i == 6)
		proportionalLabel(-0.25, 0.5, expression(paste('x = -1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
	if(i == 11)
		proportionalLabel(-0.25, 0.5, expression(paste('x = 0',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
	if(i == 16)
		proportionalLabel(-0.25, 0.5, expression(paste('x = 1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
	if(i == 21)
		proportionalLabel(-0.25, 0.5, expression(paste('x = 2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
}
dev.off()

#########################################
##  Adjusting for Run-specific variation
#########################################
m3.summ$X1
m3a.summ$X1

##  Calculate Predicted Lines
m3Fast5   <-  inv_logit(m3.coef[1] + m3.coef[2] * nSperm_z)
m3Fast55  <-  inv_logit((m3.coef[1] + m3.coef[4]) + (m3.coef[2] + m3.coef[6]) * nSperm_z)
m3Slow5   <-  inv_logit((m3.coef[1] + m3.coef[3]) + (m3.coef[2] + m3.coef[5]) * nSperm_z)
m3Slow55  <-  inv_logit((m3.coef[1] + m3.coef[3] + m3.coef[7]) + (m3.coef[2] + m3.coef[5] + m3.coef[8]) * nSperm_z)
m3Fast    <-  inv_logit((m3.coef[1] + (m3.coef[4])/2) + (m3.coef[2] + (m3.coef[6])/2) * nSperm_z)
m3Slow    <-  inv_logit((m3.coef[1] + m3.coef[3] + (0.5*(m3.coef[7]))) + (m3.coef[2] + m3.coef[5] + (0.5*(m3.coef[8]))) * nSperm_z)

m3aFast5   <-  inv_logit(m3a.coef[1] + m3a.coef[2] * nSperm_z)
m3aFast55  <-  inv_logit((m3a.coef[1] + m3a.coef[4]) + (m3a.coef[2] + m3a.coef[6]) * nSperm_z)
m3aSlow5   <-  inv_logit((m3a.coef[1] + m3a.coef[3]) + (m3a.coef[2] + m3a.coef[5]) * nSperm_z)
m3aSlow55  <-  inv_logit((m3a.coef[1] + m3a.coef[3] + m3a.coef[7]) + (m3a.coef[2] + m3a.coef[5] + m3a.coef[8]) * nSperm_z)
m3aFast    <-  inv_logit((m3a.coef[1] + (m3a.coef[4])/2) + (m3a.coef[2] + (m3a.coef[6])/2) * nSperm_z)
m3aSlow    <-  inv_logit((m3a.coef[1] + m3a.coef[3] + (0.5*(m3a.coef[7]))) + (m3a.coef[2] + m3a.coef[5] + (0.5*(m3a.coef[8]))) * nSperm_z)


##  Calculate overall intercept for Model 3
coefs <- as.matrix(m3.mcmc)[,621:628]
dim(coefs)
head(coefs[,-c(2,5,6,8)])
Xmatr <- as.matrix(ddply(data.frame(X), ~interaction(data$nSperm), colwise(mean))[, -1])
Xmatr[1,-c(2,5,6,8)]
Xmatr  <-  matrix(Xmatr[1,-c(2,5,6,8)], nrow=nrow(coefs), ncol=4, byrow=TRUE)
head(Xmatr)
Xnames
int  <-  rowSums(coefs[,-c(2,5,6,8)] * Xmatr)
hist(inv_logit(int))
mean(int)

slope  <-  rowSums(coefs[,c(2,5,6,8)] * Xmatr)
hist(inv_logit(slope))
mean(slope)


##  Calculate overall intercept for Model 3a
coefs.a <- as.matrix(m3a.mcmc)[,1:8]
dim(coefs.a)
head(coefs.a[,-c(2,5,6,8)])
Xmatr <- as.matrix(ddply(data.frame(X), ~interaction(data$nSperm), colwise(mean))[, -1])
Xmatr[1,-c(2,5,6,8)]
Xmatr  <-  matrix(Xmatr[1,-c(2,5,6,8)], nrow=nrow(coefs.a), ncol=4, byrow=TRUE)
head(Xmatr)
Xnames
int  <-  rowSums(coefs.a[,-c(2,5,6,8)] * Xmatr)
hist(inv_logit(int))
mean(int)

slope  <-  rowSums(coefs[,c(2,5,6,8)] * Xmatr)
hist(inv_logit(slope))
mean(slope)


Xmatr <- as.matrix(ddply(data.frame(X), ~interaction(as.factor(data$Rate),as.factor(data$EggPos)), colwise(mean))[,-1])
cellmeans.mcmc <- coefs %*% t(Xmatr)
Xnames
colnames(cellmeans.mcmc) <- sort(unique(interaction(data$Rate, data$EggPos))) 
#colnames(cellmeans.mcmc) <- c("F.N","M.N","F.Y","M.Y")
head(cellmeans.mcmc)
(ixn.means <- adply(cellmeans.mcmc,2,MCMCsum))
aggregate(data$delta.vo2max,list(data$sex,data$spawn), mean)





m3.grand  <-  inv_logit(
						( m3.coef[1] + (0.5*m3.coef[3]) + (0.5*m3.coef[4]) + (0.25*m3.coef[7])) + 
						((m3.coef[2] + (0.5*m3.coef[5]) + (0.5*m3.coef[6]) + (0.25*m3.coef[8])) * nSperm_z))
m3.grand  <-  inv_logit(mean(int) + mean(slope) * nSperm_z)

m3a.grand  <-  inv_logit(
						( m3a.coef[1] + (0.5*m3a.coef[3]) + (0.5*m3a.coef[4]) + (0.25*m3a.coef[7])) + 
						((m3a.coef[2] + (0.5*m3a.coef[5]) + (0.5*m3a.coef[6]) + (0.25*m3a.coef[8])) * nSperm_z))
m3a.grand  <-  inv_logit(mean(int) + mean(slope) * nSperm_z)


##############################
##  Add Run-specific effects
##############################
print(m3a, c("gamma0", "gamma1", "sigma_gamma0", "sigma_gamma1"), 
	  probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
print(m3, c("gamma[1,1]","gamma[2,2]","gamma[3,3]","gamma[4,4]","gamma[5,5]","gamma[6,6]","gamma[7,7]","gamma[8,8]","gamma[9,9]","gamma[10,10]","gamma[1,11]","gamma[2,12]","gamma[3,13]","gamma[4,14]","gamma[5,15]","gamma[6,16]","gamma[7,17]","gamma[8,18]","gamma[9,19]","gamma[10,20]"), probs=c(0.05, 0.25, 0.5, 0.75, 0.95));

m3.ranef    <-  m3.summ$Mean[c(1550,1561,1572,1583,1594,1605,1616,1627,1638,1649,
								1650,1661,1672,1683,1694,1705,1716,1727,1738,1749)]
m3a.ranef   <-  m3a.summ$Mean[9:28]

Runs  <-  list(
"m3.r1"  =  inv_logit((mean(int) + m3.ranef[1]) + (mean(slope) + m3.ranef[11]) * nSperm_z),
"m3.r2"  =  inv_logit((mean(int) + m3.ranef[2]) + (mean(slope) + m3.ranef[12]) * nSperm_z),
"m3.r3"  =  inv_logit((mean(int) + m3.ranef[3]) + (mean(slope) + m3.ranef[13]) * nSperm_z),
"m3.r4"  =  inv_logit((mean(int) + m3.ranef[4]) + (mean(slope) + m3.ranef[14]) * nSperm_z),
"m3.r5"  =  inv_logit((mean(int) + m3.ranef[5]) + (mean(slope) + m3.ranef[15]) * nSperm_z),
"m3.r6"  =  inv_logit((mean(int) + m3.ranef[6]) + (mean(slope) + m3.ranef[16]) * nSperm_z),
"m3.r7"  =  inv_logit((mean(int) + m3.ranef[7]) + (mean(slope) + m3.ranef[17]) * nSperm_z),
"m3.r8"  =  inv_logit((mean(int) + m3.ranef[8]) + (mean(slope) + m3.ranef[18]) * nSperm_z),
"m3.r9"  =  inv_logit((mean(int) + m3.ranef[9]) + (mean(slope) + m3.ranef[19]) * nSperm_z),
"m3.r10" =  inv_logit((mean(int) + m3.ranef[10]) + (mean(slope) + m3.ranef[20]) * nSperm_z)
)

Runs.a  <-  list(
"m3a.r1"  =  inv_logit((mean(int) + m3a.ranef[1]) + (mean(slope) + m3a.ranef[11]) * nSperm_z),
"m3a.r2"  =  inv_logit((mean(int) + m3a.ranef[2]) + (mean(slope) + m3a.ranef[12]) * nSperm_z),
"m3a.r3"  =  inv_logit((mean(int) + m3a.ranef[3]) + (mean(slope) + m3a.ranef[13]) * nSperm_z),
"m3a.r4"  =  inv_logit((mean(int) + m3a.ranef[4]) + (mean(slope) + m3a.ranef[14]) * nSperm_z),
"m3a.r5"  =  inv_logit((mean(int) + m3a.ranef[5]) + (mean(slope) + m3a.ranef[15]) * nSperm_z),
"m3a.r6"  =  inv_logit((mean(int) + m3a.ranef[6]) + (mean(slope) + m3a.ranef[16]) * nSperm_z),
"m3a.r7"  =  inv_logit((mean(int) + m3a.ranef[7]) + (mean(slope) + m3a.ranef[17]) * nSperm_z),
"m3a.r8"  =  inv_logit((mean(int) + m3a.ranef[8]) + (mean(slope) + m3a.ranef[18]) * nSperm_z),
"m3a.r9"  =  inv_logit((mean(int) + m3a.ranef[9]) + (mean(slope) + m3a.ranef[19]) * nSperm_z),
"m3a.r10" =  inv_logit((mean(int) + m3a.ranef[10]) + (mean(slope) + m3a.ranef[20]) * nSperm_z)
)


##  Calculate adjusted y-values
yhat_adj  <-  c()
for(i in 1:nrow(data)) {
	yhat_adj[i]  <-  Runs[[data$Run[i]]][i]
}

y_adj  <-  ((data$nFert - data$nControlFert)/data$nEggs) - yhat_adj

yhat_adjA  <-  c()
for(i in 1:nrow(data)) {
	yhat_adjA[i]  <-  Runs.a[[data$Run[i]]][i]
}

y_adjA  <-  ((data$nFert - data$nControlFert)/data$nEggs) - yhat_adjA



#####################################################
##  Run-adjusted plots for m3

# plot of fertilization rate ~ sperm, grouped by run
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm, data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
for (i in 1:max(as.numeric(data$Run))){
points(((data$nFert[data$Run==i] - data$nControlFert[data$Run==i])/data$nEggs[data$Run==i]) ~ data$nSperm[data$Run==i], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.7), cex=1.1)
lines(Runs[[i]][data$Run==i][order(nSperm_z[data$Run==i])] ~ data$nSperm[data$Run==i][order(nSperm_z[data$Run==i])],
                  col='grey60', lwd=3)
}
points(yhat_adj[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)])
axis(1)
axis(2, las=1)



##  Adjusting y-values to account for Run effect.
##  Plotting nSperm x Rate
par(omi=rep(0.3, 4))
plot(y_adj ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
lines(m3Fast[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(m3Slow[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
# points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
#         bg=transparentColor('dodgerblue1', 0.7),
#         col=transparentColor('dodgerblue4', 0.9), cex=1.1)
# points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
#         bg=transparentColor('orangered1', 0.7),
#         col=transparentColor('orangered4', 0.9), cex=1.1)
points((m3Fast + y_adj)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((m3Slow + y_adj)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
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




##  Adjusting y-values to account for Run effect.
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
lines(m3Fast5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3Fast55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3Slow5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3Slow55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points((m3Fast5 + y_adj)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((m3Fast55 + y_adj)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((m3Slow5 + y_adj)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points((m3Slow55 + y_adj)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
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
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )





#####################################################
##  Run-adjusted plots for m3a

# plot of fertilization rate ~ sperm, grouped by run
par(omi=rep(0.3, 4))
plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm, data=data, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
for (i in 1:max(as.numeric(data$Run))){
points(((data$nFert[data$Run==i] - data$nControlFert[data$Run==i])/data$nEggs[data$Run==i]) ~ data$nSperm[data$Run==i], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.7), cex=1.1)
lines(Runs.a[[i]][data$Run==i][order(nSperm_z[data$Run==i])] ~ data$nSperm[data$Run==i][order(nSperm_z[data$Run==i])],
                  col='grey60', lwd=3)
}
points(yhat_adj[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)])
axis(1)
axis(2, las=1)



##  Adjusting y-values to account for Run effect.
##  Plotting nSperm x Rate
par(omi=rep(0.3, 4))
plot(y_adjA ~ nSperm_z, 
    xlab='Sperm released', ylab=substitute('Fertilization rate'), 
    type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
whiteGrid()
box()
lines(m3aFast[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue3', lwd=3)
lines(m3aSlow[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered3', lwd=3)
# points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
#         bg=transparentColor('dodgerblue1', 0.7),
#         col=transparentColor('dodgerblue4', 0.9), cex=1.1)
# points(((data$nFert-data$nControlFert)/data$nEggs)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
#         bg=transparentColor('orangered1', 0.7),
#         col=transparentColor('orangered4', 0.9), cex=1.1)
points((m3aFast + y_adjA)[data$Rate == "Fast"] ~ data$nSperm[data$Rate == "Fast"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((m3aSlow + y_adjA)[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
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




##  Adjusting y-values to account for Run effect.
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
lines(m3aFast5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lwd=3)
lines(m3aFast55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='dodgerblue1', lty=2, lwd=3)
lines(m3aSlow5[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lwd=3)
lines(m3aSlow55[order(nSperm_z)] ~ data$nSperm[order(nSperm_z)],
                  col='orangered1', lty=2, lwd=3)
points((m3aFast5 + y_adjA)[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.7),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((m3aFast55 + y_adjA)[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('dodgerblue1', 0.2),
        col=transparentColor('dodgerblue4', 0.9), cex=1.1)
points((m3aSlow5 + y_adjA)[data$Rate == "Slow" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
        bg=transparentColor('orangered1', 0.7),
        col=transparentColor('orangered4', 0.9), cex=1.1)
points((m3aSlow55 + y_adjA)[data$Rate == "Slow" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Slow" & data$EggPos == "5"], pch=21, 
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
          pch     =  c(21,21,21,21),
          pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7),transparentColor('orangered1',0.2)),
          col     =  c('dodgerblue4','dodgerblue4','orangered4','orangered4'),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
