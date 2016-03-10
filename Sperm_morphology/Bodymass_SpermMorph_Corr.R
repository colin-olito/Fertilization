#*********************************************
#*********************************************
# Colin Olito. Created 21/10/2015.
#
# BODY MASS X SPERM MORPHOLOGY CORRELATION
#
# NOTES: 
#          
#
#*********************************************
#*********************************************

#******************
# GLOBAL OPTIONS
options("menu.graphics"=FALSE)

#******************
# DEPENDENCIES
#install.packages("pbkrtest", loc="C:/R_home/R-3.1.2/library")
#install.packages("AICcmodavg", loc="C:/R_home/R-3.1.2/library")
#install.packages("LaplacesDemon", loc="C:/R_home/R-3.1.2/library")

library(car, lib.loc="/usr/lib/R/library")
library(ggplot2, lib.loc="/usr/lib/R/library")
library(gridExtra, lib.loc="/usr/lib/R/library")
library(contrast, lib.loc="/usr/lib/R/library")
library(lubridate, lib.loc="/usr/lib/R/library")
library(MASS, lib.loc="/usr/lib/R/library")
library(plyr, lib.loc="/usr/lib/R/library")
#library(nlme, lib.loc="/usr/lib/R/library")
library(lme4, lib.loc="/usr/lib/R/library") # remember to detatch("package:nlme") because of conflicts
library(boot, lib.loc="/usr/lib/R/library")
library(multcomp, lib.loc="/usr/lib/R/library")
library(MuMIn, lib.loc="/usr/lib/R/library")
library(AICcmodavg, lib.loc="/usr/lib/R/library"
library(rstan, lib.loc="/usr/lib/R/library")
library(MCMCpack, lib.loc="/usr/lib/R/library")
library(MCMCglmm, lib.loc="/usr/lib/R/library")
library(lubridate, lib.loc="/usr/lib/R/library")
library(pbkrtest, lib.loc="/usr/lib/R/library")


#*******************
# Import Data
data <- read.csv("GaCa_BodyMass_SpermMorph_Corr.csv", header=TRUE, strip.white=TRUE)
data <- data.frame(data)
head(data)
data$Block <- factor(data$Block)
data$ind <- factor(data$ind)
data$Date <- as.Date(data$Date, format='%d/%m/%Y')

# Make Date useable as factors vs. continuous covariates
data$Date = factor(data$Date)

str(data)







#************************************
# A FEW EXPLORATORY PLOTS
#************************************

HL.plot <- ggplot(data=data, aes(x=WBM, y=HL, colour=ind)) +
              geom_point() +
              labs(title=expression(paste("Body Mass ~ Head Length"))) +
              scale_y_continuous(expression(paste("Head Length (", mu,"m)"))) +
              scale_x_log10(expression(paste("Body Mass (",g,")", sep=""))) +
              geom_smooth(method='lm') +
              theme_bw() 
HL.plot





HW.plot <- ggplot(data=data, aes(x=WBM, y=HW, colour=Block)) +
              geom_point() +
              labs(title=expression(paste("Body Mass ~ Head Length"))) +
              scale_y_continuous(expression(paste("Head Length (", mu,"m)"))) +
              scale_x_log10(expression(paste("Body Mass (",g,")", sep=""))) +
              geom_smooth(method='lm') +
              theme_bw() 
HW.plot





HL_HW.plot <- ggplot(data=data, aes(x=HL, y=HW, colour=ind)) +
              geom_point() +
              labs(title=expression(paste("Head Width ~ Head Length"))) +
              scale_y_continuous(expression(paste("Head Width (", mu,"m)"))) +
              scale_x_log10(expression(paste("Head Length (", mu,"m)", sep=""))) +
              geom_smooth(method='lm', stderr=FALSE) +
              theme_bw() 
HL_HW.plot





####
# EVERYTHING AFTER THIS IS CODE FROM MY HELIOCIDARIS ANALYSES....
# CAN BE IGNORED, AND EVENTUALLY EXAPTED FOR USE IN THIS ANALYSIS





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
