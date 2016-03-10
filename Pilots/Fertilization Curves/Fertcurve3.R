data=read.csv("C:\\Users\\colito\\Documents\\Projects\\PhD\\Fertilization\\Pilots\\Fertilization Curves\\Fertcurve3.csv", header=TRUE)

ebar = function(a,b,c,d,e){
arrows(a,b,a,b+c,angle=90,code=2,length=e)
arrows(a,b,a,b-d,angle=90,code=2,length=e)
}

fertrate=data$fert/data$trials
polyrate=data$poly/data$trials
failed=(data$trials-data$fert)/data$trials


# PLOTS OF REPLICATE FERTILIZATION CURVES
plot(fertrate[data$series==1] ~ log10(data$conc[data$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(fertrate[data$series==2] ~ log10(data$conc[data$series==2]), type='l', lwd=2, col=2)
points(fertrate[data$series==3] ~ log10(data$conc[data$series==3]), type='l', lwd=2, col=3)
box()
axis(1, at=log10(data$conc[data$series==1]), label=format(data$conc[data$series==1], scientific=TRUE))
axis(2)

# PLOTS OF RATE OF POLYSPERMY
plot(polyrate[data$series==1] ~ log10(data$conc[data$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(polyrate[data$series==2] ~ log10(data$conc[data$series==2]), type='l', lwd=2, col=2)
points(polyrate[data$series==3] ~ log10(data$conc[data$series==3]), type='l', lwd=2, col=3)
box()
axis(1, at=log10(data$conc[data$series==1]), label=format(data$conc[data$series==1], scientific=TRUE))
axis(2)

# CALCULATING MEAN FERTILIZATION
conc=unique(data$conclev)
mus=tapply(fertrate,data$conclev,mean)
mu.poly=tapply(polyrate,data$conclev,mean)
se=tapply(fertrate,data$conclev,function(x){sd(x)/sqrt(length(x))})
se.poly=tapply(polyrate,data$conclev,function(x){sd(x)/sqrt(length(x))})


# PLOT OF THE MEAN OF ALL 3 FERTILIZATION CURVES
plot(mus~log10(data$conc[data$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Average Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(data$conc[data$series==1]), mus, se, se, 0.05)
box()
axis(1, at=log10(data$conc[data$series==1]), label=format(data$conc[data$series==1], scientific=TRUE))
axis(2)



# PLOT OF THE MEAN OF ALL 3 POLYSPERMY RATE CURVES
plot(mu.poly~log10(data$conc[data$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Average Polyspermy Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(data$conc[data$series==1]), mu.poly, se.poly, se.poly, 0.05)
box()
axis(1, at=log10(data$conc[data$series==1]), label=format(data$conc[data$series==1], scientific=TRUE))
axis(2)
