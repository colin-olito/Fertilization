data=read.csv("C:\\Users\\colito\\Documents\\Projects\\PhD\\Fertilization\\Pilots\\Fertilization Curves\\Fertage1_contam.csv", header=TRUE)

ebar = function(a,b,c,d,e){
arrows(a,b,a,b+c,angle=90,code=2,length=e)
arrows(a,b,a,b-d,angle=90,code=2,length=e)
}

fertrate=data$fert/data$trials
polyrate=data$poly/data$trials
failed=(data$trials-data$fert)/data$trials
data=cbind(data,fertrate,polyrate,failed)

new=subset(data,data$treat == "nFnM")
med=subset(data,data$treat == "nFmM")
old=subset(data,data$treat == "nFoM")

# PLOTS OF REPLICATE FERTILIZATION CURVES
par(mfrow=c(3,1))
plot(new$fertrate[new$series==1] ~ log10(new$conc[new$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(new$fertrate[new$series==2] ~ log10(new$conc[new$series==2]), type='l', lwd=2, col=2)
box()
axis(1, at=log10(new$conc[new$series==1]), label=format(new$conc[new$series==1], scientific=TRUE))
axis(2)

plot(med$fertrate[med$series==1] ~ log10(med$conc[new$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(med$fertrate[med$series==2] ~ log10(med$conc[new$series==2]), type='l', lwd=2, col=2)
box()
axis(1, at=log10(med$conc[med$series==1]), label=format(med$conc[med$series==1], scientific=TRUE))
axis(2)

plot(old$fertrate[new$series==1] ~ log10(old$conc[old$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), ylab="Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(old$fertrate[old$series==2] ~ log10(old$conc[old$series==2]), type='l', lwd=2, col=2)
box()
axis(1, at=log10(old$conc[old$series==1]), label=format(old$conc[old$series==1], scientific=TRUE))
axis(2)




# PLOTS OF POLYSPERMY RATES
plot(new$polyrate[new$series==1] ~ log10(new$conc[new$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="New Sperm", ylab="Polyspermy Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(new$polyrate[new$series==2] ~ log10(new$conc[new$series==2]), type='l', lwd=2, col=2)
box()
axis(1, at=log10(new$conc[new$series==1]), label=format(new$conc[new$series==1], scientific=TRUE))
axis(2)

plot(med$polyrate[med$series==1] ~ log10(med$conc[new$series==1]),  ylim=c(0,1), xlim=log10(c(263, 26300000)), main="Medium Sperm", ylab="PolyspermyRate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(med$polyrate[med$series==2] ~ log10(med$conc[new$series==2]), type='l', lwd=2, col=2)
box()
axis(1, at=log10(med$conc[med$series==1]), label=format(med$conc[med$series==1], scientific=TRUE))
axis(2)

plot(old$polyrate[new$series==1] ~ log10(old$conc[old$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main = "Old Sperm", ylab="PolyspermyRate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
points(old$polyrate[old$series==2] ~ log10(old$conc[old$series==2]), type='l', lwd=2, col=2)
box()
axis(1, at=log10(old$conc[old$series==1]), label=format(old$conc[old$series==1], scientific=TRUE))
axis(2)



# CALCULATING MEAN FERTILIZATION
conc=unique(new$conclev)
newmu=tapply(new$fertrate,new$conclev,mean)
newmu.poly=tapply(new$polyrate,new$conclev,mean)
newse=tapply(new$fertrate,new$conclev,function(x){sd(x)/sqrt(length(x))})
newse.poly=tapply(new$polyrate,new$conclev,function(x){sd(x)/sqrt(length(x))})

medmu=tapply(med$fertrate,med$conclev,mean)
medmu.poly=tapply(med$polyrate,med$conclev,mean)
medse=tapply(med$fertrate,med$conclev,function(x){sd(x)/sqrt(length(x))})
medse.poly=tapply(med$polyrate,med$conclev,function(x){sd(x)/sqrt(length(x))})


oldmu=tapply(old$fertrate,old$conclev,mean)
oldmu.poly=tapply(old$polyrate,old$conclev,mean)
oldse=tapply(old$fertrate,old$conclev,function(x){sd(x)/sqrt(length(x))})
oldse.poly=tapply(old$polyrate,old$conclev,function(x){sd(x)/sqrt(length(x))})



# PLOT OF THE MEAN OF ALL 3 FERTILIZATION CURVES
plot(newmu~log10(new$conc[new$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="New Sperm", ylab="Average Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(new$conc[new$series==1]), newmu, newse, newse, 0.05)
box()
axis(1, at=log10(new$conc[new$series==1]), label=format(new$conc[new$series==1], scientific=TRUE))
axis(2)

plot(medmu~log10(med$conc[med$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="Med Sperm", ylab="Average Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(med$conc[med$series==1]), medmu, medse, medse, 0.05)
box()
axis(1, at=log10(med$conc[med$series==1]), label=format(med$conc[med$series==1], scientific=TRUE))
axis(2)

plot(oldmu~log10(old$conc[old$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="Old Sperm", ylab="Average Fertilization Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(old$conc[old$series==1]), oldmu, oldse, oldse, 0.05)
box()
axis(1, at=log10(old$conc[old$series==1]), label=format(old$conc[old$series==1], scientific=TRUE))
axis(2)




# PLOT OF THE MEAN OF ALL 3 POLYSPERMY CURVES
plot(newmu.poly~log10(new$conc[new$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="New Sperm", ylab="Average Polyspermy Rate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(new$conc[new$series==1]), newmu.poly, newse.poly, newse.poly, 0.05)
box()
axis(1, at=log10(new$conc[new$series==1]), label=format(new$conc[new$series==1], scientific=TRUE))
axis(2)

plot(medmu.poly~log10(med$conc[med$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="Med Sperm", ylab="Average PolyspermyRate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(med$conc[med$series==1]), medmu.poly, medse.poly, medse.poly, 0.05)
box()
axis(1, at=log10(med$conc[med$series==1]), label=format(med$conc[med$series==1], scientific=TRUE))
axis(2)

plot(oldmu.poly~log10(old$conc[old$series==1]), ylim=c(0,1), xlim=log10(c(263, 26300000)), main="Old Sperm", ylab="Average PolyspermyRate", xlab="Sperm Concentration per mL", axes=FALSE, type='l',lwd=2)
ebar(log10(old$conc[old$series==1]), oldmu.poly, oldse.poly, oldse.poly, 0.05)
box()
axis(1, at=log10(old$conc[old$series==1]), label=format(old$conc[old$series==1], scientific=TRUE))
axis(2)

