######### Asignment 2, Ozone

rm(list=ls())
setwd("C:/Users/ander/Documents/10. Semester/Advanced data analysis/Assignment 2")
library("ggplot2")
library("dplyr")
library("GGally")
library("car")
library("MASS")
library("ellipse")
library(gclus)
library("LaplacesDemon")
data(ozone)
head(ozone)
#ggpairs(ozone)

# log-Likelihood function, probably NOT poisson
pois <- function(lamb){
  sum(dpois(ozone$Ozone,lambda=lamb,log=TRUE))}
optpois <- optimize(pois,c(0,20),maximum=TRUE)
optpois$maximum


poisd=dpois(seq(0,40,by=1),lambda=optpois$maximum)
gammad = dgamma(seq(0,40,by=1),shape=2,rate=1/6)
lnormd = dlnorm(seq(0,40,by=1),mean=mean(log(ozone$Ozone)),sd=sd(log(ozone$Ozone)))
invgaussd = dinvgaussian(seq(0.1,40),mu=mean(ozone$Ozone),lambda=1/(1/330*sum(1/ozone$Ozone-1/mean(ozone$Ozone))))


par(mfrow=c(1,1))
a=hist(ozone$Ozone,breaks=100,freq=FALSE, xlim=c(0,40),ylim = c(0, 0.3),main="Density of ozone")
lines(poisd/sum(poisd)*2,col='red',lwd=2)
lines(gammad/sum(gammad)*2,col='green',lwd=2)
lines(lnormd/sum(lnormd)*2,col='blue',lwd=2)
lines(invgaussd/sum(invgaussd)*2,col='yellow',lwd=2)
legend("topright", legend=c("Poisson", "Gamma","log-normal","Inverse Gaussian"),col=c("red", "green","blue","yellow"), lty=1:2, cex=0.8)


### Simple presentation of the data
par(mfrow=c(2,4))
plot(ozone$Ozone~ozone$Temp, xlab='Temperature', ylab='Ozone')
plot(ozone$Ozone~ozone$InvHt, xlab='Inversion base height', ylab='Ozone',type="p")
plot(ozone$Ozone~ozone$Pres, xlab='Pressure', ylab='Ozone')
plot(ozone$Ozone~ozone$Vis, xlab='Visibility', ylab='Ozone')
plot(ozone$Ozone~ozone$Hgt, xlab='Vandenburg Height', ylab='Ozone')
plot(ozone$Ozone~ozone$Hum, xlab='Humidity', ylab='Ozone')
plot(ozone$Ozone~ozone$InvTmp, xlab='Inversion temperature', ylab='Ozone')
plot(ozone$Ozone~ozone$Wind, xlab='Wind speed', ylab='Ozone')

#ggpairs(ozone)
## Fit a general linear model, without 2nd degree terms
formula1 = Ozone ~ Temp+InvHt+Pres+Vis+Hgt+Hum+InvTmp+Wind
fitlm1 <- lm(formula = formula1, data = ozone)
summary(fitlm1)
par(mfrow=c(2,2))
plot(fitlm1)

drop1(fitlm1,test="F")
fitlm2 <- update(fitlm1,.~. -Wind)
summary(fitlm2)
plot(fitlm2)

drop1(fitlm2,test="F")
fitlm3 <-update(fitlm2,.~. -Pres)
summary(fitlm3)
plot(fitlm3)

drop1(fitlm3,test="F")
fitlm4 <- update(fitlm3,.~. -Hgt)
summary(fitlm4)
plot(fitlm4)

drop1(fitlm4,test="F")
fitlm5 <- update(fitlm4,.~. -InvTmp)
summary(fitlm5)
plot(fitlm5)

drop1(fitlm5,test="F")
fitlm6 <- update(fitlm5,.~. -Vis)
summary(fitlm6)
plot(fitlm6)

(sqrt(sum((ozone$Ozone-fitted(fitlm6))^2)))



############## Do some residual plots
par(mfrow=c(3,2))
plot(residuals(fitlm6)) #Residuals versus obs number 
plot(fitted(fitlm6),residuals(fitlm6)) #residuals versus fitted values
plot(ozone$Ozone,residuals(fitlm6)) #residuals versus response, this should always be linear!
plot(ozone$Temp,residuals(fitlm6)) #residuals versus Temperature
plot(ozone$InvHt,residuals(fitlm6)) #residuals versus Inverse height
plot(ozone$Hum,residuals(fitlm6)) #residuals versus Humidity

par(mfrow=c(1,3))
plot(fitted(fitlm6),residuals(fitlm6),xlab="Fitted values",ylab="Residuals",main="Fitted vs residuals")
plot(fitted(fitlm6),ozone$Ozone,xlab="Fitted values",ylab="observations",main="fitted vs observations")
lines(0:30,0:30,col="red")
qqPlot(fitlm6,simulate=FALSE,xlab="t Quantiles",ylab="Studentized residuals",main="QQplot",lwd=1)



## Fit a general linear model, with 2nd degree terms
formula2 = Ozone ~ Temp+InvHt+Pres+Vis+Hgt+Hum+InvTmp+Wind+I(Temp^2)+I(InvHt^2)+I(Pres^2)+I(Vis^2)+I(Hgt^2)+I(Hum^2)+I(InvTmp^2)+I(Wind^2)
fitlm21 <- lm(formula = formula2, data = ozone)
summary(fitlm21)
par(mfrow=c(2,2))
plot(fitlm1)

drop1(fitlm21,test="F")
fitlm22 <- update(fitlm21,.~. -I(Hum^2))
summary(fitlm22)
plot(fitlm22)

drop1(fitlm22,test="F")
fitlm23 <- update(fitlm22,.~. -InvHt)
summary(fitlm23)
plot(fitlm23)

drop1(fitlm23,test="F")
fitlm24 <- update(fitlm23,.~. -I(Hgt^2))
summary(fitlm24)
plot(fitlm24)

drop1(fitlm24,test="F")
fitlm25 <- update(fitlm24,.~. -Hgt)
summary(fitlm25)
plot(fitlm25)

drop1(fitlm25,test="F")
fitlm26 <- update(fitlm25,.~. -InvTmp)
summary(fitlm26)
plot(fitlm26)

drop1(fitlm26,test="F")
fitlm27 <- update(fitlm26,.~. -Wind)
summary(fitlm27)
plot(fitlm27)

drop1(fitlm27,test="F")
fitlm28 <- update(fitlm27,.~. -I(Wind^2))
summary(fitlm28)
plot(fitlm28)

drop1(fitlm28,test="F")
fitlm29 <- update(fitlm28,.~. -I(Vis^2))
summary(fitlm29)
plot(fitlm29)

drop1(fitlm29,test="F")
fitlm210 <- update(fitlm29,.~. -I(InvHt^2))
summary(fitlm210)
plot(fitlm210)

par(mfrow=c(3,3))
plot(residuals(fitlm210)) #Residuals versus obs number 
plot(fitted(fitlm210),residuals(fitlm210)) #residuals versus fitted values
plot(ozone$Ozone,residuals(fitlm210)) #residuals versus response, this should always be linear!
plot(ozone$Temp,residuals(fitlm210)) #residuals versus Temperature
plot(ozone$Pres,residuals(fitlm210)) #residuals versus Pressure
plot(ozone$Vis,residuals(fitlm210)) #residuals versus Visibility
plot(ozone$Hum,residuals(fitlm210)) #residuals versus Humidity

par(mfrow=c(1,3))
plot(fitted(fitlm210),residuals(fitlm210),xlab="Fitted values",ylab="Residuals",main="Fitted vs residuals")
plot(fitted(fitlm210),ozone$Ozone,xlab="Fitted values",ylab="observations",main="fitted vs observations")
lines(0:30,0:30,col="red")
qqPlot(fitlm210,simulate=FALSE,xlab="t Quantiles",ylab="Studentized residuals",main="QQplot",lwd=1)

(rmse = sqrt(sum(residuals(fitlm210)^2)))


## Fit a general linear model, with 2nd degree terms but with log
formula3 = log(Ozone) ~ Temp+InvHt+Pres+Vis+Hgt+Hum+InvTmp+Wind+I(Temp^2)+I(InvHt^2)+I(Pres^2)+I(Vis^2)+I(Hgt^2)+I(Hum^2)+I(InvTmp^2)+I(Wind^2)
fitlm31 <- lm(formula3, data = ozone)

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -I(Hgt^2))

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -I(InvTmp^2))

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -I(Temp^2))

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -Hgt)

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -InvTmp)

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -I(Wind^2))

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -Wind)

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -InvHt)

drop1(fitlm31,test="F")
fitlm31 <- update(fitlm31,.~. -I(Hum^2))

drop1(fitlm31,test="F")
fitlm36 <- update(fitlm31,.~. -Hum)
summary(fitlm36)

par(mfrow=c(1,3))
plot(exp(fitted(fitlm36)),residuals(fitlm36),xlab="Fitted values",ylab="Residuals",main="Fitted vs residuals")
plot(exp(fitted(fitlm36)),ozone$Ozone,xlab="Fitted values",ylab="observations",main="fitted vs observations")
lines(0:30,0:30,col="red")
qqPlot(fitlm36,simulate=FALSE,xlab="t Quantiles",ylab="Studentized residuals",main="QQplot",lwd=1)


par(mfrow=c(2,3))
plot(residuals(fitlm36)) #Residuals versus obs number
plot(fitted(fitlm36),residuals(fitlm36)) #residuals versus fitted values
plot(log(ozone$Ozone),residuals(fitlm36)) #residuals versus response, this should always be linear!
plot(ozone$Temp,residuals(fitlm36)) #residuals versus Temperature
plot(ozone$InvHt,residuals(fitlm36)) #residuals versus Inverse Height
plot(ozone$Vis,residuals(fitlm36)) #residuals versus Humidity
par(mfrow=c(1,1))
hist(residuals(fitlm36))

par(mfrow=c(1,3))
plot(fitted(fitlm36),residuals(fitlm36))
plot(fitted(fitlm36),rstandard(fitlm36))
plot(fitted(fitlm36),rstudent(fitlm36))

par(mfrow=c(1,1))
qqPlot(fitlm36,simulate=FALSE)

(sqrt(sum((ozone$Ozone-exp(fitted(fitlm36)))^2)))



########### GLM part

fitIG1= glm(formula = formula2,family=inverse.gaussian(link="log"),data=ozone)
fitIG2= glm(formula = formula2,family=inverse.gaussian(link="identity"),data=ozone)
fitIG3= glm(formula = formula2,family=inverse.gaussian(link="inverse"),data=ozone)
fitG1 = glm(formula = formula2,family=Gamma(link="log"),data=ozone)
fitG2 = glm(formula = formula2,family=Gamma(link="identity"),data=ozone)
fitG3 = glm(formula = formula2,family=Gamma(link="inverse"),data=ozone)
sum(residuals(fitIG1,type="deviance"))
sum(residuals(fitIG2,type="deviance"))
sum(residuals(fitIG3,type="deviance"))
sum(residuals(fitG1,type="deviance"))
sum(residuals(fitG2,type="deviance"))
sum(residuals(fitG3,type="deviance"))
sum(residuals(fitIG1,type="pearson"))
sum(residuals(fitIG2,type="pearson"))
sum(residuals(fitIG3,type="pearson"))
sum(residuals(fitG1,type="pearson"))
sum(residuals(fitG2,type="pearson"))
sum(residuals(fitG3,type="pearson"))
sum((ozone$Ozone-fitted(fitIG1))^2)
sum((ozone$Ozone-fitted(fitIG2))^2)
sum((ozone$Ozone-fitted(fitIG3))^2)
sum((ozone$Ozone-fitted(fitG1))^2)
sum((ozone$Ozone-fitted(fitG2))^2)
sum((ozone$Ozone-fitted(fitG3))^2)
AIC(fitIG1)
AIC(fitIG2)
AIC(fitIG3)
AIC(fitG1)
AIC(fitG2)
AIC(fitG3)
(sqrt(sum((ozone$Ozone-fitted(fitIG1))^2)))
(sqrt(sum((ozone$Ozone-fitted(fitIG2))^2)))
(sqrt(sum((ozone$Ozone-fitted(fitIG3))^2)))
(sqrt(sum((ozone$Ozone-fitted(fitG1))^2)))
(sqrt(sum((ozone$Ozone-fitted(fitG2))^2)))
(sqrt(sum((ozone$Ozone-fitted(fitG3))^2)))


fitG = glm(formula = formula2,family=Gamma(link="inverse"),data=ozone)
drop1(fitG,test="F")
fitG <- update(fitG,.~. -I(Vis^2))
drop1(fitG,test="F")
fitG <- update(fitG,.~. -I(Hgt^2))
drop1(fitG,test="F")
fitG <- update(fitG,.~. -I(InvTmp^2))
drop1(fitG,test="F")
fitG <- update(fitG,.~. -I(Wind^2))
drop1(fitG,test="F")
fitG <- update(fitG,.~. -Hgt)
drop1(fitG,test="F")
fitG <- update(fitG,.~. -Wind)
summary(fitG)
(sqrt(sum((ozone$Ozone-fitted(fitG))^2)))

par(mfrow=c(1,3))
plot(fitted(fitG),residuals(fitG),xlab="Fitted values",ylab="Residuals",main="Fitted vs residuals")
plot(fitted(fitG),ozone$Ozone,xlab="Fitted values",ylab="observations",main="fitted vs observations")
lines(0:30,0:30,col="red")
qqPlot(fitG,simulate=FALSE,xlab="t Quantiles",ylab="Studentized residuals",main="QQplot",lwd=1)


############# 2.6
mu <- predict(fitG,type="response")
Vmu <- I(mu^2)
w2 <- 1/summary(fitG)$dispersion
w1 = 1
(W1 <- diag(w1*Vmu))
(W2 <- diag(w2*Vmu))

(Sigma1 = inv(t(model.matrix(fitG))%*%W1%*%model.matrix(fitG)))
(Sigma2 = inv(t(model.matrix(fitG))%*%W2%*%model.matrix(fitG)))

S2 = summary(fitG)$cov.scaled
S1 = summary(fitG)$cov.unscaled


formatC(S2, format = "e", digits = 2)
formatC(Sigma2, format = "e", digits = 2) #These are the same
formatC(S1, format = "e", digits = 2)
formatC(Sigma1, format = "e", digits = 2) #These are the same


############## 2.7+8
formula4 = Ozone ~ Temp+InvHt+Pres+Vis+Hgt+Hum+InvTmp+Wind+
  I(Temp^2)+I(InvHt^2)+I(Pres^2)+I(Vis^2)+I(Hgt^2)+I(Hum^2)+I(InvTmp^2)+I(Wind^2)+
  I(Temp^3)+I(InvHt^3)+I(Pres^3)+I(Vis^3)+I(Hgt^3)+I(Hum^3)+I(InvTmp^3)+I(Wind^3)+
  Temp:InvHt+Temp:Pres+Temp:Vis+Temp:Hgt+Temp:Hum+Temp:InvTmp+Temp:Wind+
  InvHt:Pres+InvHt:Vis+InvHt:Hgt+InvHt:Hum+InvHt:InvTmp+InvHt:Wind+
  Pres:Vis+Pres:Hgt+Pres:Hum+Pres:InvTmp+Pres:Wind+
  Vis:Hgt+Vis:Hum+Vis:InvTmp+Vis:Wind+
  Hgt:Hum+Hgt:InvTmp+Hgt:Wind+
  Hum:InvTmp+Hum:Wind+InvTmp:Wind

fitG1 = glm(formula = formula4,family=Gamma(link="inverse"),data=ozone)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -InvTmp:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Pres:Hum)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Hgt:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -InvHt:InvTmp)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:InvTmp)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Vis:Hum)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Hgt:InvTmp)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:Hgt)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -InvHt:Hgt)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:Pres)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Vis:Hgt)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:Hum)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -InvHt:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Pres:Vis)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Vis:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -InvHt:Pres)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -InvHt:Vis)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:Vis)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Vis:InvTmp)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Hum:InvTmp)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Pres:InvTmp)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Temp:InvHt)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Pres:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Hum:Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Hgt:Hum)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(InvTmp^3))#######Nåede hertil
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Temp^3))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Wind^3))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(InvHt^3))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Pres^3))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Hgt^3))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Hgt^2))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Wind)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Wind^2))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Vis)
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Vis^3))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -I(Vis^2))
drop1(fitG1,test="F")
fitG1 <- update(fitG1,.~. -Hgt)
summary(fitG1)
par(mfrow=c(2,2))
plot(fitG1)

(sqrt(sum((ozone$Ozone-fitted(fitG1))^2)))

x=1:13
c = coef(fitG1)
c = c[2:14]
x2 = names(c)
conf = confint(fitG1)
conf = conf[2:14,]
par(mfrow=c(1,1))
plot(x,c)
arrows(x, conf[,1], x, conf[,2], length=0.05, angle=90, code=3)
axis(1, at=1:13, labels=x2)

library(coefplot)

coefplot(fitG1, horizontal = TRUE, innerCI = 0, pointSize =2,numberAngle = -90,
         coefficients=c("Temp","InvHt","Pres","Hum","InvTmp","(Temp^2)","I(InvHt^2)","I(Pres^2)","I(Hum^2)","I(InvTmp^2)","I(Hum^3)","InvHt:Hum","Pres:Hgt"))

ress = (ozone$Ozone-fitted(fitG1))
par(mfrow=c(1,2))
plot(fitted(fitG1),ress,xlab="Fitted values",ylab="Residuals",main="Fitted vs residuals")
plot(fitted(fitG1),ozone$Ozone,xlab="Fitted values",ylab="observations",main="fitted vs observations")
lines(0:30,0:30,col="red")
