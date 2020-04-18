### Clear variables and load data 
# clothingSum3 is sorted female first and male second
rm(list=ls())
setwd("C:/Users/ander/Documents/10. Semester/Advanced data analysis/Assignment 1")
CS <- read.csv(file = "clothingSum3.csv") #Sorted according to sex

library("ggplot2")
library("dplyr")
library("GGally")
library("car")
library("MASS")
library("ellipse")

## Remove subjId and day, and do exploratory analysis
CS1 <- subset(CS, select = -c(subjId,day) )
ggpairs(CS1)
corrplot(cor(CS1))

### Do backwards model selection
fitA1 <- lm(clo ~ tOut*tInOp * sex, data = CS1)
summary(fitA1)

drop1(fitA1,test="F")
fitA1 <- update(fitA1,.~. -tOut:tInOp:sex)
summary(fitA1)

drop1(fitA1,test="F")
fitA1 <- update(fitA1,.~. -tOut:tInOp)
summary(fitA1)

drop1(fitA1,test="F")
fitA1 <- update(fitA1,.~. -tOut:sex)
summary(fitA1)
####### Done with backwards model selection. 

############## Do some residual plots
par(mfrow=c(3,2))
plot(residuals(fitA1),col=CS1$sex) #Residuals versus obs number (not subj id)
plot(fitted(fitA1),residuals(fitA1),col=CS1$sex) #residuals versus fitted values
plot(CS1$clo,residuals(fitA1),col=CS1$sex) #residuals versus response, this will always be linear!
plot(CS1$tOut,residuals(fitA1),col=CS1$sex) #residuals versus tOut
plot(CS1$tInOp,residuals(fitA1),col=CS1$sex) #residuals versus tInOp
plot(CS1$sex,residuals(fitA1)) #residuals versus sex
##### The residuals plotted against sex suggests a small difference in variance. 
##### The residuals plotted against obs.nr. suggests the same (need color coding)
##### The boxplot in the initial analysis suggested the same. 
# Jan: find variance of residuals for women and men
var_women = var(residuals(fitA1)[1:70])
var_men = var(residuals(fitA1)[71:136])
var_ratio = var_women/var_men


### plot for report

par(mfrow=c(1,2))
plot(residuals(fitA1),,col=CS1$sex) #Residuals versus obs number (not subj id)
plot(CS1$sex,residuals(fitA1)) #residuals versus sex


####### Check for outliers using qqplot
par(mfrow=c(1,1))
qqPlot(fitA1,simulate=FALSE)
range(rstudent(fitA1))
range(rstandard(fitA1))
## which one is the outlier, if any?
which(abs(rstudent(fitA1))==max(abs(rstudent(fitA1))))
par(mfrow=c(2,2))
plot(fitA1)
#### We conclude that there are no outliers here




################ Weighted analysis

ll = 1
vseq = seq(0,10,by=0.01) #Possible weights
for (v in 1:length(vseq)) {
w = c(rep(1,70),rep(vseq[v],66))
fitw <- lm(clo ~ tOut + tInOp + sex + tInOp:sex, weights = w, data = CS1)
#summary(fitw)
ll[v] <- logLik(fitw) #collect the log-likelihood for each weight
}

par(mfrow=c(1,1)) #plot the weights and log-likelihood
plot(vseq,ll,xlab = "variance ratio",ylab="log-likelihood")

# Find the best weight
index = which.max(ll)
optimalweight = vseq[index]
w = c(rep(1,70),rep(optimalweight,66))

# Do a new fit with this weight
fitA2 <- lm(clo ~ tOut + tInOp + sex + tInOp:sex, weights = w, data = CS1)
summary(fitA2)
par(mfrow=c(2,2))
plot(fitA2)
########### p-values get higher, but:

par(mfrow=c(3,2))
plot(residuals(fitA2)) #Residuals versus obs number (not subj id)
plot(fitted(fitA2),residuals(fitA2)) #residuals versus fitted values
plot(CS1$clo,residuals(fitA2)) #residuals versus response, this will always be linear!
plot(CS1$tOut,residuals(fitA2)) #residuals versus tOut
plot(CS1$tInOp,residuals(fitA2)) #residuals versus tInOp
plot(CS1$sex,residuals(fitA2)) #residuals versus sex
# Jan: find pearson residuals, maybe option in residuals(fit)
# Otherwise these plots don't make sense


####### Check for outliers
par(mfrow=c(1,1))
qqPlot(fitA2,simulate=FALSE)
range(rstudent(fitA2))
range(rstandard(fitA2))
## which one is the outlier, if any?
which(abs(rstudent(fitA2))==max(abs(rstudent(fitA2))))
par(mfrow=c(2,2))
plot(fitA2)
#### We don't conclude that 79 is an outlier (do we?)

#OPtIONAL  Do a new fit without obs 79
CS2 = CS1[-79,]
w2 = w[-79]
fitA3 <- lm(clo ~ tOut + tInOp + sex + tInOp:sex, weights = w2, data = CS2)
summary(fitA3)
par(mfrow=c(1,1))
qqPlot(fitA3,simulate=FALSE)
par(mfrow=c(2,2))
plot(fitA3)
# Now looks dramatically different (and beter)

# Plot prediction interval
dev.off()
par(mfrow=c(1,2))
tinop <- mean(CS1$tInOp[1:70]) # For given tInOp!
#Sex <- 1 # For given Height!
newdat=data.frame(tOut=CS1$tOut[1:70],tInOp=tinop,sex=CS1$sex[1:70])
pred <- predict(fitA2,se=TRUE, newdata=newdat,interval="prediction")
conf <- predict(fitA2,se=TRUE, newdata=newdat,interval="confidence")

#plot(CS1$tOut[1:69],CS1$clo[1:69])
matplot(CS1$tOut[1:70],pred$fit,type="l",col=c(1,3,3),
        lty=c(1,2,2),main="Female sex",xlab="tOut",ylab="Prediction")


matlines(CS1$tOut[1:70],conf$fit,type="l",col=c(1,3,3),
         lty=c(1,2,2))
points(CS1$tOut[1:70],CS1$clo[1:70])
############# Nr 2
tinop <- mean(CS1$tInOp[71:136]) # For given tInOp!
newdat2=data.frame(tOut=CS1$tOut[71:136],tInOp=tinop,sex=CS1$sex[71:136])
pred2 <- predict(fitA2,se=TRUE, newdata=newdat2,interval="prediction",weights = rep(2.93,136,1))
conf2 <- predict(fitA2,se=TRUE, newdata=newdat2,interval="confidence")

matplot(CS1$tOut[71:136],pred2$fit,type="l",col=c(1,3,3),
        lty=c(1,2,2),main="Male sex",xlab="tOut",ylab="Prediction")


matlines(CS1$tOut[71:136],conf2$fit,type="l",col=c(1,3,3),
         lty=c(1,2,2))
points(CS1$tOut[71:136],CS1$clo[71:136])


############## See if subject ID can be excluded. 
# Load a new subset including subjId
CS3 <- subset(CS, select = -c(day) )
CS3 = CS3[-79,]
CS3$subjId=factor(CS3$subjId)
# Then we look at residuals versus subject ID
par(mfrow=c(1,1))
plot(sort(CS3$subjId),residuals(fitA3),xlab="Subject ID",ylab="Residuals") #residuals versus subjectId
### We can see that the within-subject variance is small, but between-subject is high
# Could be nice to include subject ID