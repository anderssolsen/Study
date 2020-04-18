## Assignment 2 part B

rm(list=ls())
setwd("C:/Users/ander/Documents/10. Semester/Advanced data analysis/Assignment 2")
library("ggplot2")
library("dplyr")
library("GGally")
library("car")
library("MASS")
library("ellipse")

dat = read.csv("dat_count.csv",sep=";")
dat = dat[order(dat$sex),] #sort by sex
#ggpairs(dat)

datsub = subset(dat,select = -c(day,subjId))

# Short presentation
par(mfrow=c(2,4))
plot(datsub$clo~datsub$time, xlab='time', ylab='clothing')
plot(datsub$clo~datsub$sex, xlab='time', ylab='clothing')
plot(datsub$clo~datsub$tOut, xlab='time', ylab='clothing')
plot(datsub$clo~datsub$tInOp, xlab='time', ylab='clothing')
plot(datsub$clo/datsub$nobs~datsub$time, xlab='time', ylab='response')
plot(datsub$clo/datsub$nobs~datsub$sex, xlab='time', ylab='response')
plot(datsub$clo/datsub$nobs~datsub$tOut, xlab='time', ylab='response')
plot(datsub$clo/datsub$nobs~datsub$tInOp, xlab='time', ylab='response')



################## B.1: Do a GLM based on binomial ####################

# The binomial will be built using clo/nobs as the response variable.
# Because nobs and time represent almost the same thing, we could argue
#   that time should be left out of the model. 
# Below is formula (full model) and formula 2 (no 3-way and 4-way interactions)
# Because we want to do an initial "goodness of fit" test, we fit the full model

# There are 5 liink functions for the binomial. All should be tested

datsub$resp<-cbind(dat$clo,datsub$nobs-dat$clo)
formula = resp ~ time*sex*tOut*tInOp
formula2 = resp ~ time*sex*tOut*tInOp-time:sex:tOut:tInOp-time:sex:tOut-time:sex:tInOp-sex:tOut:tInOp-time:tOut:tInOp

fitbin1 = glm(formula=formula,family=binomial(link="logit"),data=datsub)
summary(fitbin1)#AIC 290
1-pchisq(160.02,120)
fitbin2 = glm(formula=formula,family=binomial(link="probit"),data=datsub)
summary(fitbin2)#AIC 289
1-pchisq(159.82,120)
fitbin3 = glm(formula=formula,family=binomial(link="cauchit"),data=datsub)
summary(fitbin3)#AIC 289
1-pchisq(159.32,120)
fitbin4 = glm(formula=formula,family=binomial(link="log"),data=datsub)
summary(fitbin4)#AIC 290
1-pchisq(160.32,120)
fitbin5 = glm(formula=formula,family=binomial(link="cloglog"),data=datsub)
summary(fitbin5)#AIC 290
1-pchisq(160.18,120)

# None of the models are a good fit, since all are below alfa=5%
# We have many observations, so the Chisq is a good approximation. 
# Thus, we need to check if there are other problems (w7, slide 22)
# Possible reasons are:
#   Incorrect linear predictor
#   Incorrect link function (we tried all)
#   Outliers
#   Influential observations
#   Incorrect choice of distributions (we were asked for binomial)
# Thus we need to check the residuals!

resDev <- residuals(fitbin2,type="deviance")
par(mfrow=c(2,2))
plot(datsub$time,resDev)
plot(datsub$sex,resDev)
plot(datsub$tOut,resDev)
plot(datsub$tInOp,resDev)

resPears <- residuals(fitbin2,type="pearson")
par(mfrow=c(2,2))
plot(datsub$time,resPears)
plot(datsub$sex,resPears)
plot(datsub$tOut,resPears)
plot(datsub$tInOp,resPears)

# We could argue for removing the two observations with time=1hr and 2hr?
# Maybe we should do a leverage plot to confirm this?

datsub2 = datsub[order(datsub$time),] #sort by time
datsub2 = datsub2[-c(1,2),] #remove first 2

# Compute goodness of fit again

fitbin21 = glm(formula=formula,family=binomial(link="logit"),data=datsub2)
summary(fitbin21)#AIC 288
1-pchisq(158.75,118)
fitbin22 = glm(formula=formula,family=binomial(link="probit"),data=datsub2)
summary(fitbin22)#AIC 288
1-pchisq(158.62,118)
fitbin23 = glm(formula=formula,family=binomial(link="cauchit"),data=datsub2)
summary(fitbin23)#AIC 287
1-pchisq(157.28,118)
fitbin24 = glm(formula=formula,family=binomial(link="log"),data=datsub2)
summary(fitbin24)#AIC 289
1-pchisq(158.99,118)
fitbin25 = glm(formula=formula,family=binomial(link="cloglog"),data=datsub2)
summary(fitbin25)#AIC 288
1-pchisq(158.88,118)

# This did not have an effect!
# So the reason must be overdispersion. 
# Jan: We should NOT do a weighted analysis as in assignment 1
#   However, we should be able to say WHY there is overdispersion. 
# Anders: I tried to do a weighted analysis, but it didn't work.

###### Using overdispersion and the original data set
# And the chosen link function: cauchit
# Beware: the fit is overwritten for each update

fitbin31 = glm(formula=formula,family=quasibinomial(link="cauchit"),data=datsub)
summary(fitbin31) 

drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:sex:tOut:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:tOut:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:sex:tOut)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:sex:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -sex:tOut:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -tOut:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -sex:tOut)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:tOut)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time:sex)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -sex:tInOp)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -tOut)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -time)
drop1(fitbin31, test = "F")
fitbin31 = update(fitbin31,.~. -tInOp)

# So the final model contains ONLY sex as predictor. 

# What needs to be done with Binomial: 
#   Present the model with plots and confidence interval
#   Compute confidence interval and prediction interval as in last assignment. 
#   Compute odds ratio (if you're a female, there x% higher chance you'll change clothes)



################## B.2: Do a GLM based on Poisson ####################
# In binomial, we used 'nobs' as the offset. 
# In poisson, it makes more sense to talk about rates.
# So we should use time as offset and discard nobs. (Jan's words)

formulapois = clo ~ offset(log(time))+sex*tOut*tInOp


# Poisson has 3 link functions, and we test the goodness of fit again.
# However, the identity and inverse link functions say that they need starting values.
# We could supply these with start=rep(0.5,8) but still doesn't quite work
fitpois1 = glm(formula = formulapois,family=poisson(link="log"),data=datsub)
summary(fitpois1) #AIC 277
1-pchisq(143.63,128)
#fitpois2 = glm(formula = formulapois,family=poisson(link="identity"),data=datsub,start=rep(0.5,8))
#summary(fitpois2)
#1-pchisq(158.88,118)
#fitpois3 = glm(formula = formulapois,family=poisson(link="inverse"),data=datsub)
#summary(fitpois3)
#1-pchisq(158.88,118)

# The poisson model with the "log" link function passed the test!


drop1(fitpois1,test="F")
fitpois1 = update(fitpois1,.~. -sex:tOut:tInOp)
summary(fitpois1)

drop1(fitpois1,test="F")
fitpois1 = update(fitpois1,.~. -tOut:tInOp)
summary(fitpois1)

drop1(fitpois1,test="F")
fitpois1 = update(fitpois1,.~. -sex:tOut)
summary(fitpois1)

drop1(fitpois1,test="F")
fitpois1 = update(fitpois1,.~. -sex:tInOp)
summary(fitpois1)

drop1(fitpois1,test="F")
fitpois1 = update(fitpois1,.~. -tOut)
summary(fitpois1)

drop1(fitpois1,test="F")
fitpois1 = update(fitpois1,.~. -tInOp)
summary(fitpois1)

# We end up with the same model as for the binomial!
# Now, we still need to provide some residual plots and present the model as before
# Odds ratio could also be nice.
# Beware of comparing AIC directly


par(mfrow=c(2,2))
plot(fitpois1)

Rd<-residuals(fitpois1,type='deviance')
Rp<-residuals(fitpois1,type='pearson')
par(mfrow=c(2,2))
plot(datsub$sex,Rd, xlab='sex', ylab='Deviance residuals')
plot(fitted(fitpois1),Rd, xlab='fitted', ylab='Deviance residuals')
plot(datsub$sex,Rp, xlab='sex', ylab='Deviance residuals')
plot(fitted(fitpois1),Rp, xlab='fitted', ylab='Deviance residuals')






