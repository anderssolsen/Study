rm(list=ls())
##################################################
# Data
##################################################
par(mfrow=c(1,1))
V<-c(1065,1071,1075,1083,1089,1094,1100,1107,1111,
     1120,1128,1135)
BD<-c(2,3,5,11,10,21,29,48,56,88,98,99)
tr<-rep(100,12)
data=list(V=V,BD=BD,tr=tr)
plot(BD/tr~V,data=data)
abline(lm((BD/tr~V)))
##################################################

##################################################
# Link-func.
##################################################

par(mfrow=c(2,2))
# Identity
x <- seq(-5,5,by=0.001)
plot(x,x,lwd=2,type="l",col=2,xlab=expression(mu),
     ylab=expression(g(mu)),
     main=expression(bold(R) %->% bold(R)))
legend(x=min(x),y=max(x),legend="Identity",lty=1,
       col=2,lwd=2)


# log
x <- seq(0,10,by=0.001)
plot(x,log(x),col=2,lwd=2,type="l",lty=1,
     xlab=expression(mu),
     ylab=expression(g(mu)),
     main=expression(bold(R)["+"] %->% bold(R)))
legend(x=10,y=-5,legend="log",lty=1,
       col=2,lwd=2,xjust=1)





x <- seq(0,10,by=0.001)
matplot(x,cbind(sqrt(x),1/x,1/x^2),col=1:3,lwd=2,
        type="l",lty=1,ylim=c(0,4),
        xlab=expression(mu),
        ylab=expression(g(mu)),
        main=expression(bold(R)["+"] %->% bold(R)["+"]))
legend(x=2,y=4,legend=c(expression(sqrt(mu)),
                        expression(1/mu),
                        expression(1/mu^2)),
       lty=1,col=1:3,lwd=2,xjust=0)



x <- seq(0.001,0.999,by=0.001)
matplot(x,cbind(log(x/(1-x)),log(-log(x)),
                log(-log(1-x)),qnorm(x),qcauchy(x)),
        col=1:5,lwd=2,type="l",lty=1,ylim=c(-4,5),
        xlab=expression(mu),
        ylab=expression(g(mu)),
        main=expression(group("[",list(0, 1),"]") %->%
                          bold(R)))
#legend(x=0.2,y=5,legend=c("logit","log-log","clog-log",
#                         "probit","Cauchit"),lty=1,col=1:5,
#     lwd=2,xjust=0)
##################################################



##################################################
###  Canonical link (Spark example)
par(mfrow=c(1,1))
plot(log((BD/tr)/(1-BD/tr))~V,data=data)


##################################################
#  Specification in R
model<-cbind(BD,tr-BD)~V
logis.glm<-glm(model,family=binomial,
               data=data)

summary(logis.glm)

par(mfrow=c(1,1))
plot(log((BD/tr)/(1-BD/tr))~V,data=data)
lines(data$V,predict(logis.glm,type="link"))


##################################################
## Properties of ML estimator (canonical link)
mu <- logis.glm$fitted.values
W <- diag((mu*(1-mu))*data$tr)
diag(W)

X <- model.matrix(logis.glm)

## Parameter cavariance
Vbeta <- solve(t(X)%*%W%*%X)
Vbeta
summary(logis.glm)$cov.unscaled


##################################################
## Dianostics
par(mfrow=c(2,2),mar=c(4,4,4,4))
plot(logis.glm)
## Not perfect 
##################################################

##################################################
## Linear predictor
par(mfrow=c(1,1))
Vplot <- seq(1050,1150)
Xplot <- cbind(1,Vplot)
eta <- Xplot%*%coef(logis.glm)

plot(Vplot,eta,type="l")
p <- BD/tr
points(V,log(p/(1-p)),pch=19)
D <- Xplot %*% Vbeta %*% t(Xplot)
matlines(Vplot,cbind(eta-2*sqrt(diag(D)),
                     eta+2*sqrt(diag(D))),lty=2,col=2)

## Comparing with R
points(V,predict(logis.glm,type="link",se=TRUE)$fit,pch=19,cex=0.5)
points(V,predict(logis.glm,type="link",se=TRUE)$fit+
      2*predict(logis.glm,type="link",se=TRUE)$se.fit,col=2,pch=19,cex=0.5)
points(V,predict(logis.glm,type="link",se=TRUE)$fit-
      2*predict(logis.glm,type="link",se=TRUE)$se.fit,col=2,pch=19,cex=0.5)
##################################################

##################################################
## Fitted values
mu <- as.numeric(exp(eta)/(1+exp(eta)))
plot(Vplot,mu,type="l")
p <- BD/tr
points(V,p,pch=19)
Dmu <- diag((mu*(1-mu)))%*%D%*%diag((mu*(1-mu)))
matlines(Vplot,cbind(mu-2*sqrt(diag(Dmu)),
                     mu+2*sqrt(diag(Dmu))),lty=2,col=2)

## Comparing with R
mu <- predict(logis.glm,type="response",se=TRUE)$fit
se <- predict(logis.glm, type="response", se=TRUE)$se.fit
points(V, mu, pch=19,cex=0.5)
points(V,mu + 2 * se, col=2, pch=19, cex=0.5)
points(V,mu - 2 * se, col=2, pch=19, cex=0.5)
##################################################

##################################################
## Residuals
## Leverage
H <- sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W)

par(mfrow=c(2,2),mar=c(4,4,4,4))
rD <- residuals(logis.glm,type="deviance")
d <- 2*(p*log(p/mu)+(1-p)*log((1-p)/(1-mu))) ## unit deviance (binomial)
sign(p-mu)*sqrt(100*d)-rD ## diff it works..

plot(logis.glm,which=1)## Deviance residals
points(predict(logis.glm,type="link"),rD,pch=19)

stdDr <- rD * 10/sqrt(1-diag(H)) ## Studentlized deviance residuals
n <- length(rD)
plot(logis.glm,which=2)
points(qnorm((1:n-0.5)/n),sort(stdDr), pch=19) ## Deviance residuals

plot(logis.glm,which=3)
points(predict(logis.glm,type="link"),sqrt(abs(stdDr)), pch=19) ## std. Deviance res

# Studentlized Pearson residuals
rPstd <- residuals(logis.glm,type="pearson")/sqrt(1-diag(H)) 
plot(logis.glm,which=5)
points(diag(H),residuals(logis.glm,type="pearson")/sqrt(1-diag(H)), pch=19)
##################################################


##################################################
## Skip (Cooks distance)
## plot(logis.glm,which=4)
## cookD <- rPstd^2/2*diag(H)/(1-diag(H))
## points(1:12,cookD, pch=1)
## 
## plot(logis.glm,which=6)
## points(diag(H)/(1-diag(H)),cookD, pch=19)
##################################################

##################################################
## A closer look
par(mfrow=c(1,1))
plot(data$V,rD,pch=19) ## deviance res vs. Volt

## Second order model 
model2<-cbind(BD,tr-BD)~V+I(V^2)
logis.glm2<-glm(model2,family=binomial,
               data=data)

par(mfrow=c(2,2))
plot(logis.glm2)
## seems like a reasonable model..

##################################################
# Analysis of deviance table
anova(logis.glm,test="Chisq") ## Start by simple model
100*sum(d)
## null model
mu.null <- sum(BD)/sum(tr)
## unit deviance (null model)
d.null <- 2*(p*log(p/mu.null)+(1-p)*log((1-p)/(1-mu.null))) 
100*sum(d.null)

## Test for sufficiency 
1 - pchisq(21.02,df=10) ## I.e. fail the goodness of fit test
## Test for parameter
pchisq(762.1,df=1,lower.tail=FALSE) ## I.e. Volt significant
                   ##(but we fail the goodness of fit test)

##################################################
## More complicated model
anova(logis.glm2,test="Chisq")
## Test for sufficiency
1 - pchisq(6.16,df=9) ## I.e. pass the goodness of fit test
## Test for parameter
pchisq(14.86,df=1,lower.tail=FALSE) ## I.e square dependece significant

## Compare models
anova(logis.glm,logis.glm2,test="Chisq") ## Same conclusion

summary(logis.glm2)
plot(logis.glm2)

##################################################
## Overdispersion (using the simple model)
anova(logis.glm)
logis.glm.OD <- glm(model,family=quasibinomial,data=data)
summary(logis.glm.OD)
sigsq.dev <- 21.02/10
(sigsq.pear <- sum(residuals(logis.glm,type="pearson")^2)/10)
(sigsq.devR <- sum(residuals(logis.glm,type="deviance")^2)/10)



## analysis of deviance table
anova(logis.glm.OD,test="F")
762.1/2.014502
1-pf(378.31,df1=1,df2=10)

##################################################
## Comparing 3 models
Vplot <- data.frame(V=Vplot)
mu1 <- predict(logis.glm,newdata=Vplot,type="response",se=TRUE)$fit
se1 <- predict(logis.glm,newdata=Vplot, type="response", se=TRUE)$se.fit
mu2 <- predict(logis.glm2,newdata=Vplot,type="response",se=TRUE)$fit
se2 <- predict(logis.glm2,newdata=Vplot, type="response", se=TRUE)$se.fit
muOD <- predict(logis.glm.OD,newdata=Vplot,type="response",se=TRUE)$fit
seOD <- predict(logis.glm.OD,newdata=Vplot, type="response", se=TRUE)$se.fit

par(mfrow=c(1,1))
matplot(Vplot,cbind(mu1,mu2,muOD),type="l")
polygon(c(Vplot$V,rev(Vplot$V)),c(muOD+2*seOD,rev(muOD-2*seOD)),border=FALSE,
        col=gray(0.75))
polygon(c(Vplot$V,rev(Vplot$V)),c(mu1+2*se1,rev(mu1-2*se1)),border=FALSE,
        col=gray(0.5))
polygon(c(Vplot$V,rev(Vplot$V)),c(mu2+2*se2,rev(mu2-2*se2)),border=FALSE,
        col=gray(0.25))
matlines(Vplot,cbind(mu1,mu2,muOD),lty=1,lwd=2)
points(V,p,pch=19,col="blue")


##################################################
# Parameter estimates
summary(logis.glm)
summary(logis.glm.OD) ## note "large" p-values (overdispersion parameter)
summary(logis.glm2)
##################################################


##################################################
## Other link functions
anova(glm(model,family=binomial(probit)))
anova(glm(model,family=binomial(cauchit)))
anova(glm(model,family=binomial(cloglog)))
## Seems like cloglog would be appropriate

cloglog.glm<-glm(model,data=data,
                 family=binomial(cloglog))
par(mfrow=c(2,2))
plot(cloglog.glm)
summary(cloglog.glm)
## x11()
par(mfrow=c(2,2))
plot(logis.glm2)

cloglog.glm2<-glm(model2,data=data,
                 family=binomial(cloglog))
anova(cloglog.glm2,test="Chisq")
## Second order term not sigificant

par(mfrow=c(1,1))
matplot(V,cbind(residuals(logis.glm2),residuals(cloglog.glm)))
lines(c(0,1300),c(0,0),col=2,lty=2,lwd=2)



cloglog.glm$aic
cloglog.glm2$aic
logis.glm2$aic
## So we should choose cloglog.glm


##################################################
## Data domain (compare the models)
V <- 1020:1150
par(mfrow=c(1,2))
pred<-predict(cloglog.glm,newdata=data.frame(V=V),
              type="response",se=TRUE)
pred2<-predict(logis.glm2,newdata=data.frame(V=V),
               type="response",se=TRUE)
pred.logis<-predict(logis.glm,newdata=data.frame(V=V),
              type="response",se=TRUE)
matplot(V,pred$fit+cbind(0,2*pred$se.fit,-2*pred$se.fit),
        type="l",lty=c(1,2,2),col=c(1,2,2),ylab="p(V)",xlab="V")
points(data$V,BD/tr,pch=19,cex=0.5)
lines(c(0,1300),c(0,0),col=1,lty=1,lwd=1)
lines(c(0,1300),c(1,1),col=1,lty=1,lwd=1)

matplot(V,pred2$fit+cbind(0,2*pred2$se.fit,-2*pred2$se.fit),
        type="l",lty=c(1,2,2),col=c(3,4,4),ylab="p(V)",xlab="V",
        ylim=c(0,1))
points(data$V,BD/tr,pch=19,cex=0.5)
lines(c(0,1300),c(0,0),col=1,lty=1,lwd=1)
lines(c(0,1300),c(1,1),col=1,lty=1,lwd=1)


# Linear domain
pred<-predict(cloglog.glm,newdata=data.frame(V=V),
              type="link",se=TRUE)
matplot(V,pred$fit+cbind(0,2*pred$se.fit,-2*pred$se.fit),
        type="l",lty=c(1,2,2),col=c(1,2,2),ylab="p(V)",xlab="V")
points(data$V,log(-log(1-BD/tr)),pch=19)

pred<-predict(logis.glm2,newdata=data.frame(V=V),
              type="link",se=TRUE)
matplot(V,pred$fit+cbind(0,2*pred$se.fit,-2*pred$se.fit),
        type="l",lty=c(1,2,2),col=c(1,2,2),ylab="p(V)",xlab="V")
points(data$V,log(BD/tr/(1-BD/tr)),pch=19)
#########################################################


##################################################
## Prediction interval
V<-c(1065,1071,1075,1083,1089,1094,1100,1107,1111,
     1120,1128,1135)
BD<-c(2,3,5,11,10,21,29,48,56,88,98,99)
tr<-rep(100,12)

pred <- predict(cloglog.glm,type="response",se=TRUE)
pred2 <- predict(cloglog.glm,type="link",se=TRUE)

## Direct
confint <- cbind(pred$fit,pred$fit+1.96*pred$se.fit,
                 pred$fit-1.96*pred$se.fit)
confint[confint>1]<-1 ## Might have values above 1
confint[confint<0]<-0 ## Might have values below 0
n <-100

## By inverse link
confint2 <- 1-  exp(-exp(cbind(pred2$fit,pred2$fit+1.96*pred2$se.fit,
                 pred2$fit-1.96*pred2$se.fit)))

par(mfrow=c(1,1))
matplot(V,cbind(pred$fit,
                pred$fit+2*pred$se.fit,
                pred$fit-2*pred$se.fit),type="l",col="red",
        lty=c(1,2,2))
points(V,BD/tr,pch=19)


## Bootstrap part
set.seed(1255)
k <- 1000
Sim <- matrix(ncol=12,nrow=k)
Sim2 <- Sim
Sim3 <- Sim
for(i in 1:k){
    r <- rnorm(12, pred$fit,pred$se.fit)
    r[r>1]<-1; r[r<0]<-0
    Sim[i,] <- rbinom(12,prob=r,size=n)

    r <- rnorm(12, pred2$fit,pred2$se.fit)
    r <- 1-exp(-exp(r))
    Sim3[i,] <- rbinom(12,prob=r,size=n)

    
    r <- rbinom(12,prob=pred$fit,size=100)
    fit <- glm(cbind(r, tr - r) ~ V,
               family=binomial(cloglog))
    pred.tmp <- predict(fit,type="response")
    Sim2[i,] <- rbinom(12,prob=pred.tmp,size=n)
}


Pred.Int1 <- t(apply(Sim,2,quantile,probs=c(0.025,0.975))/n)
Pred.Int2 <- t(apply(Sim2, 2, quantile,
                     probs=c(0.025,0.975))/n)
Pred.Int3 <- t(apply(Sim3, 2, quantile,
                     probs=c(0.025,0.975))/n)

matlines(V,Pred.Int1,col="blue",lty=2)

## Compare:

## Lower bound
round(cbind(Pred.Int1[ ,1], Pred.Int2[ ,1],
            Pred.Int3[ ,1]),
      digits=3)

## Upper bound
round(cbind(Pred.Int1[ ,2], Pred.Int2[ ,2],
            Pred.Int3[ ,2]),
      digits=3)
## Very close


############################################################
## Odds ratio (change in odds ratio for one unit of change in Volt)
par(mfrow=c(1,1))
p1 <- predict(cloglog.glm,newdata=data.frame(V=seq(1060,1140)),type="response")
p2 <- predict(cloglog.glm,newdata=data.frame(V=seq(1061,1141)),type="response")
OR <- p1/(1-p1)/(p2/(1-p2))
plot(seq(1061,1141),OR,type="l",ylim=c(0,1))

##################################################
## Canonical link
p1 <- predict(logis.glm,newdata=data.frame(V=seq(1060,1140)),type="response")
p2 <- predict(logis.glm,newdata=data.frame(V=seq(1061,1141)),type="response")
OR <- p1/(1-p1)/(p2/(1-p2))
lines(seq(1061,1141),OR,type="l",ylim=c(0,1),col=2)
##################################################


##################################################
## Canonical link (non-linear model)
p1 <- predict(logis.glm2,newdata=data.frame(V=seq(1060,1140)),type="response")
p2 <- predict(logis.glm2,newdata=data.frame(V=seq(1061,1141)),type="response")
OR <- p1/(1-p1)/(p2/(1-p2))
lines(seq(1061,1141),OR,type="l",ylim=c(0,1),col=3)

##################################################
##################################################


##############################
## Appendix: Leverage for non-canonical link
mu <- cloglog.glm$fitted.values
W <- diag(100* (log(1-mu)^2*(1-mu))/mu) ## The change is here
X <- model.matrix(logis.glm)

H <- sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W)
diag(H)
sum(diag(H))



par(mfrow=c(2,2),mar=c(4,4,4,4))
rD <- residuals(cloglog.glm,type="deviance")
plot(cloglog.glm,which=1)
points(predict(cloglog.glm,type="link"),rD,pch=19)

stdDr <- rD * 10/sqrt(1-diag(H))
plot(cloglog.glm,which=2)
points(qnorm((1:12-0.5)/12),sort(stdDr), pch=19)

plot(cloglog.glm,which=3)
points(predict(cloglog.glm,type="link"),sqrt(abs(stdDr)), pch=19)

rPstd <- residuals(cloglog.glm,type="pearson")/sqrt(1-diag(H))
plot(cloglog.glm,which=5)
points(diag(H),residuals(cloglog.glm,type="pearson")/sqrt(1-diag(H)), pch=19)


##################################################
## End
##################################################
