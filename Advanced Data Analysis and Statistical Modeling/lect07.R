setwd("C:/Users/ander/Documents/10. Semester/Advanced data analysis/Week07")
##################################################
## 1: Seeds (binomial) Slide 21.
##################################################
rm(list=ls())
dat <- read.table("seed.dat", sep = ";",  header = TRUE)
head(dat)
str(dat)

dat$variety <- as.factor(dat$variety)
dat$root <- as.factor(dat$root)
dat$resp <- cbind(dat$y, (dat$n - dat$y))
fit1 <- glm(resp ~ variety * root,
            family = binomial,
            data = dat)
fit1

par(mfrow=c(2,2))
plot(fit1)
## No concerns here

## just to illustrate calculations....
w <- dat$n
X <- model.matrix(fit1)
eta <- predict(fit1,type="link")
mu <- predict(fit1,type="response")
Vmu <- mu*(1-mu)
W <- diag(1/(1+exp(-eta))^4*exp(-2*eta)*1/Vmu)*w
H <- sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W)
h <- diag(H)
stdDevRes <- residuals(fit1,type="deviance")*sqrt(w)/(sqrt(1-diag(H)))

par(mfrow=c(2,2)) 
plot(fit1,which=1)
points(predict(fit1,type="link"),residuals(fit1,type="deviance"),
       pch=19)

plot(fit1,which=2)
points(qnorm((1:21-0.5)/21),sort(stdDevRes),
       pch=19)

plot(fit1,which=3)
points(predict(fit1,type="link"),sqrt(abs(stdDevRes)),
       pch=19)

plot(fit1,which=5)
rps <- (dat$y/dat$n-mu)/sqrt(mu*(1-mu)*(1-h)/dat$n)
points(h,rps,pch=19)

# Goodness of fit test
summary(fit1)
(pval <- 1 - pchisq(33.28, 17))
## So the assumption is rejected, Hence overdispersion

par(mfrow=c(1,2))
resDev <- residuals(fit1,type="deviance")
plot(jitter(as.numeric(dat$variety), amount=0.1), resDev,
     xlab="Variety", ylab="Deviance residuals", cex=0.6,
     axes=FALSE)
box()
axis(1,label=c("O.a. 75", "O.a. 73" ),at=c(1,2))
axis(2)
plot(jitter(as.numeric(dat$root),
            amount=0.1), resDev, xlab= "Root" ,
ylab="Deviance residuals", cex=0.6, axes=FALSE)
box()
axis(1,label=c("Bean","Cucumber"),at=c(1,2))
axis(2)
## Doesn't seems to be systematic effects here..

# Including overdispersion
fit2 <- glm(resp ~ variety * root,
          family = quasibinomial, data = dat)
summary(fit2)

# JUST TO COMPARE THIS MODEL IS CONSIDERED WRONG HERE
confint(fit1)
confint(fit2) # Wider CI's

## Model reduction
drop1(fit2, test = "F")

fit3 <- glm(resp ~ variety + root,
            family = quasibinomial,
            data = dat)
drop1(fit3, test = "F")

fit4<-glm(resp ~ root,
          family = quasibinomial, data = dat)
drop1(fit4, test="F")
summary(fit4)

## Interpretation of parameters
par <- coef(fit4)
par ## Germination more likeli on root2

std<-sqrt(diag(vcov(fit4)))
std

## CI of parameters (Using t-dist)
par+std%o%c(lower=-1,upper=1)*qt(0.975,19)

## CI of parameters (Using normal dist)
confint.default(fit4)
# or with profile likelihiid
confint(fit4) ## profile likelihood based

# Odds ratio
exp(coef(fit4)[2])
exp(confint(fit4)[2, ]) ## 95% confidence interval



##################################################
## 2: Accident rates (Posison) Slide 29
##################################################
dat <- data.frame(sex = c("F", "M"),
                  years = c(17.3,21.4),
                  y = c(175, 320))

fit1<-glm(y~offset(log(years))+sex,family=poisson,data=dat)
## Note the included offset.
anova(fit1,test= "Chisq")

summary(fit1)
## More likely that elderly males are involved in accidents

## Odds ratio:
exp(coef(fit1)[2])
## or
(320/21.4)/(175/17.3)

####################################################################
## 3: Challenger (Binomial) Slide 32
####################################################################
## slightly more detailed dataset
dat<-read.table("challenger_data.txt",header=TRUE)
dat
par(mfrow=c(1,2))
plot(dat$temp, dat$failed, xlab= ' Temperature ' ,
     ylab= ' No damaged (out of 6) ' )
plot(dat$pres, dat$failed, xlab= ' Pressure ' ,
     ylab= ' No damaged (out of 6) ' )

## Response variable
dat$resp<-cbind(dat$failed,dat$n-dat$failed)

## Model 
fit0<-glm(formula = resp ~ temp+pres,
          family = binomial(link = logit),
          data = dat)
summary(fit0)

## Different parametrization
fit02<-glm(formula = I(failed/n) ~ temp+pres,
          family = binomial(link = logit),
          data = dat,weights=dat$n)

## Same results
summary(fit02)
summary(fit0)


1-pchisq(9.4,df=20) ## note the chisq approximation is bad
par(mfrow=c(2,2))
plot(fit0)
## Not that many zeros and small samples imply "bad" dianostic plots

drop1(fit0, test= "Chisq" )


fit1<-glm(formula = resp ~ temp,
          family = binomial(link = logit),
          data = dat)
drop1(fit1, test= "Chisq" )

summary(fit1)
1-pchisq(9.5,df=21) ## accept, but the chi^2 app. is not good in this case

## CI direct use of Wald
par(mfrow=c(1,2))
tmp <- 31:85
pred<-predict(fit1, type= "response" , newdata=data.frame(temp=tmp),se=TRUE)
plot(tmp, pred$fit, type= "l" , lwd=3, col= "red" , xlab= "Temperature" , ylab= "P(damage)" )
points(dat$temp,dat$fail/dat$n)
lines(tmp,pred$fit+2*pred$se.fit)
lines(tmp,pred$fit-2*pred$se.fit)
coef(fit1)

## CI wald + inverse transformation
logit <- function(x){  exp(x)/(1+exp(x))    }
pred<-predict(fit1, type= "link" , newdata=data.frame(temp=tmp),se=TRUE)
plot(tmp, logit(pred$fit), type= "l" , lwd=3, col= "red" , xlab= "Temperature" , ylab= "P(damage)" )
points(dat$temp,dat$fail/dat$n)
lines(tmp,logit(pred$fit+2*pred$se.fit))
lines(tmp,logit(pred$fit-2*pred$se.fit))
coef(fit1)
##################################################
##
##################################################



####################################################
## 4: Example 4.12 (Multinomial) Slide 37
####################################################
rm(list=ls())
vdiss <- c( 234,  41,  42, 35)
diss  <- c( 559, 100,  76, 48)
neu   <- c(1157, 145,  89, 39)
sat   <- c(5826, 602, 254, 95)
vsat  <- c(2553, 237,  72, 27)

## Setting up the table
(tab <- cbind(vdiss, diss, neu, sat, vsat))


## Cumulative table
tot <- rowSums(tab)
tabp <- tab/tot
(accp <- t(apply(tabp,1,cumsum)))
t <- c(0,2,5,7)

## Plot observations
par(mfrow=c(1,1))
matplot(t, accp[ ,-5], pch = 19, type = "b",col=1,ylim=c(0,1))
for(i in 1:length(t)){
    for(j in 1:4){
        ph <- sum(tab[i,1:j])/sum(tab[i, ])
        CI <- ph +c(-1,1)*2*sqrt(ph*(1-ph)/sum(tab[i, ]))
        matlines(t[i]*c(1,1),CI)
    }
}


## Data frame for modelling
df.bus <- data.frame(Freq = c(vdiss,diss,neu,sat,vsat),
                     sat = factor(rep(c("vdiss","diss","neu","sat","vsat"),
                                      each=4),
                                  ordered=TRUE,
                                  levels=c("vdiss","diss","neu","sat","vsat")),
                     delay = rep(c(0,2,5,7),5))
df.bus$sat ## Ordered

## Proportional odds model 
library(MASS)
m1.p <- polr( sat ~ factor(delay), data = df.bus,weights=Freq)
m2.p <- polr( sat ~ delay, data = df.bus,weights=Freq)
anova(m1.p,m2.p) ## So simple model 

m3.p <- polr( sat ~ delay + I(delay^2), data = df.bus,weights=Freq)
anova(m2.p,m3.p) ## Still simple model 

## Other link functions
m4.p <- polr( sat ~ delay, data = df.bus,weights=Freq,method="loglog")
m5.p <- polr( sat ~ delay, data = df.bus,weights=Freq,method="cloglog")
AIC(m4.p)
AIC(m2.p) ## So we should use loglog

## Model not including the ordering (i.e. not proportional odd)
library(nnet)
m1.m <- multinom( sat ~ factor(delay), data = df.bus,weights=Freq)
m2.m <- multinom( sat ~ delay, data = df.bus,weights=Freq)
anova(m1.m,m2.m) ## Still the simple one

## Result
m2.m

## Plot the result
delay <- seq(0,7)
pred.m2 <- t(apply(predict(m2.m,,newdata=data.frame(delay=delay), type = "p"),1,
                   cumsum))
pred.p2 <- t(apply(predict(m2.p,newdata=data.frame(delay=delay), type = "p"),1,
                   cumsum))

## Example of calculation
pred.p2[1, ]
1/(1+exp(-m2.p$zeta))
pred.p2[2, ]
1/(1+exp(-m2.p$zeta+delay[2]*coef(m2.p)))

## using multinom
diff(pred.m2[1, ])
exp(coef(m2.m)[ ,1])/(1+sum(exp(coef(m2.m)[ ,1])))

diff(pred.m2[2, ])
exp(coef(m2.m)[ ,1]+coef(m2.m)[ ,2]*delay[2])/
    (1+sum(exp(coef(m2.m)[ ,1]+coef(m2.m)[ ,2]*delay[2])))
##################################################

matlines(delay,pred.m2[ ,-5],type="l",col=2)
matlines(delay,pred.p2[ ,-5],type="l",col=3)

## log-odds ration
log.odds <- function(p1,p2){
    log(p1/(1-p1)/(p2/(1-p2)))
}

## Proportional odds (Constant)
log.odds(pred.p2[1 ,-5],pred.p2[2 ,-5])
log.odds(pred.p2[2 ,-5],pred.p2[3 ,-5])

## Multinomial model log-odds ratio (Not constant)
log.odds(pred.m2[1 ,-5],pred.m2[2 ,-5])
log.odds(pred.m2[2 ,-5],pred.m2[3 ,-5])


## Plot in logit domain
par(mfrow=c(1,1))
delay=c(0,2,5,7)
## Model 
matplot(delay,-delay%o%c(1,1,1,1)*coef(m2.p)+
              matrix(m2.p$zeta,ncol=4,nrow=length(delay),byrow=T),
        type="l",col=1,xlab="delay",ylab="eta")
## Observations
matpoints(delay,log(accp[ ,-5]/(1-accp[ ,-5])),type="p",pch=19,col=1,cex=0.5)
## uncertainty
for(i in 1:length(t)){
    for(j in 1:4){
        ph <- sum(tab[i,1:j])/sum(tab[i, ])
        CI <- ph +c(-1,1)*2*sqrt(ph*(1-ph)/sum(tab[i, ]))
        matlines(t[i]*c(1,1),log(CI/(1-CI)))
    }
}



##################################################
## 5: Example 4.13 (Gamma) Slide 41
##################################################
var <- c(1.29,1.223,0.8248)
nsheets <- c(3,5,7)
wgt <- c(2,4,6)/2

fit0 <- glm(var~-1+nsheets,family=Gamma(link=inverse),weights=wgt)
summary(fit0)

par(mfrow=c(2,2))
plot(fit0) ## useless with 3 obs


## Predictions
nshet <- 1:10
pred.link <- predict(fit0, type= "link" ,
                     newdata=data.frame(nsheets=nshet),se=TRUE)
pred.resp <-predict(fit0, type= "response" ,
                    newdata=data.frame(nsheets=nshet),se=TRUE)

pred.link
pred.resp

##################################################
## end
##################################################
