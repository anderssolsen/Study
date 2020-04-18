####################################
y <- rnorm(10)
y

## estimate (assume sigma=1)
nl <- function(mu,y){
    - sum(dnorm(y,mean=mu,sd=1,log=TRUE))
}

nlminb(0,nl,y=y)

## estimate (sigma unknown)
nl <- function(theta,y){
    - sum(dnorm(y,mean=theta[1],sd=theta[2],log=TRUE))
}



nlminb(c(0,1),nl,y=y,lower=c(-Inf,0))
